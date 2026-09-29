#!/usr/bin/env python3

import argparse
import csv
import gzip
import os
import subprocess
import sys
import tempfile
from collections import defaultdict

FLAG_COLS = {"is_cycle_skip", "is_forward"}  # removed unused is_duplicate

ALLOWED_CONTIGS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

NEW_INFO_DEFS = [
    '##INFO=<ID=triN,Number=1,Type=String,Description="Normalized trinucleotide context (MutationalPatterns-style), e.g., T[T>C]C">',
    '##INFO=<ID=prev_3bp,Number=1,Type=String,Description="3bp context on the 5-prime side of the variant">',
    '##INFO=<ID=next_3bp,Number=1,Type=String,Description="3bp context on the 3-prime side of the variant">',
    '##INFO=<ID=best_threshold,Number=1,Type=Float,Description="Best trinucleotide-denoising threshold selected for this record">',
    '##INFO=<ID=effective_thr,Number=1,Type=Float,Description="Effective threshold applied for trinucleotide-denoising for this record">',
    '##INFO=<ID=DUP_COUNT_FILTERED,Number=1,Type=Integer,Description="Number of TSV rows for the same CHROM:POS:REF:ALT before dedup (after trinuc denoising)">',
    '##INFO=<ID=DUP_COUNT_RAW_MULTI_ALLELE,Number=1,Type=Integer,Description="In outMap.vcf.gz: number of records sharing the same CHROM:POS (REF/ALT ignored)">',
    '##INFO=<ID=DUP_COUNT_RAW,Number=1,Type=Integer,Description="In outMap.vcf.gz: number of records sharing the same CHROM:POS:REF:ALT">',
    '##INFO=<ID=DUP_COUNT_SNVQ40,Number=1,Type=Integer,Description="In outMap.vcf.gz: among DUP_COUNT_RAW, number of records with FILTER==PASS">',
]


def run_cmd(cmd: list) -> None:
    """Run external command, raise on failure."""
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        sys.stderr.write(f"[ERROR] Command failed: {' '.join(cmd)}\n{p.stderr}\n")
        raise SystemExit(p.returncode)


def read_vcf_header(vcf_gz: str) -> list:
    """Read VCF header lines up to and including #CHROM line."""
    lines = []
    with gzip.open(vcf_gz, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                lines.append(line)
                if line.startswith("#CHROM"):
                    break
            else:
                break
    if not lines or not lines[-1].startswith("#CHROM"):
        raise RuntimeError(f"Could not find #CHROM header in {vcf_gz}")
    return lines


def inject_new_info_defs(header_lines: list) -> list:
    """Inject new INFO definitions before #CHROM line if not already present."""
    existing = {ln.split("##INFO=<ID=", 1)[1].split(",", 1)[0]
                for ln in header_lines if ln.startswith("##INFO=<ID=")}
    out = []
    for ln in header_lines:
        if ln.startswith("#CHROM"):
            for new_def in NEW_INFO_DEFS:
                new_id = new_def.split("##INFO=<ID=", 1)[1].split(",", 1)[0]
                if new_id not in existing:
                    out.append(new_def)
        out.append(ln)
    return out


def build_info_from_row(row: dict) -> str:
    """Convert TSV row dict to VCF INFO string."""
    exclude = {"chrom", "pos", "ref", "alt", "qual", "filt", "sample_id"}
    parts = []
    for k, v in row.items():
        if k in exclude or v is None:
            continue
        v = str(v).strip()
        if v in ("", "NA"):
            continue
        if k in FLAG_COLS:
            if v in ("1", "true", "t", "yes"):
                parts.append(k)
        else:
            parts.append(f"{k}={v}")
    return ";".join(parts)


def bgzip_and_tabix(vcf_path: str, vcf_gz_path: str) -> None:
    """Compress VCF with bgzip and index with tabix."""
    with open(vcf_gz_path, "wb") as out_f:
        p = subprocess.run(["bgzip", "-c", vcf_path], stdout=out_f, stderr=subprocess.PIPE)
    if p.returncode != 0:
        raise SystemExit(f"[ERROR] bgzip failed: {p.stderr.decode()}")
    run_cmd(["tabix", "-f", "-p", "vcf", vcf_gz_path])


def query_outmap_counts(
    outmap_vcfgz: str,
    variant_keys: list,
    threads: int,
    progress_every: int,
) -> tuple:
    """Query outMap VCF for DUP_COUNT annotations."""
    # ensure tabix index exists
    if not (os.path.exists(outmap_vcfgz + ".tbi") or os.path.exists(outmap_vcfgz + ".csi")):
        sys.stderr.write(f"[INFO] indexing: {outmap_vcfgz}\n")
        run_cmd(["bcftools", "index", "-t", outmap_vcfgz])

    want_key, want_pos, seen_pos, bed_rows = set(), set(), set(), []
    for chrom, pos, ref, alt in variant_keys:
        poskey, key = f"{chrom}:{pos}", f"{chrom}:{pos}:{ref}:{alt}"
        want_pos.add(poskey)
        want_key.add(key)
        if poskey not in seen_pos:
            seen_pos.add(poskey)
            bed_rows.append((chrom, pos - 1, pos))

    raw_multi, raw, pas = defaultdict(int), defaultdict(int), defaultdict(int)
    if not bed_rows:
        return raw_multi, raw, pas

    with tempfile.TemporaryDirectory() as td:
        bed_path = os.path.join(td, "pos.bed")
        with open(bed_path, "wt") as bed:
            bed.writelines(f"{c}\t{s}\t{e}\n" for c, s, e in bed_rows)

        sys.stderr.write(f"[INFO] outMap query: unique_positions={len(bed_rows):,}, unique_keys={len(want_key):,}\n")

        cmd = ["bcftools", "view", "-H", "--threads", str(threads), "-R", bed_path, outmap_vcfgz]
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)

        for n_lines, line in enumerate(p.stdout, 1):
            cols = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt, filt = cols[0], cols[1], cols[3], cols[4], cols[6]
            poskey, key = f"{chrom}:{pos}", f"{chrom}:{pos}:{ref}:{alt}"
            if poskey in want_pos:
                raw_multi[poskey] += 1
            if key in want_key:
                raw[key] += 1
                if filt == "PASS":
                    pas[key] += 1
            if progress_every > 0 and n_lines % progress_every == 0:
                sys.stderr.write(f"[PROGRESS] outMap_lines={n_lines:,}\n")

        p.stdout.close()
        if p.wait() != 0:
            raise SystemExit("[ERROR] bcftools view -R outMap failed\n" + p.stderr.read())

    return raw_multi, raw, pas


def write_outputs(header_lines, records, raw_multi, raw, pas,
                  out_vcf_gz, out_multi_gz, out_single_gz,
                  out_true_single_gz, out_true_multi_gz, true_multi_raw_max):
    """Write all output VCF files."""
    with tempfile.TemporaryDirectory() as td:
        paths = {
            "all":         os.path.join(td, "out.vcf"),
            "multi":       os.path.join(td, "out.multi.vcf")        if out_multi_gz else None,
            "single":      os.path.join(td, "out.single.vcf")       if out_single_gz else None,
            "true_single": os.path.join(td, "out.true_single.vcf")  if out_true_single_gz else None,
            "true_multi":  os.path.join(td, "out.true_multi.vcf")   if out_true_multi_gz else None,
        }

        handles = {k: open(v, "wt") for k, v in paths.items() if v}

        try:
            for f in handles.values():
                for ln in header_lines:
                    f.write(ln + "\n")

            for chrom, pos, ref, alt, qual, info0, dupcount in records:
                poskey, key = f"{chrom}:{pos}", f"{chrom}:{pos}:{ref}:{alt}"
                dc_raw_multi = int(raw_multi.get(poskey, 0))
                dc_raw       = int(raw.get(key, 0))
                dc_pass      = int(pas.get(key, 0))

                info = (f"{info0};" if info0 else "") + \
                       f"DUP_COUNT_FILTERED={dupcount};DUP_COUNT_RAW_MULTI_ALLELE={dc_raw_multi};DUP_COUNT_RAW={dc_raw};DUP_COUNT_SNVQ40={dc_pass}"
                line = f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{qual:.2f}\tPASS\t{info}\n"

                handles["all"].write(line)
                if "multi" in handles and dupcount >= 2:
                    handles["multi"].write(line)
                if "single" in handles and dupcount == 1:
                    handles["single"].write(line)
                if "true_single" in handles and dupcount == 1 and dc_raw_multi == 1:
                    handles["true_single"].write(line)
                if "true_multi" in handles and dupcount >= 2 and dc_raw <= int(true_multi_raw_max):
                    handles["true_multi"].write(line)

        finally:
            for f in handles.values():
                f.close()

        bgzip_and_tabix(paths["all"], out_vcf_gz)
        for src, dst in [(paths["multi"], out_multi_gz), (paths["single"], out_single_gz),
                         (paths["true_single"], out_true_single_gz), (paths["true_multi"], out_true_multi_gz)]:
            if src and dst:
                bgzip_and_tabix(src, dst)


def main():
    ap = argparse.ArgumentParser(
        description="TSV -> VCF (dedup by best QUAL) + outMap dupcount annotations + singleton/multi outputs"
    )
    ap.add_argument("-p", "--previous_vcf",       required=True,  help="Previous VCF (.vcf.gz) to copy header from")
    ap.add_argument("-i", "--input_tsv",           required=True,  help="Input trinuc_denoised TSV")
    ap.add_argument("-o", "--output_vcf",          required=True,  help="Output VCF (.vcf.gz)")
    ap.add_argument("-m", "--output_multiread",    default=None,   help="Output VCF for DUP_COUNT_FILTERED>=2")
    ap.add_argument("-s", "--output_singleton",    default=None,   help="Output VCF for DUP_COUNT_FILTERED==1")
    ap.add_argument("-t", "--output_true_singleton",default=None,  help="Output VCF for DUP_COUNT_FILTERED==1 AND DUP_COUNT_RAW_MULTI_ALLELE==1")
    ap.add_argument("-u", "--output_true_multiread",default=None,  help="Output VCF for DUP_COUNT_FILTERED>=2 AND DUP_COUNT_RAW<=true_multi_raw_max")
    ap.add_argument("--true_multi_raw_max",        type=int, default=None)
    ap.add_argument("--outmap",                    required=True,  help="outMap VCF (.vcf.gz) for dupcount annotation")
    ap.add_argument("--threads",                   type=int, default=1)
    ap.add_argument("--progress_every",            type=int, default=1_000)
    args = ap.parse_args()

    if args.output_true_multiread and args.true_multi_raw_max is None:
        raise SystemExit("[ERROR] -u requires --true_multi_raw_max")

    header = inject_new_info_defs(read_vcf_header(args.previous_vcf))

    best, dup, order_keys = {}, {}, []
    alts_per_site = defaultdict(set)
    n_rows = n_skip = 0

    sys.stderr.write(f"[INFO] reading TSV: {args.input_tsv}\n")
    with open(args.input_tsv, "rt", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            n_rows += 1
            chrom = row["chrom"].strip()
            if chrom not in ALLOWED_CONTIGS:
                n_skip += 1
                continue
            pos, ref, alt, qual = int(row["pos"]), row["ref"].strip(), row["alt"].strip(), float(row["qual"])
            key = (chrom, pos, ref, alt)
            alts_per_site[(chrom, pos, ref)].add(alt)
            dup[key] = dup.get(key, 0) + 1
            if key not in best:
                order_keys.append(key)
            info0 = build_info_from_row(row)
            if key not in best or qual > best[key][0]:
                best[key] = (qual, info0)
            if n_rows % 1_000_000 == 0:
                sys.stderr.write(f"[INFO] TSV rows processed={n_rows:,}\n")

    multi_sites = {site for site, alts in alts_per_site.items() if len(alts) >= 2}
    records, n_drop = [], 0
    for chrom, pos, ref, alt in order_keys:
        if (chrom, pos, ref) in multi_sites:
            n_drop += 1
            continue
        qual, info0 = best[(chrom, pos, ref, alt)]
        records.append((chrom, pos, ref, alt, qual, info0, dup[(chrom, pos, ref, alt)]))

    sys.stderr.write(f"[INFO] input_rows={n_rows:,}, skipped={n_skip:,}, multi_allelic_dropped={n_drop:,}, variants={len(records):,}\n")

    variant_keys = [(c, p, r, a) for c, p, r, a, _, __, ___ in records]
    raw_multi, raw, pas = query_outmap_counts(
        args.outmap, variant_keys, args.threads, args.progress_every
    )

    write_outputs(
        header, records, raw_multi, raw, pas,
        args.output_vcf, args.output_multiread, args.output_singleton,
        args.output_true_singleton, args.output_true_multiread, args.true_multi_raw_max
    )

    # Summary
    n_multi  = sum(1 for r in records if r[-1] >= 2)
    n_single = sum(1 for r in records if r[-1] == 1)
    n_true_single = sum(1 for c, p, r, a, _, __, dc in records
                        if dc == 1 and int(raw_multi.get(f"{c}:{p}", 0)) == 1)

    sys.stderr.write(f"[OK] wrote: {args.output_vcf}\n")
    sys.stderr.write(f"[INFO] n_variants={len(records):,} (singleton={n_single:,}, multiread={n_multi:,}, true_singleton={n_true_single:,})\n")


if __name__ == "__main__":
    main()
