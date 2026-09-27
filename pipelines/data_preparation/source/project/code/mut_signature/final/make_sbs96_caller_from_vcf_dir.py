#!/usr/bin/env python3
import argparse
import gzip
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

BASES = set("ACGT")


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def make_sbs96_index():
    subs = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]
    bases = ["A", "C", "G", "T"]
    return [
        f"{left}[{ref}>{alt}]{right}"
        for ref, alt in (s.split(">") for s in subs)
        for left in bases
        for right in bases
    ]


SBS96_INDEX = make_sbs96_index()
SBS96_SET = set(SBS96_INDEX)


def default_preferred_order():
    pta = [f"PTA_single_neuron_{i:02d}" for i in range(1, 17)]
    return [
        "Illumina",
        "ppmSeq_multiread",
        "ppmSeq_singleton_HC",
        # "ppmSeq_putative_multiread",
        *pta,
        "HiDEFseq_double_stranded",
        "HiDEFseq_single_stranded",
        "UDSeq",
        "NanoSeq",
        "AmpliSeq",
    ]


def normalize_chrom(name, fasta_refs):
    name = str(name)
    if name in fasta_refs:
        return name
    if name.startswith("chr") and name[3:] in fasta_refs:
        return name[3:]
    if (not name.startswith("chr")) and ("chr" + name) in fasta_refs:
        return "chr" + name
    return None


def snv_to_sbs96(chrom, pos1, ref, alt, fasta, fasta_refs, require_ref_match=True):
    ref = str(ref).upper()
    alt = str(alt).upper()
    if len(ref) != 1 or len(alt) != 1:
        return None, "not_snv"
    if ref not in BASES or alt not in BASES or ref == alt:
        return None, "non_acgt_or_ref_eq_alt"

    chrom2 = normalize_chrom(chrom, fasta_refs)
    if chrom2 is None:
        return None, "chrom_not_in_fasta"

    start0 = int(pos1) - 2
    end0 = int(pos1) + 1
    if start0 < 0:
        return None, "near_contig_start"

    try:
        tri = fasta.fetch(chrom2, start0, end0).upper()
    except Exception:
        return None, "fasta_fetch_failed"

    if len(tri) != 3 or any(b not in "ACGT" for b in tri):
        return None, "bad_trinuc"

    left, mid, right = tri[0], tri[1], tri[2]
    if require_ref_match and mid != ref:
        return None, "ref_mismatch_to_fasta"

    # Convert to pyrimidine-centered representation used by SBS96.
    if ref in {"A", "G"}:
        tri = revcomp(tri)
        left, mid, right = tri[0], tri[1], tri[2]
        ref = revcomp(ref)
        alt = revcomp(alt)

    key = f"{left}[{ref}>{alt}]{right}"
    if key not in SBS96_SET:
        return None, "not_in_sbs96"
    return key, None


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def find_vcfs(vcf_dir: Path, suffix: str):
    files = []
    if suffix:
        patterns = [f"*{suffix}", f"*{suffix}.gz"]
    else:
        patterns = ["*.vcf", "*.vcf.gz"]
    for pat in patterns:
        files.extend(sorted(p for p in vcf_dir.glob(pat) if p.is_file()))
    # Deduplicate and keep deterministic ordering.
    return sorted(set(files), key=lambda p: str(p))


def sample_name_from_file(path: Path, suffix: str):
    name = path.name
    if name.endswith(".gz"):
        name = name[:-3]
    if suffix and name.endswith(suffix):
        return name[: -len(suffix)]
    if name.endswith(".vcf"):
        return name[:-4]
    return name


def decide_sample_order(sample_names, explicit_order):
    if explicit_order:
        base_order = [x.strip() for x in explicit_order.split(",") if x.strip()]
    else:
        base_order = default_preferred_order()
    in_set = set(sample_names)
    ordered = [s for s in base_order if s in in_set]
    ordered.extend(sorted(s for s in sample_names if s not in set(ordered)))
    return ordered


def parse_one_vcf(
    vcf_path,
    sample_name,
    fasta,
    fasta_refs,
    pass_only=True,
    require_ref_match=True,
    dedup_within_file=True,
    collect_annotations=False,
):
    counts = Counter()
    reasons = Counter()
    seen = set() if dedup_within_file else None
    ann_rows = [] if collect_annotations else None

    with open_text(vcf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                reasons["malformed_vcf_line"] += 1
                continue

            chrom, pos_str, _id, ref, alt_field, qual, filt, info = parts[:8]
            if pass_only and filt not in {"PASS", "."}:
                reasons["filtered_non_pass"] += 1
                continue

            try:
                pos1 = int(pos_str)
            except ValueError:
                reasons["bad_pos"] += 1
                continue

            ref = ref.upper()
            for alt in alt_field.split(","):
                alt = alt.upper().strip()
                if not alt:
                    reasons["empty_alt"] += 1
                    continue

                key = (chrom, pos1, ref, alt)
                if dedup_within_file and key in seen:
                    reasons["duplicate_variant_in_file"] += 1
                    continue
                if dedup_within_file:
                    seen.add(key)

                sbs96, why = snv_to_sbs96(
                    chrom,
                    pos1,
                    ref,
                    alt,
                    fasta,
                    fasta_refs,
                    require_ref_match=require_ref_match,
                )
                if sbs96 is None:
                    reasons[why] += 1
                    continue

                counts[sbs96] += 1
                if collect_annotations:
                    ann_rows.append((sample_name, chrom, pos1, ref, alt, sbs96))

    return counts, reasons, ann_rows


def main():
    ap = argparse.ArgumentParser(
        description=(
            "Create a SigProfilerExtractor SBS96 matrix from per-platform VCF files "
            "(one VCF per platform/sample)."
        )
    )
    ap.add_argument("--vcf-dir", required=True, help="Directory containing platform VCF files")
    ap.add_argument("--fasta", required=True, help="Reference FASTA (with .fai)")
    ap.add_argument("-o", "--output", default="counts.SBS96.platforms.tsv", help="Output SBS96 matrix TSV")
    ap.add_argument(
        "--vcf-suffix",
        default=".vcf",
        help="Suffix that identifies input VCFs (default: .vcf). "
        "Use empty string to include all *.vcf/*.vcf.gz files.",
    )
    ap.add_argument(
        "--sample-order",
        default=None,
        help="Optional comma-separated sample order. Missing names are ignored; remaining samples are appended sorted.",
    )
    ap.add_argument("--allow-ref-mismatch", action="store_true", help="Do not require FASTA reference base to match REF")
    ap.add_argument(
        "--pass-only",
        action="store_true",
        default=True,
        help="Only count VCF records with FILTER PASS or . (default: on)",
    )
    ap.add_argument("--no-pass-only", dest="pass_only", action="store_false", help="Do not filter by FILTER column")
    ap.add_argument(
        "--dedup-within-file",
        action="store_true",
        default=True,
        help="Deduplicate duplicate (CHROM,POS,REF,ALT) within each input file (default: on)",
    )
    ap.add_argument("--no-dedup-within-file", dest="dedup_within_file", action="store_false")
    ap.add_argument(
        "--annotate-variants",
        default=None,
        help="Optional TSV path for long-form annotations: sample,chr,pos,ref,alt,SBS96",
    )
    args = ap.parse_args()

    vcf_dir = Path(args.vcf_dir)
    if not vcf_dir.exists() or not vcf_dir.is_dir():
        raise ValueError(f"--vcf-dir is not a directory: {vcf_dir}")

    vcf_files = find_vcfs(vcf_dir, args.vcf_suffix)
    if not vcf_files:
        raise ValueError(f"No VCF files found in {vcf_dir} using suffix '{args.vcf_suffix}'")

    sample_to_vcf = {}
    for path in vcf_files:
        sample = sample_name_from_file(path, args.vcf_suffix)
        if sample in sample_to_vcf:
            raise ValueError(
                f"Duplicate sample name '{sample}' from files:\n"
                f"  {sample_to_vcf[sample]}\n"
                f"  {path}\n"
                "Adjust --vcf-suffix or file names to keep sample names unique."
            )
        sample_to_vcf[sample] = path

    samples = decide_sample_order(list(sample_to_vcf.keys()), args.sample_order)

    fasta = pysam.FastaFile(args.fasta)
    fasta_refs = set(fasta.references)

    all_counts = {}
    all_reasons = {}
    ann_rows = []

    for sample in samples:
        vcf_path = sample_to_vcf[sample]
        counts, reasons, this_ann = parse_one_vcf(
            vcf_path,
            sample,
            fasta=fasta,
            fasta_refs=fasta_refs,
            pass_only=args.pass_only,
            require_ref_match=(not args.allow_ref_mismatch),
            dedup_within_file=args.dedup_within_file,
            collect_annotations=bool(args.annotate_variants),
        )
        all_counts[sample] = counts
        all_reasons[sample] = reasons
        if this_ann:
            ann_rows.extend(this_ann)

    out = pd.DataFrame({"MutationType": SBS96_INDEX})
    for sample in samples:
        cts = all_counts[sample]
        out[sample] = np.array([cts[k] for k in SBS96_INDEX], dtype=int)
    out.to_csv(args.output, sep="\t", index=False)

    if args.annotate_variants:
        ann = pd.DataFrame(ann_rows, columns=["sample", "chr", "pos", "ref", "alt", "SBS96"])
        ann.to_csv(args.annotate_variants, sep="\t", index=False)

    totals = {s: int(out[s].sum()) for s in samples}
    print(f"Wrote {args.output}", file=sys.stderr)
    print(f"Platforms/samples: {len(samples)}", file=sys.stderr)
    print(f"Totals per sample: {totals}", file=sys.stderr)
    for s in samples:
        if all_reasons[s]:
            print(f"Skipped reasons for {s}: {dict(all_reasons[s])}", file=sys.stderr)
    if args.annotate_variants:
        print(f"Wrote {args.annotate_variants}", file=sys.stderr)


if __name__ == "__main__":
    main()


# Example:
# python ./work/project/code/mut_signature/final/make_sbs96_caller_from_vcf_dir.py \
#   --vcf-dir ./work/regional_controls/08.final_tsv/01.filtered_vcf \
#   --fasta ./resources/hg38/Homo_sapiens_assembly38.fasta \
#   --output ./work/project/data/mutational_signatures/sigprofiler/caller.counts.SBS96.platforms.from_vcfs.tsv \
#   --annotate-variants ./work/project/data/mutational_signatures/sigprofiler/caller.variants.with_SBS96.tsv

# python -c "import pandas as pd; p='./work/project/data/mutational_signatures/sigprofiler/caller.counts.SBS96.platforms.from_vcfs.tsv'; o='./work/project/data/mutational_signatures/sigprofiler/caller.counts.SBS96.platforms.from_vcfs.noAmpliSeq.tsv'; df=pd.read_csv(p,sep='\t'); assert 'AmpliSeq' in df.columns, 'AmpliSeq column not found'; df=df.drop(columns=['AmpliSeq']); df.to_csv(o,sep='\t',index=False); print('Wrote:',o,'shape=',df.shape)"
