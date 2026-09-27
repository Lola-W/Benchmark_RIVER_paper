#!/usr/bin/env python3
import argparse
import pandas as pd
import numpy as np
import pyBigWig
import pysam

# ---------- bigWigs ----------
DEFAULT_BW = {
  "segdup": ("./resources/hg38/bw_mask/hg38.segdup.bw", 1,  "gte0.1"),
  "gaps":   ("./resources/hg38/bw_mask/hg38_UCSC_gaps.bw", 1, "gte0.1"),
  "sat":    ("./resources/hg38/bw_mask/satellites_rDNA.bw", 1, "gte0.1"),
  "umap":   ("./resources/hg38/k50.Umap.MultiTrackMappability.bw", 20, "lt0.4"),
  "encode": ("./resources/hg38/bw_mask/ENCODE_hg38-blacklist.v2.bw", 20, "gte0.1"),
}

DEFAULT_GNOMAD_VCF = "./resources/hg38/af-only-gnomad.hg38.vcf.gz"
DEFAULT_KNOWN_INDEL_GATK = "./resources/hg38/Homo_sapiens_assembly38.known_indels.vcf.gz"
DEFAULT_KNOWN_INDEL_MILLS = "./resources/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

BIT_ORDER = ["segdup", "gaps", "satellites_rDNA", "low_mappability", "encode_blacklist"]

def parse_threshold(s: str):
    s = s.strip().lower()
    for op in ("lte", "gte", "lt", "gt"):
        if s.startswith(op):
            val = float(s[len(op):])
            if op == "lt":  return lambda x: x <  val
            if op == "lte": return lambda x: x <= val
            if op == "gt":  return lambda x: x >  val
            if op == "gte": return lambda x: x >= val
    raise ValueError(f"Bad threshold format: {s} (expected lt0.4 / gte0.1 / etc)")

def normalize_chrom(chrom: str, chromset):
    if chrom in chromset:
        return chrom
    if chrom.startswith("chr"):
        c2 = chrom[3:]
        if c2 in chromset:
            return c2
    else:
        c2 = "chr" + chrom
        if c2 in chromset:
            return c2
    return None  # <-- IMPORTANT: signal “not found”

def ensure_locus_columns(df: pd.DataFrame) -> pd.DataFrame:
    # Accept either:
    #  (A) chr,pos,ref,alt
    #  (B) variant_id like "1-10022248-C-T" or "chr1-10022248-C-T"
    cols = set(df.columns)

    has_locus = {"chr","pos","ref","alt"}.issubset(cols)
    if has_locus:
        df["chr"] = df["chr"].astype(str)
        df["pos"] = df["pos"].astype(int)
        df["ref"] = df["ref"].astype(str)
        df["alt"] = df["alt"].astype(str)
        return df

    if "variant_id" not in cols:
        raise ValueError("Input must contain either chr/pos/ref/alt OR variant_id")

    parts = df["variant_id"].astype(str).str.split("-", n=3, expand=True)
    if parts.shape[1] < 4:
        raise ValueError("variant_id must look like '1-123-A-G' (4 fields split by '-')")

    df["chr"] = parts[0].astype(str)
    df.loc[~df["chr"].str.startswith("chr"), "chr"] = "chr" + df.loc[~df["chr"].str.startswith("chr"), "chr"]
    df["pos"] = parts[1].astype(int)
    df["ref"] = parts[2].astype(str)
    df["alt"] = parts[3].astype(str)
    return df


def bw_mean(bw, chrom: str, pos1: int, binsize: int):
    if chrom is None:
        return 0.0

    chroms = bw.chroms()
    if chrom not in chroms:
        return 0.0

    clen = chroms[chrom]  # chromosome length in bigWig
    pos0 = pos1 - 1

    # out of bounds position
    if pos0 < 0 or pos0 >= clen:
        return 0.0

    if binsize <= 1:
        start, end = pos0, pos0 + 1
    else:
        half = binsize // 2
        start = pos0 - half
        end = start + binsize

    # clamp to [0, clen]
    start = max(0, start)
    end = min(clen, end)

    # ensure valid interval
    if end <= start:
        return 0.0

    try:
        v = bw.stats(chrom, start, end, type="mean")[0]
    except RuntimeError:
        return 0.0

    return 0.0 if v is None or np.isnan(v) else float(v)


def vcf_match_flag(vcf: pysam.VariantFile, chrom: str, pos1: int, ref: str, alt: str):
    if chrom is None:
        return 0
    try:
        recs = vcf.fetch(chrom, pos1 - 1, pos1)
    except Exception:
        return 0
    for r in recs:
        if r.pos == pos1 and r.ref == ref and alt in list(r.alts or []):
            return 1
    return 0


def gnomad_af_at_site(vcf: pysam.VariantFile, chrom: str, pos1: int, ref: str, alt: str):
    if chrom is None:
        return 0, 0.0
    try:
        recs = vcf.fetch(chrom, pos1 - 1, pos1)
    except Exception:
        return 0, 0.0
    for r in recs:
        if r.pos != pos1 or r.ref != ref:
            continue
        alts = list(r.alts or [])
        if alt not in alts:
            continue
        af = r.info.get("AF", None)
        if af is None:
            return 1, 0.0
        if isinstance(af, (list, tuple)):
            return 1, float(af[alts.index(alt)] or 0.0)
        return 1, float(af)
    return 0, 0.0


def repeatmasker_lookup(tb: pysam.TabixFile, chrom: str, pos1: int):
    if chrom is None:
        return 0, "", "", ""
    start, end = pos1 - 1, pos1
    try:
        hits = tb.fetch(chrom, start, end)
    except Exception:
        return 0, "", "", ""
    rep_class, rep_family, rep_name = [], [], []
    for line in hits:
        f = line.rstrip("\n").split("\t")
        if len(f) >= 6:
            rep_class.append(f[3])
            rep_family.append(f[4])
            rep_name.append(f[5])
    if not rep_class:
        return 0, "", "", ""
    def uniq_join(x): return ";".join(pd.unique(pd.Series(x)))
    return 1, uniq_join(rep_class), uniq_join(rep_family), uniq_join(rep_name)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i","--input", required=True, help="Input TSV with chr,pos,ref,alt")
    ap.add_argument("-o","--output", required=True)
    ap.add_argument("--gnomad-vcf", default=DEFAULT_GNOMAD_VCF)
    ap.add_argument("--known-indel-gatk", default=DEFAULT_KNOWN_INDEL_GATK)
    ap.add_argument("--known-indel-mills", default=DEFAULT_KNOWN_INDEL_MILLS)
    ap.add_argument("--germline-vcf", required=True, help="Your matched-normal germline VCF.gz")
    ap.add_argument("--repeatmasker-bed-gz", required=True, help="Tabix BED.gz: chrom start end repClass repFamily repName")
    ap.add_argument("--chunksize", type=int, default=200000)
    args = ap.parse_args()

    # open bigWigs
    bws, bw_chromsets, thr, bins = {}, {}, {}, {}
    for key, (path, binsize, threshold) in DEFAULT_BW.items():
        bw = pyBigWig.open(path)
        bws[key] = bw
        bw_chromsets[key] = set(bw.chroms().keys())
        thr[key] = parse_threshold(threshold)
        bins[key] = int(binsize)

    # open VCFs
    vcf_gnomad = pysam.VariantFile(args.gnomad_vcf)
    vcf_gatk   = pysam.VariantFile(args.known_indel_gatk)
    vcf_mills  = pysam.VariantFile(args.known_indel_mills)
    vcf_germ   = pysam.VariantFile(args.germline_vcf)

    chromset_gnomad = set(vcf_gnomad.header.contigs)
    chromset_gatk   = set(vcf_gatk.header.contigs)
    chromset_mills  = set(vcf_mills.header.contigs)
    chromset_germ   = set(vcf_germ.header.contigs)

    # repeatmasker tabix
    tb_rmsk = pysam.TabixFile(args.repeatmasker_bed_gz)
    # tabix contigs may be in header
    chromset_rmsk = set(tb_rmsk.contigs)

    first = True
    with open(args.output, "w") as out_fh:
        for chunk in pd.read_csv(args.input, sep="\t", dtype=str, chunksize=args.chunksize):
            chunk = ensure_locus_columns(chunk)
            
            n = len(chunk)
            in_segdup = np.zeros(n, np.int8)
            in_gaps   = np.zeros(n, np.int8)
            in_sat    = np.zeros(n, np.int8)
            low_map   = np.zeros(n, np.int8)
            in_encode = np.zeros(n, np.int8)

            g_present = np.zeros(n, np.int8)
            g_af      = np.zeros(n, float)

            indel_gatk  = np.zeros(n, np.int8)
            indel_mills = np.zeros(n, np.int8)

            is_germline = np.zeros(n, np.int8)

            in_rmsk = np.zeros(n, np.int8)
            rmsk_class  = np.array([""]*n, dtype=object)
            rmsk_family = np.array([""]*n, dtype=object)
            rmsk_name   = np.array([""]*n, dtype=object)

            for i, r in enumerate(chunk.itertuples(index=False)):
                chrom_in = getattr(r, "chr")
                pos1 = int(getattr(r, "pos"))
                ref = getattr(r, "ref")
                alt = getattr(r, "alt")

                # ---------- blacklist bigWigs ----------
                c_seg = normalize_chrom(chrom_in, bw_chromsets["segdup"])
                c_gap = normalize_chrom(chrom_in, bw_chromsets["gaps"])
                c_sat = normalize_chrom(chrom_in, bw_chromsets["sat"])
                c_uma = normalize_chrom(chrom_in, bw_chromsets["umap"])
                c_enc = normalize_chrom(chrom_in, bw_chromsets["encode"])

                segv = bw_mean(bws["segdup"], c_seg, pos1, bins["segdup"])
                gapv = bw_mean(bws["gaps"],   c_gap, pos1, bins["gaps"])
                satv = bw_mean(bws["sat"],    c_sat, pos1, bins["sat"])
                umav = bw_mean(bws["umap"],   c_uma, pos1, bins["umap"])
                encv = bw_mean(bws["encode"], c_enc, pos1, bins["encode"])

                in_segdup[i] = 1 if thr["segdup"](segv) else 0
                in_gaps[i]   = 1 if thr["gaps"](gapv) else 0
                in_sat[i]    = 1 if thr["sat"](satv) else 0
                low_map[i]   = 1 if thr["umap"](umav) else 0
                in_encode[i] = 1 if thr["encode"](encv) else 0

                # ---------- gnomAD AF ----------
                c_gnom = normalize_chrom(chrom_in, chromset_gnomad)
                present, afv = gnomad_af_at_site(vcf_gnomad, c_gnom, pos1, ref, alt)
                g_present[i] = present
                g_af[i] = afv

                # ---------- known indel sites ----------
                c_gatk  = normalize_chrom(chrom_in, chromset_gatk)
                c_mills = normalize_chrom(chrom_in, chromset_mills)
                indel_gatk[i]  = vcf_match_flag(vcf_gatk,  c_gatk,  pos1, ref, alt)
                indel_mills[i] = vcf_match_flag(vcf_mills, c_mills, pos1, ref, alt)

                # ---------- germline ----------
                c_germ = normalize_chrom(chrom_in, chromset_germ)
                is_germline[i] = vcf_match_flag(vcf_germ, c_germ, pos1, ref, alt)

                # ---------- repeatmasker ----------
                c_rmsk = normalize_chrom(chrom_in, chromset_rmsk)
                hit, cl, fa, na = repeatmasker_lookup(tb_rmsk, c_rmsk, pos1)
                in_rmsk[i] = hit
                rmsk_class[i] = cl
                rmsk_family[i] = fa
                rmsk_name[i] = na

            # attach columns
            chunk["in_segdup"] = in_segdup
            chunk["in_ucsc_gaps"] = in_gaps
            chunk["in_satellites_rDNA"] = in_sat
            chunk["low_mappability"] = low_map
            chunk["in_encode_blacklist"] = in_encode

            # Use numpy arrays
            b0 = chunk["in_segdup"].to_numpy(dtype=np.int64)
            b1 = chunk["in_ucsc_gaps"].to_numpy(dtype=np.int64)
            b2 = chunk["in_satellites_rDNA"].to_numpy(dtype=np.int64)
            b3 = chunk["low_mappability"].to_numpy(dtype=np.int64)
            b4 = chunk["in_encode_blacklist"].to_numpy(dtype=np.int64)

            bitmask = (b0 * 1) | (b1 * 2) | (b2 * 4) | (b3 * 8) | (b4 * 16)

            chunk["blacklist_any"] = (bitmask > 0).astype(np.int8)
            chunk["blacklist_bitmask"] = bitmask.astype(np.int64)
            chunk["blacklist_binary"] = [format(int(x), "05b") for x in bitmask]


            chunk["gnomad_present"] = g_present
            chunk["gnomad_af"] = g_af

            chunk["in_known_indel_gatk"] = indel_gatk
            chunk["in_known_indel_mills1000G"] = indel_mills
            chunk["in_known_indel_any"] = ((indel_gatk + indel_mills) > 0).astype(np.int8)

            chunk["is_germline"] = is_germline

            chunk["in_repeatmasker"] = in_rmsk
            chunk["repeatmasker_class"] = rmsk_class
            chunk["repeatmasker_family"] = rmsk_family
            chunk["repeatmasker_name"] = rmsk_name

            chunk.to_csv(out_fh, sep="\t", index=False, header=first)
            first = False

    # close handles
    for bw in bws.values():
        bw.close()
    vcf_gnomad.close(); vcf_gatk.close(); vcf_mills.close(); vcf_germ.close()
    tb_rmsk.close()

if __name__ == "__main__":
    main()

# python ./work/project/code/annotation/annotate_blacklist_and_gnomad.py   -i ./work/project/data/upset/all_ids_df.tsv   -o ./work/project/data/annotation/all_ids_df.annotated.full.tsv   --germline-vcf ./work/project/data/germline_HC_VQSR/41411195_F_illumina_matchedNormal.haplotypecaller.SNV_INDEL.VQSR_INFERRED.vcf.gz   --repeatmasker-bed-gz ./resources/hg38/hg38.rmsk.class.bed.gz