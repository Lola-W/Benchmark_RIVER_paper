#!/usr/bin/env python3
import argparse, re, sys
from collections import Counter
import numpy as np
import pandas as pd
import pysam
import os

BASES = set("ACGT")
DEPTHS = ["30x", "100x", "200x"]
DEPTH_ANY_RE = re.compile(r"(?:^|_)(30x|100x|200x)(?:_|$)", flags=re.IGNORECASE)

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

def make_sbs96_index():
    subs = ["C>A","C>G","C>T","T>A","T>C","T>G"]
    bases = ["A","C","G","T"]
    return [f"{l}[{r}>{a}]{rr}" for r,a in (s.split(">") for s in subs) for l in bases for rr in bases]

SBS96_INDEX = make_sbs96_index()
SBS96_SET = set(SBS96_INDEX)

def normalize_chrom(name, fasta_refs):
    name = str(name)
    if name in fasta_refs: return name
    if name.startswith("chr") and name[3:] in fasta_refs: return name[3:]
    if (not name.startswith("chr")) and ("chr"+name) in fasta_refs: return "chr"+name
    return None

def snv_to_sbs96(chrom, pos1, ref, alt, fasta, fasta_refs, require_ref_match=True):
    ref, alt = str(ref).upper(), str(alt).upper()
    if len(ref) != 1 or len(alt) != 1: return None, "not_snv"
    if ref not in BASES or alt not in BASES or ref == alt: return None, "non_acgt_or_ref_eq_alt"

    chrom2 = normalize_chrom(chrom, fasta_refs)
    if chrom2 is None: return None, "chrom_not_in_fasta"

    start0, end0 = int(pos1) - 2, int(pos1) + 1
    if start0 < 0: return None, "near_contig_start"

    try:
        tri = fasta.fetch(chrom2, start0, end0).upper()
    except Exception:
        return None, "fasta_fetch_failed"

    if len(tri) != 3 or any(b not in "ACGT" for b in tri): return None, "bad_trinuc"
    left, mid, right = tri[0], tri[1], tri[2]
    if require_ref_match and mid != ref: return None, "ref_mismatch_to_fasta"

    # pyrimidine-center representation
    if ref in {"A","G"}:
        tri = revcomp(tri)
        left, mid, right = tri[0], tri[1], tri[2]
        ref, alt = revcomp(ref), revcomp(alt)

    key = f"{left}[{ref}>{alt}]{right}"
    if key not in SBS96_SET: return None, "not_in_sbs96"
    return key, None

def coerce01(s):
    return pd.to_numeric(s, errors="coerce").fillna(0).astype(int).clip(0, 1)

def has_depth_any(s): return bool(DEPTH_ANY_RE.search(s))
def has_depth(s, d):  return bool(re.search(rf"(?:^|_){re.escape(d)}(?:_|$)", s, flags=re.IGNORECASE))

def base_key(s):
    out = re.sub(r"_(30x|100x|200x)_", "_", s, flags=re.IGNORECASE)
    out = re.sub(r"_(30x|100x|200x)$", "", out, flags=re.IGNORECASE)
    out = re.sub(r"__+", "_", out)
    return out

def pick_preferred(cands, d):
    def score(x):
        # prefer exact suffix _{d} (case-insensitive), then shorter names for determinism
        return (0 if re.search(rf"_{re.escape(d)}$", x, flags=re.IGNORECASE) else 1, len(x), x)
    return sorted(cands, key=score)[0] if cands else None

def derive_sample_name_from_col(col):
    # strip common prefixes; keep the rest as sample name
    for pref in ["is_existing_lower_ci_0001_", "is_existing_"]:
        if col.startswith(pref):
            col = col[len(pref):]
            break
    # normalize X casing: _30X -> _30x
    col = re.sub(r"_(30x|100x|200x)$", lambda m: "_" + m.group(1).lower(), col, flags=re.IGNORECASE)
    col = re.sub(r"_(30x|100x|200x)_", lambda m: "_" + m.group(1).lower() + "_", col, flags=re.IGNORECASE)
    return col

def build_selected(samples_set):
    samples = set(samples_set)
    bases = sorted({base_key(s) for s in samples})

    sel = {}
    sel["no_depth"] = {s for s in samples if not has_depth_any(s)}

    for d in DEPTHS:
        chosen = set()
        for b in bases:
            cands_d = [s for s in samples if base_key(s) == b and has_depth(s, d)]
            if cands_d:
                chosen.add(pick_preferred(cands_d, d))
            else:
                if b in samples:
                    chosen.add(b)
                else:
                    cands_nd = [s for s in samples if base_key(s) == b and not has_depth_any(s)]
                    if cands_nd:
                        chosen.add(sorted(cands_nd)[0])
        sel[d] = chosen
    return sel

def main():
    ap = argparse.ArgumentParser(
        description="Build SBS96 matrices from variant TSV using is_existing* binary flags; outputs no_depth/30x/100x/200x."
    )
    ap.add_argument("--in-tsv", required=True)
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out-prefix", required=True,
                    help="Prefix for outputs (directory is created). Writes: <prefix>.no_depth.tsv/.30x.tsv/.100x.tsv/.200x.tsv")

    ap.add_argument("--chr-col", default="chr")
    ap.add_argument("--pos-col", default="pos")
    ap.add_argument("--ref-col", default="ref")
    ap.add_argument("--alt-col", default="alt")

    ap.add_argument("--flag-regex", default=r"^is_existing",
                    help="Regex to select binary flag columns (default: '^is_existing').")
    ap.add_argument("--exclude-flag-cols", default="",
                    help="Comma-separated flag columns to exclude even if they match --flag-regex.")
    ap.add_argument("--germline-cols", default="is_germline",
                    help="Comma-separated germline flag columns; rows with any==1 are removed (if present).")
    ap.add_argument("--allow-ref-mismatch", action="store_true",
                    help="Do not require FASTA base == REF at POS")
    ap.add_argument("--dedup", action="store_true", default=True,
                    help="Deduplicate on (chr,pos,ref,alt) by taking max over binary flags (default: on)")
    ap.add_argument("--no-dedup", dest="dedup", action="store_false")

    ap.add_argument("--annotate-variants-prefix", default=None,
                    help="Optional prefix; writes <prefix>.variants.with_sbs96.tsv with _sbs96 + chosen flag cols")
    ap.add_argument("--include-flag-cols", default="",
                help="Comma-separated extra binary flag columns to include even if they don't match --flag-regex.")
    args = ap.parse_args()

    df = pd.read_csv(args.in_tsv, sep="\t", dtype=str, low_memory=False)

    for c in [args.chr_col, args.pos_col, args.ref_col, args.alt_col]:
        if c not in df.columns:
            raise ValueError(f"Missing required column {c}. Present: {list(df.columns)[:50]} ...")

    df[args.pos_col] = pd.to_numeric(df[args.pos_col], errors="coerce").astype("Int64")
    df[args.ref_col] = df[args.ref_col].astype(str).str.upper()
    df[args.alt_col] = df[args.alt_col].astype(str).str.upper()

    # choose flag columns
    # choose flag columns by regex + explicitly included columns
    flag_re = re.compile(args.flag_regex)
    exclude = {c.strip() for c in args.exclude_flag_cols.split(",") if c.strip()}
    flag_cols = [c for c in df.columns if flag_re.search(c) and c not in exclude]

    include = [c.strip() for c in args.include_flag_cols.split(",") if c.strip()]
    missing_inc = [c for c in include if c not in df.columns]
    if missing_inc:
        raise ValueError(f"--include-flag-cols columns not found: {missing_inc}")

    for c in include:
        if c not in exclude and c not in flag_cols:
            flag_cols.append(c)

    if not flag_cols:
        raise ValueError(f"No columns matched --flag-regex '{args.flag_regex}' (and no included cols).")

    # germline filter (if columns exist)
    requested_germ = [c.strip() for c in args.germline_cols.split(",") if c.strip()]
    germ_cols = [c for c in requested_germ if c in df.columns]
    before = len(df)
    if germ_cols:
        for gc in germ_cols:
            df[gc] = coerce01(df[gc])
        df = df[~(df[germ_cols].max(axis=1) == 1)].copy()
    after = len(df)

    # coerce binary flags
    for c in flag_cols:
        df[c] = coerce01(df[c])

    # dedup
    key_cols = [args.chr_col, args.pos_col, args.ref_col, args.alt_col]
    if args.dedup:
        df = df[key_cols + flag_cols + germ_cols].groupby(key_cols, dropna=False, as_index=False).max()

    # SBS96 label
    fasta = pysam.FastaFile(args.fasta)
    fasta_refs = set(fasta.references)

    reasons = Counter()
    labels = []
    for chrom, pos, ref, alt in zip(df[args.chr_col], df[args.pos_col], df[args.ref_col], df[args.alt_col]):
        if pd.isna(pos):
            labels.append(None); reasons["pos_na"] += 1; continue
        lab, why = snv_to_sbs96(chrom, int(pos), ref, alt, fasta, fasta_refs,
                                require_ref_match=(not args.allow_ref_mismatch))
        labels.append(lab)
        if why: reasons[why] += 1
    df["_sbs96"] = labels
    df_snv = df[df["_sbs96"].notna()].copy()

    # map: sample_name -> column_name
    sample_to_col = {derive_sample_name_from_col(c): c for c in flag_cols}
    samples = list(sample_to_col.keys())
    sel = build_selected(samples)

    os_dir = os.path.dirname(args.out_prefix)
    os.makedirs(os_dir, exist_ok=True)

    def write_matrix(out_path, chosen_samples):
        out = pd.DataFrame({"MutationType": SBS96_INDEX})
        for s in chosen_samples:
            col = sample_to_col[s]
            sub = df_snv[df_snv[col] == 1]
            cts = Counter(sub["_sbs96"].tolist())
            out[s] = np.array([cts[k] for k in SBS96_INDEX], dtype=int)
        out.to_csv(out_path, sep="\t", index=False)
        return {s: int(out[s].sum()) for s in chosen_samples}

    totals_nod = write_matrix(args.out_prefix + ".no_depth.tsv", sorted(sel["no_depth"]))
    totals_30  = write_matrix(args.out_prefix + ".30x.tsv",      sorted(sel["30x"]))
    totals_100 = write_matrix(args.out_prefix + ".100x.tsv",     sorted(sel["100x"]))
    totals_200 = write_matrix(args.out_prefix + ".200x.tsv",     sorted(sel["200x"]))

    if args.annotate_variants_prefix:
        keep = key_cols + ["_sbs96"] + flag_cols
        df[keep].to_csv(f"{args.annotate_variants_prefix}.variants.with_sbs96.tsv", sep="\t", index=False)

    print(f"[OK] Input rows: {before} | after germline filter: {after} | unique variants used: {len(df_snv)}", file=sys.stderr)
    print(f"[OK] Wrote: {args.out_prefix}.no_depth.tsv  Totals: {totals_nod}", file=sys.stderr)
    print(f"[OK] Wrote: {args.out_prefix}.30x.tsv       Totals: {totals_30}", file=sys.stderr)
    print(f"[OK] Wrote: {args.out_prefix}.100x.tsv      Totals: {totals_100}", file=sys.stderr)
    print(f"[OK] Wrote: {args.out_prefix}.200x.tsv      Totals: {totals_200}", file=sys.stderr)
    if reasons:
        print(f"[INFO] SBS96 labeling skipped reasons (post-dedup): {dict(reasons)}", file=sys.stderr)

if __name__ == "__main__":
    main()

# python ./work/project/code/mut_signature/final/make_sbs96_caller_from_final_tsv.py \
#   --in-tsv ./work/regional_controls/08.final_tsv/no_unrelated.no_germline.no_blacklist.duplex_pileup.tsv \
#   --fasta ./resources/hg38/Homo_sapiens_assembly38.fasta \
#   --out-prefix ./work/project/data/mutational_signatures/sigprofiler/bam_split/existing_SBS96 \
#   --flag-regex '^is_existing_lower_ci_0001_' \
#   --include-flag-cols 'existing_any_PTA1_16' \
#   --exclude-flag-cols 'existing_any_PTA1_16_singleton,existing_any_PTA1_16_multiple,existing_any_PTA1_16_number,existing_any_PTA1_16_which'