#!/usr/bin/env python3
import argparse
import csv
from collections import Counter

def read_header(path):
    with open(path, "r", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        return next(r)

def col_idx_map(header):
    return {c:i for i,c in enumerate(header)}

def to_int(x):
    if x is None:
        return None
    x = x.strip()
    if x == "" or x == ".":
        return None
    try:
        return int(float(x))
    except Exception:
        return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="raw_duplex_pileup.tsv")
    ap.add_argument("--out1", default="no_blacklist.duplex_pileup.tsv")
    ap.add_argument("--out2", default="no_germline.no_blacklist.duplex_pileup.tsv")
    ap.add_argument("--out3", default="no_unrelated.no_germline.no_blacklist.duplex_pileup.tsv")
    args = ap.parse_args()

    blacklist_cols = [
        "in_ucsc_gaps",
        "in_encode_blacklist",
        "low_mappability",
        "in_segdup",
        "in_satellites_rDNA",
        "in_known_indel_any",
        "gnomad_common",
        "in_simple_repeats_ucsc",
        "in_hmer7",
    ]

    germline_col = "is_germline"

    germline_like_cols = [
        "is_germline_like_upper_ci_052_ILLUMINA",
        "is_germline_like_upper_ci_052_ILLUMINA_30x",
        "is_germline_like_upper_ci_052_ILLUMINA_100x",
        "is_germline_like_upper_ci_052_ILLUMINA_200x",
        "is_germline_like_upper_ci_052_DUPLEX_PPMSEQ",
        "is_germline_like_upper_ci_052_DUPLEX_PPMSEQ_30X",
        "is_germline_like_upper_ci_052_DUPLEX_PPMSEQ_100X",
        "is_germline_like_upper_ci_052_DUPLEX_PPMSEQ_200X",
        "is_germline_like_upper_ci_052_DUPLEX_UDSEQ",
        "is_germline_like_upper_ci_052_DUPLEX_UDSEQ_30X",
        "is_germline_like_upper_ci_052_DUPLEX_NANOSEQ",
        "is_germline_like_upper_ci_052_DUPLEX_HIDEFSEQ",
    ]

    unrelated_col = "existing_in_unrelated"

    header = read_header(args.input)
    idx = col_idx_map(header)

    # check required columns
    for c in blacklist_cols + [germline_col] + germline_like_cols + [unrelated_col]:
        if c not in idx:
            raise SystemExit(f"ERROR: missing column: {c}")

    # Counters
    n_total = 0

    # step1: count blacklist ones in raw, and filter -> out1
    blacklist_ones = Counter()
    n_pass1 = 0

    # step2: from pass1, count germline==1, count germline_like_any==1, filter -> out2
    n_germline1_in_pass1 = 0
    n_germline_like_any1_in_pass1 = 0
    n_pass2 = 0

    # step3: from pass2, count unrelated==1, filter -> out3
    n_unrelated1_in_pass2 = 0
    n_pass3 = 0

    with open(args.input, "r", newline="") as fin, \
         open(args.out1, "w", newline="") as f1, \
         open(args.out2, "w", newline="") as f2, \
         open(args.out3, "w", newline="") as f3:

        r = csv.reader(fin, delimiter="\t")
        w1 = csv.writer(f1, delimiter="\t", lineterminator="\n")
        w2 = csv.writer(f2, delimiter="\t", lineterminator="\n")
        w3 = csv.writer(f3, delimiter="\t", lineterminator="\n")

        hdr_in = next(r)
        # keep original header order
        w1.writerow(hdr_in)
        w2.writerow(hdr_in)
        w3.writerow(hdr_in)

        for row in r:
            n_total += 1
            if len(row) < len(hdr_in):
                row += [""] * (len(hdr_in) - len(row))

            # ---- step1 stats + filter ----
            bl_vals = {}
            for c in blacklist_cols:
                v = to_int(row[idx[c]])
                bl_vals[c] = v
                if v == 1:
                    blacklist_ones[c] += 1

            pass1 = all(bl_vals[c] == 0 for c in blacklist_cols)
            if not pass1:
                continue

            n_pass1 += 1
            w1.writerow(row)

            # ---- step2 stats + filter (computed on pass1 stream) ----
            g = to_int(row[idx[germline_col]])
            if g == 1:
                n_germline1_in_pass1 += 1

            gl_any = 0
            for c in germline_like_cols:
                v = to_int(row[idx[c]])
                if v == 1:
                    gl_any = 1
                    break
            if gl_any == 1:
                n_germline_like_any1_in_pass1 += 1

            pass2 = (g == 0) and (gl_any == 0)
            if not pass2:
                continue

            n_pass2 += 1
            w2.writerow(row)

            # ---- step3 stats + filter ----
            u = to_int(row[idx[unrelated_col]])
            if u == 1:
                n_unrelated1_in_pass2 += 1

            pass3 = (u == 0)
            if not pass3:
                continue

            n_pass3 += 1
            w3.writerow(row)

    # ---- report ----
    print("=== INPUT ===")
    print(f"total_rows\t{n_total}")
    print()

    print("=== STEP 1: no_blacklist ===")
    for c in blacklist_cols:
        print(f"{c}_n_eq_1\t{blacklist_ones[c]}")
    print(f"step1_output_rows\t{n_pass1}")
    print(f"step1_output\t{args.out1}")
    print()

    print("=== STEP 2: no_germline.no_blacklist ===")
    print(f"{germline_col}_n_eq_1_in_step1\t{n_germline1_in_pass1}")
    print(f"germline_like_any_n_eq_1_in_step1\t{n_germline_like_any1_in_pass1}")
    print(f"step2_output_rows\t{n_pass2}")
    print(f"step2_output\t{args.out2}")
    print()

    print("=== STEP 3: no_unrelated.no_germline.no_blacklist ===")
    print(f"{unrelated_col}_n_eq_1_in_step2\t{n_unrelated1_in_pass2}")
    print(f"step3_output_rows\t{n_pass3}")
    print(f"step3_output\t{args.out3}")

if __name__ == "__main__":
    main()
