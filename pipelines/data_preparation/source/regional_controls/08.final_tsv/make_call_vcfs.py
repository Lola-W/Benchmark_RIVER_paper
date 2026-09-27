#!/usr/bin/env python3
import sys, csv, os, random, re

MAX_N = 3000
SEED = 42

raw_tsv = sys.argv[1]
fil_tsv = sys.argv[2]

ALLOWED_CHR = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CHR_ORDER = {c:i for i,c in enumerate(ALLOWED_CHR)}

def variant_key(v):
    ch, pos, ref, alt = v
    return (CHR_ORDER.get(ch, 10**9), int(pos), ref, alt)

def sanitize(name):
    return re.sub(r'[^A-Za-z0-9._+-]+', '_', name)

def process(tsv_path, outdir):
    rng = random.Random(SEED)
    os.makedirs(outdir, exist_ok=True)

    with open(tsv_path, "r", newline="") as f:
        r = csv.reader(f, delimiter="\t")
        header = next(r)
        idx = {c:i for i,c in enumerate(header)}

        # call_* + AmpliSeq, exclude call_PTA_multiple
        cols = [c for c in header if c.startswith("call_") and c != "call_PTA_multiple"]
        if "AmpliSeq" in idx:
            cols.append("AmpliSeq")

        for need in ["chr", "pos", "ref", "alt"]:
            if need not in idx:
                raise SystemExit(f"ERROR: missing required column: {need}")

        counts = {c: 0 for c in cols}
        reservoirs = {c: [] for c in cols}

        for row in r:
            if len(row) < len(header):
                row += [""] * (len(header) - len(row))

            ch = row[idx["chr"]]
            if ch not in CHR_ORDER:
                continue  # keep only chr1-22,chrX

            v = (ch, row[idx["pos"]], row[idx["ref"]], row[idx["alt"]])

            for c in cols:
                if row[idx[c]].strip() == "1":
                    counts[c] += 1
                    # reservoir sampling
                    if len(reservoirs[c]) < MAX_N:
                        reservoirs[c].append(v)
                    else:
                        j = rng.randint(0, counts[c] - 1)
                        if j < MAX_N:
                            reservoirs[c][j] = v

        # write vcf per column
        for c in cols:
            n = counts[c]
            k = min(n, MAX_N)
            sample = reservoirs[c][:k]

            # de-dup within sampled set, then sort
            sample = sorted(set(sample), key=variant_key)

            out = os.path.join(outdir, sanitize(c) + ".vcf")
            with open(out, "w", newline="") as w:
                w.write("##fileformat=VCFv4.2\n")
                w.write(f"##source=make_call_vcfs_seed{SEED}_max{MAX_N}\n")
                w.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
                for ch, pos, ref, alt in sample:
                    w.write(f"{ch}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t.\n")

        return cols, counts

def write_report(cols, counts_raw, counts_fil):
    cols_all = sorted(set(cols))
    with open("call_sampling_report.tsv", "w", newline="") as w:
        w.write("col\traw_n_eq_1\traw_written\tfiltered_n_eq_1\tfiltered_written\n")
        for c in cols_all:
            rn = counts_raw.get(c, 0)
            fn = counts_fil.get(c, 0)
            w.write(f"{c}\t{rn}\t{min(rn,MAX_N)}\t{fn}\t{min(fn,MAX_N)}\n")

cols_raw, counts_raw = process(raw_tsv, "00.raw_vcf")
cols_fil, counts_fil = process(fil_tsv, "01.filtered_vcf")
cols = sorted(set(cols_raw) | set(cols_fil))
write_report(cols, counts_raw, counts_fil)
