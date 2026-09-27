#!/usr/bin/env python3
import argparse, csv, os, re
DEPTHS=["30x","100x","200x"]
DEPTH_ANY_RE=re.compile(r"(?:^|_)(30x|100x|200x)(?:_|$)")
def has_depth_any(s): return bool(DEPTH_ANY_RE.search(s))
def has_depth(s,d): return bool(re.search(rf"(?:^|_){re.escape(d)}(?:_|$)", s))
def base_key(s):
    out=s
    for d in DEPTHS:
        out=out.replace(f"_{d}_","_")
        if out.endswith(f"_{d}"): out=out[:-(len(d)+1)]
    return out
def pick_preferred(cands,d):
    cands=sorted(cands, key=lambda x:(0 if x.endswith(f"_{d}") else 1, len(x), x))
    return cands[0] if cands else None

def build_selected(samples_set, exclude="AmpliSeq"):
    samples={s for s in samples_set if s!=exclude}
    bases=sorted({base_key(s) for s in samples})
    sel={}
    sel["base"]={s for s in samples if not has_depth_any(s)}
    for d in DEPTHS:
        chosen=set()
        for b in bases:
            cands_d=[s for s in samples if base_key(s)==b and has_depth(s,d)]
            if cands_d: chosen.add(pick_preferred(cands_d,d))
            else:
                if b in samples: chosen.add(b)
                else:
                    cands_nd=[s for s in samples if base_key(s)==b and not has_depth_any(s)]
                    if cands_nd: chosen.add(sorted(cands_nd)[0])
        sel[d]=chosen
    return sel

def write_long(infile, out_prefix, exclude):
    samples=set()
    with open(infile, newline="") as f:
        r=csv.DictReader(f, delimiter="\t")
        if "sample" not in (r.fieldnames or []): raise SystemExit("Missing required column: sample")
        for row in r: samples.add(row["sample"])
    sel=build_selected(samples, exclude)

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    outs={"base":out_prefix+".no_depth.tsv","30x":out_prefix+".30x.tsv","100x":out_prefix+".100x.tsv","200x":out_prefix+".200x.tsv"}
    fhs={k:open(v,"w",newline="") for k,v in outs.items()}
    try:
        with open(infile, newline="") as f:
            r=csv.DictReader(f, delimiter="\t")
            ws={k:csv.DictWriter(fhs[k], fieldnames=r.fieldnames, delimiter="\t", lineterminator="\n") for k in fhs}
            for w in ws.values(): w.writeheader()
            for row in r:
                s=row["sample"]
                if s==exclude: continue
                if s in sel["base"]: ws["base"].writerow(row)
                if s in sel["30x"]:  ws["30x"].writerow(row)
                if s in sel["100x"]: ws["100x"].writerow(row)
                if s in sel["200x"]: ws["200x"].writerow(row)
    finally:
        for fh in fhs.values(): fh.close()

def write_wide(infile, out_prefix, exclude):
    with open(infile, newline="") as f:
        r=csv.reader(f, delimiter="\t")
        header=next(r)
        if len(header)<2: raise SystemExit("Wide TSV must have >=2 columns (key + samples).")
        keycol=header[0]
        sample_cols=header[1:]
        sel=build_selected(sample_cols, exclude)

    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    outs={"base":out_prefix+".no_depth.tsv","30x":out_prefix+".30x.tsv","100x":out_prefix+".100x.tsv","200x":out_prefix+".200x.tsv"}

    # Choose columns to write for each output (keep keycol + chosen sample columns)
    cols_to_keep={}
    for k in ["base","30x","100x","200x"]:
        chosen=[c for c in sample_cols if c in sel["base"]] if k=="base" else [c for c in sample_cols if c in sel[k]]
        # Keep deterministic ordering like original header
        cols_to_keep[k]=[0]+[1+i for i,c in enumerate(sample_cols) if c in set(chosen)]

    fhs={k:open(v,"w",newline="") for k,v in outs.items()}
    try:
        # write headers
        for k,fh in fhs.items():
            w=csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow([header[i] for i in cols_to_keep[k]])
        # write rows
        with open(infile, newline="") as f:
            r=csv.reader(f, delimiter="\t")
            next(r)  # skip header
            for row in r:
                for k,fh in fhs.items():
                    w=csv.writer(fh, delimiter="\t", lineterminator="\n")
                    w.writerow([row[i] for i in cols_to_keep[k]])
    finally:
        for fh in fhs.values(): fh.close()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("-i","--infile", required=True)
    ap.add_argument("-o","--out_prefix", required=True, help="prefix path for outputs")
    ap.add_argument("--exclude_sample", default="AmpliSeq")
    args=ap.parse_args()

    # Detect format from header
    with open(args.infile, newline="") as f:
        header=f.readline().rstrip("\n").split("\t")
    if "sample" in header:
        write_long(args.infile, args.out_prefix, args.exclude_sample)
    else:
        write_wide(args.infile, args.out_prefix, args.exclude_sample)

if __name__=="__main__":
    main()
    
    
# python ./work/project/code/mut_signature/final/split_by_depth.py -i ./work/project/data/mutational_signatures/sigprofiler/caller.counts.SBS96.platforms.from_vcfs.tsv -o ./work/project/data/mutational_signatures/sigprofiler/caller_split/caller.counts.SBS96.platforms