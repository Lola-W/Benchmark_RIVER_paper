#!/usr/bin/env python3
import argparse, os, sys, inspect
import pandas as pd
import numpy as np
if not hasattr(np, "mat"):
    np.mat = np.asmatrix  # compat shim for nimfa on NumPy>=2.0


def make_sbs96_index():
    subs = ["C>A","C>G","C>T","T>A","T>C","T>G"]
    bases = ["A","C","G","T"]
    return [f"{l}[{r}>{a}]{rr}" for r,a in (s.split(">") for s in subs) for l in bases for rr in bases]
SBS96_INDEX = make_sbs96_index()
SBS96_SET = set(SBS96_INDEX)

def read_matrix(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()
    sep = "\t" if ext in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep)
    if "MutationType" not in df.columns: raise ValueError("Need 'MutationType' column (SBS96 labels).")
    sample_cols = [c for c in df.columns if c != "MutationType"]
    if not sample_cols: raise ValueError("Need >=1 sample column.")
    df["MutationType"] = df["MutationType"].astype(str)
    df[sample_cols] = df[sample_cols].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    return df

def set_thread_env(n: int):
    for k in ("OMP_NUM_THREADS","MKL_NUM_THREADS","OPENBLAS_NUM_THREADS",
              "NUMEXPR_NUM_THREADS","VECLIB_MAXIMUM_THREADS"):
        os.environ[k] = str(n)


def canonicalize_sbs96(df: pd.DataFrame) -> pd.DataFrame:
    sample_cols = [c for c in df.columns if c != "MutationType"]
    bad = df.loc[~df["MutationType"].isin(SBS96_SET), "MutationType"].unique().tolist()
    if bad: raise ValueError(f"Non-SBS96 MutationType labels (first 10): {bad[:10]}")
    df = df.groupby("MutationType", as_index=False)[sample_cols].sum()
    out = pd.DataFrame({"MutationType": SBS96_INDEX}).merge(df, on="MutationType", how="left")
    out[sample_cols] = out[sample_cols].fillna(0).astype(int)
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i","--input", required=True)
    ap.add_argument("-o","--output_dir", required=True)
    ap.add_argument("--project", required=True, help="Output folder name created in CWD.")
    ap.add_argument("--min_sigs", type=int, default=1)
    ap.add_argument("--max_sigs", type=int, default=8)
    ap.add_argument("--nmf_replicates", type=int, default=200)
    ap.add_argument("--resample", action="store_true", default=True)
    ap.add_argument("--no_resample", dest="resample", action="store_false")
    ap.add_argument("--cpu", type=int, default=None)
    ap.add_argument("--assignment_cpu", type=int, default=None)

    ap.add_argument("--context_type", default="96")
    ap.add_argument("--matrix_normalization", default="gmm")
    ap.add_argument("--nmf_init", default="random")
    ap.add_argument("--seeds", default="random")
    ap.add_argument("--precision", default="single", choices=["single","double"])

    ap.add_argument("--min_nmf_iterations", type=int, default=None)
    ap.add_argument("--max_nmf_iterations", type=int, default=None)
    ap.add_argument("--nmf_test_conv", type=int, default=None)
    ap.add_argument("--nmf_tolerance", type=float, default=None)

    ap.add_argument("--stability", type=float, default=None)
    ap.add_argument("--min_stability", type=float, default=None)
    ap.add_argument("--combined_stability", type=float, default=None)
    ap.add_argument("--allow_stability_drop", action="store_true", default=None)

    ap.add_argument("--reference_genome", default="GRCh38")
    ap.add_argument("--opportunity_genome", default="GRCh38")
    ap.add_argument("--exome", action="store_true", default=False)

    ap.add_argument("--cosmic_version", type=float, default=None)
    ap.add_argument("--make_decomposition_plots", action="store_true", default=None)
    ap.add_argument("--collapse_to_SBS96", action="store_true", default=None)

    ap.add_argument("--export_probabilities", action="store_true", default=None)
    ap.add_argument("--get_all_matrices", action="store_true", default=False)
    ap.add_argument("--volume", default=None)

    args = ap.parse_args()
    args.cpu = args.cpu or int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
    set_thread_env(args.cpu)
    
    outdir = os.path.abspath(args.output_dir)
    os.makedirs(outdir, exist_ok=True)

    df = canonicalize_sbs96(read_matrix(args.input))
    cleaned = os.path.abspath(os.path.join(outdir, f"{args.project}.cleaned.SBS96.tsv"))
    df.to_csv(cleaned, sep="\t", index=False)

    # Run inside outdir because SigProfilerExtractor writes output folder to CWD (README)
    cwd = os.getcwd()
    os.chdir(outdir)

    try:
        from SigProfilerExtractor import sigpro as sig
        ver = getattr(__import__("SigProfilerExtractor"), "__version__", "unknown")
    except Exception as e:
        os.chdir(cwd)
        raise RuntimeError("Failed to import SigProfilerExtractor in this env.") from e

    print(f"[INFO] SigProfilerExtractor version: {ver}", file=sys.stderr)
    print(f"[INFO] Input matrix (absolute): {cleaned}", file=sys.stderr)
    print(f"[INFO] Output folder: {os.path.join(outdir, args.project)}", file=sys.stderr)

    # Build kwargs and drop ones not supported by the installed SigProfilerExtractor
    kwargs = dict(
        reference_genome=args.reference_genome,          # applicable only for VCF per README
        opportunity_genome=args.opportunity_genome,      # key for GRCh38 in matrix workflows
        context_type=args.context_type,
        exome=args.exome,
        minimum_signatures=args.min_sigs,
        maximum_signatures=args.max_sigs,
        nmf_replicates=args.nmf_replicates,
        resample=args.resample,
        cpu=args.cpu,
        matrix_normalization=args.matrix_normalization,
        nmf_init=args.nmf_init,
        seeds=args.seeds,
        precision=args.precision,
        get_all_signature_matrices=args.get_all_matrices,
        volume=args.volume,
    )

    # Optional knobs (only pass if user set them)
    opt = {
        "assignment_cpu": args.assignment_cpu,
        "min_nmf_iterations": args.min_nmf_iterations,
        "max_nmf_iterations": args.max_nmf_iterations,
        "nmf_test_conv": args.nmf_test_conv,
        "nmf_tolerance": args.nmf_tolerance,
        "stability": args.stability,
        "min_stability": args.min_stability,
        "combined_stability": args.combined_stability,
        "allow_stability_drop": args.allow_stability_drop,
        "cosmic_version": args.cosmic_version,
        "make_decomposition_plots": args.make_decomposition_plots,
        "collapse_to_SBS96": args.collapse_to_SBS96,
        "export_probabilities": args.export_probabilities,
    }
    for k,v in opt.items():
        if v is not None: kwargs[k] = v

    sig_params = set(inspect.signature(sig.sigProfilerExtractor).parameters.keys())
    dropped = [k for k in kwargs if k not in sig_params]
    kwargs = {k:v for k,v in kwargs.items() if k in sig_params}
    if dropped:
        print(f"[WARN] Dropping unsupported params for your installed version: {dropped}", file=sys.stderr)

    sig.sigProfilerExtractor("matrix", args.project, cleaned, **kwargs)

    os.chdir(cwd)
    print("[OK] Done.", file=sys.stderr)

if __name__ == "__main__":
    main()



# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/counts.SBS96.single_cell_combined.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/single_cell_combined \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 8 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots

# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/counts.SBS96.single_neuron_separate.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/single_cell_seperate \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 8 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots


# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/counts.SBS96.single_cell_combined.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/single_cell_combined \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 8 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots

# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/counts.SBS96.with_30x100x200x_and_PTA.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/single_cell_seperate \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 8 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots


# final
# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/caller_split_noputative_final/caller.counts.SBS96.platforms.no_depth.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/caller_from_vcfs_final \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 12 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots

# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/caller_split_noputative_final/caller.counts.SBS96.platforms.30x.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/caller_30x_from_vcfs_final \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 8 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots

# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/caller_split_noputative_final/caller.counts.SBS96.platforms.100x.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/caller_100x_from_vcfs_final \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 12 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots

# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/caller_split_noputative_final/caller.counts.SBS96.platforms.200x.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/caller_200x_from_vcfs_final \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 12 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots

# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/bam_split/existing_SBS96.no_depth.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/bam.no_depth_final \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 8 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots

# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/bam_split/existing_SBS96.30x.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/bam.30x_final \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 8 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots




# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/bam_split/existing_SBS96.100x.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/bam.100x_final \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 8 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots

# python ./work/project/code/mut_signature/final/denovo_sigprofilerextractor.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/bam_split/existing_SBS96.200x.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/bam.200x_final \
#   --project 7614 \
#   --context_type 96 \
#   --opportunity_genome GRCh38 \
#   --min_sigs 1 --max_sigs 8 \
#   --nmf_replicates 200 \
#   --cpu 8 \
#   --cosmic_version 3.4 \
#   --make_decomposition_plots