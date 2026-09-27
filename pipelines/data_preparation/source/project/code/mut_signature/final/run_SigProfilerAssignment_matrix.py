#!/usr/bin/env python3
import argparse, os
import pandas as pd
from SigProfilerAssignment import Analyzer as Analyze

def read_matrix(path):
    ext = os.path.splitext(path)[1].lower()
    sep = "\t" if ext in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep)

    # Require first column to be MutationType (SigProfiler-style)
    if "MutationType" not in df.columns:
        raise ValueError("Input must contain a 'MutationType' column (SBS96 labels).")

    # Coerce sample columns to numeric ints
    sample_cols = [c for c in df.columns if c != "MutationType"]
    df[sample_cols] = (df[sample_cols]
                       .apply(pd.to_numeric, errors="coerce")
                       .fillna(0)
                       .astype(int))
    return df

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i","--input", required=True, help="SBS96 matrix (MutationType + sample cols)")
    ap.add_argument("-o","--output", required=True, help="Output directory")
    ap.add_argument("--project", default="SigProfilerAssignment", help="Label used in plot titles")
    ap.add_argument("--cosmic_version", type=float, default=3.4)
    ap.add_argument("--context_type", default="96")
    ap.add_argument("--genome_build", default="GRCh38")

    ap.add_argument("--make_plots", action="store_true", help="SigProfilerAssignment summary plots")
    ap.add_argument("--sample_reconstruction_plots", default="pdf",
                    choices=["pdf","png","both","none"],
                    help="Per-sample reconstruction plots format")

    ap.add_argument("--plot_sbs96", action="store_true",
                    help="Make publication SBS96 spectra using SigProfilerPlotting")
    ap.add_argument("--plot_percentage", action="store_true",
                    help="Plot spectra as percentages (better for cross-platform comparison)")
    ap.add_argument("--plot_format", default="pdf", choices=["pdf","png"])
    ap.add_argument("--plot_dpi", type=int, default=300)

    args = ap.parse_args()
    os.makedirs(args.output, exist_ok=True)

    df = read_matrix(args.input)
    tsv_path = os.path.join(args.output, "counts.SBS96.cleaned.tsv")
    df.to_csv(tsv_path, sep="\t", index=False)

    Analyze.cosmic_fit(
        samples=tsv_path,
        output=args.output,
        input_type="matrix",
        context_type=args.context_type,
        genome_build=args.genome_build,
        cosmic_version=args.cosmic_version,
        make_plots=args.make_plots,
        sample_reconstruction_plots=args.sample_reconstruction_plots,
        verbose=True,
    )  # make_plots + sample_reconstruction_plots are the knobs you want :contentReference[oaicite:4]{index=4}

    if args.plot_sbs96:
        import sigProfilerPlotting as sigPlt
        plot_dir = os.path.join(args.output, "Publication_Plots")
        os.makedirs(plot_dir, exist_ok=True)

        # sigProfilerPlotting.plotSBS(matrix_path, output_path, project, plot_type, percentage, ..., savefig_format, ..., dpi)
        sigPlt.plotSBS(tsv_path, plot_dir, args.project, args.context_type,
                       args.plot_percentage, None, None, None,
                       args.plot_format, None, args.plot_dpi)  # :contentReference[oaicite:5]{index=5}

        print("Publication SBS96 plots ->", plot_dir)

    print("Done. Search for figures with:")
    print(f"  find {args.output} -type f \\( -iname \"*.pdf\" -o -iname \"*.png\" \\) | sort")

if __name__ == "__main__":
    main()






# python ./work/project/code/mut_signature/final/run_SigProfilerAssignment_matrix.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/counts.SBS96.single_cell_combined.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/single_cell_combined \
#   --project 7614 \
#   --cosmic_version 3.4 \
#   --context_type 96 \
#   --genome_build GRCh38 \
#   --make_plots \
#   --sample_reconstruction_plots pdf \
#   --plot_sbs96 \
#   --plot_percentage \
#   --plot_format pdf \
#   --plot_dpi 300



# python ./work/project/code/mut_signature/final/run_SigProfilerAssignment_matrix.py \
#   -i ./work/project/data/mutational_signatures/sigprofiler/counts.SBS96.single_neuron_separate.tsv \
#   -o ./work/project/data/mutational_signatures/result/COSMIC/single_cell_seperate \
#   --project 7614 \
#   --cosmic_version 3.4 \
#   --context_type 96 \
#   --genome_build GRCh38 \
#   --make_plots \
#   --sample_reconstruction_plots pdf \
#   --plot_sbs96 \
#   --plot_percentage \
#   --plot_format pdf \
#   --plot_dpi 300