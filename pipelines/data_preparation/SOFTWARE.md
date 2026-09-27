# Software and resources

Preparation utilities `prepare_workspace.py` and
`assemble_trinucleotide_matrix.py` use Python 3.9+ and the standard library.
The preserved source collection requires different environments by stage:

- Bash 4+ (associative arrays), Slurm for job arrays, and standard shell tools.
- SAMtools, BCFtools, bgzip/tabix, BEDTools, mosdepth, GNU parallel,
  UCSC bigWigToBedGraph and liftOver where referenced by the chosen stage.
- Python: NumPy, pandas, SciPy, pysam, pyBigWig; signature stages additionally
  use SigProfilerMatrixGenerator, SigProfilerExtractor, SigProfilerAssignment
  and sigProfilerPlotting. Consult each script's imports and CLI arguments.
- R packages used by the preparation notebook include readr, dplyr, tidyr,
  stringr and plotting packages named in its individual chunks. R Markdown
  and its rendering dependencies are needed only for rendering the notebook.

This collection does not provide a complete version-pinned environment. The
historical environment names were site-specific and were removed from launcher
setup. Install/activate the appropriate software and put executables on PATH
before running a stage; do not assume one environment satisfies every stage.

Reference FASTAs, indexes, sequence dictionaries, coordinate-conversion chains,
annotation VCFs with indexes, mappability/repeat/blacklist tracks, and signature
reference data must be obtained separately. Preserve each stage's assembly and
contig naming. Input and resource paths are specified in the selected scripts,
including references constructed from variables and local filenames.
Alignment locations include additional regional/unrelated controls and cleaned
or downsampled intermediates beyond the sibling calling package's BAM inventory.
No raw data or reference downloads occur during workspace preparation.
