# Manuscript figure scripts

Run R notebooks from this directory and evaluate chunks in order; later chunks
reuse objects created earlier. R Markdown normally uses the notebook directory
as its working directory. Inputs are in `../raw_data/`, shared styling is in
`../R/`, and explicit exports go into `../output/`.

- [Figure 2b](fig2b_plot_mosdepth_method_coverage.py): coverage distributions.
- [Figure 2c](fig2c_trinucleotide_minimal.Rmd): trinucleotide composition.
- [Figures 3b–3c](fig3b_3c_filtering_caller_overlap.Rmd): filtering and caller overlap.
- [Figures 3e–3f](fig3e_3f_caller_signatures.Rmd): caller mutation signatures.
- [Figures 4b–4d](fig4b_4c_4d_read_support.Rmd): read support and allele fractions.
- [Figures 4e, 4g–4h and Extended Data 4b](fig4e_4g_4h_extData4b_correlations_signatures.Rmd): correlations and read-supported signatures.
- [Figure 5g](fig5g_lolliplot_geoclone.Rmd): four representative geoclones across 11 tissues.

The filenames use manuscript numbering; original chunk names are retained.
For Figure 5g use the `fig2a-style-four-pattern-export` chunk. Other exploratory
exports in that notebook are not the Figure 5g panel. These scripts cover the
listed panels, not every panel in the manuscript.

## Python coverage plot

```bash
python3 -m pip install -r fig2b_requirements.txt
python3 fig2b_plot_mosdepth_method_coverage.py --input-dir ../raw_data/coverage --output ../output/fig2b_coverage.pdf
```

The 25 CDF tables produce ten curves, including the mean across 16 PTA neurons.
Default paths are located from the script itself. Only CDF inputs are supplied;
`--kind pdf` selects probability-density input tables, not the output file format.
Raw curves are dashed and filtered curves solid in the retained code. Depth
labels reflect the original nominal sampling settings, not necessarily an
observed mean of exactly 30× for every method.

## R notebooks

```bash
Rscript -e 'rmarkdown::render("fig2c_trinucleotide_minimal.Rmd")'
```

Use Python 3.10+ for the pinned matplotlib 3.10.1 coverage-plot dependency.
For R, install readr, dplyr, tidyr, tibble, stringr, ggplot2, scales, UpSetR,
circlize, gtable, ragg, knitr, forcats, ggdendro, ggrepel, patchwork, and purrr
from CRAN, plus ComplexHeatmap from Bioconductor. HTML rendering also requires
rmarkdown and Pandoc; the shared style uses Arial. A tested, version-pinned R
environment is not yet supplied.

All input tables for these seven scripts are included.

For read-support methodology, cite [RIVER](../README.md#related-software-and-citation).
The supplied figure tables can be plotted without installing RIVER or rerunning
variant callers. Extended Data 2 (RIVER input-control/FDR analysis) is not
included in this seven-script set.
