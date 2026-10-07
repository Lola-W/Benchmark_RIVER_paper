# Benchmark_RIVER_paper

Code and processed input tables for the mosaicism benchmark study, organized
into analysis pipelines and figure scripts

- [Figure scripts](figure/README.md): eight scripts for the listed manuscript panels.
- [Variant calling](https://github.com/Lola-W/Benchmark_RIVER_paper/blob/main/pipelines/variant_calling/README.md): caller workflows, configurations, and alignment placeholders.
   - [ppmSeq](https://github.com/Lola-W/Benchmark_RIVER_paper/blob/main/pipelines/variant_calling/ppmseq/README.md): SRSNV scoring and post-hoc filtering of ppmSeq somatic SNV calls.
- [Data preparation](pipelines/data_preparation/README.md): annotation, filtering, read evidence, coverage, and signature analysis.
- `raw_data/`: 15 processed analysis tables and 25 coverage CDF tables used by the figure scripts, not raw sequencing reads.
- `R/nature_biotech_ggplot.R`: shared figure styling.
- `output/`: destination for generated figures and exported plot tables.

## Related software and citation

**RIVER:** RIVER contributors (2026). *RIVER: Read-level Integration for Variant
Evidence Recovery*. Version 0.1.0, commit
[`6bc202e`](https://github.com/Lola-W/RIVER/tree/6bc202eed7c345a6876124b24829b877bc7726a3).
[Source repository](https://github.com/Lola-W/RIVER). MIT License.

**ppmSeq:** [ppmSeq](pipelines/variant_calling/ppmseq/README.md) is a two-stage Snakemake workflow that runs the Ultima SRSNV pipeline on ppmSeq CRAMs, then removes germline variants, recurrent artifacts and sequence-context-specific noise to produce multiread, putative multiread and high-confidence singleton somatic SNV call sets.

## Reproduce figures

From this directory:

```bash
cd figure
python3 -m pip install -r fig2b_requirements.txt
python3 fig2b_plot_mosdepth_method_coverage.py
Rscript -e 'rmarkdown::render("fig2c_trinucleotide_minimal.Rmd")'
```

Run other R notebooks from `figure/`, evaluating chunks in order. Data paths are
relative to this package; see [figure requirements and usage](figure/README.md).
The supplied tables support the listed plotting scripts. This is not a complete
collection of every manuscript panel.

## Large figure inputs

Three pileup tables are configured for Git LFS in `.gitattributes`; the largest
is approximately 124 MiB and exceeds GitHub's 100 MiB regular-file limit.
[Install Git LFS](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage)
before the first `git add` when publishing this repository. After cloning:

```bash
git lfs install
git lfs pull
```

Use the complete TSV files when plotting, not Git LFS pointer files. This local
package contains the complete tables. A separate DOI-backed data deposit can
also distribute them, preserving their filenames and relative locations.

## Analysis scope and publication status

The pipeline READMEs describe software, references, input locations, and cluster
execution. BAM/CRAM files, reference genomes, containers, and large intermediate
datasets must be supplied separately. Input placeholders contain filenames only.

Existing analysis parameters and plotting code are preserved.

Before public release, complete the following:

- Replace the ppmSeq repository placeholder and provide reviewer/reader access to RIVER.
- Supply sequencing-data accessions and stable sources for custom masks, models, and other derived resources listed in the pipeline documentation.
- Specify the license for the original benchmark code and reuse terms for the processed tables. RIVER's MIT license applies to the separate RIVER package. Preserve attribution and license terms for reused third-party code, including the GPL-3.0 lab workflow source cited in the calling README.

## Contact info

- Jiaming Weng: jmweng@ucsd.edu
- Do Hyeon Cha, MD, MS: eric1@kaist.ac.kr
