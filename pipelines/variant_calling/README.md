# Variant calling

Study caller workflows and settings. Install the third-party callers listed in
[SOFTWARE.md](SOFTWARE.md) and supply the configured inputs.

## Included workflows

| Directory | Input | Calling output |
|---|---|---|
| `illumina/paired` | Brain and matched-normal BAMs | `results/mutect_filter_vcf/*_somatic.vcf.gz` and `results/strelka_vcf/*_somatic.vcf.gz` |
| `illumina/mutect2_single` | Brain BAM and public panel of normals | `results/scratch/passed_single_mode/*_somatic.vcf.gz` |
| `illumina/mosaicforecast` | Brain BAM and PASS tumor-only Mutect2 VCF | `results/prediction_results/*_prediction.txt`; phasing files |
| `illumina/deepmosaic` | Brain BAM and PASS tumor-only Mutect2 VCF | `features.txt`, `output.txt` |
| `illumina/mosaichunter` | Brain BAM | `results/{sample}/{chunk}/final.passed.tsv` |
| `illumina/deepsomatic` | Brain and matched-normal BAMs | Somatic VCF and gVCF |
| `dupcaller` | NanoSeq or UDSeq brain/control BAMs, germline VCF, noise mask | DupCaller call outputs at the configured prefix |
| `hidefseq` | CCS reads, PacBio germline BAM and Clair3/DeepVariant germline VCFs | HiDEF-seq `finalCalls/` tables, including `1-22X` SBS calls |
| `pta` | Single-neuron and bulk BAMs | SCAN2 `somatic_genotypes.rda` and `callable_regions.rda` per neuron |

**Scope:** caller outputs are not yet the final benchmark union or the input
tables in `../../raw_data/`. Cross-caller consensus, cross-sample filtering,
amplicon validation, plotting-table generation, and coordinate conversion belong
to the subsequent processing package. ppmSeq calling is external; see the
[separate repository placeholder](ppmseq/README.md). This package includes no
ppmSeq calling scripts or caller configurations.

## Setup

Use Linux with Slurm and a shared filesystem. Install the callers and their
dependencies listed in [SOFTWARE.md](SOFTWARE.md). These Snakemake launchers use
the **Snakemake 7** `--cluster` interface. Load/activate the appropriate software
environment before submission; shell scripts do not activate personal Conda
environments. Singularity must be available for container-based workflows.

Populate the paths listed in [inputs/BAM_PLACEHOLDERS.tsv](inputs/BAM_PLACEHOLDERS.tsv).
This is a text inventory only: no real or empty BAM files are included. Place the
BAMs and indexes in `inputs/`, except the PTA files in `pta/input/`. Existing names
and sample identifiers are retained; edit sample sheets if your filenames differ.
Place reference resources in `resources/`, callers in `software/`, and container
images in `containers/` as indicated by each configuration. Files may be linked,
but container jobs must also bind the target location if it is outside this
package. Required non-BAM inputs are listed in
[external_inputs.tsv](external_inputs.tsv). These are placeholders, not downloads.

Host paths in the distributed configurations are relative. Shell configurations
resolve their own location at runtime. `prepare_paths.py` converts relative YAML
and sample-sheet paths to temporary `.runtime/` files before Snakemake, HiDEF-seq,
and DeepMosaic start. This keeps the original Snakefiles/calling commands intact
and prevents relative BAM symlink targets from breaking in output directories.
Only values starting `./` or `../` are resolved; use those prefixes when editing
paths. Source configurations are never overwritten. Use the supplied launchers
so path preparation is applied. Do not commit generated `.runtime/` files.
The fixed `/hidef/...` paths inside the HiDEF-seq container are retained because
they name software within that image, not locations on the author's computer.

Install Python 3 and PyYAML 6.x in the launcher environment, in addition to the
caller dependencies. Provide BAM indexes, FASTA `.fai`/sequence dictionaries,
and VCF indexes; MosaicForecast expects `<file>.bam.bai`. The CCS-read BAM is a
PacBio read input, not necessarily a coordinate-sorted alignment BAM.

Illumina uses `GRCh38_full_analysis_set_plus_decoy_hla.fa`;
DupCaller and HiDEF-seq use
`GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`. Do not interchange these
references simply because both are GRCh38. Resource filenames and the study's
configuration values are retained in the configurations and
[external input inventory](external_inputs.tsv).

PTA uses the same SCAN2 wrapper for either assembly. Its publication configuration
selects hg38 resources as requested by the study author. Supply hg38 BAMs, dbSNP,
and a **SHAPEIT-format hg38 reference panel**. The legacy configuration also
listed an Eagle panel, but the retained wrapper calls `--shapeit-refpanel`;
these panel formats are not interchangeable. This hg38 configuration was not
executed or independently validated during packaging. See
[pta/input/README.md](pta/input/README.md) for BAM naming.

Recorded sex values are retained in each method's configuration. They differ
between some original pipelines; they have not been silently harmonized here.

## Submit

Commands below are examples for future use, not commands run during packaging.
Set the account, partition and QoS for your cluster. On TSCC, submit computational
work through Slurm; do not run a workflow on the login node.

```bash
# From variant_calling/; activate the appropriate environment first.
export SBATCH_ACCOUNT=your_account
export SBATCH_PARTITION=your_partition
export SBATCH_QOS=your_qos

sbatch run_snakemake.sbatch illumina/paired
sbatch run_snakemake.sbatch illumina/mutect2_single
# After tumor-only Mutect2 finishes:
sbatch run_snakemake.sbatch illumina/mosaicforecast
sbatch run_snakemake.sbatch illumina/mosaichunter

sbatch illumina/deepmosaic/run.sbatch illumina/deepmosaic/config.sh
sbatch illumina/deepsomatic/run.sbatch illumina/deepsomatic/config.sh
sbatch dupcaller/run.sbatch dupcaller/nano.config.sh
sbatch dupcaller/run.sbatch dupcaller/ud.config.sh

# Submit from hidefseq/ after editing its YAML and Slurm/container configuration:
cd hidefseq
sbatch run.sbatch
```

For PTA, recreate the recorded environment from
`pta/workflow/envs/SCAN2.yaml` on a compute node, activate it, populate `pta/input`,
and submit `sbatch run.sbatch` from `pta/`. Its launcher processes neurons
sequentially with the original 24-thread SCAN2 setting.

The shared launcher defaults to 40 jobs, 12 CPUs and 96 GB per worker. Override
`JOBS`, `JOB_CPUS`, `JOB_MEM`, and `JOB_TIME` for your cluster.
These workflows have not been executed end to end in this submission layout.


## Attribution

Prior lab WGS workflows were developed by Xin Xu and Xiaoxu Yang with Martin
Breuss, based on work by Renee D. George. See the
[reference repository](https://github.com/shishenyxx/Adult_brain_somatic_mosaicism)
for attribution and its GPL-3.0 license. Third-party callers retain their own
licenses.
