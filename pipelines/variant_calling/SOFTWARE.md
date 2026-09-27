# Software and resources

Recorded software versions and upstream installation sources:

| Component | Recorded version or identifier | Upstream |
|---|---|---|
| GATK, Illumina workflows | 4.6.2.0 | [GATK](https://github.com/broadinstitute/gatk/tree/4.6.2.0) |
| Strelka2 | 2.9.10; Python 2 required | [Strelka](https://github.com/Illumina/strelka/releases/tag/v2.9.10) |
| DeepSomatic | 1.9.0 CPU image; WGS model | [DeepSomatic](https://github.com/google/deepsomatic/tree/v1.9.0) |
| DeepMosaic | Exact release unresolved; `efficientnet-b4_epoch_6.pt`, batch 12 | [DeepMosaic](https://github.com/XiaoxuYangLab/DeepMosaic) |
| MosaicForecast | Exact source release unresolved; `200xRFmodel_addRMSK_Refine.rds` | [MosaicForecast](https://github.com/parklab/MosaicForecast) |
| MosaicHunter | Parent commit `acb1922591e08be8f52e8127faeddf40ac9e9e68`; local BLAT-thread patch | [MosaicHunter](https://github.com/zzhang526/MosaicHunter/tree/acb1922591e08be8f52e8127faeddf40ac9e9e68) |
| DupCaller | Reported 1.0.5; commit `918572d4a9e5d3c14a7c01ffb513e8cbe69f352f` | [DupCaller](https://github.com/AlexandrovLab/DupCaller/tree/918572d4a9e5d3c14a7c01ffb513e8cbe69f352f) |
| HiDEF-seq workflow | 3.1, commit `15c522072cbc9ea04e521c75ac047029c53c2d23`; Nextflow 25.10.0 | [HiDEF-seq](https://github.com/evronylab/HiDEF-seq/tree/15c522072cbc9ea04e521c75ac047029c53c2d23) |
| HiDEF-seq container | Study YAML names `hidef-seq_3.0.sif`; distinct from workflow release | See populated `hidefseq/params.yaml` |
| SCAN2 | Conda `scan2=2.0=2`, `r-scan2=1.0`; full historical environment included | [SCAN2](https://github.com/parklab/SCAN2) |

The outer Snakemake launchers require the version-7 command-line interface.
The archived SCAN2 environment specifically records Snakemake 7.6.1, Python
3.8.13, GATK 4.2.6.1, SHAPEIT 2.r904 and SAMtools 1.15.1. Do not replace its
GATK with the Illumina workflow's version. Other caller environments should be
installed separately following upstream dependency instructions. The relative-path launcher helper additionally requires Python 3 and PyYAML 6.x.
Host tools
include Java, bgzip/tabix, SAMtools/BCFtools, Python with pandas/NumPy, and R for
MosaicForecast. DeepMosaic additionally needs ANNOVAR and its hg38 refGene data.

For MosaicHunter, obtain the parent source, apply `illumina/mosaichunter/12threads.patch`
from its repository root, and build using its upstream instructions. The local
source change adds `-threads=12` to BLAT; use a BLAT implementation supporting
that argument. The supplied patch excludes unrelated file-mode differences.

Custom resources still need a public archive/accession for a complete release:
MosaicHunter's lifted/combined masks, HiDEF-seq's custom BSgenome and masks,
DupCaller noise and germline inputs,
and any exact model/container artifact without an immutable public source.
Large resources are not copied into this code package. SRA accessions for reads
do not by themselves provide those derived resources.

ppmSeq calling software and its study version must be cited from the
[separate ppmSeq repository](ppmseq/README.md), whose URL is currently a
placeholder. Calling scripts and caller configurations are not bundled here.
Final cross-caller consensus and plotting-input generation belong to subsequent
benchmark processing.
