Place indexed, coordinate-sorted hg38 BAMs here (symlinks are sufficient):

- `bulk.bam` and `bulk.bam.bai`: matched bulk control.
- `7614_single_neuron_1.sc.bam` and `.bam.bai`, continuing through neuron 16.

Run from `pta/`; the Snakefile discovers `input/{sample}.sc.bam`.
The hg38 resource selection is a submission adaptation of the shared SCAN2
wrapper, as confirmed by the study author. No hg38 SCAN2 run was performed
or independently validated during packaging. Its `--shapeit-refpanel` option
requires a build-matched SHAPEIT panel; the old configuration's separate
`eagle_1000g_panel` entry is not interchangeable with it.
