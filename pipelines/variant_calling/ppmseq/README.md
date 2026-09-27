# ppmSeq variant calling — external repository

ppmSeq variant calling is maintained separately from `Benchmark_RIVER_paper`.

**Repository placeholder:** `https://github.com/OWNER/PPMSEQ_VARIANT_CALLING`

This URL is intentionally a placeholder, not an available or verified repository.
Replace it with the actual repository URL and specify the study version/commit,
authors, and archived release DOI before publication. The external repository
should document its software, models, resources, input requirements, and the
final denoising steps used to produce the benchmark callsets.

No ppmSeq calling scripts, Snakefiles, caller configurations, or calling sample
sheets are included in this package. Obtain finalized ppmSeq callsets from the
external workflow before running downstream benchmark preparation. Downstream
consumers expect the multiread and putative-multiread VCFs described in their
input paths, including the 30×/100×/200× and full-depth comparisons.

Read filtering, coverage calculation, consumption of existing ppmSeq calls, and
processed plotting tables remain in this repository. They are downstream
benchmark analyses, not a replacement ppmSeq caller.
