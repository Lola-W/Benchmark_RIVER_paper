# Alignment placeholders

ALIGNMENT_PLACEHOLDERS.tsv lists literal BAM/CRAM paths plus the 16 PTA BAMs
constructed in loops. Paths are relative to data_preparation/. No alignments or
empty stand-in files are included. Supply aligned regional/unrelated controls
with the assembly required by their script, and BAM/CRAM indexes.

Paths under work/ may be inputs to a stage or outputs of earlier stages.
Additional cleaned/intermediate BAM names are constructed from variables in
scripts; this inventory does not replace checking the selected stage's inputs.
Existing alignments may be linked; do not rename or copy large sequencing files
merely to duplicate the sibling variant_calling input placeholders.
