
samtools view -@ 12 -s 0.3 -b ./inputs/alignments/7614-Kidney-L_D.bam > ./work/project/data/donwsampling/udseq/30x_7614-Kidney-L_D.bam
samtools index -@ 12 ./work/project/data/donwsampling/udseq/30x_7614-Kidney-L_D.bam
# samtools view -@ 12 -s 42.3 \
#   -T ./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa \
#   -O cram,version=3.0 \
#   -o ./work/project/data/donwsampling/ppmseq/100x_7614.cram \
#   ./inputs/alignments/422077-25-5853-DNA-1-ppm0022-CACAACATATCAGAT.cram


# samtools index -@ 12 ./work/project/data/donwsampling/ppmseq/100x_7614.cram

