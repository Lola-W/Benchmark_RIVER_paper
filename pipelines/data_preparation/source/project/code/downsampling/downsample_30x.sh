
samtools view -@ 12 -s 0.3 -b ./inputs/alignments/7614-Cortex-L-T_D.bam > ./work/project/data/donwsampling/udseq/30x_7614-Cortex-L-T_D.bam
samtools index -@ 12 ./work/project/data/donwsampling/udseq/30x_7614-Cortex-L-T_D.bam

# samtools view -@ 12 -s 42.3 -T ./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa \
#  -C ./inputs/alignments/422073-25-5852-DNA-1-ppm0021-CGCATCCTCACAGAT.cram > ./work/project/data/donwsampling/ppmseq/100x_7614.cram

# samtools index -@ 12 ./work/project/data/donwsampling/ppmseq/100x_7614.cram

samtools view -@ 12 -s 42.6 \
  -T ./resources/GRCh38_full_analysis_set_plus_decoy_hla.fa \
  -O cram,version=3.0 \
  -o ./work/project/data/donwsampling/ppmseq/200x_7614.cram \
  ./inputs/alignments/422077-25-5853-DNA-1-ppm0022-CACAACATATCAGAT.cram
samtools index -@ 12 ./work/project/data/donwsampling/ppmseq/200x_7614.cram



samtools view -@ 24 -s 4.6 -b ./work/project/data/illumina_hg38_bam/7614_illumina_brain_somatic.bam > ./work/project/data/illumina_hg38_bam/illumina_200x.bam
samtools index -@ 24 ./work/project/data/illumina_hg38_bam/illumina_200x.bam