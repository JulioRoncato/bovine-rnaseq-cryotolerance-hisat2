#!/bin/bash

# RNA-seq based variant calling - Bovine embryos

# Build index (once)
hisat2-build bovine_genome.fa bovine_index

# Alignment
hisat2 --dta -p 8 \
  -x bovine_index \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -S sample.sam

# Process alignments
samtools view -bS sample.sam | samtools sort -o sample.sorted.bam
samtools index sample.sorted.bam

# Variant calling
bcftools mpileup -Ou -f bovine_genome.fa sample.sorted.bam | \
bcftools call -mv -Ov -o sample_variants.vcf
