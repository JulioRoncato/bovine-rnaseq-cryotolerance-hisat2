#!/bin/bash

# ============================================
# RNA-seq pipeline - Bovine embryos
# Alignment with HISAT2
# ============================================

# Build index (once)
hisat2-build bovine_genome.fa bovine_index

# Alignment
hisat2 -p 8 \
  -x bovine_index \
  -1 sample_R1.fastq.gz \
  -2 sample_R2.fastq.gz \
  -S sample.sam

# Convert to BAM, sort and index
samtools view -bS sample.sam | samtools sort -o sample.sorted.bam
samtools index sample.sorted.bam

# Quantification
featureCounts -T 8 \
  -a bovine_annotation.gtf \
  -o counts.txt \
  sample.sorted.bam
