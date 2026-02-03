# 🧬 RNA-seq Analysis of Bovine Embryos with Differential Cryotolerance

This repository contains a complete RNA-seq analysis pipeline developed during my Master's research, focusing on transcriptional and functional differences between bovine embryos with high and low cryotolerance.

---

## 🔬 Biological Background

Embryo cryotolerance is a critical factor affecting the success of cryopreservation in bovine assisted reproduction. Molecular mechanisms related to cellular stress response, energy metabolism, membrane integrity, and mitochondrial function play key roles in embryo survival after freezing and thawing.

This project investigates transcriptomic differences associated with **high** and **low cryotolerance phenotypes** in bovine embryos using RNA sequencing and functional enrichment analysis.

---

## 🧪 Experimental Design

RNA-seq data were generated from bovine embryos classified into two experimental groups:

| Phenotype | Description |
|---------|-------------|
| High Cryotolerance | Embryos with high post-thaw survival |
| Low Cryotolerance | Embryos with low post-thaw survival |

---

## ⚙️ Bioinformatics Workflow

### 1️⃣ RNA-seq Alignment and Quantification (Linux / HISAT2)

RNA-seq preprocessing was performed in a Linux environment and included:

- Quality control of raw FASTQ files  
- Alignment to the bovine reference genome using **HISAT2**  
- Sorting and indexing of BAM files  
- Gene-level read quantification using `featureCounts`  

---

### 2️⃣ Differential Expression Analysis (R / DESeq2)

Downstream analyses were conducted in **R** using **DESeq2**, including:

- Import of raw gene count matrices
- Data normalization (VST)
- Exploratory analysis (PCA)
- Differential gene expression analysis between cryotolerance groups

---

## 🧠 Functional Enrichment Analysis

Differentially expressed genes were subjected to functional enrichment analysis using **GO** and **KEGG** databases to identify biological processes and pathways associated with embryo cryotolerance.

---

## 🔐 Data Availability

Due to data confidentiality, the datasets provided in this repository are simulated or derived from publicly available bovine RNA-seq data. All analytical steps reproduce the original workflow.

---

## 🛠️ Tools and Packages

- HISAT2
- featureCounts
- R
- DESeq2
- clusterProfiler
- ggplot2
- org.Bt.eg.db

---

## 📌 Key Takeaways

This project demonstrates a complete RNA-seq pipeline using HISAT2, integrating transcriptomic analysis with biological interpretation of pathways associated with embryo cryotolerance.
