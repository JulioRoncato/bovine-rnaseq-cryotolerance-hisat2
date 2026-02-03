# 🧬 Variant Calling Analysis in Bovine Embryos with Differential Cryotolerance

This repository presents a genomic variant analysis pipeline developed during my Master's research, focusing on the identification and functional interpretation of genetic variants associated with cryotolerance in bovine embryos.

---

## 🔬 Biological Background

Cryotolerance is a key factor influencing the success of bovine embryo cryopreservation in assisted reproduction programs. Genetic variation can affect cellular mechanisms involved in stress response, membrane integrity, and energy metabolism, ultimately impacting embryo survival after freezing and thawing.

This project investigates genomic variants associated with **high** and **low cryotolerance phenotypes** in bovine embryos and explores the biological pathways potentially impacted by these variants.

---

## 🧪 Experimental Design

Whole-genome sequencing data were obtained from bovine embryos classified into two phenotypic groups:

| Phenotype | Description |
|---------|-------------|
| High Cryotolerance | Embryos with high post-thaw survival |
| Low Cryotolerance | Embryos with low post-thaw survival |

Variants were identified and compared between groups to detect candidate genes and pathways associated with cryotolerance.

---

## ⚙️ Bioinformatics Workflow

### 1️⃣ Variant Calling Pipeline (Linux)

Variant discovery was performed in a Linux environment and included:

- Read alignment to the bovine reference genome
- Sorting and indexing of BAM files
- Variant calling
- Variant filtering and quality control
- Generation of final VCF files

---

### 2️⃣ Variant Annotation and Gene Mapping

Identified variants were annotated to genes and genomic features, enabling the identification of genes harboring potentially impactful variants.

---

### 3️⃣ Functional Enrichment Analysis (R)

Genes associated with variants were subjected to functional enrichment analysis in **R**, including:

- Gene Ontology (GO) enrichment
- Pathway analysis (e.g., KEGG, Reactome)
- Visualization of enriched biological processes

---

## 📊 Data Visualization

- Variant distribution plots
- Summary statistics of variant types
- Functional enrichment plots
- Pathway-level visualizations

---

## 🔐 Data Availability and Reproducibility

Due to data confidentiality, the datasets provided in this repository are simulated or derived from publicly available bovine genomic data. The analytical workflow, scripts, and parameters faithfully reproduce the original analysis.

---

## 🛠️ Tools and Packages

- Linux-based variant calling tools
- R
- clusterProfiler
- ggplot2
- biomaRt

---

## 📌 Key Takeaways

This project demonstrates:
- A complete genomic variant analysis pipeline
- Integration of Linux-based pipelines with R-based functional analysis
- Application of genomics to reproductive biology
- Biological interpretation of variant-driven gene and pathway effects
