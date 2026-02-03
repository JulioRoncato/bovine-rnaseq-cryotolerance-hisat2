library(clusterProfiler)
library(org.Bt.eg.db)
library(ggplot2)

set.seed(123)

variant_genes <- sample(
  keys(org.Bt.eg.db, keytype = "ENTREZID"),
  300
)

kegg_enrich <- enrichKEGG(
  gene = variant_genes,
  organism = "bta",
  pvalueCutoff = 0.3
)

if (nrow(kegg_enrich@result) > 0) {

  png("figures/pathway_enrichment.png", width = 900, height = 700)

  print(
    dotplot(kegg_enrich, showCategory = 10) +
      ggtitle("Pathway Enrichment of Genes Harboring Variants")
  )

  dev.off()
}
