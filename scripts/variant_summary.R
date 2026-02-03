library(VariantAnnotation)
library(ggplot2)

vcf <- readVcf("data/simulated_variants.vcf")

variant_types <- table(info(vcf)$TYPE)

df <- as.data.frame(variant_types)
colnames(df) <- c("VariantType", "Count")

p <- ggplot(df, aes(x = VariantType, y = Count)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Distribution of Variant Types")

ggsave("figures/variant_distribution.png", p, width = 7, height = 5)
