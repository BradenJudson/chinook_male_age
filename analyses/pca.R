setwd("~/chinook_male_age/analyses")

library(tidyverse); library(vcfR); library(ggplot2); library(adegenet)
library(SeqArray); library(poppr); library(SNPRelate); library(GENESIS);
library(GGally); library(SeqVarTools)

nc <- 12 # Number of cores to use in parallel functions.

# Read in genetic data and convert to genlight object.
R.utils::gunzip("../data/global_vcf/snps_maf005_singletons.vcf.gz")
rad <- read.vcfR("../data/global_vcf/snps_maf001_singletons.vcf")
(radgl <- vcfR2genlight(rad))

# Site info
samples <- as.factor(c(radgl$ind.names))

info <- read.csv("../data/chRADseq_sampleinfo.csv") %>% 
  filter(Sample %in% c(samples)) %>% 
  arrange(factor(Sample, levels = samples)) %>% 
  mutate(Population = tools::toTitleCase(tolower(Population)))

radgl@pop <- as.factor(info$Population)


# PCA --------------------------------------------------------------------------

rad_pca <- glPca(radgl, nf = 3, parallel = T, n.cores = 14)
write.csv(rad_pca$scores, "rad_pca_scores.csv", row.names = TRUE)

pca_scores <- read.csv("rad_pca_scores.csv")

# Isolate PC scores.
pc_scores <- as.data.frame(pca_scores) %>% 
  rename(Sample = X) %>% 
  merge(., info, by = "Sample")

# Isolate eigenvalues (% variation explained for each PC axis).
(pc_var <- c(rad_pca$eig/sum(rad_pca$eig)*100)[1:5]); barplot(rad_pca$eig)

(fullpca <- ggplot(data = pc_scores, aes(x = PC1, y = PC2)) +
    geom_point(aes(fill = Population), colour = "black",
               shape = 21, size = 3/2, alpha = 3/5) +
    labs(x = paste0("PC1 (", round(pc_var[1], 1), "%)"),
         y = paste0("PC2 (", round(pc_var[2], 1), "%)")) +
    theme_bw() +
    theme(legend.title = element_blank(),
          legend.position = "top"))

ggsave("../plots/rad_pca.tiff", dpi = 300, width = 8, height = 6)

