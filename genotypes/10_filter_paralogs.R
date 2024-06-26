setwd("~/chinook_male_age/genotypes")

library(tidyverse); library(vcfR); set.seed(121)

# Load HDplot function from https://github.com/gjmckinney/HDplot
# Also see McKinney et al., 2017 (DOI: 10.1111/1755-0998.12613).
source("./HDplot/HDplot.R")

input <- read.vcfR("../data/global_vcf/snps_maf005.vcf.gz")

# Tresholds in Chinook described in McKinney et al., 2017.
Hmax <- 0.6
Dmax <- 7.0

HDplotResults <- HDplot(input) %>% 
  mutate(SNP = as.factor(case_when(H >  Hmax | abs(D) >  Dmax ~ "Duplicate",
                                   H <= Hmax | abs(D) <= Dmax ~ "Singleton")))

summary(HDplotResults$SNP); head(HDplotResults)

(HD_HD <- ggplot(data = HDplotResults) +
    geom_vline(xintercept = Hmax, linewidth = 1,
               color = "gray70") +
    geom_hline(yintercept = c(Dmax, Dmax*-1), 
               linewidth = 1, color = "gray70") +
    geom_point(aes(x = H, y = D, fill = SNP, 
                   color = SNP), 
               alpha = 2/5) +
    scale_color_manual(values = c("red3", "blue3")) +
    labs(x = "H", y = "D (read-ratio deviation)") +
    theme_bw() + theme(legend.position = "none",
                       panel.grid = element_blank()))

ggsave("../plots/hdplot_fig.tiff", dpi = 300,
       width = 10, height = 8)

# Write putative paralog positions to directory. 
# Exclude using vcftools --exclude-positions
write.table(HDplotResults[HDplotResults$SNP == "Duplicate", c("CHROM", "POS")],
            "./stats/duplicateSNP_IDs_maf005.txt", sep = "\t", quote = F, 
            row.names = F, col.names = F)
