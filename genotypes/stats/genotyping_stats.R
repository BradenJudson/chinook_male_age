setwd("~/chinook_male_age/genotyping/stats")

library(tidyverse); library(ggExtra); library(cowplot); set.seed(2140)

reads <- read.table("../../data/sample_reads.txt", col.names = c("bam", "reads")) %>% 
  mutate(fishID = gsub(".1.sorted.bam", "", bam))

samples <- read.csv("../../data/chRADseq_samples.csv") %>% 
  mutate(fishID = gsub("[[:space:]]Chelex", "", fishID))

dat <- merge(reads, samples[,c(1,5)]) %>% 
  mutate(Population = tools::toTitleCase(tolower(gsub("\\_.*", "", pop))))

(bps <- ggplot(data = dat, aes(x = Population, y = reads/1e6)) +
  geom_boxplot(alpha = 1/10, outlier.colour = NA, fill = "grey80") +
  geom_jitter(size  = 3/2, shape = 21, color = "black",
             fill = 'grey80', width = 0.08) +
  labs(x = NULL, y = "Reads (Millions)") +
  theme_bw())

(mds <- ggMarginal(bps, type = "density", linewidth = 1,
                   fill = "gray90", alpha = 0.7))

save_plot("../../plots/read_dist.tiff", mds, ncol = 2, 
          base_height = 6, base_asp = 4/5)

# Global summary of reads.
summary(dat$reads)

# Local summary of reads.
dat %>% 
  group_by(Population) %>% 
  summarise(mean = mean(reads),
            sd   = sd(reads),
            min  = min(reads),
            max  = max(reads)) %>% 
  as.data.frame(row.names = F)


# Unfiltered stats -------------------------------------------------------------

# Individual missingness -------------------------------------------------------

imiss <- read.table("unfiltered_stats/unfiltered_vcf.imiss", header = T)

hist(imiss$F_MISS); summary(imiss$F_MISS)

pass <- imiss[imiss$F_MISS < 1/4, ]


# Individual depth -------------------------------------------------------------

idepth <- read.table("unfiltered_stats/unfiltered_vcf.idepth", header = T)

hist(idepth$MEAN_DEPTH); summary(idepth$MEAN_DEPTH)

ind_stat <- merge(imiss[,c(1,5)], idepth[,c(1,3)]) %>% 
  pivot_longer(cols = c(F_MISS, MEAN_DEPTH))

ggplot(data = ind_stat, aes(x = value)) +
  geom_histogram(colour = "black", fill = "grey90") + 
  facet_wrap(~name, scales = "free") + theme_bw() +
  labs(x = NULL, y = "Count")

ggsave("../../plots/indv_stats_unfiltered.tiff", dpi = 300,
       width = 14, height = 7)


# SNP depth --------------------------------------------------------------------

ldepth <- read.table("unfiltered_stats/unfiltered_vcf.ldepth.mean", header = T)

hist(ldepth$MEAN_DEPTH); summary(ldepth$MEAN_DEPTH)

# SNP Missingness --------------------------------------------------------------

lmiss <- read.table("unfiltered_stats/unfiltered_vcf.lmiss", header = T)

hist(lmiss$F_MISS); summary(lmiss$F_MISS)

snpdat <- merge(lmiss[,c(1:2,6)], ldepth[,c(1:3)], by = c(1:2)) %>% 
  pivot_longer(cols = c(F_MISS, MEAN_DEPTH))

ggplot(data = snpdat, aes(x = value)) +
  geom_histogram(colour = "black", fill = "grey90") + 
  facet_wrap(~name, scales = "free") + theme_bw() +
  labs(x = NULL, y = "Count")

ggsave("../../plots/snp_stats_unfiltered.tiff", dpi = 300,
       width = 14, height = 7)



