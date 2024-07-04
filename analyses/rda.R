setwd("~/chinook_male_age/analyses")

library(tidyverse); library(vegan); library(cowplot)

source("../scripts/outliers.R")

# -------------------------------------------------------------------------
# need to correct for population structure in RDA (Z  =  PC1)
# Use BEAGLE to impute missing genotypes, not most common
# use PLINK to incorporate map file to understand genomic positions of outliers
# Manhattan plots for RDA1
# Global RDA with either PC1 or Q values as conditioning
# Only use imputed datasets where needed (e.g., RDA and RFs)
# RDAs with only the outlier loci
# -------------------------------------------------------------------------



# Organize data ----------------------------------------------------------------

# Convert VCF to raw using PLINK and read in here.
geno <- read.table("../data/global_vcf/imputed/snps_maf005_raw.raw", 
                   header = T, sep = " ", row.names = 1, check.names = F) %>% 
  select(!c(1:5))

# Number of missing genotypes vs. non-missing genotypes. 
sum(is.na(geno)); sum(!is.na(geno))

# For keeping sample order consistent throughout.
samples <- as.factor(rownames(geno))

# Read in and format phenotypic information.
info <- read.csv("../data/chRADseq_samples.csv") %>%
  select(c("fishID", "pop", "age", "year")) %>% 
  filter(fishID %in% samples) %>% 
  mutate(pop = as.factor(sub("\\_.*", "", 
                         tools::toTitleCase(tolower(pop)))),
         jack = as.factor(case_when(age == 2 ~ 1,
                          age >  2 ~ 0))) %>% 
  arrange(factor(fishID, levels = samples))

# Read in individual PC scores representing population structure.
indpc <- read.csv("../data/rad_pca_scores.csv") %>% 
  filter(X %in% samples) %>% 
  arrange(factor(X, levels = samples)) %>% 
  column_to_rownames("X")

################################################################################
# Global RDAs ------------------------------------------------------------------
################################################################################

# Global RDA 1: Jacks ----------------------------------------------------------

jack.rda <- rda(geno ~ jack, Z = indpc$PC1, data = info, scale = T)
RsquareAdj(jack.rda)
plot(jack.rda, scaling = 3)

j.rdaload <- scores(jack.rda, choices = 1, display = "species")

# F = 1.0271, Pr(>F) = 0.255, Var = 44.7.
(jack_sig <- anova.cca(jack.rda, parallel=getOption("mc.cores")))

hist(j.rdaload[,1], main = NULL)


# Global RDA 2: Age ------------------------------------------------------------

age.rda <- rda(geno ~ age, Z = indpc$PC1, data = info, scale = T)
RsquareAdj(age.rda); screeplot(age.rda)
plot(age.rda, scaling = 3)

# F = 1.0273, Pr(>F) = 0.285, Var = 44.7.
(age_sig <- anova.cca(age.rda, parallel=getOption("mc.cores")))

a.rdaload <- scores(age.rda, choices = 1, display = "species")
hist(a.rdaload[,1], main = NULL)

hist(scores(age.rda, choices = 1, display = "species")[,1], main = NULL)



################################################################################
# Local RDAs -------------------------------------------------------------------
################################################################################

# Function for performing RDAs within each population.
pop_rda <- function(population, phenotype, maf) {
  
  # Read in raw genotype information for each population. 
  # No missing data because these were imputed w/ Beagle.
  geno <- read.table(paste0("../data/pop_vcfs/raw/", population,
                            "_imp", sprintf("%03d", maf), "_raw.raw"),
                     header = T, sep = " ", row.names = 1,
                     check.names = FALSE) %>% select(!c(1:5)) 
                     # Discard unnecessary columns.

  # Samples within each population.
  samples <- as.factor(rownames(geno))

  # Print overview of each dataset going in.
  print(paste0("Population: ", population, ", n = ", length(samples)), quote = F)

  # Isolate phenotypes for each indv. within each pop.
  info <- read.csv("../data/chRADseq_samples.csv") %>%
    select(c("fishID", "pop", "age", "year")) %>%
    filter(fishID %in% samples) %>%
    mutate(pop = as.factor(sub("\\_.*", "", tools::toTitleCase(tolower(pop)))),
           jack = as.factor(case_when(age == 2 ~ 1, age >  2 ~ 0))) %>%
    arrange(factor(fishID, levels = samples))

  # Run the RDA. 
  model <- rda(geno ~ info[,phenotype],
               data = info, scale = T)
  
  # Isolate model R2. 
  (rsq  <- RsquareAdj(model)[["r.squared"]])
  
  # Significance of the RDA.
  (sig <- anova.cca(model, parallel=getOption("mc.cores")))
  
  # Specify model outputs and names and return as a list. 
  output <- list(
    model = model, R2  = rsq,
    anova = sig,   maf = maf,
    phenotype = phenotype,
    samples = nrow(info)
  ); return(output)
  
}


# Write a function to extract model-specific information from above.
pop_rda_summ <- function(rda_list) {
  modsumm <- data.frame(
    phenotype = unlist(lapply(rda_list, function(x) x$phenotype)), 
    maf       = unlist(lapply(rda_list, function(x) x$maf)),
    R2        = unlist(lapply(rda_list, function(x) x$R2))) %>% 
    bind_cols(bind_rows(lapply(rda_list, function(x) broom::tidy(x$anova)[1,])))
}

# Run all Qualicum models.
q.age.1  <- pop_rda(population = "Qualicum", phenotype = "age",  maf = 1)
q.jack.1 <- pop_rda(population = "Qualicum", phenotype = "jack", maf = 1)
q.age.5  <- pop_rda(population = "Qualicum", phenotype = "age",  maf = 5)
q.jack.5 <- pop_rda(population = "Qualicum", phenotype = "jack", maf = 5)

# Summarize Qualicum models.
(qualicum_models <- pop_rda_summ(rda_list = list(q.age.1, q.jack.1, q.age.5, q.jack.5)))

# Run all Puntledge models.
p.age.1  <- pop_rda(population = "Puntledge", phenotype = "age",  maf = 1)
p.jack.1 <- pop_rda(population = "Puntledge", phenotype = "jack", maf = 1)
p.age.5  <- pop_rda(population = "Puntledge", phenotype = "age",  maf = 5)
p.jack.5 <- pop_rda(population = "Puntledge", phenotype = "jack", maf = 5)

# Summarize Puntledge models.
(puntledge_models <- pop_rda_summ(rda_list = list(p.age.1, p.jack.1, p.age.5, p.jack.5)))

# Run all Chilliwack models.
c.age.1  <- pop_rda(population = "Chilliwack", phenotype = "age",  maf = 1)
c.jack.1 <- pop_rda(population = "Chilliwack", phenotype = "jack", maf = 1)
c.age.5  <- pop_rda(population = "Chilliwack", phenotype = "age",  maf = 5)
c.jack.5 <- pop_rda(population = "Chilliwack", phenotype = "jack", maf = 5)

# Summarize Puntledge models.
(chilliwack_models <- pop_rda_summ(rda_list = list(c.age.1, c.jack.1, c.age.5, c.jack.5)))

# Population and model-specific summary stats. 
(local_rda_summary <- bind_rows(
  qualicum_models %>% mutate(pop = "Qualicum"),
  puntledge_models %>% mutate(pop = "Puntledge"),
  chilliwack_models %>% mutate(pop = "Chilliwack")) %>% 
    select(c(pop, phenotype, maf, Variance, statistic, p.value)))

write.csv(local_rda_summary, "../data/local_rda_summary.csv", row.names = F)


# SNP Loadings -----------------------------------------------------------------

# SNP identities and locations from the parent VCF (w/  maf = 1% to not miss any).
map <- read.table("../data/global_vcf/imputed/global_maf001.map")[,c(1,2,4)] %>% 
  `colnames<-`(., c("chr", "SNP", "pos"))

# Expand on Forester's outlier function to build a dataframe of outlier SNPs
# and their respective RDA loadings per model per population. Also use the
# parent map file above to assign genomic positions to each locus.
find_outliers <- \(rda_list, z) {
  
  # z = number of standard deviations away from the mean loading.
  outlier_snps <- rda_list %>% 
    lapply(., FUN = \(x) scores(x[["model"]], choices = 1)[[1]] %>% 
             `names<-`(., rownames(scores(x[["model"]], choices = 1)[[1]])) %>% 
             outliers(., z) %>% as.data.frame() %>% rownames_to_column() %>% 
             `colnames<-`(., c("SNP", "loading")) %>% 
             mutate(SNP = str_sub(SNP, end = -3)) %>% 
             merge(., map, by = "SNP"))
  
  } 

# Isolate outlier loci for Qualicum.
q.outliers <- find_outliers(rda_list = tibble::lst(q.age.1, q.jack.1, q.age.5, q.jack.5), z = 3)

# Isolate outlier loci for Puntledge.
p.outliers <- find_outliers(rda_list = tibble::lst(p.age.1, p.jack.1, p.age.5, p.jack.5), z = 3)

# Isolate outlier loci for Chilliwack.
c.outliers <- find_outliers(rda_list = tibble::lst(c.age.1, c.jack.1, c.age.5, c.jack.5), z = 3)

# Tabulate SNP outliers for each model and add to summary table.
(local_rda_summary <- local_rda_summary %>% 
  mutate(n.outliers = c(unlist(lapply(q.outliers, nrow)),
                        unlist(lapply(p.outliers, nrow)),
                        unlist(lapply(c.outliers, nrow)))))


(outliers_full <- bind_rows(c(q.outliers, p.outliers, c.outliers),
                            .id = "dataset") %>% group_by(chr) %>%
                            tally())


plot_outliers <- \(model) {
  
  model_name <- rlang::as_label(eval(parse(text=enquo(model)))[[2]])
  print(model_name)
  
  outlier_list <- paste0(substr(model_name, 0, 1), ".outliers")
  
  j <- get(x = outlier_list, envir = .GlobalEnv)
  
  m <- get(x = model_name,   envir = .GlobalEnv)
  
  outlier_snps <- as.vector(j[[model_name]]$SNP)
  print(outlier_snps)
  
  df <- data.frame(m$model$CCA$v) %>% 
     mutate(PC1 = model$model$CA$v[,1],
            SNP = str_sub(rownames(.), end = -3),
            out = case_when(
               SNP %in% outlier_snps ~ "Y",
              !SNP %in% outlier_snps ~ "N"
            ))

   (rda_plot <- ggplot() +
       geom_point(data = df[df$out == "N", ],
                  aes(y = PC1, x = RDA1), 
                  shape = 21, size = 2, fill = 'gray95', 
                  colour = 'gray60', alpha = 4/5) +
     geom_point(data = df[df$out == "Y", ],
                aes(y = PC1, x = RDA1), 
                shape = 21, size = 2, fill = '#1f78b4') +
       # scale_fill_manual(values = c('#f1eef6', '#1f78b4')) +
       annotate("segment", xend = min(df$RDA1), x = 0,
                yend = 0, y = 0, size = 1, colour = "black",
                arrow = arrow(type = "open", length = unit(0.02, "npc"))) +
       annotate("text", x = min(df$RDA1), y = 0.002, 
                label = tools::toTitleCase(c(m$phenotype))) +
       theme_bw() + theme(legend.position = "none") +
       ggtitle(paste0("maf = ", model$maf/100)))

  return(rda_plot)
  
}

(qualicum_rdas <- cowplot::plot_grid(plotlist = list(
  plot_outliers(q.age.1), plot_outliers(q.jack.1),
  plot_outliers(q.age.5), plot_outliers(q.jack.5)),
  ncol = 2, nrow = 2))

ggsave("../plots/qualicum_rdas.tiff", dpi = 300,
       height = 10, width = 12)

(puntledge_rdas <- cowplot::plot_grid(plotlist = list(
  plot_outliers(p.age.1), plot_outliers(p.jack.1),
  plot_outliers(p.age.5), plot_outliers(p.jack.5)),
  ncol = 2, nrow = 2))

ggsave("../plots/puntledge_rdas.tiff", dpi = 300,
       height = 10, width = 12)

(chilliwack_rdas <- cowplot::plot_grid(plotlist = list(
  plot_outliers(c.age.1), plot_outliers(c.jack.1),
  plot_outliers(c.age.5), plot_outliers(c.jack.5)),
  ncol = 2, nrow = 2))

ggsave("../plots/chilliwack_rdas.tiff", dpi = 300,
       height = 10, width = 12)



# RDA Manhattan plots ----------------------------------------------------------

# Illustrate where RDA outliers occur throughout the genome.
rda_manhattan <- \(model) {
  
  model_name <- rlang::as_label(eval(parse(text=enquo(model)))[[2]])
  print(model_name)
  
  outlier_list <- paste0(substr(model_name, 0, 1), ".outliers")
  
  j <- get(x = outlier_list, envir = .GlobalEnv)
  
  m <- get(x = model_name,   envir = .GlobalEnv)
  
  outlier_snps <- as.vector(j[[model_name]]$SNP)
  print(outlier_snps)
  
  snps <- data.frame(m$model$CCA$v) %>% 
    mutate(SNP = str_sub(rownames(.), end = -3),
           out = case_when(
             SNP %in% outlier_snps ~ "Y",
             !SNP %in% outlier_snps ~ "N"
           ))  %>% merge(., map, by = "SNP")
  
  df <- data.frame(m$model$CCA$v) %>% 
    mutate(SNP = str_sub(rownames(.), end = -3)) %>% 
    merge(., map, by = "SNP") %>% 
    group_by(chr) %>% 
    summarise(chr_len = max(pos)) %>% 
    mutate(tot = cumsum(chr_len) - chr_len,
           chr = as.factor(chr)) %>%
    select(-chr_len) %>% 
    left_join(snps, ., by = c("chr" = "chr")) %>% 
    arrange(chr, pos) %>% 
    mutate(bpcum = pos + tot)
  
  axisdf <- df %>%
    mutate(LG = paste0("Ots", 
                sprintf("%02d",as.numeric(as.factor(chr))))) %>% 
    group_by(LG) %>% 
    summarize(m = min(bpcum))
  
  zc_manhattan <- ggplot() +
    geom_point(data = df[df$out == "N",],
               aes(x = bpcum, y = RDA1,
                   colour = as.factor(chr))) +
    scale_colour_manual(values = rep(c("gray80", "gray30"), 17)) +
    geom_point(data = df[df$out == "Y",],
               aes(x = bpcum, y = RDA1), colour = "red2") +
    theme_classic() + 
    geom_hline(yintercept = 0, alpha = 1/3) +
    scale_x_continuous(breaks = axisdf$m,labels = axisdf$LG) +
    labs(x = NULL) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 90, vjust = 0.5))

  return(zc_manhattan)
  
}

# Qualicum RDA Manhattan:
cowplot::plot_grid(plotlist = list(rda_manhattan(q.age.5), rda_manhattan(q.jack.5)),
                   ncol = 1, labels = c("a", "b"))
ggsave("../plots/qualicum_rda_manhattan.tiff", dpi = 300, width  = 12, height = 12)

# Puntledge RDA Manhattan:
cowplot::plot_grid(plotlist = list(rda_manhattan(p.age.5), rda_manhattan(p.jack.5)),
                   ncol = 1, labels = c("a", "b"))
ggsave("../plots/puntledge_rda_manhattan.tiff", dpi = 300, width  = 12, height = 12)

# Chilliwack RDA Manhattan:
cowplot::plot_grid(plotlist = list(rda_manhattan(c.age.5), rda_manhattan(c.jack.5)),
                   ncol = 1, labels = c("a", "b"))
ggsave("../plots/chilliwack_rda_manhattan.tiff", dpi = 300, width  = 12, height = 12)


# RDAs Part II: Outliers only --------------------------------------------------

# Print outlier SNP positions. 






