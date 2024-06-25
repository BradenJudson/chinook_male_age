setwd("~/chinook_male_age/analyses")

library(tidyverse); library(vegan)

source("../scripts/outliers.R")

# -------------------------------------------------------------------------
# need to correct for population structure in RDA (Z  =  PC1)
# Use BEAGLE to impute missing genotypes, not most common
# use PLINK to incorporate map file to understand genomic positions of outliers
# -------------------------------------------------------------------------



# Organize data ----------------------------------------------------------------

# Convert VCF to raw using PLINK and read in here.
geno <- read.table("../data/global_vcf/imputed/snps_maf005_raw.raw", 
                   header = T, sep = " ", row.names = 1, check.names = F) %>% 
  select(!c(1:5))

# geno <- read.table("../data/pop_vcfs/raw/Chilliwack_imp005_raw.raw", 
#                    header = T, sep = " ", row.names = 1) %>% 
#   select(starts_with("X"))


# Number of missing genotypes vs. non-missing genotypes. 
sum(is.na(geno)); sum(!is.na(geno))

# Imputation: replace missing genotypes with the most common allele at that position.
#gen_imp <- apply(geno, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

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
q.age.1  <- pop_rda(population = "Qualicum", phenotype = "age", maf = 1)
q.jack.1 <- pop_rda(population = "Qualicum", phenotype = "jack", maf = 1)
q.age.5  <- pop_rda(population = "Qualicum", phenotype = "age", maf = 5)
q.jack.5 <- pop_rda(population = "Qualicum", phenotype = "jack", maf = 5)

# Summarize Qualicum models.
(qualicum_models <- pop_rda_summ(rda_list = list(q.age.1, q.jack.1, q.age.5, q.jack.5)))

# Run all Puntledge models.
p.age.1  <- pop_rda(population = "Puntledge", phenotype = "age", maf = 1)
p.jack.1 <- pop_rda(population = "Puntledge", phenotype = "jack", maf = 1)
p.age.5  <- pop_rda(population = "Puntledge", phenotype = "age", maf = 5)
p.jack.5 <- pop_rda(population = "Puntledge", phenotype = "jack", maf = 5)

# Summarize Puntledge models.
(puntledge_models <- pop_rda_summ(rda_list = list(p.age.1, p.jack.1, p.age.5, p.jack.5)))

# Run all Chilliwack models.
c.age.1  <- pop_rda(population = "Chilliwack", phenotype = "age", maf = 1)
c.jack.1 <- pop_rda(population = "Chilliwack", phenotype = "jack", maf = 1)
c.age.5  <- pop_rda(population = "Chilliwack", phenotype = "age", maf = 5)
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


find_outliers <- function(rda_list, z) {
  
  outlier_snps <- rda_list %>% 
    lapply(., FUN = function(x) scores(x[["model"]], choices = 1)[[1]] %>% 
             `names<-`(., rownames(scores(x[["model"]], choices = 1)[[1]])) %>% 
             outliers(., 3))
  
  }

q.outliers <- find_outliers(rda_list = list(q.age.1, q.jack.1, q.age.5, q.jack.5), z = 3)

p.outliers <- find_outliers(rda_list = list(p.age.1, p.jack.1, p.age.5, p.jack.5), z = 3)

c.outliers <- find_outliers(rda_list = list(c.age.1, c.jack.1, c.age.5, c.jack.5), z = 3)


(local_rda_summary <- local_rda_summary %>% 
  mutate(n.outliers = c(unlist(lapply(q.outliers, length)),
                        unlist(lapply(p.outliers, length)),
                        unlist(lapply(c.outliers, length)))))


length(unique(names(q.outliers %>% unlist())))
length(unique(names(p.outliers %>% unlist())))
length(unique(names(c.outliers %>% unlist())))


















