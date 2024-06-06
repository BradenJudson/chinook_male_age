setwd("~/chinook_male_age/analyses")

library(tidyverse); library(vcfR); library(ggplot2); library(adegenet)
library(SeqArray); library(purrr); library(SNPRelate); library(GENESIS);
library(GGally); library(SeqVarTools); library(GWASTools)

nc <- 12 # Number of cores to use in parallel functions.

# Read in genetic data and convert to genlight object.
# gunzip("../data/snps_maf001_singletons.vcf.gz")
rad <- read.vcfR("../data/snps_maf001_singletons_sub.recode.vcf")
radgl <- vcfR2genlight(rad)

# Site info
samples <- as.factor(c(radgl$ind.names))

info <- read.csv("../data/chRADseq_sampleinfo.csv") %>% 
  filter(Sample %in% c(samples)) %>% 
  arrange(factor(Sample, levels = samples)) %>% 
  mutate(Population = tools::toTitleCase(tolower(Population)))

radgl@pop <- as.factor(info$Population)

radgl <- radgl[!indNames(radgl) %in% c("25300", "223831", "221434")]

subsamples <- as.factor(c(radgl@ind.names))

# write.table(x = info[,c("Population", "Sample")], "../data/sample_pops_n707.txt",
#             row.names = FALSE, col.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------------------

# https://uw-gac.github.io/SISG_2019/pc-relate.html
# https://bioconductor.org/packages/devel/bioc/vignettes/GENESIS/inst/doc/assoc_test.html

# Kinship ----------------------------------------------------------------------

# To use KING algorithm, convert VCF to GDS object.
rad_gds <- snpgdsVCF2GDS(vcf.fn = "../data/snps_maf001_singletons_sub.recode.vcf",
                         out.fn = "../data/snps_maf001_sub.gds")

g <- snpgdsOpen("../data/snps_maf001_sub.gds", allow.duplicate = T)

# Conduct KING analysis using file generated above. 
# Must specify autosomes.only = F otherwise all SNPs are discarded.
    # All SNPs already have known positions on nuclear autosomes. 
king <- snpgdsIBDKING(gdsobj = g, autosome.only = F,
                      sample.id = indNames(radgl))


# Conduct a PCA on the kinship matrix and add population-level information.
kpcadf <- as.data.frame(dudi.pca(king$kinship, nf = 2, scannf = F, 
                               center = T, scale = T)$li) %>% 
  mutate(Sample = indNames(radgl)) %>% 
  merge(., info, by = "Sample")

# Visualize the PCA of the kinship matrix.
(kinshippca <- ggplot(data = kpcadf, 
       aes(x = Axis1, y = Axis2, 
           colour = Population)) + 
  geom_point(size = 3/2, alpha = 2/3) +
  theme_bw() +
    labs(x = "PC1", y = "PC2"))

# Isolate kinship matrix and label individuals.
kmat <- king$kinship
rownames(kmat) <- colnames(kmat) <- indNames(radgl)

kpc <- pcairPartition(kinobj = kmat, divobj = kmat)
names(kpc); sapply(kpc, length)

# Run a PCA on the unrelated set of individuals and use those data
# to set project values for the related set of individuals. 
unrel_pca <- snpgdsPCA(gdsobj = g, sample.id = kpc$unrels, autosome.only = FALSE)
snp_load  <- snpgdsPCASNPLoading(unrel_pca, gdsobj = g)
samp_load <- snpgdsPCASampLoading(loadobj = snp_load, gdsobj = g, sample.id = kpc$rels)

# Combined unrelated and related PCs, order samples and reformat. 
pcs <- rbind(unrel_pca$eigenvect, samp_load$eigenvect) %>% 
  `rownames<-`(c(unrel_pca$sample.id, samp_load$sample.id)) %>% 
  merge(., info[,c(4,3)], by.x = 0, by.y = "Sample") %>% 
  rename("Sample" = Row.names) %>% relocate("Population", 1) %>% 
  `colnames<-`(., gsub("V", "PC", colnames(.))) %>% 
  column_to_rownames(var = "Sample")

# Plot data from above using GGally.
# Looks like good separation along PC1, but separation along PCs > 2 is poor.
# Only plot first 12 PCs as an example.
(pclp <- ggparcoord(pcs, columns = 2:14, 
           groupColumn = "Population", 
           scale = "uniminmax") +
  theme_bw() + labs(x = NULL, y = "PC Score")) +
  theme(legend.position = "top", legend.title = element_blank())

ggsave("../plots/pc_air_lineplot.tiff", dpi = 300, width = 12, height = 8)


# In preperation of using PC-Relate some file (re-)formatting is required.
gdsfmt::showfile.gds(closeall = T)
SeqArray::seqSNP2GDS(gds.fn = "../data/snps_maf001.gds",
                     out.fn = "../data/snps_seqmaf001.gds")
seqDat <- SeqVarData("../data/snps_seqmaf001.gds")
PC1 <- as.matrix(pcs[,2]); rownames(PC1) <- rownames(pcs)
iterator <- SeqVarBlockIterator(seqData = seqDat)

# Run PCRelate using the most informative PC identified above, specifying to use all samples but
# only unrelated samples are used to "train" the calculation of PC beta coefficients.
pcrel <- pcrelate(iterator, pcs = PC1, sample.include=row.names(pcs), training.set = kpc$unrels)

# Convert the above information into a pairwise matrix (i.e., a `pcrelate` object).
# Reset filtering of base file for next steps.
pcrelMat <- pcrelateToMatrix(pcrel); seqResetFilter(seqDat, verbose = FALSE)

# Perform the PCAir analysis finally.
# Include relatedness and divergence matrices, include individual IDs.
kinrel <- pcair(seqDat, 
                kinobj = pcrelMat,
                divobj = kmat,
                sample.include = rownames(pcs),
                autosome.only = FALSE)

# Extract data from above.
pcair_dat <- as.data.frame(kinrel$vectors) %>% 
  `colnames<-`(., c(paste0("PC", 1:ncol(.)))) %>% 
  merge(., info[,c(4,3)], by.x = 0, by.y = "Sample") %>% 
  rename("Sample" = Row.names)

# Plot the data.First PC separates samples.
(pcair_vis <- ggplot(data = pcair_dat, 
       aes(x = PC1, colour = Population, 
           fill = Population)) + 
  geom_density(alpha = 1/2) + theme_bw() +
  labs(y = "Density"))

# Combine i) Relatedness PCs, and ii) Relatedness PCs corrected for population structure.
cowplot::plot_grid(kinshippca + theme(legend.position = "top",
                                      legend.title = element_blank()) + xlab(NULL), 
                   pcair_vis  + theme(legend.position = "none"), 
                   ncol = 1, rel_heights = c(1.1, 1), align = "VH",
                   labels = c("A", "B"))

ggsave("../plots/relatedness_pcs.tiff", dpi = 300, height = 10, width = 8)
write.csv(pcair_dat[,c(1:2)], "../data/pcair1.csv", row.names = F)

write.csv(as.matrix(pcrelMat), "../data/pcrel_mat.csv", row.names = T)


# GWAS I: Global ---------------------------------------------------------------

gdsfmt::showfile.gds(closeall = T)

# Read in and reformat individual IDs, locations and phenotypes.
info <- read.csv("../data/chRADseq_samples.csv") %>%
  select(c("fishID", "pop", "age", "year")) %>% 
  mutate(pop = as.factor(sub("\\_.*", "", 
               tools::toTitleCase(tolower(pop)))),
         jack = case_when(age == 2 ~ 1,
                          age >  2 ~ 0)) %>% 
  rename("Sample" = fishID) 

# Read in PC1 values representing relatedness after accounting for population structure.
pcair <- read.csv("../data/pcair1.csv")

# Join data together and format for GWAS. 
(ind_dat <- ScanAnnotationDataFrame(merge(info, pcair) %>% rename("scanID" = Sample)))

# Write a function for conducting a global, multi-population GWAS that allows 
# easy changing of response variables and distributions.
global_gwas <- function(phenotype, distribution) {
  
  # Make GDS file readable/accessible.
  gdsfmt::showfile.gds(closeall = T)
  
  # Specify null model.
  ageNull <- fitNullModel(ind_dat, 
                          outcome = phenotype, 
                          covars = "PC1",
                          family = distribution, 
                          cov.mat = as.matrix(pcrelMat))
  
  # Specify genotypic data.
  geno <- GenotypeData(GdsGenotypeReader(file = "../data/snps_maf001_sub.gds"))
  
  # snpBlock = 100 from McKinney 2021: 10.1111/mec.15712
  genoIterator <- GenotypeBlockIterator(geno, snpBlock = 100)
  
  # Run the actual GWAS model in serial.
  age_assoc <- assocTestSingle(genoIterator, null.model = ageNull, test = "Score", 
                               BPPARAM = BiocParallel::SerialParam())
  age_assoc$chr <- as.vector(radgl@chromosome) # Return chromosome information that is lost in the GDS conversion.
  age_assoc$fdrp <- p.adjust(age_assoc$Score.pval, method = "fdr") # Compute corrected p-values.
  age_assoc$phenotype <- phenotype # Make a uniform column specifying target phenotype.
  return(age_assoc) # And save the corresponding dataframe to the environment.
  
}

jack_gwas_global <- global_gwas(phenotype = "jack", distribution = "binomial")
age_gwas_global  <- global_gwas(phenotype = "age", distribution = "gaussian")

# Plotting function for GWAS Manhattan plots.
singlepop_manhattan <- function(df) {
  
  sigval <- 0.05/nrow(df)
  
  ggplot(data = df, aes(x = variant.id, y = -log(Score.pval))) +
  geom_point(aes(colour = chr)) + theme_bw() +
  scale_color_manual(values = c(rep(c("gray20", "gray70"), 17))) +
  geom_point(data = age_assoc[age_assoc$outlier == "Outlier",],
             aes(x = variant.id, y = -log(Score.pval)),
             colour = "red") +
  theme(legend.title = element_blank(),
        legend.position = "none") +
  labs(x = NULL) +
  geom_hline(yintercept = -log(sigval)) }

(jack_manh <- singlepop_manhattan(df = age_gwas_global))
ggsave("../plots/global_manhattan_jacksbinom.tiff", dpi = 300, width  = 10, height = 5)


(age_manh <- singlepop_manhattan(df = jack_gwas_global))
ggsave("../plots/global_manhattan_agegauss.tiff", dpi = 300, width  = 10, height = 5)


# GWAS II: Local ----- ---------------------------------------------------------

# Read in GRM if starting fresh from this point.
pcrelMat <- read.csv("../data/pcrel_mat.csv", row.names = 1) %>% 
  `colnames<-`(., gsub("\\.", "-" , sub("^X", "", colnames(.))))

# Read in sample information and rearrange slightly.
ind_info <- read.csv("../data/chRADseq_samples.csv") %>% 
  select(c("fishID", "age")) %>% 
  mutate(jack = case_when(age == 2 ~ 1,
                          age >  2 ~ 0)) %>% 
  filter(fishID %in% subsamples)

gwaslist <- list(); gdsfmt::showfile.gds(closeall = T)

local_gwas <- function(phenotype, distribution) {
  
  gwaslist <- list(); gdsfmt::showfile.gds(closeall = T)
  
  for (pop in c(unique(info$pop))) {
  
  # Isolate individuals for each population.
  indvs <- as.vector(info[info$pop == pop, "Sample"])
  
  # Subset GRM for each population.
  local.relMat <- as.matrix(pcrelMat[rownames(pcrelMat) %in% indvs,
                                     colnames(pcrelMat) %in% indvs])
  
  # Technically only need to run this part once to obtain the population-level gds files. 
  # Global VCF was split into population-level VCFs and all monomorphic loci were removed.
  snpgdsVCF2GDS(vcf.fn = paste0("../data/pop_vcfs/", pop, ".vcf"),
                out.fn = paste0("../data/pop_vcfs/", pop, ".gds"))
  
  # For some reason the above conversion loses chromosomal information.
  # Not efficient, but I just pull out those data from the original vcf files and add them later.
  popchrs <- as.vector(read.vcfR(paste0("../data/pop_vcfs/", pop, ".vcf"))@fix[,"CHROM"])
  
  # Read in correctly formatted genotypic data for each population.
  popgeno <- GenotypeData(GdsGenotypeReader(paste0("../data/pop_vcfs/", pop, ".gds")))
  
  # Create an annotated dataframe with all relevant individual information.
  ind_data <- ScanAnnotationDataFrame(
    data = ind_info[ind_info$fishID %in% indvs, ] %>% 
      rename("scanID" = fishID))
  
  # Specify null model.
  nullmodel <- fitNullModel(ind_data,
                            outcome = phenotype,
                            family = distribution,
                            cov.mat = local.relMat)
  
  # Specify genotype block iterator.
  gi <- GenotypeBlockIterator(popgeno, snpBlock = 100)
  
  # Run the GWAS model. Output is a dataframe.    
  agemod <- assocTestSingle(gi, test = "Score",
                            null.model = nullmodel,
                            BPPARAM = BiocParallel::SerialParam())
  agemod$chr <- popchrs # Re-add chromosome info.
  
  # Determine adjusted significance threshold and identify SNPs that
  # are above/below that value as outliers or not. Add those values to 
  # the original dataframe for plotting purposes later on.
  agemod$sigval  <- popFDR <- 0.05/nrow(agemod)
  agemod$outlier <- case_when(agemod$Score.pval <  popFDR ~ "Outlier",
                              agemod$Score.pval >= popFDR ~ "Non-outlier")
  
  # Assign GWAS outputs to named list.
  gwaslist[[length(gwaslist)+1]] <- agemod
  
}

# Name list elements and convert to long-form dataframe.
names(gwaslist) <- c(unique(info$pop))
popGWAS <- map_df(gwaslist, ~as.data.frame(.x), .id = "pop") 
return(popGWAS)
}

jack_gwas_local <- local_gwas(phenotype = "jack", distribution = "binomial")
age_gwas_local  <- local_gwas(phenotype = "age",  distribution = "gaussian")

multi_manhattan <- function(df) {
  sigval <- 0.05/nrow(df)
  ggplot(data = df, 
         aes(x = variant.id, 
             y = -log(Score.pval))) +
    geom_point(aes(colour = chr)) + theme_bw() +
    scale_color_manual(values = c(rep(c("gray20", "gray70"), 17))) +
    theme(legend.title = element_blank(),
          legend.position = "none") +
    labs(x = NULL) +
    facet_wrap(~pop, ncol = 1) +
    geom_hline(yintercept = -log(sigval)) 
}

multi_manhattan(jack_gwas_local)
ggsave("../plots/jack_gwas_local_manhattan.tiff", dpi = 300, width = 12, height = 8)

multi_manhattan(age_gwas_local)
ggsave("../plots/age_gwas_local_manhattan.tiff", dpi = 300, width = 12, height = 8)
