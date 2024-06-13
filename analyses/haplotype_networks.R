setwd("~/chinook_male_age/analyses")

library(tidyverse); library(cowplot); library(igraph); library(vcfR)
library(gplots)


# To-do ------------------------------------------------------------------------

# Consider checking LD patterns by population.
# Label cowplot panels.
# for heatmap2 plots: replace the labels with a number line for genomic position with landmarks that roughly correspond to the major blocks
# add chromosome label to LD plots
# maybe do boxplot of LD patterns by chromosome?
# Number/label chromosomes consistently so they are not coloured/labelled as a series.
# filter vcf to those in haplotype groups - need to modify bash with output txt file
# check filtering of groups wit hn = 4, 6 and 8 SNPs
# make sure all individuals have at least one unclassified haplotype
# tidy figures and tables
# do statistical tests for age ~ hap


# LD Plots ---------------------------------------------------------------------

# For each ld file, read it in and order it appropriately, then rename the list elements.
chromLD <- lapply(X = paste0("../data/chr_LDs/", 
                  list.files(path = "../data/chr_LDs", pattern = "*.ld")),
           FUN = function(x) read.delim(x, sep = "")) %>% 
           lapply(., \(x) x[order(x$R2), ]) %>% 
           `names<-`(.,  c(paste0("Ots", sprintf("%02d", seq(1, 34, 1)))))

# Make a plot for each chromosome analysed.
LDplots <- lapply(chromLD, FUN = function(x)
  ggplot(data = x, aes(x = BP_A/1e6,
                       y = BP_B/1e6,
                       colour = R2)) +
    geom_point() + theme_classic() +
    scale_color_gradient(low  = "gray99", high="black",
                         name = expression(R^{2})) +
    labs(x = "Position (Mbp)", y = "Position (Mbp)") +
    theme(legend.position = c(.95, .2)))

# Title each LD plot with chromosome name.
for (i in 1:34) { LDplots[[i]] <- LDplots[[i]] + 
    ggtitle(paste(names(chromLD)[i])) }

# Merge plots together in a single image and save.
(genomeLD <- cowplot::plot_grid(plotlist = LDplots, ncol = 6))
ggsave("../plots/ots_full.jpg",width = 6000, height = 6000, units = "px")


# Make a plot of just Ots17, 18 and 30.
cowplot::plot_grid(plotlist = LDplots[c("Ots17", "Ots18", "Ots30")], ncol = 2,
                   labels = c("Ots17", "Ots18", "Ots30"), hjust = 0.5)
ggsave("../plots/ots17_18_30LD.jpg", width = 12, height = 12)


chStats <- lapply(chromLD, FUN = function(x) x %>% 
                  summarise(med = median(R2), 
                            ymin = min(R2),
                            ymax = max(R2),
                            lowQ = as.numeric(quantile(R2, probs = 1/4)),
                            hiQ  = as.numeric(quantile(R2, probs = 3/4)))) %>% 
  bind_rows() %>%
  mutate(chrom = c(paste0("Ots", sprintf("%02d", seq(1, 34, 1)))))



ggplot(data = chStats, aes(x = as.factor(chrom))) +
  geom_boxplot(aes(
    lower = lowQ, 
    upper = hiQ, 
    middle = med,
    ymin = lowQ,    #### ADD MIN/MAX But scale y-axis limits and indicate w/arrow?
    ymax = hiQ),
    stat = "identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(y = "R2", x = NULL) 
  


# Network analysis I: Global, Ots17 --------------------------------------------

# Isolate SNPs in high LD on Ots17.
Ots17HL <- chromLD[["Ots17"]] %>% 
  select(c("SNP_A", "SNP_B", "R2")) %>% 
  filter(R2 > 0.3)

# Read in SNP map from plink (anchor SNP IDs to places in the genome).
snpmap <- read.delim("../data/Ots17.map", header = FALSE)[,c(1:2,4)] %>% 
  `colnames<-`(., c("chr", "SNP", "position"))

# Get frequencies for each SNP in the high LD dataset.
nodes <- data.frame(table(c(Ots17HL$SNP_A, Ots17HL$SNP_B))) %>% 
  `colnames<-`(., c("SNP", "Counts"))

net17 <- graph_from_data_frame(Ots17HL, nodes, directed = T)
plot(net17,vertex.label=NA,layout=layout.fruchterman.reingold)

sym17 <- as.undirected(net17, mode = "collapse",
                       edge.attr.comb = list(weight = "sum", "ignore"))

ceb17 <- cluster_edge_betweenness(sym17); plot_dendrogram(ceb17, "hclust")

mem17 <- as.matrix(membership(ceb17)) %>% 
  data.frame() %>% rownames_to_column() %>% 
  `colnames<-`(., c("SNP", "Group"))

hist(mem17$Group, main = NULL, xlab = "Group", breaks = nrow(mem17))

(memGrp <- mem17 %>% group_by(Group) %>% 
  tally() %>% filter(n > 4) %>% 
  arrange(desc(n)))

(topGrpSnps <- mem17 %>%
  filter(Group %in% unlist(memGrp[, 1])) %>%
  merge(., snpmap, by = "SNP") %>% 
  group_by(Group) %>%
  summarise(minPos   = min(position),
            maxPos   = max(position),
            numSnps  = n()) %>%
  mutate(snpDistance = maxPos - minPos,
         chromosome  = "NC_056445.1"))


haploSnps <- mem17[mem17$Group %in% unlist(memGrp$Group),] %>% 
  merge(., snpmap, by = "SNP") %>% select(c(3,4)) %>% 
  mutate(position2 = position, rid = paste0("id_", rownames(.)))

write.table(haploSnps, quote = F, sep = "\t",
            "../data/ots17_haploSnpsBCF.txt",
            row.names = F, col.names = F)

# SNPs in high LD (R2 > .3) used in haplotype analysis.
vcf17 <- read.vcfR("../data/Ots17haplotypes.vcf")

# Read in phased haplotype information.
hap <- read.delim("../data/Ots17phased_hapguess_switch.out", 
                  header = FALSE, stringsAsFactors = FALSE,
                  skip = 21) %>% 
  filter(V1 != "END GENOTYPES") 

# Organize haplotype data output from fastPHASE.
hapTable <- data.frame(
  id  = gsub("^# ID ", "", hap[seq(1, nrow(hap), 3), 1],),
  fwd = hap[seq(2, nrow(hap), 3), 1],
  rev = hap[seq(3, nrow(hap), 3), 1] ) %>% 
  pivot_longer(cols = c("fwd", "rev"), 
               values_to = "haplotype") %>% 
  select(c(1,3)) %>% 
  separate_wider_delim(haplotype, delim = " ", names = c(vcf17@fix[,3])) %>% 
  mutate(id = make.unique(id, sep = "_")) %>% 
  column_to_rownames(var = "id")

# From McKinney et al., 2021.
source("phasedBase_012.R")

ots17hap <- apply(hapTable, 2 , phasedBase_012)
rownames(ots17hap) <- rownames(hapTable)
dim(ots17hap) == dim(hapTable)

(phased_heatmap <- heatmap.2(ots17hap, trace = "none",
                            key = FALSE, labRow = FALSE, labCol = FALSE,
                            hclustfun = function(x) hclust(x, method = "ward.D")))


plot(phased_heatmap$rowDendrogram, nodePar=list(lab.cex = .2))
abline(h=180, col = 'red2', lty = 'dashed')
abline(v= (ncol(vcf17@gt)/2)+0.5, col = 'red2')
phasedDendro <- as.hclust(phased_heatmap$rowDendrogram)
dendroGroups <- cutree(phasedDendro, h = 180)

hGrpDF <- data.frame(hapGrp = dendroGroups) %>% 
  rownames_to_column(var = "id")

dat <- read.csv("../data/chRADseq_samples.csv")[,c(1,5,6,7)]

hGrpDF <- data.frame(Hap1 = dendroGroups[seq(1, length(dendroGroups), 2)],
                     Hap2 = dendroGroups[seq(2, length(dendroGroups), 2)]) %>% 
  rownames_to_column(var = "fishID") %>% 
  merge(., dat, by = "fishID") %>% 
  pivot_longer(cols = c("Hap1", "Hap2"), values_to = "haplotype") %>% 
  mutate(haplotype = as.factor(haplotype),
         pop = as.factor(pop))

write.csv(hGrpDF, "../data/haplotype_groups.csv", row.names = F)

# Haplotype associations -------------------------------------------------------

hapCounts <- hGrpDF %>% 
  group_by(pop, haplotype, age) %>% tally() %>% 
  ungroup() %>% 
  group_by(age) %>% 
  mutate(propHap = n/sum(n))

(hapC <- ggplot(data = hT, aes(x = age, y = n, fill = haplotype)) +
  geom_bar(position = "dodge", stat = "identity",
           colour = "gray50") +
  theme_bw() + labs(x = "Age", y = "Samples") +
  theme(legend.position = "top") + labs(x = NULL) +
  guides(fill = guide_legend(nrow = 1)))

(hapP <- ggplot(data = hapCounts[hapCounts$age < 5,], 
                aes(x = age, y = propHap, fill = haplotype)) +
  geom_bar(position = "dodge", stat = "identity",
           colour = "gray50") +
  theme_bw() + labs(x = "Age", y = "Proportion of all samples") +
  theme(legend.position = "none"))

cowplot::plot_grid(hapC, hapP, ncol = 1)
ggsave("../plots/haplotype_ages.tiff", dpi = 300, width = 8, height = 8)

chisq.test(hGrpDF$haplotype, hGrpDF$age)
chisq.test(table(hGrpDF$haplotype, hGrpDF$pop))

(p <- lm(data = hGrpDF, age ~ haplotype*pop)); summary(p)
anova(p)


# cat Ots17.ld | cut -f2,5 | tr "\t" "\n" | sort | uniq | wc -l

popHap<- hGrpDF %>% 
  group_by(pop, haplotype, age) %>% tally() %>% 
  ungroup() %>% 
  group_by(pop) %>% 
  mutate(propHap = n/max(n))


(hapPop <- ggplot(data = popHap[popHap$age < 5,], 
                aes(x = pop, y = propHap, fill = haplotype)) +
    geom_bar(position = "dodge", stat = "identity",
             colour = "gray50") +
    theme_bw() + labs(x = "Age", y = "Proportion of all samples") +
    theme(legend.position = "none"))


# Network analysis II: Local, Ots17 --------------------------------------------

# Read in LD information for each population and format into a list.
pop17LD <- lapply(X = paste0("../data/pop_haplotypes/",
                             list.files(path = "../data/pop_haplotypes/",
                                        pattern = "*.ld")),
                  FUN = function(x) read.delim(x, sep = "")) %>% 
                        lapply(., \(x) x[order(x$R2), ]) %>% 
                  `names<-`(., c("Chilliwack", "Puntledge", "Qualicum"))
                  # Files are alphabetical by default - specify that here.

# Read in sample information for later.
dat <- read.csv("../data/chRADseq_samples.csv")[,c(1,5,6,7)]

# Read in SNP map from plink (anchor SNP IDs to places in the genome).
snpmap <- read.delim("../data/Ots17.map", header = FALSE)[,c(1:2,4)] %>% 
  `colnames<-`(., c("chr", "SNP", "position"))


# Write function to perform haplotype analyses for each population without
# having to copy script repeatedly. 
hapFunc <- function(population, memLimit, dendSplit) {
  
  # Isolate population of interest from list specified above.
  x <- pop17LD[[population]]
  
  # Isolate high LD SNPs.
  highLD <- x %>% select(c("SNP_A", "SNP_B", "R2")) %>% filter(R2 > 0.3) 
  
  # Occurrences of high LD SNPs comprise graph nodes.
  nodes  <- data.frame(table(c(highLD$SNP_A, highLD$SNP_B))) %>% 
    `colnames<-`(., c("SNP", "Counts"))

  # Visualize relationships between nodes.
  net <- graph_from_data_frame(highLD, nodes, directed = T)
  nodePlot <- plot(net, vertex.label = NA, layout = layout.fruchterman.reingold)
  
  # Collapse graph into clusters based on nodes/edge layout.
  sym <- as.undirected(net, mode = "collapse",
         edge.attr.comb = list(weight = "sum", "ignore"))
  ceb <- cluster_edge_betweenness(sym)
  
  # Determine SNP contributions to group memberships.
  membership <- data.frame(as.matrix(membership(ceb))) %>% 
    rownames_to_column() %>% 
    `colnames<-`(., c("SNP", "Group")) 
  
  # Number of SNPs contributing to clusters.
  memGrp <- membership %>% group_by(Group) %>% 
    tally() %>% filter(n > memLimit) %>% 
    arrange(desc(n))
  
  GrpSnps <- membership %>% 
    filter(Group %in% unlist(memGrp[,1])) %>% 
    merge(., snpmap, by = "SNP") %>% 
    group_by(Group) %>% 
    summarise(minPos   = min(position),
              maxPos   = max(position),
              numSnps  = n()) %>% 
    mutate(snpDistance = maxPos - minPos,
           chromosome  = "NC_056445.1")
  
  # Read in original VCF used for the haplotype/phasing analyses.
  PopVcf <- read.vcfR(paste0("../data/pop_vcfs/HLD/", 
                      population, "_HLD.vcf"), verbose = FALSE)

  # Read in phased haplotype information.
  hap    <- read.delim(paste0("../data/pop_haplotypes/", population,
                              "_hapguess_switch.out"), 
                       header = FALSE, stringsAsFactors = FALSE,
                       skip = 21) %>% filter(V1 != "END GENOTYPES") 
  
  # Orient haplotype information by individual and haplotype 1 or 2.
  hapTab <- data.frame( id  = gsub("^# ID ", "", hap[seq(1, nrow(hap), 3), 1],),
                        fwd = hap[seq(2, nrow(hap), 3), 1],
                        rev = hap[seq(3, nrow(hap), 3), 1] ) %>% 
    pivot_longer(cols = c("fwd", "rev"), values_to = "haplotype") %>% 
    select(c(1,3)) %>% 
    separate_wider_delim(haplotype, delim = " ", names = c(PopVcf@fix[,2])) %>% 
    mutate(id = make.unique(id, sep = "_")) %>% 
    column_to_rownames(var = "id")
  
  ots17hap <- apply(hapTab, 2, phasedBase_012)
  rownames(ots17hap) <- rownames(hapTab)
  
  # Below saves the heatmap to the plot directory but also saves the data
  # that comprise the heatmap to the resulting list of features for each population.
  png(width = 2500, height = 1500, units = "px", 
      filename = paste0("../plots/", population, "_heatmap2_GP.png"))
  PHM <- heatmap.2(ots17hap, trace = "none", cexCol = 1,
                   Colv = FALSE,
                   key = FALSE, labRow = FALSE, margins = c(8, 5),
                   hclustfun = function(x) hclust(x, method = "ward.D"))
  print(paste0("Heatmap printed to ../plots/", population, "_heatmap2_GP.png")); dev.off()
  
  
  # Dendrogram - individual-based with line designating haplotype clusters.
  png(width = 3000, height = 1000, units = "px", # Plot and save dendrogram with cutoff value.
      filename = paste0("../plots/", population, "_rowDendro.png"))
  plot(PHM$rowDendrogram,  nodePar=list(lab.cex = 4/5, pch = c(NA,NA))) 
  abline(h=as.numeric(dendSplit), col = 'red2', lty = 'dashed'); 
  abline(v=(ncol(PopVcf@gt)), col = 'red2', lty = 'dashed'); dev.off()
  print(paste0("Dendrogram printed to ../plots/", population, "_rowDendro.png"))
  
  # Isolate haplotype groups per individual.  
  hapGrp <- data.frame(dGrp =  cutree(as.hclust(PHM$rowDendrogram), 
            h = dendSplit)) %>% rownames_to_column(var = "id")
  
  # Make a dataframe of which haplotypes correspond to which individuals.
  hapDF <- data.frame(fishID = hapGrp[seq(1,nrow(hapGrp),2), "id"],
                      Hap1 = hapGrp[seq(1,nrow(hapGrp),2), "dGrp"],
                      Hap2 = hapGrp[seq(2,nrow(hapGrp),2), "dGrp"]) %>% 
    merge(., dat, by = "fishID") %>% 
    pivot_longer(cols = c("Hap1", "Hap2"), values_to = "haplotype") %>% 
    mutate(haplotype = as.factor(haplotype),
           pop = as.factor(pop)) # pop factor for merging later on.
  
  # Summarise haplotype distributions.
  hapSummary <- hapDF %>% 
    group_by(haplotype, age) %>% 
    tally() %>% ungroup() %>% 
    group_by(age) %>% 
    mutate(pHap = n/sum(n))
  
  # Re-assign everything relevant back to a list and populate.
  outputList <- list(
    grpHist     = hist(membership$Group, main = NULL, breaks = nrow(membership)),
    membership  = memGrp, 
    groupSnps   = GrpSnps,
    haplotypes  = hapTab,
    heatmapdata = ots17hap,
    hapDatFull  = hapDF,
    haploSumm   = hapSummary,
    q = membership
  )

  return(outputList)
  
  }


ChHaps <- hapFunc(population = "Chilliwack", memLimit = 4, dendSplit = 40)
PuHaps <- hapFunc(population = "Puntledge",  memLimit = 4, dendSplit = 80)
QuHaps <- hapFunc(population = "Qualicum",   memLimit = 4, dendSplit = 60)


all_haps <- bind_rows(list(ChHaps$haploSumm %>% mutate(pop = "Chilliwack"),
                    PuHaps$haploSumm %>% mutate(pop = "Puntledge"),
                    QuHaps$haploSumm %>% mutate(pop = "Qualicum")))

(age_hap <- ggplot(data = all_haps, aes(x = age, y = n, fill = haplotype)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_bw() + labs(x = "Age", y = "Samples") +
  facet_wrap(. ~ pop) +
  theme(legend.position = "top"))

(age_per <- ggplot(data = all_haps[all_haps$age < 5,],
                   aes(x = age, y = pHap, fill = haplotype)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme_bw() + labs(x = "Age", y = "Proportion of Samples") +
    facet_wrap(. ~ pop) +
    theme(legend.position = "none"))

cowplot::plot_grid(age_hap, age_per, ncol = 1, rel_heights = c(1.15, 1))
ggsave("../plots/withinpop_agehaplo.tiff", dpi = 300, width = 10, height = 8)


## Use broom or something here to do lm and/or anova/chisq by population.


hapTest <- bind_rows(list(
  ChHaps$hapDatFull %>% mutate(pop = "Chilliwack"),
  PuHaps$hapDatFull %>% mutate(pop = "Puntledge"),
  QuHaps$hapDatFull %>% mutate(pop = "Qualicum")
))













