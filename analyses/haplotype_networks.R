setwd("~/chinook_male_age/analyses")

library(tidyverse); library(cowplot); library(igraph); library(vcfR)
library(gplots); library(broom)

source("../scripts/phasedBase_012.R") # From McKinney et al., 2021: 10.1111/mec.15712

# To-do ------------------------------------------------------------------------

# Number/label chromosomes consistently so they are not coloured/labelled as a series.

# LD Plots ---------------------------------------------------------------------

# For each ld file, read it in and order it appropriately, then rename the list elements.
chromLD <- lapply(X = paste0("../data/chr_LDs/maf005/", 
                  list.files(path = "../data/chr_LDs/maf005", pattern = "*.ld")),
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
    theme(legend.position = c(.85, .35)))

# Title each LD plot with chromosome name.
for (i in 1:34) { LDplots[[i]] <- LDplots[[i]] + 
    ggtitle(paste(names(chromLD)[i])) }

# Merge plots together in a single image and save.
(genomeLD <- cowplot::plot_grid(plotlist = LDplots, ncol = 6))
ggsave("../plots/ots_full_maf005.jpg",width = 6000, height = 6000, units = "px")


# Make a plot of just Ots17, 18 and 30.
cowplot::plot_grid(plotlist = LDplots[c("Ots17", "Ots18", "Ots30")], ncol = 2,
                   labels = c("Ots17", "Ots18", "Ots30"), hjust = 0.5)
ggsave("../plots/ots17_18_30LD.jpg", width = 12, height = 12)

# Get some summary stats of linkage per chromosome. 
chStats <- lapply(chromLD, FUN = function(x) x %>% 
                  summarise(med = median(R2), ymin = min(R2), ymax = max(R2),
                            lowQ = as.numeric(quantile(R2, probs = 1/4)),
                            hiQ  = as.numeric(quantile(R2, probs = 3/4)))) %>% 
  bind_rows() %>% mutate(chrom = c(paste0("Ots", sprintf("%02d", seq(1, 34, 1)))))

ggplot(data = chStats, aes(x = as.factor(chrom))) +
  geom_segment(aes(y = 0, yend = 0.012), 
               arrow = arrow(length = unit(0.1, "cm")), 
               colour = "gray80") +
  geom_boxplot(aes(lower = lowQ, upper = hiQ, 
    middle = med, ymin = lowQ, ymax = hiQ),
    stat = "identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(limits = c(0, 0.012), 
                     breaks = seq(0, 0.012, 0.005)) +
  labs(y = bquote(~R^2), x = NULL) 

ggsave("../plots/chromLDpatterns.tiff", dpi = 300, width = 8, height = 8)
  

# Network analysis I: Global, Ots17 --------------------------------------------

# # Isolate SNPs in high LD on Ots17.
# Ots17HL <- chromLD[["Ots17"]] %>% 
#   select(c("SNP_A", "SNP_B", "R2")) %>% 
#   filter(R2 > 0.3)
# 
# # Read in SNP map from plink (anchor SNP IDs to places in the genome).
# snpmap <- read.delim("../data/global_vcf/global_maf005.map", header = FALSE)[,c(1:2,4)] %>% 
#   `colnames<-`(., c("chr", "SNP", "position"))
# 
# # Get frequencies for each SNP in the high LD dataset.
# nodes <- data.frame(table(c(Ots17HL$SNP_A, Ots17HL$SNP_B))) %>% 
#   `colnames<-`(., c("SNP", "Counts"))
# 
# net17 <- graph_from_data_frame(Ots17HL, nodes, directed = T)
# plot(net17,vertex.label=NA,layout=layout.fruchterman.reingold)
# 
# sym17 <- as.undirected(net17, mode = "collapse",
#                        edge.attr.comb = list(weight = "sum", "ignore"))
# 
# ceb17 <- cluster_edge_betweenness(sym17); plot_dendrogram(ceb17, "hclust")
# 
# mem17 <- as.matrix(membership(ceb17)) %>% 
#   data.frame() %>% rownames_to_column() %>% 
#   `colnames<-`(., c("SNP", "Group"))
# 
# hist(mem17$Group, main = NULL, xlab = "Group", breaks = nrow(mem17))
# 
# (memGrp <- mem17 %>% group_by(Group) %>% 
#   tally() %>% filter(n > 4) %>% 
#   arrange(desc(n)))
# 
# (topGrpSnps <- mem17 %>%
#     filter(Group %in% unlist(memGrp[, 1])) %>%
#     merge(., snpmap, by = "SNP"))
# 
# haploSnps <- mem17[mem17$Group %in% unlist(memGrp$Group),] %>% 
#   merge(., snpmap, by = "SNP") %>% select(c(3,4)) %>% 
#   mutate(position2 = position, rid = paste0("id_", rownames(.)))
# 
# # Write this file and use to filter vcf w/ bcftools.
# write.table(haploSnps[,c(1:2)], quote = F, sep = "\t",
#             "../data/ots17_haploSnpsBCF.txt",
#             row.names = F, col.names = F)
# 
# # SNPs in high LD (R2 > .3) used in haplotype analysis.
# # Filtered to SNPs identified above.
# vcf17 <- read.vcfR("../data/Ots17haplotypes.vcf")
# 
# q <- rownames(extract.gt(vcf17)) %in% topGrpSnps$SNP
# 
# # Read in phased haplotype information.
# hap <- read.delim("../data/ots17maf005_hapguess_switch.out", 
#                   header = FALSE, stringsAsFactors = FALSE,
#                   skip = 21) %>% 
#   filter(V1 != "END GENOTYPES") 
# 
# # Organize haplotype data output from fastPHASE.
# hapTable <- data.frame(
#   id  = gsub("^# ID ", "", hap[seq(1, nrow(hap), 3), 1],),
#   fwd = hap[seq(2, nrow(hap), 3), 1],
#   rev = hap[seq(3, nrow(hap), 3), 1] ) %>% 
#   pivot_longer(cols = c("fwd", "rev"), 
#                values_to = "haplotype") %>% 
#   select(c(1,3)) %>% 
#   separate_wider_delim(haplotype, delim = " ", names_sep = ".") %>%
#   .[,c(TRUE, q)] %>% 
#   mutate(id = make.unique(id, sep = "_")) %>% 
#   column_to_rownames(var = "id") %>% 
#   `colnames<-`(., c(topGrpSnps$position))
# 
# ots17hap <- apply(hapTable, 2 , phasedBase_012)
# rownames(ots17hap) <- rownames(hapTable)
# dim(ots17hap) == dim(hapTable)
# 
# png(width = 2500, height = 1500, units = "px", "../plots/global_heatmap2_maf001.png")
# (phased_heatmap <- heatmap.2(ots17hap, trace = "none",Colv = FALSE,
#                             key = FALSE, labRow = FALSE, cexCol = 4/5,
#                             hclustfun = function(x) hclust(x, method = "ward.D")))
# dev.off()
# 
# plot(phased_heatmap$rowDendrogram, nodePar=list(lab.cex = .2))
# abline(h = 180, col = 'red2', lty = 'dashed')
# abline(v = (nrow(ots17hap)/2)+0.5, col = 'red2')
# phasedDendro <- as.hclust(phased_heatmap$rowDendrogram)
# dendroGroups <- cutree(phasedDendro, h = 180)
# 
# hGrpDF <- data.frame(hapGrp = dendroGroups) %>% 
#   rownames_to_column(var = "id")
# 
# dat <- read.csv("../data/chRADseq_samples.csv")[,c(1,5,6,7)]
# 
# hGrpDF <- data.frame(Hap1 = dendroGroups[seq(1, length(dendroGroups), 2)],
#                      Hap2 = dendroGroups[seq(2, length(dendroGroups), 2)]) %>% 
#   rownames_to_column(var = "fishID") %>% 
#   merge(., dat, by = "fishID") %>% 
#   pivot_longer(cols = c("Hap1", "Hap2"), values_to = "haplotype") %>% 
#   mutate(haplotype = as.factor(haplotype),
#          pop = as.factor(pop))
# 
# write.csv(hGrpDF, "../data/haplotype_groups.csv", row.names = F)
# 
# # Haplotype associations -------------------------------------------------------
# 
# hapCounts <- hGrpDF %>% 
#   filter(age != 5) %>% 
#   group_by(pop, haplotype, age) %>% 
#   tally() %>% 
#   ungroup() %>% 
#   group_by(age) %>% 
#   mutate(propHap = n/sum(n))
# 
# (hapC <- ggplot(data = hapCounts, 
#                 aes(x = age, y = n, fill = haplotype)) +
#   geom_bar(position = "dodge", stat = "identity",
#            colour = "gray50") +
#   theme_bw() + labs(x = "Age", y = "Samples") +
#   theme(legend.position = "top") + labs(x = NULL) +
#   guides(fill = guide_legend(nrow = 1)))
# 
# (hapP <- ggplot(data = hapCounts[hapCounts$age < 5,], 
#                 aes(x = age, y = propHap, fill = haplotype)) +
#   geom_bar(position = "dodge", stat = "identity",
#            colour = "gray50") +
#   theme_bw() + labs(x = "Age", y = "Proportion of all samples") +
#   theme(legend.position = "none"))
# 
# cowplot::plot_grid(hapC, hapP, ncol = 1)
# ggsave("../plots/haplotype_ages.tiff", dpi = 300, width = 8, height = 8)
# 
# chisq.test(hGrpDF$haplotype, hGrpDF$age)
# chisq.test(table(hGrpDF$haplotype, hGrpDF$pop))
# 
# (p <- lm(data = hGrpDF, age ~ haplotype*pop)); summary(p)
# anova(p)
# 
# 
# popHap<- hGrpDF %>% 
#   group_by(pop, haplotype, age) %>% tally() %>% 
#   ungroup() %>% 
#   group_by(pop) %>% 
#   mutate(propHap = n/max(n))
# 
# 
# (hapPop <- ggplot(data = popHap[popHap$age < 5,], 
#                 aes(x = pop, y = propHap, fill = haplotype)) +
#     geom_bar(position = "dodge", stat = "identity",
#              colour = "gray50") +
#     theme_bw() + labs(x = NULL, y = "Proportion of all samples") +
#     theme(legend.position = "none"))


# Network analysis II: Local, Ots17 --------------------------------------------

# Read in LD information for each population and format into a list.
pop17LD <- lapply(X = paste0("../data/pop_haplotypes/maf001/",
                             list.files(path = "../data/pop_haplotypes/maf001/",
                                        pattern = "*.ld")),
                  FUN = function(x) read.delim(x, sep = "")) %>% 
                        lapply(., \(x) x[order(x$R2), ]) %>% 
                  `names<-`(., c("Chilliwack", "Puntledge", "Qualicum"))
                  # Files are alphabetical by default - specify that here.

# Read in sample information for later.
dat <- read.csv("../data/chRADseq_samples.csv")[,c(1,5,6,7)]

# Read in SNP map from plink (anchor SNP IDs to places in the genome).
snpmap <- read.delim("../data/global_vcf/global_maf001.map", header = FALSE)[,c(1:2,4)] %>% 
  `colnames<-`(., c("chr", "SNP", "position"))

# Write function to perform haplotype analyses for each population without
# having to copy script repeatedly. 
hapFunc <- function(population, memLimit, dendSplit, maf) {
  
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
    merge(., snpmap, by = "SNP")
  
  # Read in original VCF used for the haplotype/phasing analyses.
  # Extract positional SNP names by position in the original VCF input.
  PopVcf   <- read.vcfR(paste0("../data/pop_vcfs/HLD_maf", maf, "/",
                       population,"_maf", maf, "_HLD.vcf"), verbose = FALSE)
  keepSNPs <- rownames(extract.gt(PopVcf)) %in% GrpSnps$SNP

  # Read in phased haplotype information.
  hap      <- read.delim(paste0("../data/pop_haplotypes/maf", maf, "/", 
                                population, "_hapguess_switch.out"), 
                       header = FALSE, stringsAsFactors = FALSE,
                       skip = 21) %>% filter(V1 != "END GENOTYPES") 
  
  # Orient haplotype information by individual and haplotype 1 or 2.
  hapTab   <- data.frame( id  = gsub("^# ID ", "", hap[seq(1, nrow(hap), 3), 1],),
                        fwd = hap[seq(2, nrow(hap), 3), 1],
                        rev = hap[seq(3, nrow(hap), 3), 1] ) %>% 
    pivot_longer(cols = c("fwd", "rev"), values_to = "haplotype") %>% 
    select(c(1,3))  %>% 
    separate_wider_delim(haplotype, delim = " ", names_sep = "_") %>%
    .[,c(TRUE, keepSNPs)] %>% 
    mutate(id = make.unique(id, sep = "_")) %>% 
    column_to_rownames(var = "id") %>% 
    `colnames<-`(., c(GrpSnps[GrpSnps$SNP %in% rownames(extract.gt(PopVcf)), "position"]))
  
  ots17hap <- apply(hapTab, 2, phasedBase_012)
  rownames(ots17hap) <- rownames(hapTab)

  
  # png(width = 2500, height = 1500, units = "px", 
  #     filename = paste0("../plots/", population, "_heatmap2_GP_maf", maf, ".png"))
  PHM <- heatmap.2(ots17hap, trace = "none", cexCol = 1, Colv = FALSE, cexRow = 2/3,
                    key = FALSE, margins = c(8, 5),
                    hclustfun = function(x) hclust(x, method = "ward.D"))
  # print(paste0("Heatmap printed to ../plots/", population, "_heatmap2_GP2.png")); dev.off()
  
  # Dendrogram - individual-based with line designating haplotype clusters.
  png(width = 3000, height = 1000, units = "px", # Plot and save dendrogram with cutoff value.
      filename = paste0("../plots/", population, "_maf", maf, "_rowDendro.png"))
  plot(PHM$rowDendrogram,  nodePar = list(lab.cex = 4/5, pch = c(NA, NA))) 
  abline(h=as.numeric(dendSplit), col = 'red2', lty = 'dashed'); 
  abline(v=(ncol(PopVcf@gt)), col = 'red2', lty = 'dashed'); dev.off()
  print(paste0("Dendrogram printed to ../plots/", population, "_rowDendro.png"))
  
  # Isolate haplotype groups per individual.  
  hapGrp <- data.frame(dGrp =  cutree(as.hclust(PHM$rowDendrogram), 
            h = dendSplit)) %>% rownames_to_column(var = "id")
  
  # Assign colours to each age in the order of the row-oriented dendrogram.
  dend_lab_cols <- data.frame(fishID = rownames(ots17hap)[PHM$rowInd]) %>% 
    mutate(fishID = gsub("_1", "", fishID)) %>% 
    merge(., dat[,c("fishID", "age")]) %>% 
    mutate(colour = case_when(
      age == 2 ~ "blue2",  age == 3 ~ "forestgreen",
      age == 4 ~ "orange", age == 5 ~ "darkred"))
  
  png(width = 2500, height = 1500, units = "px", 
      filename = paste0("../plots/", population, "_heatmap2_GP_maf", maf, ".png"))
  PHM2 <- heatmap.2(ots17hap, trace = "none", cexCol = 3/2, Colv = FALSE, colRow = NULL,
                   RowSideColors = c(dend_lab_cols$colour), key = FALSE, margins = c(8, 1/5),
                   hclustfun = function(x) hclust(x, method = "ward.D"))
  print(paste0("Heatmap printed to ../plots/", population, "_heatmap2_GP2.png")); dev.off()
  
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
    heatmap     = PHM2,
    Grp         = hapGrp
  )

  return(outputList)
  
  }


# High-level haplogroups and age-associations ----------------------------------

# Runs the entire haplotype "pipeline" for each population. Maf needs to be a character string.
ChHapsHL <- hapFunc(population = "Chilliwack", memLimit = 8, dendSplit = 40, maf = "001")
PuHapsHL <- hapFunc(population = "Puntledge",  memLimit = 4, dendSplit = 50, maf = "001")
QuHapsHL <- hapFunc(population = "Qualicum",   memLimit = 1, dendSplit = 60, maf = "001")


# Build a dataframe used to plot haplotype counts and proportions by age.
all_hapsHL <- bind_rows(list(ChHapsHL$haploSumm %>% mutate(pop = "Chilliwack"),
                             PuHapsHL$haploSumm %>% mutate(pop = "Puntledge"),
                             QuHapsHL$haploSumm %>% mutate(pop = "Qualicum")))

(age_hapHL <- ggplot(data = all_hapsHL, aes(x = age, y = n, fill = haplotype)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme_bw() + labs(x = "Age", y = "Samples") +
    facet_wrap(. ~ pop) +
    theme(legend.position = "top"))

(age_perHL <- ggplot(data = all_hapsHL[all_hapsHL$age < 5,],
                   aes(x = age, y = pHap, fill = haplotype)) +
    geom_bar(position = "dodge", stat = "identity", color = 'gray75') +
    theme_bw() + labs(x = "Age", y = "Proportion of Samples") +
    facet_wrap(. ~ pop) +
    theme(legend.position = "none"))

cowplot::plot_grid(age_hapHL, age_perHL, ncol = 1, rel_heights = c(1.15, 1))
ggsave("../plots/withinpop_agehaplo.tiff", dpi = 300, width = 10, height = 8)

# Test if age is predicted by haplotpe. 
(hapTest <- bind_rows(list(
  ChHapsHL$hapDatFull %>% mutate(pop = "Chilliwack"),
  PuHapsHL$hapDatFull %>% mutate(pop = "Puntledge"),
  QuHapsHL$hapDatFull %>% mutate(pop = "Qualicum"))) %>% 
    nest(data = -pop) %>% 
    mutate(lmTest = map(data, ~ lm(age ~ haplotype, data = .x)),
           lmSumm = map(lmTest, glance)) %>% 
    unnest(lmSumm))



# Pop-check --------------------------------------------------------------------

phs <- list()
for (i in c(20, 50)) {
  hs <- hapFunc(population = "Puntledge",  memLimit = 8, dendSplit = i, maf = "001")

  phs[[length(phs)+1]] <- hs$Grp %>%
    `colnames<-`(., c("id", paste0("split_", i)))
  }
p <- reduce(phs, full_join, by = "id") %>%
  mutate(haplotype = case_when(
     split_20 %in% c(8, 10, 5, 6, 4) ~ 1,
    !split_20 %in% c(8, 10, 5, 6, 4) ~ split_50
  )) %>% mutate(sc = case_when(haplotype == 1 ~ "X",
                               haplotype != 1 ~ "Y")) %>%
  mutate(id = gsub("_1", "", id),
         hap = rep(c("hap1", "hap2"), 236)) %>%
  select(c("id", "sc", "hap")) %>%
  pivot_wider(values_from = "sc", names_from = "hap")
nrow(p[p$hap1 == "Y" & p$hap2 == "Y",])
nrow(p[p$hap1 == "X" & p$hap2 == "X",])
nrow(p[p$hap1 != p$hap2,])




pop_check <- \(hap) {
  
  pX <- hap %>% 
    group_by(haplotype) %>% 
    tally() %>% 
    arrange(desc(n)) %>% 
    .[1,1]

  print(paste0("Most common haplotype is ", pX), quote = F)
    
  j <- hap %>% 
    mutate(sc = case_when( haplotype %in% pX ~ "X",
                          !haplotype %in% pX ~ "Y")) %>% 
    select(-c(haplotype)) %>% 
    pivot_wider(values_from = sc, names_from = name)
  
  print(paste0("Number of XY samples = ", nrow(j[j$Hap1 != j$Hap2,])), quote = F)
  print(paste0("Number of YY samples = ", nrow(j[j$Hap1 == "Y" & j$Hap2 == "Y", ])), quote = F)
  print(paste0("Number of XX samples = ", nrow(j[j$Hap1 == "X" & j$Hap2 == "X", ])), quote = F)
  
}

pop_check(QuHapsHL$hapDatFull)
pop_check(PuHapsHL$hapDatFull)
pop_check(ChHapsHL$hapDatFull)


qhs <- list()
for (i in c(20, 60)) {
  hs <- hapFunc(population = "Qualicum",  memLimit = 8, dendSplit = i, maf = "001")
  
  qhs[[length(qhs)+1]] <- hs$Grp %>%
    `colnames<-`(., c("id", paste0("split_", i)))
}
q <- reduce(qhs, full_join, by = "id") %>% 
  mutate(haplotype = case_when(
    split_20 %in% c(4,5,2,1,9,7) ~ 2,
    !split_20 %in% c(4,5,2,1,9,7) ~ split_60
  )) %>% mutate(sc = case_when(haplotype == 2 ~ "X",
                               haplotype != 2 ~ "Y")) %>%
  mutate(id = gsub("_1", "", id),
         hap = rep(c("hap1", "hap2"), 233)) %>%
  select(c("id", "sc", "hap")) %>%
  pivot_wider(values_from = "sc", names_from = "hap")
nrow(q[q$hap1 == "Y" & q$hap2 == "Y",])
nrow(q[q$hap1 == "X" & q$hap2 == "X",])
nrow(q[q$hap1 != q$hap2,])


# Fine-scale haplogroups and age-associations ----------------------------------

# First, adjust puntledge.
PuHapsFS <- hapFunc(population = "Puntledge", memLimit = 8, dendSplit = 25, maf = "001")
length(unique(PuHapsFS$Grp$dGrp))

# Second, adjust Qualicum.
# Requires manual clustering.
QuFS <- list()
for (i in c(20, 25)) {
  
  QuHaps <- hapFunc(population = "Qualicum", memLimit = 8, dendSplit = i,  maf = "001")
  
  QuFS[[length(QuFS)+1]] <-  QuHaps$Grp %>%
    `colnames<-`(., c("id", paste0("split_", i)))
}

QuHapsFS <- reduce(QuFS, full_join, by = "id") %>% 
  mutate(haplotype1 = as.numeric(case_when(
     split_20 %in% c("2", "4") ~ split_25,
    !split_20 %in% c("2", "4") ~ split_20
  )), # Rename some of the factors below for consistency.
  haplotype = as.factor(case_when(
    haplotype1 < 4 ~ haplotype1,
    haplotype1 > 4 ~ haplotype1 - 1
  ))) %>% merge(., QuHapsHL$hapDatFull[,c("fishID", "age")],
               by.x = "id", by.y = "fishID") %>% 
  select(-c("haplotype1"))

length(unique(QuHapsFS$haplotype))

# Third, adjust Chilliwack.
ChFS <- list()
for (i in c(20, 50)) {
  
  ChHaps <- hapFunc(population = "Chilliwack", memLimit = 8, dendSplit = i,  maf = "001")
  
  ChFS[[length(ChFS)+1]] <- ChHaps$Grp %>%
    `colnames<-`(., c("id", paste0("split_", i)))
}

ChHapsFS <- reduce(ChFS, full_join, by = "id") %>% 
  mutate(haplotype = as.factor(case_when(
     split_20 %in% c("2", "1") ~ "6",
     split_20 %in% c("7", "9") ~ "7",
    !split_20 %in% c("2", "1", "7", "9") ~ as.factor(split_50)
  ))) %>% merge(., ChHapsHL$hapDatFull[,c("fishID", "age")],
              by.x = "id", by.y = "fishID") 

hapsum <- \(haps) { haps %>% group_by(haplotype, age) %>% 
    tally() %>% ungroup() %>% group_by(age) %>% 
    mutate(pHap = n/sum(n)) }

length(unique(ChHapsFS$haplotype))

# Join together into a single dataframe.
all_hapsFS <- bind_rows(list(PuHapsFS$haploSumm %>% mutate(pop = "Puntledge"),
                             hapsum(ChHapsFS)   %>% mutate(pop = "Chilliwack"),
                             hapsum(QuHapsFS)   %>% mutate(pop = "Qualicum")))

(age_hapFS <- ggplot(data = all_hapsFS[all_hapsFS$age < 5,],
                     aes(x = age, y = n, fill = haplotype)) +
    geom_bar(position = "dodge", stat = "identity") +
    theme_bw() + labs(x = "Age", y = "Samples") +
    facet_wrap(. ~ pop) +
    theme(legend.position = "top"))

(age_perFS <- ggplot(data = all_hapsFS[all_hapsFS$age < 5,],
                     aes(x = age, y = pHap, fill = haplotype)) +
    geom_bar(position = "dodge", stat = "identity", colour = "gray75") +
    theme_bw() + labs(x = "Age", y = "Proportion of Samples") +
    facet_wrap(. ~ pop) +
    theme(legend.position = "none"))

(hapTestFS <- bind_rows(PuHaps$hapDatFull %>% mutate(pop = "Puntledge"),
                        ChHapsFS[,c(1,4:5)] %>% mutate(pop = "Chilliwack"),
                        QuHapsFS[,c(1,4:5)] %>% mutate(pop = "Qualicum")) %>% 
    nest(data = -pop) %>% 
    mutate(lmTest = map(data, ~ lm(age ~ haplotype, data = .x)),
           lmSumm = map(lmTest, glance)) %>% 
    unnest(lmSumm))


################################################################################
################################################################################
################################################################################
