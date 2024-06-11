setwd("~/chinook_male_age/analyses")

library(tidyverse); library(cowplot); library(igraph); library(vcfR)
library(gplots)


# To-do ------------------------------------------------------------------------

# Consider checking LD patterns by population.
# Label cowplot panels.
# for heatmap2 plots: replace the labels with a number line for genomic position with landmarks that roughly correspond to the major blocks

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
    scale_color_gradient(low="gray99",high="black",
                         name = expression(R^{2})) +
    labs(x = "Position (Mbp)", y = "Position (Mbp)") +
    theme(legend.position = c(.95, .2)))

# Merge plots together in a single image and save.
(genomeLD <- cowplot::plot_grid(plotlist = LDplots, ncol = 6))
ggsave("../plots/ots_full.jpg",width = 5000, height = 5000, units = "px")


# Make a plot of just Ots17, 18 and 30.
cowplot::plot_grid(plotlist = LDplots[c("Ots17", "Ots18", "Ots30")], ncol = 2,
                   labels = c("Ots17", "Ots18", "Ots30"), hjust = 0.5)
ggsave("../plots/ots17_18_30LD.jpg", width = 12, height = 12)

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

#function to convert to 0/1 format with 0 being the most common allele
phasedBase_012<-function(bases){
  orgAlleles<-unique(bases)
  #create named vector to convert
  newAlleles<-c("0","1")
  names(newAlleles)<-orgAlleles
  #convert bases to numbers
  numBases<-str_replace_all(bases,newAlleles)
  #swap 0 and 1 if necessary so 0 is most common
  if(sum(numBases==0)>sum(numBases==1)){
    numBases<-chartr("01","10",numBases)
  }
  numBases<-as.numeric(numBases)
  return(numBases)
}

ots17hap <- apply(hapTable, 2 , phasedBase_012)
rownames(ots17hap) <- rownames(hapTable)
dim(ots17hap) == dim(hapTable)

(phased_heatmap <- heatmap.2(ots17hap, trace = "none",
                            key = FALSE, labRow = FALSE, labCol = FALSE,
                            hclustfun = function(x) hclust(x, method = "ward.D")))


plot(phased_heatmap$rowDendrogram, xaxt="n")
abline(h=180, col = 'red2', lty = 'dashed')
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







