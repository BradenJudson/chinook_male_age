setwd("~/chinook_male_age/analyses")

library(tidyverse); library(cowplot); library(igraph); library(vcfR)


# To-do ------------------------------------------------------------------------

# Consider checking LD patterns by population.
# Label cowplot panels.

# LD Plots ---------------------------------------------------------------------

# For each ld file, read it in and order it appropriately, then rename the list elements.
chromLD <- lapply(X = paste0("../data/", list.files(path = "../data/", pattern = "*.ld")),
           FUN = function(x) read.delim(x, sep = "")) %>% 
           lapply(., \(x) x[order(x$R2), ]) %>% 
           `names<-`(.,  c("Ots17", "Ots18", "Ots30"))

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
cowplot::plot_grid(plotlist = LDplots, ncol = 2)
ggsave("../plots/ots17_18_30LD.jpg", dpi = 300, width = 12, height = 12)


# Network analysis -------------------------------------------------------------


# Read in SNP map from plink (anchor SNP IDs to places in the genome).
snpmap <- read.delim("../data/Ots17.map", header = FALSE)[,c(1:2,4)] %>% 
  `colnames<-`(., c("chr", "SNP", "position"))

# Isolate SNPs in high LD on Ots17.
Ots17HL <- chromLD[["Ots17"]] %>% 
  select(c("SNP_A", "SNP_B", "R2")) %>% 
  filter(R2 > 0.3)

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


################################################################################
##### MCKINNEY chooses top 2 groups with n > 4. Why? CHECK THIS ################
################################################################################


(memGrp <- mem17 %>% group_by(Group) %>% 
  tally() %>% filter(n > 4) %>% 
  arrange(desc(n)))


################################################################################
##### WHY DO HAPLOTYPE GROUPS OVERLAP? #########################################
################################################################################


topGrpSnps <- mem17 %>%
  filter(Group %in% unlist(memGrp[1:4, 1])) %>%
  merge(., snpmap, by = "SNP") %>% 
  group_by(Group) %>%
  summarise(minPos   = min(position),
            maxPos   = max(position),
            numSnps  = n()) %>%
  mutate(snpDistance = maxPos - minPos,
         chromosome  = "NC_056445.1")


write.table(topGrpSnps[,c(6,2,3)],quote = F,
            "../data/ots17_haplotype_regs.txt",
            col.names = FALSE, row.names = F)


# Use bcftools --regions-file and above to extract vcf w/ subset of SNPs.

vcf17 <- read.vcfR()







