# This script analyses the module genes from DIAMoND
# Created and modified by Alberto Zenere, 2022-04-11

# 11_Analysis_modules ####
# 11.1 Setup
# 11.2 Load files
# 11.3 Save 
# 11.4 Load seed genes
# 11.5 MS-enrichment
# 11.6 Enrichment of Activation
# 11.7 Overlap with P4 genes

# Setup ####
rm(list=ls()) # remove all entries in the global environment 

pathlist <- .libPaths();
newpath <-  "C:/Users/albze08/.conda/envs/P3/Lib/R/library"
.libPaths(newpath)

pack_R <- c("limma","stringr","dplyr","tidyverse","scatterplot3d","pracma","ggrepel","readxl",
            "lumi", "matrixStats")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}
require("lumi")
require("pracma")
require("ggpubr")
library("ggvenn")

# Set directory structure 
main_dir <-"C:/Users/albze08/Desktop/phd/P4/methylation"
FOLDER_RDS <- "RDS_files/"
FOLDER_FIGURES <- "figures/"
setwd("C:/Users/albze08/Desktop/phd/P4/methylation")

# My functions 
source('probe_to_gene.R')
source('fisher_test.R')

# 11.2 Load files ####
#MS-associated genes
disgenet <- read.csv("data/MS_DisGeNET.txt", sep = "\t")
disgenet_genes <- disgenet$geneSymbol %>% unique()

science_genes <- read.csv(file = "C:/Users/albze08/Desktop/phd/P4/RNAseq/data/complete_MS_annotation_p10e-6.csv")
science_genes <- unique(science_genes$SYMBOL)

MS_genes <- union(disgenet_genes, science_genes)

# Load modules 
diamond_methyl_cd4 <- readRDS("RDS_files/diamond_methyl_cd4.RDS")$module_genes 
diamond_methyl_cd8 <- readRDS("RDS_files/diamond_methyl_cd8.RDS")$module_genes 
diamond_rna_cd4 <- readRDS("RDS_files/diamond_rna_cd4.RDS")$module_genes 
diamond_rna_cd8 <- readRDS("RDS_files/diamond_rna_cd8.RDS")$module_genes 

# Overlap 
g_cd4 <- intersect(diamond_methyl_cd4, diamond_rna_cd4)
g_cd8 <- intersect(diamond_methyl_cd8, diamond_rna_cd8)

global_module <- intersect(diamond_methyl_cd4,diamond_rna_cd4) %>% intersect(diamond_methyl_cd8) %>% intersect(diamond_rna_cd8)

# 11.3 Save ####
write.table(diamond_methyl_cd4, "diamond_methyl_cd4.txt", col.names = F, row.names = F)
write.table(diamond_methyl_cd8, "diamond_methyl_cd8.txt", col.names = F, row.names = F)
write.table(diamond_rna_cd4, "diamond_rna_cd4.txt", col.names = F, row.names = F)
write.table(diamond_rna_cd8, "diamond_rna_cd8.txt", col.names = F, row.names = F)

write.table(g_cd4, "diamond_common_cd4.txt", col.names = F, row.names = F)
write.table(g_cd8, "diamond_common_cd8.txt", col.names = F, row.names = F)

# 11.4 Load seed genes ####
rna_cd4_seed <- read.table("gene_rebound_common_cd4.txt", header=F, row.names = NULL)$V1
rna_cd8_seed <- read.table("gene_rebound_common_cd8.txt", header=F, row.names = NULL)$V1
methyl_cd4_seed <- read.table("cpg_rebound_common_cd4_gene.txt", header=F, row.names = NULL)$V1
methyl_cd8_seed <- read.table("cpg_rebound_common_cd8_gene_01.txt", header=F, row.names = NULL)$V1

# Load Differential analysis files ####
ProbeFeatures <- readRDS(file=paste0(FOLDER_RDS, "TH2636_ProbeFeatures_211123.RDS"))

# CD4
CD4_3rd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD4_3rd_HP_methyl_all.RDS")
CD4_3rd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD4_3rd_MS_methyl_all.RDS")

CD4_PP_3rd_HP_methyl <- readRDS("RDS_files/DMR/CD4_PP_HP_methyl_all.RDS")
CD4_PP_3rd_MS_methyl <- readRDS("RDS_files/DMR/CD4_PP_MS_methyl_all.RDS")

CD4_2nd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD4_2nd_1st_HP_all.RDS")
CD4_2nd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD4_2nd_1st_MS_all.RDS")

# CD8
CD8_3rd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD8_3rd_HP_methyl_all.RDS")
CD8_3rd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD8_3rd_MS_methyl_all.RDS")

CD8_PP_3rd_HP_methyl <- readRDS("RDS_files/DMR/CD8_PP_HP_methyl_all.RDS")
CD8_PP_3rd_MS_methyl <- readRDS("RDS_files/DMR/CD8_PP_MS_methyl_all.RDS")

CD8_2nd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD8_2nd_1st_HP_all.RDS")
CD8_2nd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD8_2nd_1st_MS_all.RDS")

# Load RNA-seq ####
#CD4 RNA-seq
CD4_3rd_1st_HP_rna <- readRDS("RDS_files/DMG/CD4_3rd_HP_rna_all.RDS"); 
CD4_3rd_1st_MS_rna <- readRDS("RDS_files/DMG/CD4_3rd_MS_rna_all.RDS"); 

CD4_PP_3rd_HP_rna <- readRDS("RDS_files/DMG/CD4_PP_HP_rna_all.RDS");
CD4_PP_3rd_MS_rna <- readRDS("RDS_files/DMG/CD4_PP_MS_rna_all.RDS"); 

CD4_2nd_1st_HP_rna <- readRDS("RDS_files/DMG/HP_2nd_1st_cd4_all.RDS")
CD4_2nd_1st_MS_rna <- readRDS("RDS_files/DMG/MS_2nd_1st_cd4_all.RDS")

#CD8 RNA-seq
CD8_3rd_1st_HP_rna <- readRDS("RDS_files/DMG/CD8_3rd_HP_rna_all.RDS"); 
CD8_3rd_1st_MS_rna <- readRDS("RDS_files/DMG/CD8_3rd_MS_rna_all.RDS"); 

CD8_PP_3rd_HP_rna <- readRDS("RDS_files/DMG/CD8_PP_HP_rna_all.RDS"); 
CD8_PP_3rd_MS_rna <- readRDS("RDS_files/DMG/CD8_PP_MS_rna_all.RDS"); 

CD8_2nd_1st_HP_rna <- readRDS("RDS_files/DMG/HP_2nd_1st_cd8_all.RDS")
CD8_2nd_1st_MS_rna <- readRDS("RDS_files/DMG/MS_2nd_1st_cd8_all.RDS")

# Universe RNA-seq 
universe_cd4 <- readRDS("RDS_files/DMG/CD4_3rd_HP_rna_all.RDS") %>% rownames()
universe_cd8 <- readRDS("RDS_files/DMG/CD8_3rd_HP_rna_all.RDS") %>% rownames()

# Universe Methylation
universe_file <- "C:/Users/albze08/Desktop/phd/P4/methylation/RDS_files/universe_CD8.RDS" #name of file containing the list of genes
universe_methyl <- readRDS(universe_file) %>% unique()

# Universe of cpgs ####
list_cpg <- readRDS("RDS_files/DMR/CD4_2nd_1st_HP_all.RDS") %>% rownames()
cpg_gene <- probe_to_gene(list_cpg)
cpg_gene$feature <- ProbeFeatures$feature[match( cpg_gene$cpg,  rownames(ProbeFeatures))] %>% as.character()
cpg_gene <- cpg_gene[grep("TSS",cpg_gene$feature),]

# 11.5 MS-enrichment ####
ppi_network <- readRDS("data/ppi_network_symbol.RDS")
ppi_genes <- unique(c(ppi_network$symbol1, ppi_network$symbol2))

#seed genes
f_rna_cd4_seed <- fisher_test(rna_cd4_seed, universe_cd4, intersect(MS_genes,universe_cd4))
f_rna_cd8_seed <- fisher_test(rna_cd8_seed, universe_cd8, intersect(MS_genes,universe_cd8))
f_methyl_cd4_seed <- fisher_test(methyl_cd4_seed, universe_methyl, intersect(MS_genes,universe_methyl))
f_methyl_cd8_seed <- fisher_test(methyl_cd8_seed, universe_methyl, intersect(MS_genes,universe_methyl))

df <- data.frame(pval= c(f_rna_cd4_seed$p_value, f_rna_cd8_seed$p_value, f_methyl_cd4_seed$p_value, f_methyl_cd8_seed$p_value),
                 odds=c(f_rna_cd4_seed$odds_ratio, f_rna_cd8_seed$odds_ratio, f_methyl_cd4_seed$odds_ratio, f_methyl_cd8_seed$odds_ratio),
                 name=c("RNA CD4", "RNA CD8", "Methyl CD4", "Methyl CD8"))
df$name <- factor(df$name, levels=c("Methyl CD4", "Methyl CD8", "RNA CD4", "RNA CD8"))

p <- ggplot(df, aes(x=name, y=odds)) + geom_bar(stat="identity") +
  theme_bw()  +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pdf("figures_manus/MSenrich_seed.pdf")
p
dev.off()

#module genes
#CD4
universe_module_cd4_rna <- union(universe_cd4,ppi_genes) 
universe_module_cd4_methyl <- union(universe_methyl,ppi_genes) 
universe_module_cd4 <- intersect(universe_module_cd4_rna,universe_module_cd4_methyl)

f_methyl_cd4_module <- fisher_test(diamond_methyl_cd4, universe_module_cd4_methyl, intersect(MS_genes,universe_module_cd4_methyl))
f_rna_cd4_module <- fisher_test(diamond_rna_cd4, universe_module_cd4_rna, intersect(MS_genes,universe_module_cd4_rna))
f_both_cd4_module <- fisher_test(g_cd4, universe_module_cd4, intersect(MS_genes,universe_module_cd4))

#CD8
universe_module_cd8_rna <- union(universe_cd8,ppi_genes) 
universe_module_cd8_methyl <- union(universe_methyl,ppi_genes) 
universe_module_cd8 <- intersect(universe_module_cd8_rna,universe_module_cd8_methyl)

f_methyl_cd8_module <- fisher_test(diamond_methyl_cd8, universe_module_cd8_methyl, intersect(MS_genes,universe_module_cd8_methyl))
f_rna_cd8_module <- fisher_test(diamond_rna_cd8, universe_module_cd8_rna, intersect(MS_genes,universe_module_cd8_rna))
f_both_cd8_module <- fisher_test(g_cd8, universe_module_cd8, intersect(MS_genes,universe_module_cd8))

df <- data.frame(pval= c(f_methyl_cd4_module$p_value, f_rna_cd4_module$p_value, f_both_cd4_module$p_value, f_methyl_cd8_module$p_value, f_rna_cd8_module$p_value, f_both_cd8_module$p_value),
                 odds= c(f_methyl_cd4_module$odds_ratio, f_rna_cd4_module$odds_ratio, f_both_cd4_module$odds_ratio, f_methyl_cd8_module$odds_ratio, f_rna_cd8_module$odds_ratio, f_both_cd8_module$odds_ratio),
                 n_genes= c(f_methyl_cd4_module$n_genes, f_rna_cd4_module$n_genes, f_both_cd4_module$n_genes, f_methyl_cd8_module$n_genes, f_rna_cd8_module$n_genes, f_both_cd8_module$n_genes),
                 name=c("Methyl CD4", "RNA CD4", "Both CD4", "Methyl CD8", "RNA CD8", "Both CD8"))
df$name <- factor(df$name, levels=c("Methyl CD4", "Methyl CD8", "RNA CD4",  "RNA CD8", "Both CD4", "Both CD8"))

p <- ggplot(df, aes(x=name, y=odds)) + geom_bar(stat="identity") +
  theme_bw()  +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pdf("figures_manus/MSenrich_modules.pdf")
p
dev.off()


# 11.6 Enrichment of Activation ####
State_cd4 <- readRDS(paste0(FOLDER_RDS, "DMG/State_cd4.RDS"))
genes_state_CD4 <- rownames(State_cd4)[State_cd4$adj.P.Val<0.05]


State_cd8 <- readRDS(paste0(FOLDER_RDS, "DMG/State_cd8.RDS"))
genes_state_CD8 <- rownames(State_cd8)[State_cd8$adj.P.Val<0.05]

#Seed
fisher_test(rna_cd4_seed, universe_cd4, genes_state_CD4)
fisher_test(rna_cd8_seed, universe_cd8, genes_state_CD8)

p1 <- ggvenn(list(a=rna_cd4_seed, b=genes_state_CD4), c("a","b"),
             show_percentage = F)
p2 <- ggvenn(list(a=rna_cd8_seed, b=genes_state_CD8), c("a","b"),
             show_percentage = F)

pdf("figures_manus/Act_enrich_common_genes.pdf")
ggarrange(p1,p2, nrow=1, ncol=2)
dev.off()


#Modules
fisher_test(g_cd4, union(universe_module_cd4,genes_state_CD4), genes_state_CD4)
fisher_test(g_cd8, union(universe_module_cd8,genes_state_CD8), genes_state_CD8)

p1 <- ggvenn(list(a=g_cd4, b=genes_state_CD4), c("a","b"),
             show_percentage = F)
p2 <- ggvenn(list(a=g_cd8, b=genes_state_CD8), c("a","b"),
             show_percentage = F)

pdf("figures_manus/Act_enrich_common.pdf")
ggarrange(p1,p2, nrow=1, ncol=2)
dev.off()

# Overlap module-module 
fisher_test(diamond_rna_cd4, universe_cd4 %>% union(universe_methyl) %>% union(ppi_genes), diamond_methyl_cd4)
fisher_test(diamond_rna_cd8, universe_cd8 %>% union(universe_methyl) %>% union(ppi_genes), diamond_methyl_cd8)

# Overlap module-seed ####
fisher_test(rna_cd4_seed, universe_cd4 %>% union(universe_methyl) %>% union(ppi_genes), diamond_methyl_cd4)
fisher_test(methyl_cd4_seed, universe_cd4 %>% union(universe_methyl) %>% union(ppi_genes), diamond_rna_cd4)

fisher_test(rna_cd8_seed, universe_cd8 %>% union(universe_methyl) %>% union(ppi_genes), diamond_methyl_cd8)
fisher_test(methyl_cd8_seed, universe_cd8 %>% union(universe_methyl) %>% union(ppi_genes), diamond_rna_cd8)

intersect(rna_cd4_seed,diamond_methyl_cd4)
intersect(methyl_cd4_seed,diamond_rna_cd4)
intersect(rna_cd8_seed,diamond_methyl_cd8)
intersect(methyl_cd8_seed,diamond_rna_cd8)


# 10.7 Overlap with P4 genes ####
P4_genes_down <- read.table("data/stimP424vsstim24_DEGs_down.txt")$SYMBOL %>% unique()

#RNA-seq CD4 seed
f_rna_cd4_seed <- fisher_test(rna_cd4_seed, union(universe_cd4,P4_genes_down), P4_genes_down)

#RNA-seq CD4 module
f_rna_cd4_diamond <- fisher_test(diamond_rna_cd4, union(universe_cd4,P4_genes_down)  %>% union(ppi_genes), P4_genes_down)

#Methylation CD4 seed
f_methyl_cd4_seed <- fisher_test(methyl_cd4_seed, union(universe_methyl,P4_genes_down), P4_genes_down)

#Methylation CD4 module
f_methyl_cd4_diamond <- fisher_test(diamond_methyl_cd4, union(universe_methyl,P4_genes_down) %>% union(ppi_genes), P4_genes_down)


#RNA-seq CD8 seed
f_rna_cd8_seed <- fisher_test(rna_cd8_seed, union(universe_cd8,P4_genes_down), P4_genes_down)

#RNA-seq CD8 module
f_rna_cd8_diamond <- fisher_test(diamond_rna_cd8, union(universe_cd8,P4_genes_down)  %>% union(ppi_genes), P4_genes_down)

#Methylation CD8 seed
f_methyl_cd8_seed <- fisher_test(methyl_cd8_seed, union(universe_methyl,P4_genes_down), P4_genes_down)

#Methylation CD8 module
f_methyl_cd8_diamond <- fisher_test(diamond_methyl_cd8, union(universe_methyl,P4_genes_down) %>% union(ppi_genes), P4_genes_down)

#CD4 module
f_cd4 <- fisher_test(g_cd4, union(universe_module_cd4,P4_genes_down) %>% union(ppi_genes), P4_genes_down)
f_cd8 <- fisher_test(g_cd8, union(universe_module_cd8,P4_genes_down) %>% union(ppi_genes), P4_genes_down)


df <- data.frame(pval = -log10(c(f_rna_cd4_seed$p_value, f_rna_cd4_diamond$p_value, f_methyl_cd4_seed$p_value, f_methyl_cd4_diamond$p_value,
                          f_rna_cd8_seed$p_value, f_rna_cd8_diamond$p_value, f_methyl_cd8_seed$p_value, f_methyl_cd8_diamond$p_value,
                          f_cd4$p_value, f_cd8$p_value)),
                 OR = c(f_rna_cd4_seed$odds_ratio, f_rna_cd4_diamond$odds_ratio, f_methyl_cd4_seed$odds_ratio, f_methyl_cd4_diamond$odds_ratio,
                        f_rna_cd8_seed$odds_ratio, f_rna_cd8_diamond$odds_ratio, f_methyl_cd8_seed$odds_ratio, f_methyl_cd8_diamond$odds_ratio,
                        f_cd4$odds_ratio, f_cd8$odds_ratio),
                 group = c("RNA CD4 seed", "RNA CD4 DIAMOnD", "Methyl CD4 seed", "Methyl CD4 DIAMOnD",
                           "RNA CD8 seed", "RNA CD8 DIAMOnD", "Methyl CD8 seed", "Methyl CD8 DIAMOnD",
                           "Common CD4", "Common CD8"))
df$group <- factor(df$group, levels = c("RNA CD4 seed", "RNA CD4 DIAMOnD", "Methyl CD4 seed", "Methyl CD4 DIAMOnD",
                                        "RNA CD8 seed", "RNA CD8 DIAMOnD", "Methyl CD8 seed", "Methyl CD8 DIAMOnD",
                                        "Common CD4", "Common CD8"))
p <- ggplot(df, aes(x=group, y=OR)) + geom_bar(stat="identity") 

pdf("figures_manus/P4enrich.pdf", width=12)
p
dev.off()



