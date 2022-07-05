# This script generates figure 5 (comparison of MS and HC)
# Created and modified by Alberto Zenere, 2021-01-9

# 09_Figure_5 ####
# 9.1 Setup
# 9.2 Load MS-associated genes
# 9.3 Load DEGs and DMPs
# 9.4 Select rebound cpgs and genes for CD8
# 9.5 Plot for CD8
# 9.6 Select rebound cpgs and genes for CD4
# 9.7 Plot for CD4
# 9.8 Enrichment of Activation
# 9.9 Create list CD8 methylation with stricter threshold
# 9.10 MS enrichment on genes that rebound in MS but not HC

# 9.1 Setup ####
rm(list=ls()) # remove all entries in the global environment 

pathlist <- .libPaths();
newpath <-  "C:/Users/albze08/.conda/envs/graMS/Lib/R/library"
.libPaths(newpath)

# Set directory structure 
main_dir <-"C:/Users/albze08/Desktop/phd/P4/methylation"
FOLDER_RDS <- "RDS_files/"
FOLDER_FIGURES <- "figures/"

setwd("C:/Users/albze08/Desktop/phd/P4/methylation")

# Load packages 

pack_R <- c("limma","stringr","dplyr","tidyverse","pracma","ggrepel")

for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

require("lumi")
require("ggpubr")
require("pheatmap")
library("VennDiagram") 
library("gridExtra")
library("ggalluvial")
library("ggvenn")

set.seed(206)

# My functions 
source('probe_to_gene.R')
source('fisher_test.R')
source('fisher_test_less.R')


# 9.2 Load MS-associated genes ####
# Load DisGeNET
disgenet <- read.csv("data/MS_DisGeNET.txt", sep = "\t")
disgenet_genes <- disgenet$geneSymbol %>% unique()

# Load Science 2019
science_genes <- read.csv(file = "C:/Users/albze08/Desktop/phd/P4/RNAseq/data/complete_MS_annotation_p10e-6.csv")
science_genes <- unique(science_genes$SYMBOL)

#combine
MS_genes <- union(disgenet_genes, science_genes)

# 9.3 Load DEGs and DMPs ####
#CD4 RNA
CD4_3rd_HP_rna <- readRDS("RDS_files/DMG/CD4_3rd_HP_rna_all.RDS")
CD4_PP_HP_rna <- readRDS("RDS_files/DMG/CD4_PP_HP_rna_all.RDS")

CD4_3rd_MS_rna <- readRDS("RDS_files/DMG/CD4_3rd_MS_rna_all.RDS")
CD4_PP_MS_rna <- readRDS("RDS_files/DMG/CD4_PP_MS_rna_all.RDS")

#CD8 RNA
CD8_3rd_HP_rna <- readRDS("RDS_files/DMG/CD8_3rd_HP_rna_all.RDS")
CD8_PP_HP_rna <- readRDS("RDS_files/DMG/CD8_PP_HP_rna_all.RDS")

CD8_3rd_MS_rna <- readRDS("RDS_files/DMG/CD8_3rd_MS_rna_all.RDS")
CD8_PP_MS_rna <- readRDS("RDS_files/DMG/CD8_PP_MS_rna_all.RDS")

#CD4 methyl
CD4_3rd_HP_methyl <- readRDS("RDS_files/DMR/CD4_3rd_HP_methyl_all.RDS")
CD4_PP_HP_methyl <- readRDS("RDS_files/DMR/CD4_PP_HP_methyl_all.RDS")

CD4_3rd_MS_methyl <- readRDS("RDS_files/DMR/CD4_3rd_MS_methyl_all.RDS")
CD4_PP_MS_methyl <- readRDS("RDS_files/DMR/CD4_PP_MS_methyl_all.RDS")

#CD8 methyl
universe_methyl <- readRDS("C:/Users/albze08/Desktop/phd/P4/methylation/RDS_files/universe_CD8.RDS")

CD8_3rd_HP_methyl <- readRDS("RDS_files/DMR/CD8_3rd_HP_methyl_all.RDS")
CD8_PP_HP_methyl <- readRDS("RDS_files/DMR/CD8_PP_HP_methyl_all.RDS")

CD8_3rd_MS_methyl <- readRDS("RDS_files/DMR/CD8_3rd_MS_methyl_all.RDS")
CD8_PP_MS_methyl <- readRDS("RDS_files/DMR/CD8_PP_MS_methyl_all.RDS")


# Re-order ####
list_cd4_rna <- rownames(CD4_3rd_HP_rna)
list_cd8_rna <- rownames(CD8_3rd_HP_rna)

list_cd4_methyl <- rownames(CD4_3rd_HP_methyl)
list_cd8_methyl <- rownames(CD8_3rd_HP_methyl)

CD4_PP_HP_rna <- CD4_PP_HP_rna[list_cd4_rna,]
CD4_3rd_MS_rna <- CD4_3rd_MS_rna[list_cd4_rna,]
CD4_PP_MS_rna <- CD4_PP_MS_rna[list_cd4_rna,]

CD8_PP_HP_rna <- CD8_PP_HP_rna[list_cd8_rna,]
CD8_3rd_MS_rna <- CD8_3rd_MS_rna[list_cd8_rna,]
CD8_PP_MS_rna <- CD8_PP_MS_rna[list_cd8_rna,]

CD4_PP_HP_methyl <- CD4_PP_HP_methyl[list_cd4_methyl,]
CD4_3rd_MS_methyl <- CD4_3rd_MS_methyl[list_cd4_methyl,]
CD4_PP_MS_methyl <- CD4_PP_MS_methyl[list_cd4_methyl,]

CD8_PP_HP_methyl <- CD8_PP_HP_methyl[list_cd8_methyl,]
CD8_3rd_MS_methyl <- CD8_3rd_MS_methyl[list_cd8_methyl,]
CD8_PP_MS_methyl <- CD8_PP_MS_methyl[list_cd8_methyl,]


# 9.4 Select rebound cpgs and genes for CD8 ####

#RNA
rebound_gene_cd8_hp <- intersect(rownames(CD8_3rd_HP_rna)[CD8_3rd_HP_rna$P.Value<0.05],
                                  rownames(CD8_PP_HP_rna)[CD8_PP_HP_rna$P.Value<0.05])

rebound_gene_cd8_ms <- intersect(rownames(CD8_3rd_MS_rna)[CD8_3rd_MS_rna$P.Value<0.05],
                                     rownames(CD8_PP_MS_rna)[CD8_PP_MS_rna$P.Value<0.05])

rebound_gene_cd8 <- union(rebound_gene_cd8_hp, rebound_gene_cd8_ms)


#Methylation
rebound_cpg_cd8_hp <- intersect(rownames(CD8_3rd_HP_methyl)[CD8_3rd_HP_methyl$P.Value<0.05 & abs(CD8_3rd_HP_methyl$delta_beta)>0.05],
                                 rownames(CD8_PP_HP_methyl)[CD8_PP_HP_methyl$P.Value<0.05 & abs(CD8_PP_HP_methyl$delta_beta)>0.05])

rebound_cpg_cd8_ms <- intersect(rownames(CD8_3rd_MS_methyl)[CD8_3rd_MS_methyl$P.Value<0.05 & abs(CD8_3rd_MS_methyl$delta_beta)>0.05],
                                 rownames(CD8_PP_MS_methyl)[CD8_PP_MS_methyl$P.Value<0.05 & abs(CD8_PP_MS_methyl$delta_beta)>0.05])

rebound_cpg_cd8 <- union(rebound_cpg_cd8_hp, rebound_cpg_cd8_ms)
rebound_cpg_cd8_gene <- probe_to_gene(rebound_cpg_cd8)$gene

# 9.5 Plot for CD8 ####

plot_logFC_vs_logFC <- function(v1, v2, gene_list=NULL){
  
  #gather
  Npoints <- min(length(v1), 2e2)
  df <- data.frame(v1=v1, v2=v2)
  
  if (!is.null(gene_list)){
    is_MS_associated <- gene_list %in% MS_genes
  } else {
      is_MS_associated <- rep(FALSE,nrow(df))
  }
  df$is_MS <- is_MS_associated
  df_subset <- df[sample(1:nrow(df), Npoints),]
  
  
  #calculate percentages
  bottom_right <- round( 100*sum(df$v1>0 & v2<0)/length(v1), digits = 1 )
  up_left <- round( 100*sum(v1<0 & v2>0)/length(v1), digits = 1 )
  
  up_right_down <- round( 100*sum(v1>v2 & v1>0 & v2>0)/length(v1), digits = 1 )
  up_right_up <- round( 100*sum(v1<v2 & v1>0 & v2>0)/length(v1), digits = 1 )
  
  down_left_down <- round( 100*sum(v1>v2 & v1<0 & v2<0)/length(v1), digits = 1 )
  down_left_up <- round( 100*sum(v1<v2 & v1<0 & v2<0)/length(v1), digits = 1 )
  
  
  #figure limits
  x_low <- min(df_subset$v1)
  x_high <- max(df_subset$v1)
  y_low <- min(df_subset$v2)
  y_high <- max(df_subset$v2)
  
  x_low <- min(x_low, y_low)
  y_low <- x_low
  x_high <- max(x_high, y_high)
  y_high <- x_high
  
  
  #plot
  p <- ggplot() +
    geom_point( data=df_subset, aes(x=v1, y=v2, col=is_MS), alpha=1, shape=19, size=.15) +
    scale_color_manual(values = c("black", "#FF3030")) +
    stat_density_2d_filled(data=df, aes(x=v1, y=v2), alpha=0.7, geom = "polygon", bins=4) +
    scale_fill_manual(values = c("white", "#87CEFF", "#1874CD", "blue4")) +
    stat_density_2d(bins=4, colour="black") +
    annotate("text",x=x_high, y=y_low, label= paste0( bottom_right, "%"), color="black") +
    annotate("text",x=x_low, y=y_high, label= paste0( up_left, "%"), color="black") +
    annotate("text",x=x_high, y=.1*y_high, label= paste0( up_right_down, "%"), color="black") +
    annotate("text",x=0.2*x_high, y=y_high, label= paste0( up_right_up, "%"), color="black") +
    annotate("text",x=0.2*x_low, y=y_low, label= paste0( down_left_down, "%"), color="black") +
    annotate("text",x=x_low, y=0.7*y_low, label= paste0( down_left_up, "%"), color="black") +
    geom_abline(intercept = 0, slope = 1) +
    geom_vline(xintercept=0) +
    geom_hline(yintercept=0) +
    theme_bw()  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    xlim(1.1*c(x_low,x_high)) + ylim(1.1*c(y_low,y_high)) +
    theme(legend.position="none", text = element_text(size=15)) +
    xlab("HC") + ylab("MS")

    
  return(p)
 
}


p_3rd_rna_CD8 <- plot_logFC_vs_logFC(CD8_3rd_HP_rna[rebound_gene_cd8, "logFC"], CD8_3rd_MS_rna[rebound_gene_cd8, "logFC"], rebound_gene_cd8) +
  ggtitle("RNA-seq: 3rd-1st") 
p_PP_rna_CD8 <- plot_logFC_vs_logFC(CD8_PP_HP_rna[rebound_gene_cd8, "logFC"], CD8_PP_MS_rna[rebound_gene_cd8, "logFC"], rebound_gene_cd8) +
  ggtitle("RNA-seq: PP-3rd") 

p_3rd_methyl_CD8 <- plot_logFC_vs_logFC(CD8_3rd_HP_methyl[rebound_cpg_cd8, "logFC"], CD8_3rd_MS_methyl[rebound_cpg_cd8, "logFC"], rebound_cpg_cd8_gene) +
  ggtitle("Methylation: 3rd-1st") 
p_PP_methyl_CD8 <- plot_logFC_vs_logFC(CD8_PP_HP_methyl[rebound_cpg_cd8, "logFC"], CD8_PP_MS_methyl[rebound_cpg_cd8, "logFC"], rebound_cpg_cd8_gene) +
  ggtitle("Methylation: PP-3rd") 


# 9.6 Select rebound cpgs and genes for CD4 ####

#RNA
rebound_gene_cd4_hp <- intersect(rownames(CD4_3rd_HP_rna)[CD4_3rd_HP_rna$P.Value<0.05],
                                 rownames(CD4_PP_HP_rna)[CD4_PP_HP_rna$P.Value<0.05])

rebound_gene_cd4_ms <- intersect(rownames(CD4_3rd_MS_rna)[CD4_3rd_MS_rna$P.Value<0.05],
                                 rownames(CD4_PP_MS_rna)[CD4_PP_MS_rna$P.Value<0.05])

rebound_gene_cd4 <- union(rebound_gene_cd4_hp, rebound_gene_cd4_ms)


#Methylation
rebound_cpg_cd4_hp <- intersect(rownames(CD4_3rd_HP_methyl)[CD4_3rd_HP_methyl$P.Value<0.05 & abs(CD4_3rd_HP_methyl$delta_beta)>0.05],
                                rownames(CD4_PP_HP_methyl)[CD4_PP_HP_methyl$P.Value<0.05 & abs(CD4_PP_HP_methyl$delta_beta)>0.05])

rebound_cpg_cd4_ms <- intersect(rownames(CD4_3rd_MS_methyl)[CD4_3rd_MS_methyl$P.Value<0.05 & abs(CD4_3rd_MS_methyl$delta_beta)>0.05],
                                rownames(CD4_PP_MS_methyl)[CD4_PP_MS_methyl$P.Value<0.05 & abs(CD4_PP_MS_methyl$delta_beta)>0.05])

rebound_cpg_cd4 <- union(rebound_cpg_cd4_hp, rebound_cpg_cd4_ms)
rebound_cpg_cd4_gene <- probe_to_gene(rebound_cpg_cd4)$gene


# 9.7 Plot for CD8 ####

p_3rd_rna_CD4 <- plot_logFC_vs_logFC(CD4_3rd_HP_rna[rebound_gene_cd4, "logFC"], CD4_3rd_MS_rna[rebound_gene_cd4, "logFC"], rebound_gene_cd4) +
  ggtitle("RNA-seq: 3rd-1st") 
p_PP_rna_CD4 <- plot_logFC_vs_logFC(CD4_PP_HP_rna[rebound_gene_cd4, "logFC"], CD4_PP_MS_rna[rebound_gene_cd4, "logFC"], rebound_gene_cd4) +
  ggtitle("RNA-seq: PP-3rd") 

p_3rd_methyl_CD4 <- plot_logFC_vs_logFC(CD4_3rd_HP_methyl[rebound_cpg_cd4, "logFC"], CD4_3rd_MS_methyl[rebound_cpg_cd4, "logFC"], rebound_cpg_cd4_gene) +
  ggtitle("Methylation: 3rd-1st") 
p_PP_methyl_CD4 <- plot_logFC_vs_logFC(CD4_PP_HP_methyl[rebound_cpg_cd4, "logFC"], CD4_PP_MS_methyl[rebound_cpg_cd4, "logFC"], rebound_cpg_cd4_gene) +
  ggtitle("Methylation: PP-3rd") 

pdf("figures_manus/MSvsHP.pdf", width=14)
ggarrange(p_3rd_methyl_CD4, p_PP_methyl_CD4, p_3rd_methyl_CD8, p_PP_methyl_CD8,
          p_3rd_rna_CD4, p_PP_rna_CD4, p_3rd_rna_CD8, p_PP_rna_CD8, nrow=2, ncol=4)
dev.off()


#Enrichment of MS-associated genes ####
#CD4 RNA
gene_rebound_common_cd4 <- list_cd4_rna[ (CD4_3rd_HP_rna$logFC*CD4_3rd_MS_rna$logFC>0) &
                                           (CD4_PP_HP_rna$logFC*CD4_PP_MS_rna$logFC>0)] %>% intersect(rebound_gene_cd4)

fisher_cd4_rna <- fisher_test(gene_rebound_common_cd4, list_cd4_rna, intersect(MS_genes,list_cd4_rna))

#CD8 RNA
gene_rebound_common_cd8 <- list_cd8_rna[ (CD8_3rd_HP_rna$logFC*CD8_3rd_MS_rna$logFC>0) &
                                           (CD8_PP_HP_rna$logFC*CD8_PP_MS_rna$logFC>0)] %>% intersect(rebound_gene_cd8)

fisher_cd8_rna <- fisher_test(gene_rebound_common_cd8, list_cd8_rna, intersect(MS_genes,list_cd8_rna))


#CD4 methyl
cpg_rebound_common_cd4 <- list_cd4_methyl[ (CD4_3rd_HP_methyl$logFC*CD4_3rd_MS_methyl$logFC>0) &
                                           (CD4_PP_HP_methyl$logFC*CD4_PP_MS_methyl$logFC>0)] %>% intersect(rebound_cpg_cd4)
cpg_rebound_common_cd4_gene <- probe_to_gene(cpg_rebound_common_cd4)
cpg_rebound_common_cd4_gene <- cpg_rebound_common_cd4_gene$gene %>% unique()
cpg_rebound_common_cd4_gene <- cpg_rebound_common_cd4_gene[cpg_rebound_common_cd4_gene!=""]

fisher_cd4_methyl <- fisher_test(cpg_rebound_common_cd4_gene, universe_methyl, intersect(MS_genes,universe_methyl))


#CD8 methyl
cpg_rebound_common_cd8 <- list_cd8_methyl[ (CD8_3rd_HP_methyl$logFC*CD8_3rd_MS_methyl$logFC>0) &
                                             (CD8_PP_HP_methyl$logFC*CD8_PP_MS_methyl$logFC>0)] %>% intersect(rebound_cpg_cd8)

cpg_rebound_common_cd8_gene <- probe_to_gene(cpg_rebound_common_cd8)
cpg_rebound_common_cd8_gene <- cpg_rebound_common_cd8_gene$gene %>% unique()
cpg_rebound_common_cd8_gene <- cpg_rebound_common_cd8_gene[cpg_rebound_common_cd8_gene!=""]

fisher_cd8_methyl <- fisher_test(cpg_rebound_common_cd8_gene, universe_methyl, intersect(MS_genes,universe_methyl))

#Plot enrichment results
df <- data.frame(group=c("CD4+ RNA-seq", "CD8+ RNA-seq", "CD4+ methylation", "CD8+ methylation"),
                 pval= -log10(c(fisher_cd4_rna$p_value, fisher_cd8_rna$p_value, fisher_cd4_methyl$p_value,
                        fisher_cd8_methyl$p_value)),
                 odds = c(fisher_cd4_rna$odds_ratio, fisher_cd8_rna$odds_ratio, fisher_cd4_methyl$odds_ratio,
                          fisher_cd8_methyl$odds_ratio))
p1 <- ggplot(df, aes(x=group, y=pval)) + geom_bar(stat="identity") + 
  geom_hline(yintercept=-log10(0.05)) +
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=-90, hjust=1)) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
p2 <- ggplot(df, aes(x=group, y=odds)) + geom_bar(stat="identity") + 
  theme_bw() + 
  theme(axis.title.x = element_blank(), axis.text.x=element_text(angle=-90, hjust=1)) +  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pdf("figures_manus/MS_enrich_common_genes.pdf")
ggarrange(p1,p2, nrow=2, ncol=1)
dev.off()

#Write lists ####
write.table(gene_rebound_common_cd4, "gene_rebound_common_cd4.txt", row.names = F, col.names = F)
write.table(gene_rebound_common_cd8, "gene_rebound_common_cd8.txt", row.names = F, col.names = F)

write.table(cpg_rebound_common_cd4_gene, "cpg_rebound_common_cd4_gene.txt", row.names = F, col.names = F)
write.table(cpg_rebound_common_cd8_gene, "cpg_rebound_common_cd8_gene.txt", row.names = F, col.names = F)



# 9.8 Enrichment of Activation ####
State_cd4 <- readRDS(paste0(FOLDER_RDS, "DMG/State_cd4.RDS"))
genes_state_CD4 <- rownames(State_cd4)[State_cd4$adj.P.Val<0.05]

fisher_test(gene_rebound_common_cd4, list_cd4_rna, intersect(genes_state_CD4,list_cd4_rna))


State_cd8 <- readRDS(paste0(FOLDER_RDS, "DMG/State_cd8.RDS"))
genes_state_CD8 <- rownames(State_cd8)[State_cd8$adj.P.Val<0.05]

fisher_test(gene_rebound_common_cd8, list_cd8_rna, intersect(genes_state_CD8,list_cd8_rna))


p1 <- ggvenn(list(a=gene_rebound_common_cd4, b=genes_state_CD4), c("a","b"),
       show_percentage = F)
p2 <- ggvenn(list(a=gene_rebound_common_cd8, b=genes_state_CD8), c("a","b"),
       show_percentage = F)

pdf("figures_manus/Act_enrich_common_genes.pdf")
ggarrange(p1,p2, nrow=1, ncol=2)
dev.off()



# 9.9 Create list CD8 methylation with stricter threshold #####
rebound_cpg_cd8_hp <- intersect(rownames(CD8_3rd_HP_methyl)[CD8_3rd_HP_methyl$P.Value<0.05 & abs(CD8_3rd_HP_methyl$delta_beta)>0.1],
                                rownames(CD8_PP_HP_methyl)[CD8_PP_HP_methyl$P.Value<0.05 & abs(CD8_PP_HP_methyl$delta_beta)>0.1])

rebound_cpg_cd8_ms <- intersect(rownames(CD8_3rd_MS_methyl)[CD8_3rd_MS_methyl$P.Value<0.05 & abs(CD8_3rd_MS_methyl$delta_beta)>0.05],
                                rownames(CD8_PP_MS_methyl)[CD8_PP_MS_methyl$P.Value<0.05 & abs(CD8_PP_MS_methyl$delta_beta)>0.05])

rebound_cpg_cd8 <- union(rebound_cpg_cd8_hp, rebound_cpg_cd8_ms)
cpg_rebound_common_cd8 <- list_cd8_methyl[ (CD8_3rd_HP_methyl$logFC*CD8_3rd_MS_methyl$logFC>0) &
                                             (CD8_PP_HP_methyl$logFC*CD8_PP_MS_methyl$logFC>0)] %>% intersect(rebound_cpg_cd8)

cpg_rebound_common_cd8_gene <- probe_to_gene(cpg_rebound_common_cd8)
cpg_rebound_common_cd8_gene <- cpg_rebound_common_cd8_gene$gene %>% unique()
cpg_rebound_common_cd8_gene <- cpg_rebound_common_cd8_gene[cpg_rebound_common_cd8_gene!=""]

write.table(cpg_rebound_common_cd8_gene, "cpg_rebound_common_cd8_gene_01.txt", row.names = F, col.names = F)



# 9.10 MS enrichment on genes that rebound in MS but not HC ####
#Methyl CD4
rebound_cpg_cd4_ms <- intersect(rownames(CD4_3rd_MS_methyl)[CD4_3rd_MS_methyl$P.Value<0.05 & abs(CD4_3rd_MS_methyl$delta_beta)>0.05],
                                rownames(CD4_PP_MS_methyl)[CD4_PP_MS_methyl$P.Value<0.05 & abs(CD4_PP_MS_methyl$delta_beta)>0.05])

cpg_rebound_ms_cd4 <- list_cd4_methyl[ (CD4_3rd_MS_methyl$logFC*CD4_PP_MS_methyl$logFC<0) &
                                             (CD4_3rd_HP_methyl$logFC*CD4_PP_HP_methyl$logFC>0)] %>% intersect(rebound_cpg_cd4_ms)

cpg_rebound_ms_cd4 <- probe_to_gene(cpg_rebound_ms_cd4)
cpg_rebound_ms_cd4 <- cpg_rebound_ms_cd4$gene %>% unique()
cpg_rebound_ms_cd4_gene <- cpg_rebound_ms_cd4[cpg_rebound_ms_cd4!=""]

fisher_test(cpg_rebound_ms_cd4_gene, universe_methyl, intersect(MS_genes, universe_methyl))

#Methyl CD8
rebound_cpg_cd8_ms <- intersect(rownames(CD8_3rd_MS_methyl)[CD8_3rd_MS_methyl$P.Value<0.05 & abs(CD8_3rd_MS_methyl$delta_beta)>0.05],
                                rownames(CD8_PP_MS_methyl)[CD8_PP_MS_methyl$P.Value<0.05 & abs(CD8_PP_MS_methyl$delta_beta)>0.05])

cpg_rebound_ms_cd8 <- list_cd8_methyl[ (CD8_3rd_MS_methyl$logFC*CD8_PP_MS_methyl$logFC<0) &
                                         (CD8_3rd_HP_methyl$logFC*CD8_PP_HP_methyl$logFC>0)] %>% intersect(rebound_cpg_cd8_ms)

cpg_rebound_ms_cd8 <- probe_to_gene(cpg_rebound_ms_cd8)
cpg_rebound_ms_cd8 <- cpg_rebound_ms_cd8$gene %>% unique()
cpg_rebound_ms_cd8_gene <- cpg_rebound_ms_cd8[cpg_rebound_ms_cd8!=""]

fisher_test(cpg_rebound_ms_cd8_gene, universe_methyl, intersect(MS_genes, universe_methyl))


#RNA CD4
gene_rebound_ms_cd4 <- intersect(rownames(CD4_3rd_MS_rna)[CD4_3rd_MS_rna$P.Value<0.05],
                                rownames(CD4_PP_MS_rna)[CD4_PP_MS_rna$P.Value<0.05])
gene_rebound_ms_cd4 <- list_cd4_rna[ (CD4_3rd_MS_rna$logFC*CD4_PP_MS_rna$logFC<0) & 
                                       (CD4_3rd_HP_rna$logFC*CD4_PP_HP_rna$logFC>0)] %>% intersect(gene_rebound_ms_cd4)

fisher_test(gene_rebound_ms_cd4, list_cd4_rna, intersect(MS_genes, list_cd4_rna))

#RNA CD8
gene_rebound_ms_cd8 <- intersect(rownames(CD8_3rd_MS_rna)[CD8_3rd_MS_rna$P.Value<0.05],
                                 rownames(CD8_PP_MS_rna)[CD8_PP_MS_rna$P.Value<0.05])
gene_rebound_ms_cd8 <- list_cd8_rna[ (CD8_3rd_MS_rna$logFC*CD8_PP_MS_rna$logFC<0) & 
                                       (CD8_3rd_HP_rna$logFC*CD8_PP_HP_rna$logFC>0)] %>% intersect(gene_rebound_ms_cd8)

fisher_test(gene_rebound_ms_cd8, list_cd8_rna, intersect(MS_genes, list_cd8_rna))

