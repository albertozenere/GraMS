
rm(list=ls()) # remove all entries in the global environment 

# Set R library ####
pathlist <- .libPaths();
newpath <-  "C:/Users/albze08/.conda/envs/P3/Lib/R/library"
.libPaths(newpath)

# Load packages ####
library("reticulate")
library("MODifieR")
library("DOSE")
library("biomaRt")
library("ggplot2")
library("ggpubr")
library("xlsx")

py_config() #this tells you which python installation reticulate (used by modifier) is using. If packages are missing then they should be installed there using pip from cmd (be in that folder)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

# Set directory structure ####
main_dir <-"C:/Users/albze08/Desktop/phd/P4/methylation"
FOLDER_RDS <- "RDS_files/"
FOLDER_FIGURES <- "figures/"
setwd("C:/Users/albze08/Desktop/phd/P4/methylation")


# Import and format input files ####
disgenet <- read.csv("data/MS_DisGeNET.txt", sep = "\t")
disgenet_genes <- disgenet$geneSymbol %>% unique()

science_genes <- read.csv(file = "C:/Users/albze08/Desktop/phd/P4/RNAseq/data/complete_MS_annotation_p10e-6.csv")
science_genes <- unique(science_genes$SYMBOL)


# My functions ####
#Do the hypergeometric test on new list
source("hyge_test_2.R")
source("fisher_test.R")

# Load PPI network ####
#ppi_network <- read.delim("data/ppi_network.txt", 
#                          stringsAsFactors = FALSE)
#ppi_genes <- union(ppi_network$entrez1, ppi_network$entrez2) %>% unique()
# symbol_entrez <- getBM(attributes <- c("entrezgene_id", "external_gene_name"), 
#                        filters= "entrezgene_id",
#                        values = unique(c(ppi_network$entrez1, ppi_network$entrez2)), 
#                        mart = mart)
# 
# ord1 <- match(ppi_network$entrez1, symbol_entrez$entrezgene_id)
# ord2 <- match(ppi_network$entrez2, symbol_entrez$entrezgene_id)
# 
# idx <- (!is.na(ord1)) & (!is.na(ord2))
# 
# 
# symbol1 <- symbol_entrez$external_gene_name[ord1[idx]]
# symbol2 <- symbol_entrez$external_gene_name[ord2[idx]]
# ppi_network <- data.frame(symbol1, symbol2)
# 
# colnames(ppi_network) <- c("symbol1", "symbol2")
# saveRDS(ppi_network, "data/ppi_network_symbol.RDS")

ppi_network <- readRDS("data/ppi_network_symbol.RDS")
ppi_genes <- unique(c(ppi_network$symbol1, ppi_network$symbol2))

#diamond function ####
# build_clique_db(ppi_network = ppi_network, #to be run only one time, on the entrez-id ppi network
#                 db_folder = "./data",
#                 db_name = "diamond_db")

diamond_genes <- function(list_genes, ppi_network){
  
  seed_input <- data.frame(seed = list_genes, pval = rep(0,length(list_genes))) # add a fake p-value to keep them in the module
  diamond_input <- create_custom_microarray_input_object(diff_genes = seed_input) #input for module methods
  
  module_diamond <- diamond(MODifieR_input = diamond_input,
                            ppi_network = ppi_network,
                            deg_cutoff = 0.05,
                            n_output_genes = 200,
                            seed_weight = 10,
                            include_seed = TRUE)
  genes <- module_diamond$module_genes %>% unique()
  genes <- genes[genes!=""]
  genes <- na.omit(genes)
  
  return(list(seed=list_genes, module_genes=genes))
}

#Module Inference ####

#seed genes
rna_cd4_seed <- read.table("gene_rebound_common_cd4.txt", header=F, row.names = NULL)$V1
rna_cd8_seed <- read.table("gene_rebound_common_cd8.txt", header=F, row.names = NULL)$V1
methyl_cd4_seed <- read.table("cpg_rebound_common_cd4_gene.txt", header=F, row.names = NULL)$V1
methyl_cd8_seed <- read.table("cpg_rebound_common_cd8_gene_01.txt", header=F, row.names = NULL)$V1

# diamond ####
diamond_methyl_cd4 <- diamond_genes(methyl_cd4_seed, ppi_network)
diamond_methyl_cd8 <- diamond_genes(methyl_cd8_seed, ppi_network)

diamond_rna_cd4 <- diamond_genes(rna_cd4_seed, ppi_network)
diamond_rna_cd8 <- diamond_genes(rna_cd8_seed, ppi_network)

#Save ####
saveRDS(diamond_methyl_cd4, "RDS_files/diamond_methyl_cd4.RDS")
saveRDS(diamond_methyl_cd8, "RDS_files/diamond_methyl_cd8.RDS")
saveRDS(diamond_rna_cd4, "RDS_files/diamond_rna_cd4.RDS")
saveRDS(diamond_rna_cd8, "RDS_files/diamond_rna_cd8.RDS")

#module lists
df <- data.frame(module_genes = diamond_methyl_cd4$module_genes, is_seed_gene = diamond_methyl_cd4$module_genes %in% diamond_methyl_cd4$seed)
write.xlsx( df, file="DIAMOnD_modules.xlsx", sheetName="Methylation CD4", row.names=FALSE)
df <- data.frame(module_genes = diamond_rna_cd4$module_genes, is_seed_gene = diamond_rna_cd4$module_genes %in% diamond_rna_cd4$seed)
write.xlsx( df, file="DIAMOnD_modules.xlsx", sheetName="RNA-seq CD4", append=TRUE, row.names=FALSE)
df <- data.frame(module_genes = diamond_methyl_cd8$module_genes, is_seed_gene = diamond_methyl_cd8$module_genes %in% diamond_methyl_cd8$seed)
write.xlsx( df, file="DIAMOnD_modules.xlsx", sheetName="Methylation CD8", append=TRUE, row.names=FALSE)
df <- data.frame(module_genes = diamond_rna_cd8$module_genes, is_seed_gene = diamond_rna_cd8$module_genes %in% diamond_rna_cd8$seed)
write.xlsx( df, file="DIAMOnD_modules.xlsx", sheetName="RNA-seq CD8", append=TRUE, row.names=FALSE)


module_cd4 <- intersect(diamond_rna_cd4$module_genes, diamond_methyl_cd4$module_genes)
df <- data.frame(module_genes = module_cd4, is_seed_gene = module_cd4 %in% union(diamond_rna_cd4$seed, diamond_methyl_cd4$seed))
write.xlsx( df, file="DIAMOnD_modules.xlsx", sheetName="CD4 combined", append=TRUE, row.names=FALSE)
module_cd8 <- intersect(diamond_rna_cd8$module_genes, diamond_methyl_cd8$module_genes)
df <- data.frame(module_genes = module_cd8, is_seed_gene = module_cd8 %in% union(diamond_rna_cd8$seed, diamond_methyl_cd8$seed))
write.xlsx( df, file="DIAMOnD_modules.xlsx", sheetName="CD8 combined", append=TRUE, row.names=FALSE)

#Read ####
diamond_methyl_cd4 <- readRDS("RDS_files/diamond_methyl_cd4.RDS")
diamond_methyl_cd8 <- readRDS("RDS_files/diamond_methyl_cd8.RDS")
diamond_rna_cd4 <- readRDS("RDS_files/diamond_rna_cd4.RDS")
diamond_rna_cd8 <- readRDS("RDS_files/diamond_rna_cd8.RDS")

# Load files ####
#universe for seed genes
universe_methyl <- readRDS("C:/Users/albze08/Desktop/phd/P4/methylation/RDS_files/universe_CD8.RDS") %>% as.character()

universe_rna_cd4 <- readRDS("C:/Users/albze08/Desktop/phd/P4/RNAseq/RDS/universe/universe_cd4.RDS")
universe_rna_cd8 <- readRDS("C:/Users/albze08/Desktop/phd/P4/RNAseq/RDS/universe/universe_cd8.RDS")

# Enrichment ####
#Universe for modules
universe_module_methyl <- union(universe_methyl, ppi_genes)

universe_module_cd4 <- union(universe_rna_cd4, ppi_genes)
universe_module_cd8 <- union(universe_rna_cd8, ppi_genes)

#MS-associated genes
MS_genes <- union(science_genes, disgenet_genes)

#methylation seed
hyge_methyl_cd4 <- fisher_test(methyl_cd4_seed, universe_methyl, intersect(MS_genes,universe_methyl))
hyge_methyl_cd8 <- fisher_test(methyl_cd8_seed, universe_methyl, intersect(MS_genes,universe_methyl))

#methylation modules
hyge_methyl_cd4_module <- fisher_test(diamond_methyl_cd4$module_genes, universe_module_methyl, intersect(MS_genes,universe_module_methyl))
hyge_methyl_cd8_module <- fisher_test(diamond_methyl_cd8$module_genes, universe_module_methyl, intersect(MS_genes,universe_module_methyl))

#rna seed
hyge_rna_cd4 <- fisher_test(rna_cd4_seed, universe_rna_cd4, intersect(MS_genes,universe_rna_cd4))
hyge_rna_cd8 <- fisher_test(rna_cd8_seed, universe_rna_cd8, intersect(MS_genes,universe_rna_cd8))

#rna modules
hyge_rna_cd4_module <- fisher_test(diamond_rna_cd4$module_genes, universe_module_cd4, intersect(MS_genes,universe_module_cd4))
hyge_rna_cd8_module <- fisher_test(diamond_rna_cd8$module_genes, universe_module_cd8, intersect(MS_genes,universe_module_cd8))


#Overlap between modules ####
fisher_test(methyl_cd4_seed, union(universe_methyl, universe_rna_cd4), rna_cd4_seed)
fisher_test(diamond_methyl_cd4$module_genes, union(universe_module_methyl, universe_module_cd4), diamond_rna_cd4$module_genes)

fisher_test(methyl_cd8_seed, union(universe_methyl, universe_rna_cd8), rna_cd8_seed)
fisher_test(diamond_methyl_cd8$module_genes, union(universe_module_methyl, universe_module_cd8), diamond_rna_cd8$module_genes)

#library("ggvenn")
# p1 <- ggvenn(list(methyl=diamond_methyl_cd4$module_genes, rna=diamond_rna_cd4$module_genes), c("methyl","rna"),
#              show_percentage = F)
# p2 <- ggvenn(list(methyl=diamond_methyl_cd8$module_genes, rna=diamond_rna_cd8$module_genes), c("methyl","rna"),
#              show_percentage = F)
# 
# pdf("figures_manus/Overlap_modules_diamond.pdf")
# ggarrange(p1,p2, nrow=1, ncol=2)
# dev.off()

# Plot enrichment ####
df <- rbind(hyge_methyl_cd4[1,1:3], hyge_methyl_cd8[1,1:3], hyge_rna_cd4[1,1:3], hyge_rna_cd8[1,1:3],
            hyge_methyl_cd4_module[1,1:3], hyge_methyl_cd8_module[1,1:3], hyge_rna_cd4_module[1,1:3], hyge_rna_cd8_module[1,1:3])

df$module <- c(rep("Seed genes",4), rep("Module genes",4))
df$p_value <- -log10(df$p_value)
df$omic <- c("methylation", "methylation", "RNA-seq", "RNA-seq",
             "methylation", "methylation", "RNA-seq", "RNA-seq")
df$cell <- rep(c("CD4+", "CD8+"), 4)

significance_legend <- "*: p-value<0.05\n **: p-value<0.01\n ***: p-value<0.001 \n ***: p-value<2.2e-16"


pdf("figures_manus/Modules_diamond.pdf")
ggplot(df, aes(y=odds_ratio, x=1:8, fill = module)) + 
  geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") +
  geom_text(aes(label=c("*", "***", "**", "*", "****", "****", "****", "****")), position=position_dodge(width=0.9), vjust=-0.25) + 
  annotate("text", x = 2, y = 5, label = significance_legend) +
  theme_bw()  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x = "") + 
  scale_y_continuous(breaks=1:16) 
dev.off()

