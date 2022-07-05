# This script generates Figure 3 (Correlation bewteen 3rd-1st and PP-3rd)
# Created and modified by Alberto Zenere, 2021-14-11

# 06.Figure_3 ####
# 6.1 Setup 
# 6.2 Load Differential analysis results
# 6.3 Correlation plots

rm(list=ls()) # remove all entries in the global environment 

# 6.1 Setup ####
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

# Set directory structure 
main_dir <-"C:/Users/albze08/Desktop/phd/P4/methylation"
FOLDER_RDS <- "RDS_files/"
FOLDER_FIGURES <- "figures/"
setwd("C:/Users/albze08/Desktop/phd/P4/methylation")


# 6.2 Load Differential analysis results ####
#CD4 methylation
CD4_PP_3rd_HP_methyl <- readRDS("RDS_files/DMR/CD4_PP_HP_methyl_all.RDS")
CD4_PP_3rd_MS_methyl <- readRDS("RDS_files/DMR/CD4_PP_MS_methyl_all.RDS")

CD4_3rd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD4_3rd_HP_methyl_all.RDS")
CD4_3rd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD4_3rd_MS_methyl_all.RDS")

#CD8 methylation
CD8_PP_3rd_HP_methyl <- readRDS("RDS_files/DMR/CD8_PP_HP_methyl_all.RDS")
CD8_PP_3rd_MS_methyl <- readRDS("RDS_files/DMR/CD8_PP_MS_methyl_all.RDS")

CD8_3rd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD8_3rd_HP_methyl_all.RDS")
CD8_3rd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD8_3rd_MS_methyl_all.RDS")

#CD4 RNA-seq
CD4_PP_3rd_HP_rna <- readRDS("RDS_files/DMG/CD4_PP_HP_rna_all.RDS")
CD4_PP_3rd_MS_rna <- readRDS("RDS_files/DMG/CD4_PP_MS_rna_all.RDS"); CD4_PP_3rd_MS_rna <- CD4_PP_3rd_MS_rna[rownames(CD4_PP_3rd_HP_rna),]

CD4_3rd_1st_HP_rna <- readRDS("RDS_files/DMG/CD4_3rd_HP_rna_all.RDS"); CD4_3rd_1st_HP_rna <- CD4_3rd_1st_HP_rna[rownames(CD4_PP_3rd_HP_rna),]
CD4_3rd_1st_MS_rna <- readRDS("RDS_files/DMG/CD4_3rd_MS_rna_all.RDS"); CD4_3rd_1st_MS_rna <- CD4_3rd_1st_MS_rna[rownames(CD4_PP_3rd_HP_rna),]

#CD8 RNA-seq
CD8_PP_3rd_HP_rna <- readRDS("RDS_files/DMG/CD8_PP_HP_rna_all.RDS")
CD8_PP_3rd_MS_rna <- readRDS("RDS_files/DMG/CD8_PP_MS_rna_all.RDS"); CD8_PP_3rd_MS_rna <- CD8_PP_3rd_MS_rna[rownames(CD8_PP_3rd_HP_rna),]

CD8_3rd_1st_HP_rna <- readRDS("RDS_files/DMG/CD8_3rd_HP_rna_all.RDS"); CD8_3rd_1st_HP_rna <- CD8_3rd_1st_HP_rna[rownames(CD8_PP_3rd_HP_rna),]
CD8_3rd_1st_MS_rna <- readRDS("RDS_files/DMG/CD8_3rd_MS_rna_all.RDS"); CD8_3rd_1st_MS_rna <- CD8_3rd_1st_MS_rna[rownames(CD8_PP_3rd_HP_rna),]


# 6.3 Correlation plots ####
pdf("figures_manus/Correlations.pdf", width=10)

par(mfrow=c(2,4))

#HP CD4
smoothScatter(CD4_PP_3rd_HP_rna$logFC ~ CD4_3rd_1st_HP_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_PP_3rd_HP_rna$logFC, CD4_3rd_1st_HP_rna$logFC))/100))


smoothScatter(CD4_PP_3rd_HP_methyl$logFC ~ CD4_3rd_1st_HP_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_PP_3rd_HP_methyl$logFC, CD4_3rd_1st_HP_methyl$logFC))/100))

#HP CD8
smoothScatter(CD8_PP_3rd_HP_rna$logFC ~ CD8_3rd_1st_HP_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_PP_3rd_HP_rna$logFC, CD8_3rd_1st_HP_rna$logFC))/100))


smoothScatter(CD8_PP_3rd_HP_methyl$logFC ~ CD8_3rd_1st_HP_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_PP_3rd_HP_methyl$logFC, CD8_3rd_1st_HP_methyl$logFC))/100))




#MS CD4
smoothScatter(CD4_PP_3rd_MS_rna$logFC ~ CD4_3rd_1st_MS_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_PP_3rd_MS_rna$logFC, CD4_3rd_1st_MS_rna$logFC))/100))


smoothScatter(CD4_PP_3rd_MS_methyl$logFC ~ CD4_3rd_1st_MS_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_PP_3rd_MS_methyl$logFC, CD4_3rd_1st_MS_methyl$logFC))/100))

#MS CD8
smoothScatter(CD8_PP_3rd_MS_rna$logFC ~ CD8_3rd_1st_MS_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_PP_3rd_MS_rna$logFC, CD8_3rd_1st_MS_rna$logFC))/100))


smoothScatter(CD8_PP_3rd_MS_methyl$logFC ~ CD8_3rd_1st_MS_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_PP_3rd_MS_methyl$logFC, CD8_3rd_1st_MS_methyl$logFC))/100))


dev.off()




