
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

# Set directory structure ####
main_dir <-"C:/Users/albze08/Desktop/phd/P4/methylation"
FOLDER_RDS <- "RDS_files/"
FOLDER_FIGURES <- "figures/"
setwd("C:/Users/albze08/Desktop/phd/P4/methylation")

# My functions ####
source('probe_to_gene.R')
source('fisher_test.R')

# Load ####
#CD4 methylation
CD4_2nd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD4_2nd_1st_HP_all.RDS")
CD4_2nd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD4_2nd_1st_MS_all.RDS")

CD4_3rd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD4_3rd_HP_methyl_all.RDS")
CD4_3rd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD4_3rd_MS_methyl_all.RDS")

CD4_3rd_2nd_HP_methyl <- readRDS("RDS_files/DMR/CD4_3rd_2nd_HP_all.RDS")
CD4_3rd_2nd_MS_methyl <- readRDS("RDS_files/DMR/CD4_3rd_2nd_MS_all.RDS")

CD4_PP_2nd_HP_methyl <- readRDS("RDS_files/DMR/CD4_PP_2nd_HP_all.RDS")
CD4_PP_2nd_MS_methyl <- readRDS("RDS_files/DMR/CD4_PP_2nd_MS_all.RDS")

#CD8 methylation
CD8_2nd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD8_2nd_1st_HP_all.RDS")
CD8_2nd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD8_2nd_1st_MS_all.RDS")

CD8_3rd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD8_3rd_HP_methyl_all.RDS")
CD8_3rd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD8_3rd_MS_methyl_all.RDS")

CD8_3rd_2nd_HP_methyl <- readRDS("RDS_files/DMR/CD8_3rd_2nd_HP_all.RDS")
CD8_3rd_2nd_MS_methyl <- readRDS("RDS_files/DMR/CD8_3rd_2nd_MS_all.RDS")

CD8_PP_2nd_HP_methyl <- readRDS("RDS_files/DMR/CD8_PP_2nd_HP_all.RDS")
CD8_PP_2nd_MS_methyl <- readRDS("RDS_files/DMR/CD8_PP_2nd_MS_all.RDS")

#CD4 RNA-seq
CD4_2nd_1st_HP_rna <- readRDS("RDS_files/DMG/HP_2nd_1st_cd4_all.RDS")
CD4_2nd_1st_MS_rna <- readRDS("RDS_files/DMG/MS_2nd_1st_cd4_all.RDS"); CD4_2nd_1st_MS_rna <- CD4_2nd_1st_MS_rna[rownames(CD4_2nd_1st_HP_rna),]

CD4_3rd_1st_HP_rna <- readRDS("RDS_files/DMG/CD4_3rd_HP_rna_all.RDS"); CD4_3rd_1st_HP_rna <- CD4_3rd_1st_HP_rna[rownames(CD4_2nd_1st_HP_rna),]
CD4_3rd_1st_MS_rna <- readRDS("RDS_files/DMG/CD4_3rd_MS_rna_all.RDS"); CD4_3rd_1st_MS_rna <- CD4_3rd_1st_MS_rna[rownames(CD4_2nd_1st_HP_rna),]

CD4_3rd_2nd_HP_rna <- readRDS("RDS_files/DMG/HP_3rd_2nd_cd4_all.RDS"); CD4_3rd_2nd_HP_rna <- CD4_3rd_2nd_HP_rna[rownames(CD4_2nd_1st_HP_rna),]
CD4_3rd_2nd_MS_rna <- readRDS("RDS_files/DMG/MS_3rd_2nd_cd4_all.RDS"); CD4_3rd_2nd_MS_rna <- CD4_3rd_2nd_MS_rna[rownames(CD4_2nd_1st_HP_rna),]

CD4_PP_2nd_HP_rna <- readRDS("RDS_files/DMG/HP_PP_2nd_cd4_all.RDS"); CD4_PP_2nd_HP_rna <- CD4_PP_2nd_HP_rna[rownames(CD4_2nd_1st_HP_rna),]
CD4_PP_2nd_MS_rna <- readRDS("RDS_files/DMG/MS_PP_2nd_cd4_all.RDS"); CD4_PP_2nd_MS_rna <- CD4_PP_2nd_MS_rna[rownames(CD4_2nd_1st_HP_rna),]

#CD8 RNA-seq
CD8_2nd_1st_HP_rna <- readRDS("RDS_files/DMG/HP_2nd_1st_cd8_all.RDS")
CD8_2nd_1st_MS_rna <- readRDS("RDS_files/DMG/MS_2nd_1st_cd8_all.RDS"); CD8_2nd_1st_MS_rna <- CD8_2nd_1st_MS_rna[rownames(CD8_2nd_1st_HP_rna),]

CD8_3rd_1st_HP_rna <- readRDS("RDS_files/DMG/CD8_3rd_HP_rna_all.RDS"); CD8_3rd_1st_HP_rna <- CD8_3rd_1st_HP_rna[rownames(CD8_2nd_1st_HP_rna),]
CD8_3rd_1st_MS_rna <- readRDS("RDS_files/DMG/CD8_3rd_MS_rna_all.RDS"); CD8_3rd_1st_MS_rna <- CD8_3rd_1st_MS_rna[rownames(CD8_2nd_1st_HP_rna),]

CD8_3rd_2nd_HP_rna <- readRDS("RDS_files/DMG/HP_3rd_2nd_cd8_all.RDS"); CD8_3rd_2nd_HP_rna <- CD8_3rd_2nd_HP_rna[rownames(CD8_2nd_1st_HP_rna),]
CD8_3rd_2nd_MS_rna <- readRDS("RDS_files/DMG/MS_3rd_2nd_cd8_all.RDS"); CD8_3rd_2nd_MS_rna <- CD8_3rd_2nd_MS_rna[rownames(CD8_2nd_1st_HP_rna),]

CD8_PP_2nd_HP_rna <- readRDS("RDS_files/DMG/HP_PP_2nd_cd8_all.RDS"); CD8_PP_2nd_HP_rna <- CD8_PP_2nd_HP_rna[rownames(CD8_2nd_1st_HP_rna),]
CD8_PP_2nd_MS_rna <- readRDS("RDS_files/DMG/MS_PP_2nd_cd8_all.RDS"); CD8_PP_2nd_MS_rna <- CD8_PP_2nd_MS_rna[rownames(CD8_2nd_1st_HP_rna),]

#Universe ####
list_cd4_methyl <- rownames(CD4_2nd_1st_HP_methyl)
list_cd8_methyl <- rownames(CD8_2nd_1st_HP_methyl)

list_cd4_rna <- rownames(CD4_2nd_1st_HP_rna)
list_cd8_rna <- rownames(CD8_2nd_1st_HP_rna)


# Correlation plots  2nd-1st vs 3rd-1st ####
pdf("figures_manus/Correlations_Suppl.pdf", width=10)

par(mfrow=c(2,4))

#HP CD4
smoothScatter(CD4_3rd_1st_HP_rna$logFC ~ CD4_2nd_1st_HP_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_3rd_1st_HP_rna$logFC, CD4_2nd_1st_HP_rna$logFC))/100))


smoothScatter(CD4_3rd_1st_HP_methyl$logFC ~ CD4_2nd_1st_HP_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_3rd_1st_HP_methyl$logFC, CD4_2nd_1st_HP_methyl$logFC))/100))

# HP CD8
smoothScatter(CD8_3rd_1st_HP_rna$logFC ~ CD8_2nd_1st_HP_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_3rd_1st_HP_rna$logFC, CD8_2nd_1st_HP_rna$logFC))/100))


smoothScatter(CD8_3rd_1st_HP_methyl$logFC ~ CD8_2nd_1st_HP_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_3rd_1st_HP_methyl$logFC, CD8_2nd_1st_HP_methyl$logFC))/100))




#MS CD4
smoothScatter(CD4_3rd_1st_MS_rna$logFC ~ CD4_2nd_1st_MS_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_3rd_1st_MS_rna$logFC, CD4_2nd_1st_MS_rna$logFC))/100))

smoothScatter(CD4_3rd_1st_MS_methyl$logFC ~ CD4_2nd_1st_MS_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_3rd_1st_MS_methyl$logFC, CD4_2nd_1st_MS_methyl$logFC))/100))


#MS CD8
smoothScatter(CD8_3rd_1st_MS_rna$logFC ~ CD8_2nd_1st_MS_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_3rd_1st_MS_rna$logFC, CD8_2nd_1st_MS_rna$logFC))/100))


smoothScatter(CD8_3rd_1st_MS_methyl$logFC ~ CD8_2nd_1st_MS_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_3rd_1st_MS_methyl$logFC, CD8_2nd_1st_MS_methyl$logFC))/100))



dev.off()


# Correlation plots  3rd-2nd vs 3rd-1st ####
pdf("figures_manus/Correlations_Suppl_2.pdf", width=12)

par(mfrow=c(2,4))

#HP CD4
smoothScatter(CD4_3rd_1st_HP_rna$logFC ~ CD4_3rd_2nd_HP_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_3rd_1st_HP_rna$logFC, CD4_3rd_2nd_HP_rna$logFC))/100))


smoothScatter(CD4_3rd_1st_HP_methyl$logFC ~ CD4_3rd_2nd_HP_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_3rd_1st_HP_methyl$logFC, CD4_3rd_2nd_HP_methyl$logFC))/100))

# HP CD8
smoothScatter(CD8_3rd_1st_HP_rna$logFC ~ CD8_3rd_2nd_HP_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_3rd_1st_HP_rna$logFC, CD8_3rd_2nd_HP_rna$logFC))/100))


smoothScatter(CD8_3rd_1st_HP_methyl$logFC ~ CD8_3rd_2nd_HP_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_3rd_1st_HP_methyl$logFC, CD8_3rd_2nd_HP_methyl$logFC))/100))




#MS CD4
smoothScatter(CD4_3rd_1st_MS_rna$logFC ~ CD4_3rd_2nd_MS_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_3rd_1st_MS_rna$logFC, CD4_3rd_2nd_MS_rna$logFC))/100))

smoothScatter(CD4_3rd_1st_MS_methyl$logFC ~ CD4_3rd_2nd_MS_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_3rd_1st_MS_methyl$logFC, CD4_3rd_2nd_MS_methyl$logFC))/100))


#MS CD8
smoothScatter(CD8_3rd_1st_MS_rna$logFC ~ CD8_3rd_2nd_MS_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_3rd_1st_MS_rna$logFC, CD8_3rd_2nd_MS_rna$logFC))/100))


smoothScatter(CD8_3rd_1st_MS_methyl$logFC ~ CD8_3rd_2nd_MS_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_3rd_1st_MS_methyl$logFC, CD8_3rd_2nd_MS_methyl$logFC))/100))



dev.off()



# Correlation plots  PP-2nd vs 2nd-1st ####
pdf("figures_manus/Correlations_Suppl_3.pdf", width=12)

par(mfrow=c(2,4))

#HP CD4
smoothScatter(CD4_PP_2nd_HP_rna$logFC ~ CD4_2nd_1st_HP_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_PP_2nd_HP_rna$logFC, CD4_2nd_1st_HP_rna$logFC))/100))


smoothScatter(CD4_PP_2nd_HP_methyl$logFC ~ CD4_2nd_1st_HP_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_PP_2nd_HP_methyl$logFC, CD4_2nd_1st_HP_methyl$logFC))/100))

# HP CD8
smoothScatter(CD8_PP_2nd_HP_rna$logFC ~ CD8_2nd_1st_HP_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_PP_2nd_HP_rna$logFC, CD8_2nd_1st_HP_rna$logFC))/100))


smoothScatter(CD8_PP_2nd_HP_methyl$logFC ~ CD8_2nd_1st_HP_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_PP_2nd_HP_methyl$logFC, CD8_2nd_1st_HP_methyl$logFC))/100))




#MS CD4
smoothScatter(CD4_PP_2nd_MS_rna$logFC ~ CD4_2nd_1st_MS_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "PP-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_PP_2nd_MS_rna$logFC, CD4_2nd_1st_MS_rna$logFC))/100))

smoothScatter(CD4_PP_2nd_MS_methyl$logFC ~ CD4_2nd_1st_MS_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "PP-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_PP_2nd_MS_methyl$logFC, CD4_2nd_1st_MS_methyl$logFC))/100))


#MS CD8
smoothScatter(CD8_PP_2nd_MS_rna$logFC ~ CD8_2nd_1st_MS_rna$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_PP_2nd_MS_rna$logFC, CD8_2nd_1st_MS_rna$logFC))/100))


smoothScatter(CD8_PP_2nd_MS_methyl$logFC ~ CD8_2nd_1st_MS_methyl$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple sclerosis",
              ylab = "PP-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_PP_2nd_MS_methyl$logFC, CD8_2nd_1st_MS_methyl$logFC))/100))



dev.off()


#Second figure: DMPs and DEGs ####

#CD4 methylation
CD4_2nd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD4_2nd_1st_HP.RDS")
CD4_2nd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD4_2nd_1st_MS.RDS")

CD4_3rd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD4_3rd_HP_methyl.RDS")
CD4_3rd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD4_3rd_MS_methyl.RDS")

CD4_3rd_2nd_HP_methyl <- readRDS("RDS_files/DMR/CD4_3rd_2nd_HP.RDS")
CD4_3rd_2nd_MS_methyl <- readRDS("RDS_files/DMR/CD4_3rd_2nd_MS.RDS")

CD4_PP_3rd_HP_methyl <- readRDS("RDS_files/DMR/CD4_PP_HP_methyl.RDS")
CD4_PP_3rd_MS_methyl <- readRDS("RDS_files/DMR/CD4_PP_MS_methyl.RDS")

CD4_PP_2nd_HP_methyl <- readRDS("RDS_files/DMR/CD4_PP_2nd_HP.RDS")
CD4_PP_2nd_MS_methyl <- readRDS("RDS_files/DMR/CD4_PP_2nd_MS.RDS")

#CD8 methylation
CD8_2nd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD8_2nd_1st_HP.RDS")
CD8_2nd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD8_2nd_1st_MS.RDS")

CD8_3rd_1st_HP_methyl <- readRDS("RDS_files/DMR/CD8_3rd_HP_methyl.RDS")
CD8_3rd_1st_MS_methyl <- readRDS("RDS_files/DMR/CD8_3rd_MS_methyl.RDS")

CD8_3rd_2nd_HP_methyl <- readRDS("RDS_files/DMR/CD8_3rd_2nd_HP.RDS")
CD8_3rd_2nd_MS_methyl <- readRDS("RDS_files/DMR/CD8_3rd_2nd_MS.RDS")

CD8_PP_3rd_HP_methyl <- readRDS("RDS_files/DMR/CD8_PP_HP_methyl.RDS")
CD8_PP_3rd_MS_methyl <- readRDS("RDS_files/DMR/CD8_PP_MS_methyl.RDS")

CD8_PP_2nd_HP_methyl <- readRDS("RDS_files/DMR/CD8_PP_2nd_HP.RDS")
CD8_PP_2nd_MS_methyl <- readRDS("RDS_files/DMR/CD8_PP_2nd_MS.RDS")

#CD4 RNA-seq
CD4_2nd_1st_HP_rna <- readRDS("RDS_files/DMG/CD4_2nd_HP_rna.RDS")
CD4_2nd_1st_MS_rna <- readRDS("RDS_files/DMG/CD4_2nd_MS_rna.RDS")

CD4_3rd_1st_HP_rna <- readRDS("RDS_files/DMG/CD4_3rd_HP_rna.RDS")
CD4_3rd_1st_MS_rna <- readRDS("RDS_files/DMG/CD4_3rd_MS_rna.RDS")

CD4_3rd_2nd_HP_rna <- readRDS("RDS_files/DMG/CD4_3rd_2nd_HP_rna.RDS")
CD4_3rd_2nd_MS_rna <- readRDS("RDS_files/DMG/CD4_3rd_2nd_MS_rna.RDS")

CD4_PP_3rd_HP_rna <- readRDS("RDS_files/DMG/CD4_PP_HP_rna.RDS")
CD4_PP_3rd_MS_rna <- readRDS("RDS_files/DMG/CD4_PP_MS_rna.RDS")

CD4_PP_2nd_HP_rna <- readRDS("RDS_files/DMG/CD4_PP_2nd_HP_rna.RDS")
CD4_PP_2nd_MS_rna <- readRDS("RDS_files/DMG/CD4_PP_2nd_MS_rna.RDS")

#CD8 RNA-seq
CD8_2nd_1st_HP_rna <- readRDS("RDS_files/DMG/CD8_2nd_HP_rna.RDS")
CD8_2nd_1st_MS_rna <- readRDS("RDS_files/DMG/CD8_2nd_MS_rna.RDS")

CD8_3rd_1st_HP_rna <- readRDS("RDS_files/DMG/CD8_3rd_HP_rna.RDS")
CD8_3rd_1st_MS_rna <- readRDS("RDS_files/DMG/CD8_3rd_MS_rna.RDS")

CD8_3rd_2nd_HP_rna <- readRDS("RDS_files/DMG/CD8_3rd_2nd_HP_rna.RDS")
CD8_3rd_2nd_MS_rna <- readRDS("RDS_files/DMG/CD8_3rd_2nd_MS_rna.RDS")

CD8_PP_3rd_HP_rna <- readRDS("RDS_files/DMG/CD8_PP_HP_rna.RDS")
CD8_PP_3rd_MS_rna <- readRDS("RDS_files/DMG/CD8_PP_MS_rna.RDS")

CD8_PP_2nd_HP_rna <- readRDS("RDS_files/DMG/CD8_PP_2nd_HP_rna.RDS")
CD8_PP_2nd_MS_rna <- readRDS("RDS_files/DMG/CD8_PP_2nd_MS_rna.RDS")

#Plot Methylation
df_methyl <- data.frame(n_DMP=c(nrow(CD4_2nd_1st_HP_methyl), nrow(CD4_2nd_1st_MS_methyl), nrow(CD8_2nd_1st_MS_methyl),
                                nrow(CD4_3rd_1st_HP_methyl), nrow(CD4_3rd_1st_MS_methyl), nrow(CD8_3rd_1st_MS_methyl),                
                                nrow(CD4_3rd_2nd_HP_methyl), nrow(CD4_3rd_2nd_MS_methyl), nrow(CD8_3rd_2nd_MS_methyl),                
                                nrow(CD4_PP_3rd_HP_methyl), nrow(CD4_PP_3rd_MS_methyl), nrow(CD8_PP_3rd_MS_methyl),                
                                nrow(CD4_PP_2nd_HP_methyl), nrow(CD4_PP_2nd_MS_methyl), nrow(CD8_PP_2nd_MS_methyl)))

df_methyl$Disease <- rep(c("HC", "MS", "MS"),5)
df_methyl$Comparison <- c( rep("2nd-1st", 3), rep("3rd-1st", 3), rep("3rd-2nd", 3), rep("PP-3rd", 3), rep("PP-2nd", 3) )
df_methyl$Comparison <- factor(df_methyl$Comparison, levels=c("2nd-1st","3rd-2nd","3rd-1st","PP-3rd", "PP-2nd"))
df_methyl$Cell <- rep(c("CD4+","CD4+", "CD8+"),5)

p1 <- ggplot(df_methyl, aes(y=n_DMP, x=interaction(Cell,Disease), fill=Comparison)) + geom_bar(stat="identity", position = "dodge") +
  theme_bw()  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


#Plot CD8 HC 
df_methyl_cd8_hp <- data.frame(n_DMP=c(nrow(CD8_2nd_1st_HP_methyl),
                                nrow(CD8_3rd_1st_HP_methyl),           
                                nrow(CD8_3rd_2nd_HP_methyl),
                                nrow(CD8_PP_3rd_HP_methyl),
                                nrow(CD8_PP_2nd_HP_methyl)))
df_methyl_cd8_hp$n_DMP <- log10(df_methyl_cd8_hp$n_DMP)
df_methyl_cd8_hp$Disease <- rep("HC",5)
df_methyl_cd8_hp$Comparison <- c( "2nd-1st", "3rd-1st", "3rd-2nd", "PP-3rd", "PP-2nd")
df_methyl_cd8_hp$Comparison <- factor(df_methyl_cd8_hp$Comparison, levels=c("2nd-1st","3rd-2nd","3rd-1st","PP-3rd","PP-2nd"))
df_methyl_cd8_hp$Cell <- rep("CD8+",5)

p2 <- ggplot(df_methyl_cd8_hp, aes(y=n_DMP, x=Comparison, fill=Comparison)) + geom_bar(stat="identity", position = "dodge") +
  theme_bw()  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

#Plot RNA-seq
df_rna <- data.frame(n_DEG=c(nrow(CD4_2nd_1st_HP_rna), nrow(CD4_2nd_1st_MS_rna), nrow(CD8_2nd_1st_HP_rna), nrow(CD8_2nd_1st_MS_rna),
                             nrow(CD4_3rd_1st_HP_rna), nrow(CD4_3rd_1st_MS_rna), nrow(CD8_3rd_1st_HP_rna), nrow(CD8_3rd_1st_MS_rna),                
                             nrow(CD4_3rd_2nd_HP_rna), nrow(CD4_3rd_2nd_MS_rna), nrow(CD8_3rd_2nd_HP_rna), nrow(CD8_3rd_2nd_MS_rna),
                             nrow(CD4_PP_3rd_HP_rna), nrow(CD4_PP_3rd_MS_rna), nrow(CD8_PP_3rd_HP_rna), nrow(CD8_PP_3rd_MS_rna),
                             nrow(CD4_PP_2nd_HP_rna), nrow(CD4_PP_2nd_MS_rna), nrow(CD8_PP_2nd_HP_rna), nrow(CD8_PP_2nd_MS_rna)))

df_rna$Disease <- rep(c("HC", "MS"),10)
df_rna$Comparison <- c( rep("2nd-1st", 4), rep("3rd-1st", 4), rep("3rd-2nd", 4), rep("PP-3rd", 4), rep("PP-2nd", 4) )
df_rna$Comparison <- factor(df_rna$Comparison, levels=c("2nd-1st","3rd-2nd","3rd-1st","PP-3rd","PP-2nd"))
df_rna$Cell <- rep(c(rep("CD4+",2), rep("CD8+",2)), 5)

p3 <- ggplot(df_rna, aes(y=n_DEG, x=interaction(Cell,Disease), fill=Comparison)) + geom_bar(stat="identity", position = "dodge") +
  theme_bw()  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

pdf("figures_manus/n_DEG_DMP.pdf", width=12)
ggarrange(p1,p2,p3, nrow=1, ncol=3)
dev.off()


# Plot overlap ####
#CD4 HC RNA
overlap_2nd_1st_rna_hp_cd4 <- fisher_test(rownames(CD4_2nd_1st_HP_rna), list_cd4_rna, rownames(CD4_PP_2nd_HP_rna))
overlap_3rd_2nd_rna_hp_cd4 <- fisher_test(rownames(CD4_3rd_2nd_HP_rna), list_cd4_rna, rownames(CD4_PP_3rd_HP_rna))
overlap_3rd_1st_rna_hp_cd4 <- fisher_test(rownames(CD4_3rd_1st_HP_rna), list_cd4_rna, rownames(CD4_PP_3rd_HP_rna))

#CD4 MS RNA
overlap_2nd_1st_rna_ms_cd4 <- fisher_test(rownames(CD4_2nd_1st_MS_rna), list_cd4_rna, rownames(CD4_PP_2nd_MS_rna))
overlap_3rd_2nd_rna_ms_cd4 <- fisher_test(rownames(CD4_3rd_2nd_MS_rna), list_cd4_rna, rownames(CD4_PP_3rd_MS_rna))
overlap_3rd_1st_rna_ms_cd4 <- fisher_test(rownames(CD4_3rd_1st_MS_rna), list_cd4_rna, rownames(CD4_PP_3rd_MS_rna))

#CD8 HC RNA
overlap_2nd_1st_rna_hp_cd8 <- fisher_test(rownames(CD8_2nd_1st_HP_rna), list_cd8_rna, rownames(CD8_PP_2nd_HP_rna))
overlap_3rd_2nd_rna_hp_cd8 <- fisher_test(rownames(CD8_3rd_2nd_HP_rna), list_cd8_rna, rownames(CD8_PP_3rd_HP_rna))
overlap_3rd_1st_rna_hp_cd8 <- fisher_test(rownames(CD8_3rd_1st_HP_rna), list_cd8_rna, rownames(CD8_PP_3rd_HP_rna))

#CD8 MS RNA
overlap_2nd_1st_rna_ms_cd8 <- fisher_test(rownames(CD8_2nd_1st_MS_rna), list_cd8_rna, rownames(CD8_PP_2nd_MS_rna))
overlap_3rd_2nd_rna_ms_cd8 <- fisher_test(rownames(CD8_3rd_2nd_MS_rna), list_cd8_rna, rownames(CD8_PP_3rd_MS_rna))
overlap_3rd_1st_rna_ms_cd8 <- fisher_test(rownames(CD8_3rd_1st_MS_rna), list_cd8_rna, rownames(CD8_PP_3rd_MS_rna))


#CD4 HC methyl
overlap_2nd_1st_methyl_hp_cd4 <- fisher_test(rownames(CD4_2nd_1st_HP_methyl), list_cd4_methyl, rownames(CD4_PP_2nd_HP_methyl))
overlap_3rd_2nd_methyl_hp_cd4 <- fisher_test(rownames(CD4_3rd_2nd_HP_methyl), list_cd4_methyl, rownames(CD4_PP_3rd_HP_methyl))
overlap_3rd_1st_methyl_hp_cd4 <- fisher_test(rownames(CD4_3rd_1st_HP_methyl), list_cd4_methyl, rownames(CD4_PP_3rd_HP_methyl))

#CD4 MS methyl
overlap_2nd_1st_methyl_ms_cd4 <- fisher_test(rownames(CD4_2nd_1st_MS_methyl), list_cd4_methyl, rownames(CD4_PP_2nd_MS_methyl))
overlap_3rd_2nd_methyl_ms_cd4 <- fisher_test(rownames(CD4_3rd_2nd_MS_methyl), list_cd4_methyl, rownames(CD4_PP_3rd_MS_methyl))
overlap_3rd_1st_methyl_ms_cd4 <- fisher_test(rownames(CD4_3rd_1st_MS_methyl), list_cd4_methyl, rownames(CD4_PP_3rd_MS_methyl))

#CD8 HC methyl
overlap_2nd_1st_methyl_hp_cd8 <- fisher_test(rownames(CD8_2nd_1st_HP_methyl), list_cd8_methyl, rownames(CD8_PP_2nd_HP_methyl))
overlap_3rd_2nd_methyl_hp_cd8 <- fisher_test(rownames(CD8_3rd_2nd_HP_methyl), list_cd8_methyl, rownames(CD8_PP_3rd_HP_methyl))
overlap_3rd_1st_methyl_hp_cd8 <- fisher_test(rownames(CD8_3rd_1st_HP_methyl), list_cd8_methyl, rownames(CD8_PP_3rd_HP_methyl))

#CD8 MS methyl
overlap_2nd_1st_methyl_ms_cd8 <- fisher_test(rownames(CD8_2nd_1st_MS_methyl), list_cd8_methyl, rownames(CD8_PP_2nd_MS_methyl))
overlap_3rd_2nd_methyl_ms_cd8 <- fisher_test(rownames(CD8_3rd_2nd_MS_methyl), list_cd8_methyl, rownames(CD8_PP_3rd_MS_methyl))
overlap_3rd_1st_methyl_ms_cd8 <- fisher_test(rownames(CD8_3rd_1st_MS_methyl), list_cd8_methyl, rownames(CD8_PP_3rd_MS_methyl))

#RNA-seq
df_overlap_rna <- data.frame(overlap=c(overlap_2nd_1st_rna_hp_cd4$n_genes, overlap_3rd_2nd_rna_hp_cd4$n_genes, overlap_3rd_1st_rna_hp_cd4$n_genes,
                           overlap_2nd_1st_rna_ms_cd4$n_genes, overlap_3rd_2nd_rna_ms_cd4$n_genes, overlap_3rd_1st_rna_ms_cd4$n_genes,
                           overlap_2nd_1st_rna_hp_cd8$n_genes, overlap_3rd_2nd_rna_hp_cd8$n_genes, overlap_3rd_1st_rna_hp_cd8$n_genes,
                           overlap_2nd_1st_rna_ms_cd8$n_genes, overlap_3rd_2nd_rna_ms_cd8$n_genes, overlap_3rd_1st_rna_ms_cd8$n_genes))
df_overlap_rna$Disease <- c(rep("HC",3), rep("MS",3), rep("HC",3), rep("MS",3))
df_overlap_rna$Comparison <- rep( c("2nd-1st", "3rd-2nd", "3rd-1st"), 4 )
df_overlap_rna$Comparison <- factor(df_overlap_rna$Comparison, levels=c("2nd-1st","3rd-2nd","3rd-1st"))
df_overlap_rna$Cell <- c(rep("CD4+",6), rep("CD8+",6))

p1 <- ggplot(df_overlap_rna, aes(y=overlap, x=interaction(Cell,Disease), fill=Comparison)) + geom_bar(stat="identity", position = "dodge") +
  theme_bw()  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 



#methylation
df_overlap_methyl <- data.frame(overlap=c(overlap_2nd_1st_methyl_hp_cd4$n_genes, overlap_3rd_2nd_methyl_hp_cd4$n_genes, overlap_3rd_1st_methyl_hp_cd4$n_genes,
                                       overlap_2nd_1st_methyl_ms_cd4$n_genes, overlap_3rd_2nd_methyl_ms_cd4$n_genes, overlap_3rd_1st_methyl_ms_cd4$n_genes,
                                       overlap_2nd_1st_methyl_ms_cd8$n_genes, overlap_3rd_2nd_methyl_ms_cd8$n_genes, overlap_3rd_1st_methyl_ms_cd8$n_genes))
df_overlap_methyl$Disease <- c(rep("HC",3), rep("MS",3), rep("MS",3))
df_overlap_methyl$Comparison <- rep( c("2nd-1st", "3rd-2nd", "3rd-1st"), 3 )
df_overlap_methyl$Comparison <- factor(df_overlap_methyl$Comparison, levels=c("2nd-1st","3rd-2nd","3rd-1st"))
df_overlap_methyl$Cell <- c(rep("CD4+",6), rep("CD8+",3))

p2 <- ggplot(df_overlap_methyl, aes(y=overlap, x=interaction(Cell,Disease), fill=Comparison)) + geom_bar(stat="identity", position = "dodge") +
  theme_bw()  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 



#methylation CD8 HC
df_overlap_methyl_cd8_hp <- data.frame(overlap=log10(c(overlap_2nd_1st_methyl_hp_cd8$n_genes, overlap_3rd_2nd_methyl_hp_cd8$n_genes, overlap_3rd_1st_methyl_hp_cd8$n_genes)))
df_overlap_methyl_cd8_hp$Disease <- rep("HC",3)
df_overlap_methyl_cd8_hp$Comparison <- c("2nd-1st", "3rd-2nd", "3rd-1st")
df_overlap_methyl_cd8_hp$Comparison <- factor(df_overlap_methyl_cd8_hp$Comparison, levels=c("2nd-1st","3rd-2nd","3rd-1st"))
df_overlap_methyl_cd8_hp$Cell <- rep("CD8+",3)

p3 <- ggplot(df_overlap_methyl_cd8_hp, aes(y=overlap, x=interaction(Cell,Disease), fill=Comparison)) + geom_bar(stat="identity", position = "dodge") +
  theme_bw()  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 


pdf("figures_manus/n_overlap.pdf", width=12)
ggarrange(p1,p2,p3, nrow=1, ncol=3)
dev.off()


