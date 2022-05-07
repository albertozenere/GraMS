
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

# Set directory structure ####
main_dir <-"C:/Users/albze08/Desktop/phd/P4/methylation/"
FOLDER_RDS <- "RDS_files/"
FOLDER_FIGURES <- "figures/"
setwd("C:/Users/albze08/Desktop/phd/P4/methylation/")


par(mar=c(1,1,1,1))
dev.off()

# My functions ####

# Load RNA-seq ####
rna_folder <- "C:/Users/albze08/Desktop/phd/P4/RNAseq/RDS/"

# CD4 ####
count_cd4 <- readRDS(file = paste0(rna_folder, "count_cd4_tmm.RDS"))

metadata_cd4 <- readRDS(file = paste0(rna_folder, "metadata_cd4.RDS"))
metadata_cd4$Sample_Group <- gsub(" ", "", metadata_cd4$Sample_Group)

# Remove non pregnant
is_P <- (!metadata_cd4$TimePoint %in% "Nonpregnant") %>% which()
count_cd4 <- count_cd4[, is_P]
metadata_cd4 <- metadata_cd4[is_P, ]

# Remove Activated
idx <- metadata_cd4$State=="Resting"
count_cd4 <- count_cd4[,idx]
metadata_cd4 <- metadata_cd4[idx,]

# CD8 ####
count_cd8 <- readRDS(file = paste0(rna_folder, "count_cd8_tmm.RDS"))

metadata_cd8 <- readRDS(file = paste0(rna_folder, "metadata_cd8.RDS"))
metadata_cd8$Sample_Group <- gsub(" ", "", metadata_cd8$Sample_Group)

# Remove non pregnant
is_P <- (!metadata_cd8$TimePoint %in% "Nonpregnant") %>% which()
count_cd8 <- count_cd8[, is_P]
metadata_cd8 <- metadata_cd8[is_P, ]

# Remove Activated
idx <- metadata_cd8$State=="Resting"
count_cd8 <- count_cd8[,idx]
metadata_cd8 <- metadata_cd8[idx,]

# Limma on CD4 ####
metadata_cd4$Sample_Group <- gsub("MS_BP", "AA", metadata_cd4$Sample_Group) 
design <- model.matrix(~ Sample_Group + Proportions_Memory + Cell_Viability + Age, metadata_cd4 )
metadata_cd4$Sample_Group <- gsub("AA", "MS_BP", metadata_cd4$Sample_Group)

count_voom <- voom(count_cd4, design, plot=F)
corfit <- duplicateCorrelation(count_voom, design, block=metadata_cd4$Individual)

contr_matrix <- makeContrasts(First = Sample_GroupMS_1st_CD4 - Sample_GroupHP_1st_CD4 ,
                              Second = Sample_GroupMS_2nd_CD4 - Sample_GroupHP_2nd_CD4,
                              Third = Sample_GroupMS_3rd_CD4 - Sample_GroupHP_3rd_CD4,
                              PP = Sample_GroupMS_PP_CD4 - Sample_GroupHP_PP_CD4,
                              Second_First = (Sample_GroupMS_2nd_CD4 - Sample_GroupMS_1st_CD4) - (Sample_GroupHP_2nd_CD4 - Sample_GroupHP_1st_CD4),
                              Third_Second = (Sample_GroupMS_3rd_CD4 - Sample_GroupMS_2nd_CD4) - (Sample_GroupHP_3rd_CD4 - Sample_GroupHP_2nd_CD4),
                              Third_First = (Sample_GroupMS_3rd_CD4 - Sample_GroupMS_1st_CD4) - (Sample_GroupHP_3rd_CD4 - Sample_GroupHP_1st_CD4),
                              PP_Third = (Sample_GroupMS_PP_CD4 - Sample_GroupMS_3rd_CD4) - (Sample_GroupHP_PP_CD4 - Sample_GroupHP_3rd_CD4),
                              PP_First = (Sample_GroupMS_PP_CD4 - Sample_GroupMS_1st_CD4) - (Sample_GroupHP_PP_CD4 - Sample_GroupHP_1st_CD4),
                              Second_MS = Sample_GroupMS_2nd_CD4 - Sample_GroupMS_1st_CD4,
                              Third_MS = Sample_GroupMS_3rd_CD4 - Sample_GroupMS_1st_CD4,
                              PP_MS = Sample_GroupMS_PP_CD4 - Sample_GroupMS_1st_CD4,
                              Second_HP = Sample_GroupHP_2nd_CD4 - Sample_GroupHP_1st_CD4,
                              Third_HP = Sample_GroupHP_3rd_CD4 - Sample_GroupHP_1st_CD4,
                              PP_HP = Sample_GroupHP_PP_CD4 - Sample_GroupHP_1st_CD4,
                              Third_Second_MS = Sample_GroupMS_3rd_CD4 - Sample_GroupMS_2nd_CD4,
                              PP_Third_MS = Sample_GroupMS_PP_CD4 - Sample_GroupMS_3rd_CD4,
                              PP_Second_MS = Sample_GroupMS_PP_CD4 - Sample_GroupMS_2nd_CD4,
                              Third_Second_HP = Sample_GroupHP_3rd_CD4 - Sample_GroupHP_2nd_CD4,
                              PP_Third_HP = Sample_GroupHP_PP_CD4 - Sample_GroupHP_3rd_CD4,
                              PP_Second_HP = Sample_GroupHP_PP_CD4 - Sample_GroupHP_2nd_CD4,
                              levels=colnames(design)
)
fit_CD4 <- lmFit(count_voom , design, block=metadata_cd4$Individual, correlation=corfit$consensus)
fit_CD4 <- contrasts.fit(fit_CD4, contrasts=contr_matrix)
fit_CD4 <- eBayes(fit_CD4)
summary(decideTests(fit_CD4))


# Gather Changes
#HP
CD4_HP_2nd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Second_HP")
CD4_HP_2nd_1st <- CD4_HP_2nd_1st[rownames(count_cd4),]

CD4_HP_3rd_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second_HP")
CD4_HP_3rd_2nd <- CD4_HP_3rd_2nd[rownames(count_cd4),]

CD4_HP_3rd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_HP")
CD4_HP_3rd_1st <- CD4_HP_3rd_1st[rownames(count_cd4),]

CD4_HP_PP_3rd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third_HP")
CD4_HP_PP_3rd <- CD4_HP_PP_3rd[rownames(count_cd4),]

#MS
CD4_MS_2nd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Second_MS")
CD4_MS_2nd_1st <- CD4_MS_2nd_1st[rownames(count_cd4),]

CD4_MS_3rd_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second_MS")
CD4_MS_3rd_2nd <- CD4_MS_3rd_2nd[rownames(count_cd4),]

CD4_MS_3rd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_MS")
CD4_MS_3rd_1st <- CD4_MS_3rd_1st[rownames(count_cd4),]

CD4_MS_PP_3rd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third_MS")
CD4_MS_PP_3rd <- CD4_MS_PP_3rd[rownames(count_cd4),]

#MSvsHP
CD4_MSvsHP_2nd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Second_First")
CD4_MSvsHP_2nd_1st <- CD4_MSvsHP_2nd_1st[rownames(count_cd4),]

CD4_MSvsHP_3rd_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second")
CD4_MSvsHP_3rd_2nd <- CD4_MSvsHP_3rd_2nd[rownames(count_cd4),]

CD4_MSvsHP_3rd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_First")
CD4_MSvsHP_3rd_1st <- CD4_MSvsHP_2nd_1st[rownames(count_cd4),]

CD4_MSvsHP_PP_3rd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third")
CD4_MSvsHP_PP_3rd <- CD4_MSvsHP_PP_3rd[rownames(count_cd4),]

CD4_MSvsHP_PP_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "PP_First")
CD4_MSvsHP_PP_1st <- CD4_MSvsHP_PP_1st[rownames(count_cd4),]


# Correlation plots CD4 ####
pdf("figures_manus/Correlations_CD4_RNAseq.pdf", width=10)

par(mfrow=c(2,5))

#HP
smoothScatter(CD4_HP_3rd_2nd$logFC ~ CD4_HP_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_HP_3rd_2nd$logFC, CD4_HP_2nd_1st$logFC))/100))


smoothScatter(CD4_HP_3rd_1st$logFC ~ CD4_HP_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_HP_3rd_1st$logFC, CD4_HP_2nd_1st$logFC))/100))


smoothScatter(CD4_HP_PP_3rd$logFC ~ CD4_HP_3rd_2nd$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-2nd", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_HP_PP_3rd$logFC, CD4_HP_3rd_2nd$logFC))/100))


smoothScatter(CD4_HP_PP_3rd$logFC ~ CD4_HP_3rd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_HP_PP_3rd$logFC, CD4_HP_3rd_1st$logFC))/100))


smoothScatter(CD4_HP_PP_3rd$logFC ~ CD4_HP_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_HP_PP_3rd$logFC, CD4_HP_2nd_1st$logFC))/100))



#MS
smoothScatter(CD4_MS_3rd_2nd$logFC ~ CD4_MS_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple Sclerosis",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_MS_3rd_2nd$logFC, CD4_MS_2nd_1st$logFC))/100))


smoothScatter(CD4_MS_3rd_1st$logFC ~ CD4_MS_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple Sclerosis",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_MS_3rd_1st$logFC, CD4_MS_2nd_1st$logFC))/100))


smoothScatter(CD4_MS_PP_3rd$logFC ~ CD4_MS_3rd_2nd$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple Sclerosis",
              ylab = "PP-3rd", xlab = "3rd-2nd", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_MS_PP_3rd$logFC, CD4_MS_3rd_2nd$logFC))/100))


smoothScatter(CD4_MS_PP_3rd$logFC ~ CD4_MS_3rd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple Sclerosis",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_MS_PP_3rd$logFC, CD4_MS_3rd_1st$logFC))/100))


smoothScatter(CD4_MS_PP_3rd$logFC ~ CD4_MS_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple Sclerosis",
              ylab = "PP-3rd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_MS_PP_3rd$logFC, CD4_MS_2nd_1st$logFC))/100))

dev.off()




# Limma on CD8 ####
metadata_cd8$Sample_Group <- gsub("MS_BP", "AA", metadata_cd8$Sample_Group) 
design <- model.matrix(~ Sample_Group + Proportions_Memory + Cell_Viability + Age, metadata_cd8 )
metadata_cd8$Sample_Group <- gsub("AA", "MS_BP", metadata_cd8$Sample_Group)

count_voom <- voom(count_cd8, design, plot=F)
corfit <- duplicateCorrelation(count_voom, design, block=metadata_cd8$Individual)

contr_matrix <- makeContrasts(First = Sample_GroupMS_1st_CD8 - Sample_GroupHP_1st_CD8 ,
                              Second = Sample_GroupMS_2nd_CD8 - Sample_GroupHP_2nd_CD8,
                              Third = Sample_GroupMS_3rd_CD8 - Sample_GroupHP_3rd_CD8,
                              PP = Sample_GroupMS_PP_CD8 - Sample_GroupHP_PP_CD8,
                              Second_First = (Sample_GroupMS_2nd_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupHP_2nd_CD8 - Sample_GroupHP_1st_CD8),
                              Third_Second = (Sample_GroupMS_3rd_CD8 - Sample_GroupMS_2nd_CD8) - (Sample_GroupHP_3rd_CD8 - Sample_GroupHP_2nd_CD8),
                              Third_First = (Sample_GroupMS_3rd_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupHP_3rd_CD8 - Sample_GroupHP_1st_CD8),
                              PP_Third = (Sample_GroupMS_PP_CD8 - Sample_GroupMS_3rd_CD8) - (Sample_GroupHP_PP_CD8 - Sample_GroupHP_3rd_CD8),
                              PP_First = (Sample_GroupMS_PP_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupHP_PP_CD8 - Sample_GroupHP_1st_CD8),
                              Second_MS = Sample_GroupMS_2nd_CD8 - Sample_GroupMS_1st_CD8,
                              Third_MS = Sample_GroupMS_3rd_CD8 - Sample_GroupMS_1st_CD8,
                              PP_MS = Sample_GroupMS_PP_CD8 - Sample_GroupMS_1st_CD8,
                              Second_HP = Sample_GroupHP_2nd_CD8 - Sample_GroupHP_1st_CD8,
                              Third_HP = Sample_GroupHP_3rd_CD8 - Sample_GroupHP_1st_CD8,
                              PP_HP = Sample_GroupHP_PP_CD8 - Sample_GroupHP_1st_CD8,
                              Third_Second_MS = Sample_GroupMS_3rd_CD8 - Sample_GroupMS_2nd_CD8,
                              PP_Third_MS = Sample_GroupMS_PP_CD8 - Sample_GroupMS_3rd_CD8,
                              PP_Second_MS = Sample_GroupMS_PP_CD8 - Sample_GroupMS_2nd_CD8,
                              Third_Second_HP = Sample_GroupHP_3rd_CD8 - Sample_GroupHP_2nd_CD8,
                              PP_Third_HP = Sample_GroupHP_PP_CD8 - Sample_GroupHP_3rd_CD8,
                              PP_Second_HP = Sample_GroupHP_PP_CD8 - Sample_GroupHP_2nd_CD8,
                              levels=colnames(design)
)
fit_CD8 <- lmFit(count_voom , design, block=metadata_cd8$Individual, correlation=corfit$consensus)
fit_CD8 <- contrasts.fit(fit_CD8, contrasts=contr_matrix)
fit_CD8 <- eBayes(fit_CD8)
summary(decideTests(fit_CD8))


# Gather Changes
#HP
CD8_HP_2nd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Second_HP")
CD8_HP_2nd_1st <- CD8_HP_2nd_1st[rownames(count_cd8),]

CD8_HP_3rd_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second_HP")
CD8_HP_3rd_2nd <- CD8_HP_3rd_2nd[rownames(count_cd8),]

CD8_HP_3rd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_HP")
CD8_HP_3rd_1st <- CD8_HP_3rd_1st[rownames(count_cd8),]

CD8_HP_PP_3rd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third_HP")
CD8_HP_PP_3rd <- CD8_HP_PP_3rd[rownames(count_cd8),]

#MS
CD8_MS_2nd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Second_MS")
CD8_MS_2nd_1st <- CD8_MS_2nd_1st[rownames(count_cd8),]

CD8_MS_3rd_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second_MS")
CD8_MS_3rd_2nd <- CD8_MS_3rd_2nd[rownames(count_cd8),]

CD8_MS_3rd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_MS")
CD8_MS_3rd_1st <- CD8_MS_3rd_1st[rownames(count_cd8),]

CD8_MS_PP_3rd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third_MS")
CD8_MS_PP_3rd <- CD8_MS_PP_3rd[rownames(count_cd8),]

#MSvsHP
CD8_MSvsHP_2nd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Second_First")
CD8_MSvsHP_2nd_1st <- CD8_MSvsHP_2nd_1st[rownames(count_cd8),]

CD8_MSvsHP_3rd_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second")
CD8_MSvsHP_3rd_2nd <- CD8_MSvsHP_3rd_2nd[rownames(count_cd8),]

CD8_MSvsHP_3rd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_First")
CD8_MSvsHP_3rd_1st <- CD8_MSvsHP_2nd_1st[rownames(count_cd8),]

CD8_MSvsHP_PP_3rd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third")
CD8_MSvsHP_PP_3rd <- CD8_MSvsHP_PP_3rd[rownames(count_cd8),]

CD8_MSvsHP_PP_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "PP_First")
CD8_MSvsHP_PP_1st <- CD8_MSvsHP_PP_1st[rownames(count_cd8),]

# Correlation plots CD8 ####
pdf("figures_manus/Correlations_CD8_RNAseq.pdf", width=10)

par(mfrow=c(2,5))

#HP
smoothScatter(CD8_HP_3rd_2nd$logFC ~ CD8_HP_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_HP_3rd_2nd$logFC, CD8_HP_2nd_1st$logFC))/100))


smoothScatter(CD8_HP_3rd_1st$logFC ~ CD8_HP_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_HP_3rd_1st$logFC, CD8_HP_2nd_1st$logFC))/100))


smoothScatter(CD8_HP_PP_3rd$logFC ~ CD8_HP_3rd_2nd$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-2nd", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_HP_PP_3rd$logFC, CD8_HP_3rd_2nd$logFC))/100))


smoothScatter(CD8_HP_PP_3rd$logFC ~ CD8_HP_3rd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_HP_PP_3rd$logFC, CD8_HP_3rd_1st$logFC))/100))


smoothScatter(CD8_HP_PP_3rd$logFC ~ CD8_HP_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Healthy Pregnant",
              ylab = "PP-3rd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_HP_PP_3rd$logFC, CD8_HP_2nd_1st$logFC))/100))


#MS
smoothScatter(CD8_MS_3rd_2nd$logFC ~ CD8_MS_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple Sclerosis",
              ylab = "3rd-2nd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_MS_3rd_2nd$logFC, CD8_MS_2nd_1st$logFC))/100))


smoothScatter(CD8_MS_3rd_1st$logFC ~ CD8_MS_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple Sclerosis",
              ylab = "3rd-1st", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_MS_3rd_1st$logFC, CD8_MS_2nd_1st$logFC))/100))


smoothScatter(CD8_MS_PP_3rd$logFC ~ CD8_MS_3rd_2nd$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple Sclerosis",
              ylab = "PP-3rd", xlab = "3rd-2nd", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_MS_PP_3rd$logFC, CD8_MS_3rd_2nd$logFC))/100))


smoothScatter(CD8_MS_PP_3rd$logFC ~ CD8_MS_3rd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple Sclerosis",
              ylab = "PP-3rd", xlab = "3rd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_MS_PP_3rd$logFC, CD8_MS_3rd_1st$logFC))/100))


smoothScatter(CD8_MS_PP_3rd$logFC ~ CD8_MS_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "Multiple Sclerosis",
              ylab = "PP-3rd", xlab = "2nd-1st", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_MS_PP_3rd$logFC, CD8_MS_2nd_1st$logFC))/100))

dev.off()




# Correlations MS vs HP ####
pdf("figures_manus/Correlations_MSvsHP_RNAseq.pdf", width=10)

par(mfrow=c(2,4))
#CD4
smoothScatter(CD4_MS_2nd_1st$logFC ~ CD4_HP_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "CD4 2nd-1st",
              ylab = "Multiple Sclerosis", xlab = "Healthy Pregnant", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_MS_2nd_1st$logFC, CD4_HP_2nd_1st$logFC))/100))

smoothScatter(CD4_MS_3rd_2nd$logFC ~ CD4_HP_3rd_2nd$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "CD4 3rd-2nd",
              ylab = "Multiple Sclerosis", xlab = "Healthy Pregnant", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_MS_3rd_2nd$logFC, CD4_HP_3rd_2nd$logFC))/100))

smoothScatter(CD4_MS_3rd_1st$logFC ~ CD4_HP_3rd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "CD4 3rd-1st",
              ylab = "Multiple Sclerosis", xlab = "Healthy Pregnant", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_MS_3rd_1st$logFC, CD4_HP_3rd_1st$logFC))/100))

smoothScatter(CD4_MS_PP_3rd$logFC ~ CD4_HP_PP_3rd$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "CD4 PP-3rd",
              ylab = "Multiple Sclerosis", xlab = "Healthy Pregnant", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD4_MS_PP_3rd$logFC, CD4_HP_PP_3rd$logFC))/100))


#CD8
smoothScatter(CD8_MS_2nd_1st$logFC ~ CD8_HP_2nd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "CD8 2nd-1st",
              ylab = "Multiple Sclerosis", xlab = "Healthy Pregnant", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_MS_2nd_1st$logFC, CD8_HP_2nd_1st$logFC))/100))

smoothScatter(CD8_MS_3rd_2nd$logFC ~ CD8_HP_3rd_2nd$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "CD8 3rd-2nd",
              ylab = "Multiple Sclerosis", xlab = "Healthy Pregnant", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_MS_3rd_2nd$logFC, CD8_HP_3rd_2nd$logFC))/100))

smoothScatter(CD8_MS_3rd_1st$logFC ~ CD8_HP_3rd_1st$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "CD8 3rd-1st",
              ylab = "Multiple Sclerosis", xlab = "Healthy Pregnant", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_MS_3rd_1st$logFC, CD8_HP_3rd_1st$logFC))/100))

smoothScatter(CD8_MS_PP_3rd$logFC ~ CD8_HP_PP_3rd$logFC, pch=16, cex=0.5, col="gray39", nrpoints=500, main = "CD8 PP-3rd",
              ylab = "Multiple Sclerosis", xlab = "Healthy Pregnant", cex.axis = 1.5, cex.lab = 1.5,
              colramp = colorRampPalette(c("white", "cornflowerblue", "moccasin", "lightsalmon", "tomato", "red", "darkred")))
legend("topright", legend = paste0("r = ", round(100*cor(CD8_MS_PP_3rd$logFC, CD8_HP_PP_3rd$logFC))/100))


dev.off()



# Plot n. DEGs ####
CD4_DEG_2nd_1st <- sum(CD4_MSvsHP_2nd_1st$P.Value<0.05)
CD4_DEG_3rd_1st <- sum(CD4_MSvsHP_3rd_1st$P.Value<0.05) 
CD4_DEG_PP_1st <- sum(CD4_MSvsHP_PP_1st$P.Value<0.05)


CD8_DEG_2nd_1st <- sum(CD8_MSvsHP_2nd_1st$P.Value<0.05)
CD8_DEG_3rd_1st <- sum(CD8_MSvsHP_3rd_1st$P.Value<0.05)
CD8_DEG_PP_1st <- sum(CD8_MSvsHP_PP_1st$P.Value<0.05)


df <- rbind( data.frame(n=c(CD4_DEG_2nd_1st, CD4_DEG_3rd_1st, CD4_DEG_PP_1st), time = 1:3, group = rep("CD4",3)),
             data.frame(n=c(CD8_DEG_2nd_1st, CD8_DEG_3rd_1st, CD8_DEG_PP_1st), time =1:3, group = rep("CD8",3)))

dev.off()
par(mar=c(1,1,1,1))

pdf("figures_manus/n_DEGs.pdf", width=5, height=5)
ggplot(df, aes(x=time, y=n, col=group)) + geom_line(size=2) + geom_point(size=5) +
  scale_x_continuous(breaks=c(1, 2, 3), labels=c("2nd-1st", "3rd-1st", "PP-1st")) + ylab("n. DEGs") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()








