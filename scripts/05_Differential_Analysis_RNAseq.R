# This script performs differential analysis on RNA-seq 
# Created and modified by Alberto Zenere, 2021-12-11

# Differential analysis on RNA-seq data ####
# 5.1 Set directory structure
# 5.2 Load packages
# 5.3 Load CD4 RNA-seq data
# 5.4 Differential analysis on Resting CD4 samples
# 5.5 Load CD8 RNA-seq data
# 5.6 Differential analysis on Resting CD8 samples
# 5.7 Save
# 5.8 Differential analysis on Activated vs Resting CD4
# 5.9 Differential analysis on Activated vs Resting CD8


rm(list=ls()) # remove all entries in the global environment 

pathlist <- .libPaths();
newpath <-  "C:/Users/albze08/.conda/envs/graMS/Lib/R/library"
.libPaths(newpath)

# 5.1 Set directory structure ####
main_dir <-"C:/Users/albze08/Desktop/phd/P4/methylation"
FOLDER_RDS <- "RDS_files/"
FOLDER_FIGURES <- "figures/"
rna_folder <- "C:/Users/albze08/Desktop/phd/P4/RNAseq/RDS/"

setwd("C:/Users/albze08/Desktop/phd/P4/methylation")

# 5.2 Load packages ####

pack_R <- c("limma","stringr","dplyr","tidyverse","pracma","ggrepel")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}

require("lumi")
require("ggpubr")
require("pheatmap")
library("VennDiagram") 
library("gridExtra")

set.seed(206)

# 5.3 Load CD4 RNA-seq data ####

# CD4 
count_cd4 <- readRDS(file = paste0(rna_folder, "count_cd4_tmm.RDS"))

metadata_cd4 <- readRDS(file = paste0(rna_folder, "metadata_cd4.RDS"))
metadata_cd4$Sample_Group <- gsub(" ", "", metadata_cd4$Sample_Group)

# Remove non pregnant
is_P <- (!metadata_cd4$TimePoint %in% "Nonpregnant") %>% which()
count_cd4 <- count_cd4[, is_P]
metadata_cd4 <- metadata_cd4[is_P, ]

# 5.4 Differential analysis on Resting CD4 ####
idx <- metadata_cd4$State=="Resting"
count_cd4_resting <- count_cd4[,idx]
metadata_cd4_resting <- metadata_cd4[idx,]

metadata_cd4_resting$Sample_Group <- gsub("MS_BP", "AA", metadata_cd4_resting$Sample_Group) 
design <- model.matrix(~ Sample_Group + Proportions_Memory + Cell_Viability + Age, metadata_cd4_resting )
metadata_cd4_resting$Sample_Group <- gsub("AA", "MS_BP", metadata_cd4_resting$Sample_Group)

count_voom <- voom(count_cd4_resting, design, plot=F)
corfit <- duplicateCorrelation(count_voom, design, block=metadata_cd4_resting$Individual)

contr_matrix <- makeContrasts(First = Sample_GroupMS_1st_CD4 - Sample_GroupHP_1st_CD4 ,
                              Second = Sample_GroupMS_2nd_CD4 - Sample_GroupHP_2nd_CD4,
                              Third = Sample_GroupMS_3rd_CD4 - Sample_GroupHP_3rd_CD4,
                              PP = Sample_GroupMS_PP_CD4 - Sample_GroupHP_PP_CD4,
                              Second_First = (Sample_GroupMS_2nd_CD4 - Sample_GroupMS_1st_CD4) - (Sample_GroupHP_2nd_CD4 - Sample_GroupHP_1st_CD4),
                              Third_First = (Sample_GroupMS_3rd_CD4 - Sample_GroupMS_1st_CD4) - (Sample_GroupHP_3rd_CD4 - Sample_GroupHP_1st_CD4),
                              PP_First = (Sample_GroupMS_PP_CD4 - Sample_GroupMS_1st_CD4) - (Sample_GroupHP_PP_CD4 - Sample_GroupHP_1st_CD4),
                              PP_Third = (Sample_GroupMS_PP_CD4 - Sample_GroupMS_3rd_CD4) - (Sample_GroupHP_PP_CD4 - Sample_GroupHP_3rd_CD4),
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
fit <- lmFit(count_voom , design, block=metadata_cd4_resting$Individual, correlation=corfit$consensus)
fit <- contrasts.fit(fit, contrasts=contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

#nominally DEG 3rd-1st
HP_3rd_1st_cd4 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Third_HP")
MS_3rd_1st_cd4 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Third_MS")

#nominally DEG PP-3rd
HP_PP_3rd_cd4 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="PP_Third_HP")
MS_PP_3rd_cd4 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="PP_Third_MS")

#nominally significant 2nd-1st
HP_2nd_1st_cd4 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Second_HP")
MS_2nd_1st_cd4 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Second_MS")

#nominally significant 3rd-2nd
HP_3rd_2nd_cd4 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Third_Second_HP")
MS_3rd_2nd_cd4 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Third_Second_MS")

#nominally significant PP-2nd
HP_PP_2nd_cd4 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="PP_Second_HP")
MS_PP_2nd_cd4 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="PP_Second_MS")

#ALL DEG 3rd-1st
HP_3rd_1st_cd4_all <- topTable(fit,adjust.method="BH", p.value = Inf, number = Inf, coef="Third_HP")
MS_3rd_1st_cd4_all <- topTable(fit,adjust.method="BH", p.value = Inf, number = Inf, coef="Third_MS")

#ALL DEG PP-3rd
HP_PP_3rd_cd4_all <- topTable(fit,adjust.method="BH", p.value = Inf, number = Inf, coef="PP_Third_HP")
MS_PP_3rd_cd4_all <- topTable(fit,adjust.method="BH", p.value = Inf, number = Inf, coef="PP_Third_MS")

#ALL 2nd-1st
HP_2nd_1st_cd4_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="Second_HP")
MS_2nd_1st_cd4_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="Second_MS")

#ALL 3rd-2nd
HP_3rd_2nd_cd4_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="Third_Second_HP")
MS_3rd_2nd_cd4_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="Third_Second_MS")

#ALL PP-2nd
HP_PP_2nd_cd4_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="PP_Second_HP")
MS_PP_2nd_cd4_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="PP_Second_MS")


# 5.5 load CD8 RNA-seq data ####
count_cd8 <- readRDS(file = paste0(rna_folder, "count_cd8_tmm.RDS"))

metadata_cd8 <- readRDS(file = paste0(rna_folder, "metadata_cd8.RDS"))
metadata_cd8$Sample_Group <- gsub(" ", "", metadata_cd8$Sample_Group)

# Remove non pregnant
is_P <- (!metadata_cd8$TimePoint %in% "Nonpregnant") %>% which()
count_cd8 <- count_cd8[, is_P]
metadata_cd8 <- metadata_cd8[is_P, ]

# 5.6 Differential analysis on CD8 Resting samples ####
idx <- metadata_cd8$State=="Resting"
count_cd8_resting <- count_cd8[,idx]
metadata_cd8_resting <- metadata_cd8[idx,]

metadata_cd8_resting$Sample_Group <- gsub("MS_BP", "AA", metadata_cd8_resting$Sample_Group) 
design <- model.matrix(~ Sample_Group + Proportions_Memory + Cell_Viability + Age, metadata_cd8_resting )
metadata_cd8_resting$Sample_Group <- gsub("AA", "MS_BP", metadata_cd8_resting$Sample_Group)

count_voom <- voom(count_cd8_resting, design, plot=F)
corfit <- duplicateCorrelation(count_voom, design, block=metadata_cd8_resting$Individual)

contr_matrix <- makeContrasts(First = Sample_GroupMS_1st_CD8 - Sample_GroupHP_1st_CD8 ,
                              Second = Sample_GroupMS_2nd_CD8 - Sample_GroupHP_2nd_CD8,
                              Third = Sample_GroupMS_3rd_CD8 - Sample_GroupHP_3rd_CD8,
                              PP = Sample_GroupMS_PP_CD8 - Sample_GroupHP_PP_CD8,
                              Second_First = (Sample_GroupMS_2nd_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupHP_2nd_CD8 - Sample_GroupHP_1st_CD8),
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
fit <- lmFit(count_voom , design, block=metadata_cd8_resting$Individual, correlation=corfit$consensus)
fit <- contrasts.fit(fit, contrasts=contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

# DEGs CD8 
#nominally DEG 3rd-1st
HP_3rd_1st_cd8 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Third_HP")
MS_3rd_1st_cd8 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Third_MS")

#nominally DEG PP-3rd
HP_PP_3rd_cd8 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="PP_Third_HP")
MS_PP_3rd_cd8 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="PP_Third_MS")

#nominally significant 2nd-1st
HP_2nd_1st_cd8 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Second_HP")
MS_2nd_1st_cd8 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Second_MS")

#nominally significant 3rd-2nd
HP_3rd_2nd_cd8 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Third_Second_HP")
MS_3rd_2nd_cd8 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Third_Second_MS")

#nominally significant PP-2nd
HP_PP_2nd_cd8 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="PP_Second_HP")
MS_PP_2nd_cd8 <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="PP_Second_MS")

#ALL 3rd-1st
HP_3rd_1st_cd8_all <- topTable(fit,adjust.method="BH", p.value = Inf, number = Inf, coef="Third_HP")
MS_3rd_1st_cd8_all <- topTable(fit,adjust.method="BH", p.value = Inf, number = Inf, coef="Third_MS")

#ALL DEG PP-3rd
HP_PP_3rd_cd8_all <- topTable(fit,adjust.method="BH", p.value = Inf, number = Inf, coef="PP_Third_HP")
MS_PP_3rd_cd8_all <- topTable(fit,adjust.method="BH", p.value = Inf, number = Inf, coef="PP_Third_MS")

#ALL 2nd-1st
HP_2nd_1st_cd8_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="Second_HP")
MS_2nd_1st_cd8_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="Second_MS")

#ALL 3rd-2nd
HP_3rd_2nd_cd8_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="Third_Second_HP")
MS_3rd_2nd_cd8_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="Third_Second_MS")

#ALL PP-2nd
HP_PP_2nd_cd8_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="PP_Second_HP")
MS_PP_2nd_cd8_all <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="PP_Second_MS")


# 5.7 Save ####
#rebound CD4
saveRDS(HP_3rd_1st_cd4, paste0(FOLDER_RDS, "DMG/CD4_3rd_HP_rna.RDS"))
saveRDS(HP_PP_3rd_cd4, paste0(FOLDER_RDS, "DMG/CD4_PP_HP_rna.RDS"))
saveRDS(MS_3rd_1st_cd4, paste0(FOLDER_RDS, "DMG/CD4_3rd_MS_rna.RDS"))
saveRDS(MS_PP_3rd_cd4, paste0(FOLDER_RDS, "DMG/CD4_PP_MS_rna.RDS"))

#all CD4
saveRDS(HP_3rd_1st_cd4_all, paste0(FOLDER_RDS, "DMG/CD4_3rd_HP_rna_all.RDS"))
saveRDS(HP_PP_3rd_cd4_all, paste0(FOLDER_RDS, "DMG/CD4_PP_HP_rna_all.RDS"))
saveRDS(MS_3rd_1st_cd4_all, paste0(FOLDER_RDS, "DMG/CD4_3rd_MS_rna_all.RDS"))
saveRDS(MS_PP_3rd_cd4_all, paste0(FOLDER_RDS, "DMG/CD4_PP_MS_rna_all.RDS"))


saveRDS(HP_2nd_1st_cd4_all, paste0(FOLDER_RDS, "DMG/HP_2nd_1st_cd4_all.RDS"))
saveRDS(MS_2nd_1st_cd4_all, paste0(FOLDER_RDS, "DMG/MS_2nd_1st_cd4_all.RDS"))
saveRDS(HP_3rd_2nd_cd4_all, paste0(FOLDER_RDS, "DMG/HP_3rd_2nd_cd4_all.RDS"))
saveRDS(MS_3rd_2nd_cd4_all, paste0(FOLDER_RDS, "DMG/MS_3rd_2nd_cd4_all.RDS"))
saveRDS(HP_PP_2nd_cd4_all, paste0(FOLDER_RDS, "DMG/HP_PP_2nd_cd4_all.RDS"))
saveRDS(MS_PP_2nd_cd4_all, paste0(FOLDER_RDS, "DMG/MS_PP_2nd_cd4_all.RDS"))


#rebound CD8
saveRDS(HP_3rd_1st_cd8, paste0(FOLDER_RDS, "DMG/CD8_3rd_HP_rna.RDS"))
saveRDS(HP_PP_3rd_cd8, paste0(FOLDER_RDS, "DMG/CD8_PP_HP_rna.RDS"))
saveRDS(MS_3rd_1st_cd8, paste0(FOLDER_RDS, "DMG/CD8_3rd_MS_rna.RDS"))
saveRDS(MS_PP_3rd_cd8, paste0(FOLDER_RDS, "DMG/CD8_PP_MS_rna.RDS"))

#rebound all CD8
saveRDS(HP_3rd_1st_cd8_all, paste0(FOLDER_RDS, "DMG/CD8_3rd_HP_rna_all.RDS"))
saveRDS(HP_PP_3rd_cd8_all, paste0(FOLDER_RDS, "DMG/CD8_PP_HP_rna_all.RDS"))
saveRDS(MS_3rd_1st_cd8_all, paste0(FOLDER_RDS, "DMG/CD8_3rd_MS_rna_all.RDS"))
saveRDS(MS_PP_3rd_cd8_all, paste0(FOLDER_RDS, "DMG/CD8_PP_MS_rna_all.RDS"))

# all
saveRDS(HP_2nd_1st_cd8_all, paste0(FOLDER_RDS, "DMG/HP_2nd_1st_cd8_all.RDS"))
saveRDS(MS_2nd_1st_cd8_all, paste0(FOLDER_RDS, "DMG/MS_2nd_1st_cd8_all.RDS"))
saveRDS(HP_3rd_2nd_cd8_all, paste0(FOLDER_RDS, "DMG/HP_3rd_2nd_cd8_all.RDS"))
saveRDS(MS_3rd_2nd_cd8_all, paste0(FOLDER_RDS, "DMG/MS_3rd_2nd_cd8_all.RDS"))
saveRDS(HP_PP_2nd_cd8_all, paste0(FOLDER_RDS, "DMG/HP_PP_2nd_cd8_all.RDS"))
saveRDS(MS_PP_2nd_cd8_all, paste0(FOLDER_RDS, "DMG/MS_PP_2nd_cd8_all.RDS"))


#2nd-1st CD4
saveRDS(HP_2nd_1st_cd4, paste0(FOLDER_RDS, "DMG/CD4_2nd_HP_rna.RDS"))
saveRDS(MS_2nd_1st_cd4, paste0(FOLDER_RDS, "DMG/CD4_2nd_MS_rna.RDS"))

#nominally significant 3rd-2nd
saveRDS(HP_3rd_2nd_cd4, paste0(FOLDER_RDS, "DMG/CD4_3rd_2nd_HP_rna.RDS"))
saveRDS(MS_3rd_2nd_cd4, paste0(FOLDER_RDS, "DMG/CD4_3rd_2nd_MS_rna.RDS"))

#nominally significant PP-2nd
saveRDS(HP_PP_2nd_cd4, paste0(FOLDER_RDS, "DMG/CD4_PP_2nd_HP_rna.RDS"))
saveRDS(MS_PP_2nd_cd4, paste0(FOLDER_RDS, "DMG/CD4_PP_2nd_MS_rna.RDS"))

#MSvsHP
saveRDS(MSvsHP_cd4_1st, paste0(FOLDER_RDS, "DMG/MSvsHP_cd4_1st.RDS"))
saveRDS(MSvsHP_cd4_2nd, paste0(FOLDER_RDS, "DMG/MSvsHP_cd4_2nd.RDS"))
saveRDS(MSvsHP_cd4_3rd, paste0(FOLDER_RDS, "DMG/MSvsHP_cd4_3rd.RDS"))
saveRDS(MSvsHP_cd4_PP, paste0(FOLDER_RDS, "DMG/MSvsHP_cd4_PP.RDS"))

saveRDS(MSvsHP_cd4_3rd_1st, paste0(FOLDER_RDS, "DMG/MSvsHP_cd4_3rd_1st.RDS"))
saveRDS(MSvsHP_cd4_PP_3rd, paste0(FOLDER_RDS, "DMG/MSvsHP_cd4_PP_3rd.RDS"))


#2nd-1st CD8
saveRDS(HP_2nd_1st_cd8, paste0(FOLDER_RDS, "DMG/CD8_2nd_HP_rna.RDS"))
saveRDS(MS_2nd_1st_cd8, paste0(FOLDER_RDS, "DMG/CD8_2nd_MS_rna.RDS"))

#nominally significant 3rd-2nd
saveRDS(HP_3rd_2nd_cd8, paste0(FOLDER_RDS, "DMG/CD8_3rd_2nd_HP_rna.RDS"))
saveRDS(MS_3rd_2nd_cd8, paste0(FOLDER_RDS, "DMG/CD8_3rd_2nd_MS_rna.RDS"))

#nominally significant PP-2nd
saveRDS(HP_PP_2nd_cd8, paste0(FOLDER_RDS, "DMG/CD8_PP_2nd_HP_rna.RDS"))
saveRDS(MS_PP_2nd_cd8, paste0(FOLDER_RDS, "DMG/CD8_PP_2nd_MS_rna.RDS"))



# 5.8 Differential analysis on Activated vs Resting CD4 ####
design <- model.matrix(~ Sample_Group + State + Activation + 
                         Proportions_Memory + Cell_Viability + Age, metadata_cd4 )

count_voom <- voom(count_cd4, design, plot=F)
corfit <- duplicateCorrelation(count_voom, design, block=metadata_cd4$Individual)

fit <- lmFit(count_voom , design)
fit <- eBayes(fit)
summary(decideTests(fit))

# DEGs CD4
State_cd4 <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="StateResting")
saveRDS(State_cd4, paste0(FOLDER_RDS, "DMG/State_cd4.RDS"))


# 5.9 Differential analysis on Activated vs Resting CD8 ####
design <- model.matrix(~ Sample_Group + State + Activation + 
                         Proportions_Memory + Cell_Viability + Age, metadata_cd8 )

count_voom <- voom(count_cd8, design, plot=F)
corfit <- duplicateCorrelation(count_voom, design, block=metadata_cd8$Individual)

fit <- lmFit(count_voom , design)
fit <- eBayes(fit)
summary(decideTests(fit))

# DEGs CD8
State_cd8 <- topTable(fit,adjust.method="none", p.value = Inf, number = Inf, coef="StateResting")
saveRDS(State_cd8, paste0(FOLDER_RDS, "DMG/State_cd8.RDS"))

