#Is there an overlap between the seed genes and the dynamic change (activated-resting) in RNAseq


rm(list=ls()) # remove all entries in the global environment 

pathlist <- .libPaths();
newpath <-  "C:/Users/albze08/.conda/envs/graMS/Lib/R/library"
.libPaths(newpath)

# Set directory structure ------------------------------------------------------
main_dir <-"C:/Users/albze08/Desktop/phd/P4/methylation"
FOLDER_RDS <- "RDS_files/"
FOLDER_FIGURES <- "figures/"

setwd("C:/Users/albze08/Desktop/phd/P4/methylation")

# Load packages ----------------------------------------------------------------

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


# Read in files
methyl_raw <- readRDS(file = paste0(FOLDER_RDS, "TH2636_methyl_raw_211123.RDS"))
methyl_combat <- readRDS(file = paste0(FOLDER_RDS, "TH2636_methyl_combat_211123.RDS"))
ProbeFeatures <- readRDS(file=paste0(FOLDER_RDS, "TH2636_ProbeFeatures_211123.RDS"))

# Replace , with . in the methyl_raw$pd
methyl_raw$pd$Proportions_Naive <- as.numeric(gsub(",", ".",methyl_raw$pd$Proportions_Naive))
methyl_raw$pd$Proportions_Memory <- as.numeric(gsub(",", ".",methyl_raw$pd$Proportions_Memory))
methyl_raw$pd$Cell_Viability <- as.numeric(gsub(",", ".",methyl_raw$pd$Cell_Viability))
methyl_raw$pd$Pregnancy_Week <- as.numeric(gsub(",", ".",methyl_raw$pd$Pregnancy_Week))
methyl_raw$pd$Time_Since_Treatment <- as.numeric(gsub(",", ".",methyl_raw$pd$Time_Since_Treatment))

# Set all variables to character
methyl_raw$pd[] <- lapply(methyl_raw$pd, as.character)

# Set class for the different variables
methyl_raw$pd$Age <- as.numeric(methyl_raw$pd$Age)
methyl_raw$pd$Cell_Viability <- as.numeric(methyl_raw$pd$Cell_Viability)
methyl_raw$pd$Delivery_Week <- as.numeric(methyl_raw$pd$Delivery_Week)
methyl_raw$pd$Pregnancy_Week <- as.numeric(methyl_raw$pd$Pregnancy_Week)
methyl_raw$pd$Proportions_Naive <- as.numeric(methyl_raw$pd$Proportions_Naive)
methyl_raw$pd$Proportions_Memory <- as.numeric(methyl_raw$pd$Proportions_Memory)
methyl_raw$pd$Time_Since_Treatment<- as.numeric(methyl_raw$pd$Time_Since_Treatment)
methyl_raw$pd$Sample_Group<- gsub(" ", "", methyl_raw$pd$Sample_Group)

check_class <- lapply(methyl_raw$pd, class) #make sure the covariates have the right class

#universe file
universe <- readRDS("C:/Users/albze08/Desktop/phd/P4/methylation/RDS_files/universe_CD8.RDS")

# Subset CD4
file <- "CD4"
cell_type <- methyl_raw$pd$Sample_Type
idx <- grep(file, cell_type)

methyl_CD4 <- methyl_combat[,idx]
pheno_CD4 <- methyl_raw$pd[idx,]

# remove HC
idx_HC <- which(is.element(pheno_CD4$Disease,"HC"))

methyl_NP_CD4 <- methyl_CD4[,idx_HC]
pheno_NP_CD4 <- pheno_CD4[idx_HC,]

pheno_CD4 <- pheno_CD4[-idx_HC,]
methyl_CD4 <- methyl_CD4[,-idx_HC]

#HP and MS
idx_HP <- pheno_CD4$Disease == "HP"
idx_MS <- pheno_CD4$Disease == "MS"

methyl_CD4_HP <- methyl_CD4[, idx_HP]
pheno_CD4_HP <- pheno_CD4[idx_HP, ]

methyl_CD4_MS <- methyl_CD4[, idx_MS]
pheno_CD4_MS <- pheno_CD4[idx_MS, ]

#group before samples
is_before <- grep("BP", pheno_CD4$Sample_Group)
methyl_BP_CD4 <- cbind(methyl_NP_CD4, methyl_CD4[,is_before])
pheno_BP_CD4 <- rbind(pheno_NP_CD4, pheno_CD4[is_before,])

# Subset CD8
file <- "CD8"
cell_type <- methyl_raw$pd$Sample_Type
idx <- grep(file, cell_type)

methyl_CD8 <- methyl_combat[,idx]
pheno_CD8 <- methyl_raw$pd[idx,] 

# remove HC
idx_HC <- which(is.element(pheno_CD8$Disease,"HC"))

methyl_NP_CD8 <- methyl_CD8[,idx_HC]
pheno_NP_CD8 <- pheno_CD8[idx_HC,]

pheno_CD8 <- pheno_CD8[-idx_HC,]
methyl_CD8 <- methyl_CD8[,-idx_HC]

#HP and MS
idx_HP <- pheno_CD8$Disease == "HP"
idx_MS <- pheno_CD8$Disease == "MS"

methyl_CD8_HP <- methyl_CD8[, idx_HP]
pheno_CD8_HP <- pheno_CD8[idx_HP, ]

methyl_CD8_MS <- methyl_CD8[, idx_MS]
pheno_CD8_MS <- pheno_CD8[idx_MS, ]

#group before samples
is_before <- grep("BP", pheno_CD8$Sample_Group)
methyl_BP_CD8 <- cbind(methyl_NP_CD8, methyl_CD8[,is_before])
pheno_BP_CD8 <- rbind(pheno_NP_CD8, pheno_CD8[is_before,])


# Gather beta and delta beta values ####
gather_mean_beta <- function(count,meta){
  #beta
  cpg_1st <- count[, meta$TimePoint=="1st"] %>% rowMeans()
  cpg_2nd <- count[, meta$TimePoint=="2nd"] %>% rowMeans()
  cpg_3rd <- count[, meta$TimePoint=="3rd"] %>% rowMeans()
  cpg_PP <- count[, meta$TimePoint=="Postpartum"] %>% rowMeans()
  
  beta <- data.frame( cbind(cpg_1st,cpg_2nd,cpg_3rd,cpg_PP) )
  colnames(beta) <- c("1st", "2nd", "3rd", "PP")
  rownames(beta) <- rownames(count)
  
  #delta beta
  delta_2nd <- cpg_2nd - cpg_1st
  delta_3rd <- cpg_3rd - cpg_2nd
  delta_PP <- cpg_PP - cpg_3rd
  delta_PP_2nd <- cpg_PP - cpg_2nd
  
  delta_beta <- data.frame( cbind(delta_2nd, delta_3rd, delta_PP, delta_PP_2nd) )
  colnames(delta_beta) <- c("2nd_1st", "3rd_2nd", "PP_3rd", "PP_2nd")
  rownames(delta_beta) <- rownames(count)
  
  #delta beta baseline
  delta_2nd <- cpg_2nd - cpg_1st
  delta_3rd <- cpg_3rd - cpg_1st
  delta_PP <- cpg_PP - cpg_1st
  
  delta_beta_baseline <- data.frame( cbind(delta_2nd, delta_3rd, delta_PP) )
  colnames(delta_beta_baseline) <- c("2nd_1st", "3rd_1st", "PP_1st")
  rownames(delta_beta_baseline) <- rownames(count)
  
  out <- list( beta, delta_beta, delta_beta_baseline )
  return(out)
}

beta_CD4_HP <- gather_mean_beta(methyl_CD4_HP, pheno_CD4_HP)[[1]]
beta_CD4_MS <- gather_mean_beta(methyl_CD4_MS, pheno_CD4_MS)[[1]]
beta_CD8_HP <- gather_mean_beta(methyl_CD8_HP, pheno_CD8_HP)[[1]]
beta_CD8_MS <- gather_mean_beta(methyl_CD8_MS, pheno_CD8_MS)[[1]]

delta_beta_CD4_HP <- gather_mean_beta(methyl_CD4_HP, pheno_CD4_HP)[[2]]
delta_beta_CD4_MS <- gather_mean_beta(methyl_CD4_MS, pheno_CD4_MS)[[2]]
delta_beta_CD8_HP <- gather_mean_beta(methyl_CD8_HP, pheno_CD8_HP)[[2]]
delta_beta_CD8_MS <- gather_mean_beta(methyl_CD8_MS, pheno_CD8_MS)[[2]]

delta_base_beta_CD4_HP <- gather_mean_beta(methyl_CD4_HP, pheno_CD4_HP)[[3]]
delta_base_beta_CD4_MS <- gather_mean_beta(methyl_CD4_MS, pheno_CD4_MS)[[3]]
delta_base_beta_CD8_HP <- gather_mean_beta(methyl_CD8_HP, pheno_CD8_HP)[[3]]
delta_base_beta_CD8_MS <- gather_mean_beta(methyl_CD8_MS, pheno_CD8_MS)[[3]]

# Limma on CD4 ####
M_CD4 <- beta2m(methyl_CD4)

pheno_CD4$Sample_Group <- gsub("MS_BP", "AA", pheno_CD4$Sample_Group) #to make sure it is discarded below
design_CD4 <- model.matrix(~ Sample_Group + Age + Proportions_Memory + Cell_Viability, pheno_CD4 )
pheno_CD4$Sample_Group <- gsub("AA", "MS_BP", pheno_CD4$Sample_Group)

# corfit <- duplicateCorrelation(M_CD4, design_CD4, block=pheno_CD4$Individual) #block the individual effect
# saveRDS(corfit, file=paste0(FOLDER_RDS, "DMG/MSvsHP/seed/CD4_corfit.RDS"))
corfit <- readRDS(paste0(FOLDER_RDS, "DMG/MSvsHP/seed/CD4_corfit.RDS"))

contr_matrix <- limma::makeContrasts(First = Sample_GroupMS_1st_CD4 - Sample_GroupHP_1st_CD4 ,
                                     Second = Sample_GroupMS_2nd_CD4 - Sample_GroupHP_2nd_CD4,
                                     Third = Sample_GroupMS_3rd_CD4 - Sample_GroupHP_3rd_CD4,
                                     PP = Sample_GroupMS_PP_CD4 - Sample_GroupHP_PP_CD4,
                                     # Before_First = (Sample_GroupMS_BP_CD4 - Sample_GroupMS_1st_CD4) - (Sample_GroupNP_CD4 - Sample_GroupHP_1st_CD4),
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
                                     PP_Second_MS = Sample_GroupMS_PP_CD4 - Sample_GroupMS_2nd_CD4,
                                     PP_Third_MS = Sample_GroupMS_PP_CD4 - Sample_GroupMS_3rd_CD4,
                                     Third_Second_HP = Sample_GroupHP_3rd_CD4 - Sample_GroupHP_2nd_CD4,
                                     PP_Second_HP = Sample_GroupHP_PP_CD4 - Sample_GroupHP_2nd_CD4,
                                     PP_Third_HP = Sample_GroupHP_PP_CD4 - Sample_GroupHP_3rd_CD4,
                                     levels=colnames(design_CD4))

fit_CD4 <- lmFit(M_CD4, design_CD4, block=pheno_CD4$Individual, correlation=corfit$consensus)
fit_CD4 <- contrasts.fit(fit_CD4, contr_matrix)
fit_CD4 <- eBayes(fit_CD4)
summary(decideTests(fit_CD4))

#Gather DMPs CD4 ####
THS <- 0.05
CD4_3rd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Third_HP")
CD4_3rd <- CD4_3rd[rownames(methyl_CD4),]
CD4_3rd_HP <- CD4_3rd[(CD4_3rd$P.Value<0.05) & abs(delta_base_beta_CD4_HP[,"3rd_1st"])>THS,] 

CD4_3rd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Third_MS")
CD4_3rd <- CD4_3rd[rownames(methyl_CD4),]
CD4_3rd_MS <- CD4_3rd[(CD4_3rd$P.Value<0.05) & abs(delta_base_beta_CD4_MS[,"3rd_1st"])>THS,] 

CD4_PP <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "PP_Third_HP")
CD4_PP <- CD4_PP[rownames(methyl_CD4),]
CD4_PP_HP <- CD4_PP[(CD4_PP$P.Value<0.05) & abs(delta_beta_CD4_HP[,"PP_3rd"])>THS,]

CD4_PP <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "PP_Third_MS")
CD4_PP <- CD4_PP[rownames(methyl_CD4),]
CD4_PP_MS <- CD4_PP[(CD4_PP$P.Value<0.05) & abs(delta_beta_CD4_MS[,"PP_3rd"])>THS,]

#2nd-1st
CD4_2nd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Second_HP")
CD4_2nd_1st <- CD4_2nd_1st[rownames(methyl_CD4),]
CD4_2nd_1st_HP <- CD4_2nd_1st[(CD4_2nd_1st$P.Value<0.05) & abs(delta_base_beta_CD4_HP[,"2nd_1st"])>THS,] 

CD4_2nd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Second_MS")
CD4_2nd_1st <- CD4_2nd_1st[rownames(methyl_CD4),]
CD4_2nd_1st_MS <- CD4_2nd_1st[(CD4_2nd_1st$P.Value<0.05) & abs(delta_base_beta_CD4_MS[,"2nd_1st"])>THS,] 

#3rd-2nd
CD4_3rd_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Third_Second_HP")
CD4_3rd_2nd <- CD4_3rd_2nd[rownames(methyl_CD4),]
CD4_3rd_2nd_HP <- CD4_3rd_2nd[(CD4_3rd_2nd$P.Value<0.05) & abs(delta_beta_CD4_HP[,"3rd_2nd"])>THS,] 

CD4_3rd_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Third_Second_MS")
CD4_3rd_2nd <- CD4_3rd_2nd[rownames(methyl_CD4),]
CD4_3rd_2nd_MS <- CD4_3rd_2nd[(CD4_3rd_2nd$P.Value<0.05) & abs(delta_beta_CD4_MS[,"3rd_2nd"])>THS,] 

#PP-2nd
CD4_PP_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "PP_Second_HP")
CD4_PP_2nd <- CD4_PP_2nd[rownames(methyl_CD4),]
CD4_PP_2nd_HP <- CD4_PP_2nd[(CD4_PP_2nd$P.Value<0.05) & abs(delta_beta_CD4_HP[,"PP_2nd"])>THS,] 

CD4_PP_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "PP_Second_MS")
CD4_PP_2nd <- CD4_PP_2nd[rownames(methyl_CD4),]
CD4_PP_2nd_MS <- CD4_PP_2nd[(CD4_PP_2nd$P.Value<0.05) & abs(delta_beta_CD4_MS[,"PP_2nd"])>THS,] 

#Gather all cpgs CD4 ####
CD4_3rd_HP_all <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Third_HP")
CD4_3rd_HP_all <- CD4_3rd_HP_all[rownames(methyl_CD4),]
CD4_3rd_HP_all$delta_beta <- delta_base_beta_CD4_HP[,"3rd_1st"]

CD4_3rd_MS_all <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Third_MS")
CD4_3rd_MS_all <- CD4_3rd_MS_all[rownames(methyl_CD4),]
CD4_3rd_MS_all$delta_beta <- delta_base_beta_CD4_MS[,"3rd_1st"]

CD4_PP_HP_all <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "PP_Third_HP")
CD4_PP_HP_all <- CD4_PP_HP_all[rownames(methyl_CD4),]
CD4_PP_HP_all$delta_beta <- delta_beta_CD4_HP[,"PP_3rd"]

CD4_PP_MS_all <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "PP_Third_MS")
CD4_PP_MS_all <- CD4_PP_MS_all[rownames(methyl_CD4),]
CD4_PP_MS_all$delta_beta <- delta_beta_CD4_MS[,"PP_3rd"]

#2nd-1st
CD4_2nd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Second_HP")
CD4_2nd_1st_HP_all <- CD4_2nd_1st[rownames(methyl_CD4),]

CD4_2nd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Second_MS")
CD4_2nd_1st_MS_all <- CD4_2nd_1st[rownames(methyl_CD4),]

#3rd-2nd
CD4_3rd_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Third_Second_HP")
CD4_3rd_2nd_HP_all <- CD4_3rd_2nd[rownames(methyl_CD4),]

CD4_3rd_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "Third_Second_MS")
CD4_3rd_2nd_MS_all <- CD4_3rd_2nd[rownames(methyl_CD4),]

#PP-2nd
CD4_PP_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "PP_Second_HP")
CD4_PP_2nd_HP_all <- CD4_PP_2nd[rownames(methyl_CD4),]

CD4_PP_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "inf", coef = "PP_Second_MS")
CD4_PP_2nd_MS_all <- CD4_PP_2nd[rownames(methyl_CD4),]



#MS vs HP
MSvsHP_cd4_1st <- topTable(fit_CD4,adjust.method="BH", p.value = Inf, number = Inf, coef="First")
MSvsHP_cd4_2nd <- topTable(fit_CD4,adjust.method="BH", p.value = Inf, number = Inf, coef="Second")
MSvsHP_cd4_3rd <- topTable(fit_CD4,adjust.method="BH", p.value = Inf, number = Inf, coef="Third")
MSvsHP_cd4_PP <- topTable(fit_CD4,adjust.method="BH", p.value = Inf, number = Inf, coef="PP")

MSvsHP_cd4_3rd_1st <- topTable(fit_CD4,adjust.method="BH", p.value = Inf, number = Inf, coef="Third_First")
MSvsHP_cd4_PP_3rd <- topTable(fit_CD4,adjust.method="BH", p.value = Inf, number = Inf, coef="PP_Third")


# Limma on CD8 ####
M_CD8 <- beta2m(methyl_CD8)

pheno_CD8$Sample_Group <- gsub("MS_BP", "AA", pheno_CD8$Sample_Group) #to make sure it is discarded below
design_CD8 <- model.matrix(~ Sample_Group + Age + Proportions_Memory + Cell_Viability, pheno_CD8 )
pheno_CD8$Sample_Group <- gsub("AA", "MS_BP", pheno_CD8$Sample_Group)

# corfit <- duplicateCorrelation(M_CD8, design_CD8, block=pheno_CD8$Individual) #block the individual effect
# saveRDS(corfit, file=paste0(FOLDER_RDS, "DMG/MSvsHP/seed/CD8_corfit.RDS"))
corfit <- readRDS(paste0(FOLDER_RDS, "DMG/MSvsHP/seed/CD8_corfit.RDS"))

contr_matrix <- limma::makeContrasts(First = Sample_GroupMS_1st_CD8 - Sample_GroupHP_1st_CD8 ,
                                     Second = Sample_GroupMS_2nd_CD8 - Sample_GroupHP_2nd_CD8,
                                     Third = Sample_GroupMS_3rd_CD8 - Sample_GroupHP_3rd_CD8,
                                     PP = Sample_GroupMS_PP_CD8 - Sample_GroupHP_PP_CD8,
                                     # Before_First = (Sample_GroupMS_BP_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupNP_CD8 - Sample_GroupHP_1st_CD8),
                                     Second_First = (Sample_GroupMS_2nd_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupHP_2nd_CD8 - Sample_GroupHP_1st_CD8),
                                     Third_First = (Sample_GroupMS_3rd_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupHP_3rd_CD8 - Sample_GroupHP_1st_CD8),
                                     PP_First = (Sample_GroupMS_PP_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupHP_PP_CD8 - Sample_GroupHP_1st_CD8),
                                     PP_Third = (Sample_GroupMS_PP_CD8 - Sample_GroupMS_3rd_CD8) - (Sample_GroupHP_PP_CD8 - Sample_GroupHP_3rd_CD8),
                                     Second_MS = Sample_GroupMS_2nd_CD8 - Sample_GroupMS_1st_CD8,
                                     Third_MS = Sample_GroupMS_3rd_CD8 - Sample_GroupMS_1st_CD8,
                                     PP_MS = Sample_GroupMS_PP_CD8 - Sample_GroupMS_1st_CD8,
                                     Second_HP = Sample_GroupHP_2nd_CD8 - Sample_GroupHP_1st_CD8,
                                     Third_HP = Sample_GroupHP_3rd_CD8 - Sample_GroupHP_1st_CD8,
                                     PP_HP = Sample_GroupHP_PP_CD8 - Sample_GroupHP_1st_CD8,
                                     Third_Second_MS = Sample_GroupMS_3rd_CD8 - Sample_GroupMS_2nd_CD8,
                                     PP_Second_MS = Sample_GroupMS_PP_CD8 - Sample_GroupMS_2nd_CD8,
                                     PP_Third_MS = Sample_GroupMS_PP_CD8 - Sample_GroupMS_3rd_CD8,
                                     Third_Second_HP = Sample_GroupHP_3rd_CD8 - Sample_GroupHP_2nd_CD8,
                                     PP_Second_HP = Sample_GroupHP_PP_CD8 - Sample_GroupHP_2nd_CD8,
                                     PP_Third_HP = Sample_GroupHP_PP_CD8 - Sample_GroupHP_3rd_CD8,
                                     levels=colnames(design_CD8))

fit_CD8 <- lmFit(M_CD8, design_CD8, block=pheno_CD8$Individual, correlation=corfit$consensus)
fit_CD8 <- contrasts.fit(fit_CD8, contr_matrix)
fit_CD8 <- eBayes(fit_CD8)
summary(decideTests(fit_CD8))

# Gather DMPs CD8 ####
CD8_3rd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Third_HP")
CD8_3rd <- CD8_3rd[rownames(methyl_CD8),]
CD8_3rd_HP <- CD8_3rd[(CD8_3rd$P.Value<0.05) & abs(delta_base_beta_CD8_HP[,"3rd_1st"])>THS,] 

CD8_3rd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Third_MS")
CD8_3rd <- CD8_3rd[rownames(methyl_CD8),]
CD8_3rd_MS <- CD8_3rd[(CD8_3rd$P.Value<0.05) & abs(delta_base_beta_CD8_MS[,"3rd_1st"])>THS,] 

CD8_PP <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "PP_Third_HP")
CD8_PP <- CD8_PP[rownames(methyl_CD8),]
CD8_PP_HP <- CD8_PP[(CD8_PP$P.Value<0.05) & abs(delta_beta_CD8_HP[,"PP_3rd"])>THS,]

CD8_PP <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "PP_Third_MS")
CD8_PP <- CD8_PP[rownames(methyl_CD8),]
CD8_PP_MS <- CD8_PP[(CD8_PP$P.Value<0.05) & abs(delta_beta_CD8_MS[,"PP_3rd"])>THS,]

#2nd-1st
CD8_2nd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Second_HP")
CD8_2nd_1st <- CD8_2nd_1st[rownames(methyl_CD8),]
CD8_2nd_1st_HP <- CD8_2nd_1st[(CD8_2nd_1st$P.Value<0.05) & abs(delta_base_beta_CD8_HP[,"2nd_1st"])>THS,] 

CD8_2nd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Second_MS")
CD8_2nd_1st <- CD8_2nd_1st[rownames(methyl_CD8),]
CD8_2nd_1st_MS <- CD8_2nd_1st[(CD8_2nd_1st$P.Value<0.05) & abs(delta_base_beta_CD8_MS[,"2nd_1st"])>THS,] 

#3rd-2nd
CD8_3rd_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Third_Second_HP")
CD8_3rd_2nd <- CD8_3rd_2nd[rownames(methyl_CD8),]
CD8_3rd_2nd_HP <- CD8_3rd_2nd[(CD8_3rd_2nd$P.Value<0.05) & abs(delta_beta_CD8_HP[,"3rd_2nd"])>THS,] 

CD8_3rd_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Third_Second_MS")
CD8_3rd_2nd <- CD8_3rd_2nd[rownames(methyl_CD8),]
CD8_3rd_2nd_MS <- CD8_3rd_2nd[(CD8_3rd_2nd$P.Value<0.05) & abs(delta_beta_CD8_MS[,"3rd_2nd"])>THS,] 

#PP-2nd
CD8_PP_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "PP_Second_HP")
CD8_PP_2nd <- CD8_PP_2nd[rownames(methyl_CD8),]
CD8_PP_2nd_HP <- CD8_PP_2nd[(CD8_PP_2nd$P.Value<0.05) & abs(delta_beta_CD8_HP[,"PP_2nd"])>THS,] 

CD8_PP_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "PP_Second_MS")
CD8_PP_2nd <- CD8_PP_2nd[rownames(methyl_CD8),]
CD8_PP_2nd_MS <- CD8_PP_2nd[(CD8_PP_2nd$P.Value<0.05) & abs(delta_beta_CD8_MS[,"PP_2nd"])>THS,] 

# Gather all cpgs CD8 ####
CD8_3rd_HP_all <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Third_HP")
CD8_3rd_HP_all <- CD8_3rd_HP_all[rownames(methyl_CD8),]
CD8_3rd_HP_all$delta_beta <- delta_base_beta_CD8_HP[,"3rd_1st"]

CD8_3rd_MS_all <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Third_MS")
CD8_3rd_MS_all <- CD8_3rd_MS_all[rownames(methyl_CD8),]
CD8_3rd_MS_all$delta_beta <- delta_base_beta_CD8_MS[,"3rd_1st"]

CD8_PP_HP_all <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "PP_Third_HP")
CD8_PP_HP_all <- CD8_PP_HP_all[rownames(methyl_CD8),]
CD8_PP_HP_all$delta_beta <- delta_beta_CD8_HP[,"PP_3rd"]

CD8_PP_MS_all <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "PP_Third_MS")
CD8_PP_MS_all <- CD8_PP_MS_all[rownames(methyl_CD8),]
CD8_PP_MS_all$delta_beta <- delta_beta_CD8_MS[,"PP_3rd"]


#2nd-1st
CD8_2nd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Second_HP")
CD8_2nd_1st_HP_all <- CD8_2nd_1st[rownames(methyl_CD8),]

CD8_2nd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Second_MS")
CD8_2nd_1st_MS_all <- CD8_2nd_1st[rownames(methyl_CD8),]

#3rd-2nd
CD8_3rd_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Third_Second_HP")
CD8_3rd_2nd_HP_all <- CD8_3rd_2nd[rownames(methyl_CD8),]

CD8_3rd_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "Third_Second_MS")
CD8_3rd_2nd_MS_all <- CD8_3rd_2nd[rownames(methyl_CD8),]

#PP-2nd
CD8_PP_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "PP_Second_HP")
CD8_PP_2nd_HP_all <- CD8_PP_2nd[rownames(methyl_CD8),]

CD8_PP_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "inf", coef = "PP_Second_MS")
CD8_PP_2nd_MS_all <- CD8_PP_2nd[rownames(methyl_CD8),]

#MSvsHP
MSvsHP_cd8_1st <- topTable(fit_CD8,adjust.method="BH", p.value = Inf, number = Inf, coef="First")
MSvsHP_cd8_2nd <- topTable(fit_CD8,adjust.method="BH", p.value = Inf, number = Inf, coef="Second")
MSvsHP_cd8_3rd <- topTable(fit_CD8,adjust.method="BH", p.value = Inf, number = Inf, coef="Third")
MSvsHP_cd8_PP <- topTable(fit_CD8,adjust.method="BH", p.value = Inf, number = Inf, coef="PP")

MSvsHP_cd8_3rd_1st <- topTable(fit_CD8,adjust.method="BH", p.value = Inf, number = Inf, coef="Third_First")
MSvsHP_cd8_PP_3rd <- topTable(fit_CD8,adjust.method="BH", p.value = Inf, number = Inf, coef="PP_Third")

# Save ####
#CD4 DMP
saveRDS(CD4_3rd_HP, paste0(FOLDER_RDS, "DMR/CD4_3rd_HP_methyl.RDS"))
saveRDS(CD4_PP_HP, paste0(FOLDER_RDS, "DMR/CD4_PP_HP_methyl.RDS"))
saveRDS(CD4_3rd_MS, paste0(FOLDER_RDS, "DMR/CD4_3rd_MS_methyl.RDS"))
saveRDS(CD4_PP_MS, paste0(FOLDER_RDS, "DMR/CD4_PP_MS_methyl.RDS"))

saveRDS(CD4_2nd_1st_HP, paste0(FOLDER_RDS, "DMR/CD4_2nd_1st_HP.RDS"))
saveRDS(CD4_2nd_1st_MS, paste0(FOLDER_RDS, "DMR/CD4_2nd_1st_MS.RDS"))
saveRDS(CD4_3rd_2nd_HP, paste0(FOLDER_RDS, "DMR/CD4_3rd_2nd_HP.RDS"))
saveRDS(CD4_3rd_2nd_MS, paste0(FOLDER_RDS, "DMR/CD4_3rd_2nd_MS.RDS"))
saveRDS(CD4_PP_2nd_HP, paste0(FOLDER_RDS, "DMR/CD4_PP_2nd_HP.RDS"))
saveRDS(CD4_PP_2nd_MS, paste0(FOLDER_RDS, "DMR/CD4_PP_2nd_MS.RDS"))

#CD4 all
saveRDS(CD4_3rd_HP_all, paste0(FOLDER_RDS, "DMR/CD4_3rd_HP_methyl_all.RDS"))
saveRDS(CD4_PP_HP_all, paste0(FOLDER_RDS, "DMR/CD4_PP_HP_methyl_all.RDS"))
saveRDS(CD4_3rd_MS_all, paste0(FOLDER_RDS, "DMR/CD4_3rd_MS_methyl_all.RDS"))
saveRDS(CD4_PP_MS_all, paste0(FOLDER_RDS, "DMR/CD4_PP_MS_methyl_all.RDS"))

saveRDS(CD4_2nd_1st_HP_all, paste0(FOLDER_RDS, "DMR/CD4_2nd_1st_HP_all.RDS"))
saveRDS(CD4_2nd_1st_MS_all, paste0(FOLDER_RDS, "DMR/CD4_2nd_1st_MS_all.RDS"))
saveRDS(CD4_3rd_2nd_HP_all, paste0(FOLDER_RDS, "DMR/CD4_3rd_2nd_HP_all.RDS"))
saveRDS(CD4_3rd_2nd_MS_all, paste0(FOLDER_RDS, "DMR/CD4_3rd_2nd_MS_all.RDS"))
saveRDS(CD4_PP_2nd_HP_all, paste0(FOLDER_RDS, "DMR/CD4_PP_2nd_HP_all.RDS"))
saveRDS(CD4_PP_2nd_MS_all, paste0(FOLDER_RDS, "DMR/CD4_PP_2nd_MS_all.RDS"))

#CD8 DMP
saveRDS(CD8_3rd_HP, paste0(FOLDER_RDS, "DMR/CD8_3rd_HP_methyl.RDS"))
saveRDS(CD8_PP_HP, paste0(FOLDER_RDS, "DMR/CD8_PP_HP_methyl.RDS"))
saveRDS(CD8_3rd_MS, paste0(FOLDER_RDS, "DMR/CD8_3rd_MS_methyl.RDS"))
saveRDS(CD8_PP_MS, paste0(FOLDER_RDS, "DMR/CD8_PP_MS_methyl.RDS"))

saveRDS(CD8_2nd_1st_HP, paste0(FOLDER_RDS, "DMR/CD8_2nd_1st_HP.RDS"))
saveRDS(CD8_2nd_1st_MS, paste0(FOLDER_RDS, "DMR/CD8_2nd_1st_MS.RDS"))
saveRDS(CD8_3rd_2nd_HP, paste0(FOLDER_RDS, "DMR/CD8_3rd_2nd_HP.RDS"))
saveRDS(CD8_3rd_2nd_MS, paste0(FOLDER_RDS, "DMR/CD8_3rd_2nd_MS.RDS"))
saveRDS(CD8_PP_2nd_HP, paste0(FOLDER_RDS, "DMR/CD8_PP_2nd_HP.RDS"))
saveRDS(CD8_PP_2nd_MS, paste0(FOLDER_RDS, "DMR/CD8_PP_2nd_MS.RDS"))

#CD8 all
saveRDS(CD8_3rd_HP_all, paste0(FOLDER_RDS, "DMR/CD8_3rd_HP_methyl_all.RDS"))
saveRDS(CD8_PP_HP_all, paste0(FOLDER_RDS, "DMR/CD8_PP_HP_methyl_all.RDS"))
saveRDS(CD8_3rd_MS_all, paste0(FOLDER_RDS, "DMR/CD8_3rd_MS_methyl_all.RDS"))
saveRDS(CD8_PP_MS_all, paste0(FOLDER_RDS, "DMR/CD8_PP_MS_methyl_all.RDS"))

saveRDS(CD8_2nd_1st_HP_all, paste0(FOLDER_RDS, "DMR/CD8_2nd_1st_HP_all.RDS"))
saveRDS(CD8_2nd_1st_MS_all, paste0(FOLDER_RDS, "DMR/CD8_2nd_1st_MS_all.RDS"))
saveRDS(CD8_3rd_2nd_HP_all, paste0(FOLDER_RDS, "DMR/CD8_3rd_2nd_HP_all.RDS"))
saveRDS(CD8_3rd_2nd_MS_all, paste0(FOLDER_RDS, "DMR/CD8_3rd_2nd_MS_all.RDS"))
saveRDS(CD8_PP_2nd_HP_all, paste0(FOLDER_RDS, "DMR/CD8_PP_2nd_HP_all.RDS"))
saveRDS(CD8_PP_2nd_MS_all, paste0(FOLDER_RDS, "DMR/CD8_PP_2nd_MS_all.RDS"))

#MS vs HP
saveRDS(MSvsHP_cd4_1st, paste0(FOLDER_RDS, "DMR/MSvsHP_CD4_1st.RDS"))
saveRDS(MSvsHP_cd4_2nd, paste0(FOLDER_RDS, "DMR/MSvsHP_CD4_2nd.RDS"))
saveRDS(MSvsHP_cd4_3rd, paste0(FOLDER_RDS, "DMR/MSvsHP_CD4_3rd.RDS"))
saveRDS(MSvsHP_cd4_PP, paste0(FOLDER_RDS, "DMR/MSvsHP_CD4_PP.RDS"))
saveRDS(MSvsHP_cd4_3rd_1st, paste0(FOLDER_RDS, "DMR/MSvsHP_CD4_3rd_1st.RDS"))
saveRDS(MSvsHP_cd4_PP_3rd, paste0(FOLDER_RDS, "DMR/MSvsHP_CD4_PP_3rd.RDS"))

saveRDS(MSvsHP_cd8_1st, paste0(FOLDER_RDS, "DMR/MSvsHP_CD8_1st.RDS"))
saveRDS(MSvsHP_cd8_2nd, paste0(FOLDER_RDS, "DMR/MSvsHP_CD8_2nd.RDS"))
saveRDS(MSvsHP_cd8_3rd, paste0(FOLDER_RDS, "DMR/MSvsHP_CD8_3rd.RDS"))
saveRDS(MSvsHP_cd8_PP, paste0(FOLDER_RDS, "DMR/MSvsHP_CD8_PP.RDS"))
saveRDS(MSvsHP_cd8_3rd_1st, paste0(FOLDER_RDS, "DMR/MSvsHP_CD8_3rd_1st.RDS"))
saveRDS(MSvsHP_cd8_PP_3rd, paste0(FOLDER_RDS, "DMR/MSvsHP_CD8_PP_3rd.RDS"))




