
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
main_dir <-"C:/Users/albze08/Desktop/phd/P4/methylation"
FOLDER_RDS <- "RDS_files/"
FOLDER_FIGURES <- "figures/"
setwd("C:/Users/albze08/Desktop/phd/P4/methylation")

# My functions ####
source('probe_to_gene.R')
source('hyge_test.R')

hyge_test_wrap <- function(subset_1, full_1){ #wrapper for cpg list
  
  genes_sub <- probe_to_gene( subset_1 )[,2] 
  genes_sub <- genes_sub[genes_sub!=""]
  
  genes_universe <- probe_to_gene( full_1 )[,2] 
  genes_universe <- genes_universe[genes_universe!=""]
  
  out <- hyge_test(genes_sub,genes_universe)
  return(out)
}

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
  
  delta_beta <- data.frame( cbind(delta_2nd, delta_3rd, delta_PP) )
  colnames(delta_beta) <- c("2nd_1st", "3rd_2nd", "PP_3rd")
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

delta_beta_CD4 <- delta_beta_CD4_MS - delta_beta_CD4_HP
delta_beta_CD8 <- delta_beta_CD8_MS - delta_beta_CD8_HP

delta_base_beta_CD4 <- delta_base_beta_CD4_MS - delta_base_beta_CD4_HP
delta_base_beta_CD8 <- delta_base_beta_CD8_MS - delta_base_beta_CD8_HP

# Gather cpgs with significant delta beta in any comparison ####
THS <- 0.05
cpg_delta_CD4_HP <- rownames(delta_beta_CD4_HP)[apply(abs(delta_beta_CD4_HP)>THS, 1, any)]
cpg_delta_CD4_MS <- rownames(delta_beta_CD4_MS)[apply(abs(delta_beta_CD4_MS)>THS, 1, any)]
cpg_delta_CD8_HP <- rownames(delta_beta_CD8_HP)[apply(abs(delta_beta_CD8_HP)>THS, 1, any)]
cpg_delta_CD8_MS <- rownames(delta_beta_CD8_MS)[apply(abs(delta_beta_CD8_MS)>THS, 1, any)]

# Limma on CD4 ####
M_CD4 <- beta2m(methyl_CD4)

pheno_CD4$Sample_Group <- gsub("MS_BP", "AA", pheno_CD4$Sample_Group) #to make sure it is discarded below
design_CD4 <- model.matrix(~ Sample_Group + Age + Proportions_Memory, pheno_CD4 )
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
                                     Third_Second = (Sample_GroupMS_3rd_CD4 - Sample_GroupMS_2nd_CD4) - (Sample_GroupHP_3rd_CD4 - Sample_GroupHP_2nd_CD4),
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
                                     Third_Second_HP = Sample_GroupHP_3rd_CD4 - Sample_GroupHP_2nd_CD4,
                                     PP_Third_HP = Sample_GroupHP_PP_CD4 - Sample_GroupHP_3rd_CD4,
                                     levels=colnames(design_CD4))

fit_CD4 <- lmFit(M_CD4, design_CD4, block=pheno_CD4$Individual, correlation=corfit$consensus)
fit_CD4 <- contrasts.fit(fit_CD4, contr_matrix)
fit_CD4 <- eBayes(fit_CD4)
summary(decideTests(fit_CD4))


# Gather Changes
#HP
CD4_HP_2nd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Second_HP")
CD4_HP_2nd_1st <- CD4_HP_2nd_1st[rownames(methyl_CD4),]

CD4_HP_3rd_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second_HP")
CD4_HP_3rd_2nd <- CD4_HP_3rd_2nd[rownames(methyl_CD4),]

CD4_HP_3rd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_HP")
CD4_HP_3rd_1st <- CD4_HP_3rd_1st[rownames(methyl_CD4),]

CD4_HP_PP_3rd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third_HP")
CD4_HP_PP_3rd <- CD4_HP_PP_3rd[rownames(methyl_CD4),]

#MS
CD4_MS_2nd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Second_MS")
CD4_MS_2nd_1st <- CD4_MS_2nd_1st[rownames(methyl_CD4),]

CD4_MS_3rd_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second_MS")
CD4_MS_3rd_2nd <- CD4_MS_3rd_2nd[rownames(methyl_CD4),]

CD4_MS_3rd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_MS")
CD4_MS_3rd_1st <- CD4_MS_3rd_1st[rownames(methyl_CD4),]

CD4_MS_PP_3rd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third_MS")
CD4_MS_PP_3rd <- CD4_MS_PP_3rd[rownames(methyl_CD4),]

#MSvsHP
CD4_MSvsHP_2nd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Second_First")
CD4_MSvsHP_2nd_1st <- CD4_MSvsHP_2nd_1st[rownames(methyl_CD4),]

CD4_MSvsHP_3rd_2nd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second")
CD4_MSvsHP_3rd_2nd <- CD4_MSvsHP_3rd_2nd[rownames(methyl_CD4),]

CD4_MSvsHP_3rd_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "Third_First")
CD4_MSvsHP_3rd_1st <- CD4_MSvsHP_2nd_1st[rownames(methyl_CD4),]

CD4_MSvsHP_PP_3rd <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third")
CD4_MSvsHP_PP_3rd <- CD4_MSvsHP_PP_3rd[rownames(methyl_CD4),]

CD4_MSvsHP_PP_1st <- topTable(fit_CD4, p.value=Inf, lfc=0, number = "Inf", coef = "PP_First")
CD4_MSvsHP_PP_1st <- CD4_MSvsHP_PP_1st[rownames(methyl_CD4),]


# Correlation plots CD4 ####
pdf("figures_manus/Correlations_CD4.pdf", width=10)

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
M_CD8 <- beta2m(methyl_CD8)

pheno_CD8$Sample_Group <- gsub("MS_BP", "AA", pheno_CD8$Sample_Group) #to make sure it is discarded below
design_CD8 <- model.matrix(~ Sample_Group + Age + Proportions_Memory, pheno_CD8 )
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
                                     Third_Second = (Sample_GroupMS_3rd_CD8 - Sample_GroupMS_2nd_CD8) - (Sample_GroupHP_3rd_CD8 - Sample_GroupHP_2nd_CD8),
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
                                     Third_Second_HP = Sample_GroupHP_3rd_CD8 - Sample_GroupHP_2nd_CD8,
                                     PP_Third_HP = Sample_GroupHP_PP_CD8 - Sample_GroupHP_3rd_CD8,
                                     levels=colnames(design_CD8))

fit_CD8 <- lmFit(M_CD8, design_CD8, block=pheno_CD8$Individual, correlation=corfit$consensus)
fit_CD8 <- contrasts.fit(fit_CD8, contr_matrix)
fit_CD8 <- eBayes(fit_CD8)
summary(decideTests(fit_CD8))


# Gather Changes
#HP
CD8_HP_2nd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Second_HP")
CD8_HP_2nd_1st <- CD8_HP_2nd_1st[rownames(methyl_CD8),]

CD8_HP_3rd_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second_HP")
CD8_HP_3rd_2nd <- CD8_HP_3rd_2nd[rownames(methyl_CD8),]

CD8_HP_3rd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_HP")
CD8_HP_3rd_1st <- CD8_HP_3rd_1st[rownames(methyl_CD8),]

CD8_HP_PP_3rd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third_HP")
CD8_HP_PP_3rd <- CD8_HP_PP_3rd[rownames(methyl_CD8),]

#MS
CD8_MS_2nd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Second_MS")
CD8_MS_2nd_1st <- CD8_MS_2nd_1st[rownames(methyl_CD8),]

CD8_MS_3rd_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second_MS")
CD8_MS_3rd_2nd <- CD8_MS_3rd_2nd[rownames(methyl_CD8),]

CD8_MS_3rd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_MS")
CD8_MS_3rd_1st <- CD8_MS_3rd_1st[rownames(methyl_CD8),]

CD8_MS_PP_3rd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third_MS")
CD8_MS_PP_3rd <- CD8_MS_PP_3rd[rownames(methyl_CD8),]

#MSvsHP
CD8_MSvsHP_2nd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Second_First")
CD8_MSvsHP_2nd_1st <- CD8_MSvsHP_2nd_1st[rownames(methyl_CD8),]

CD8_MSvsHP_3rd_2nd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_Second")
CD8_MSvsHP_3rd_2nd <- CD8_MSvsHP_3rd_2nd[rownames(methyl_CD8),]

CD8_MSvsHP_3rd_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "Third_First")
CD8_MSvsHP_3rd_1st <- CD8_MSvsHP_2nd_1st[rownames(methyl_CD8),]

CD8_MSvsHP_PP_3rd <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "PP_Third")
CD8_MSvsHP_PP_3rd <- CD8_MSvsHP_PP_3rd[rownames(methyl_CD8),]

CD8_MSvsHP_PP_1st <- topTable(fit_CD8, p.value=Inf, lfc=0, number = "Inf", coef = "PP_First")
CD8_MSvsHP_PP_1st <- CD8_MSvsHP_PP_1st[rownames(methyl_CD8),]

# Correlation plots CD8 ####
pdf("figures_manus/Correlations_CD8.pdf", width=10)

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
pdf("figures_manus/Correlations_MSvsHP.pdf", width=10)

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
CD4_DEG_2nd_1st <- sum((CD4_MSvsHP_2nd_1st$P.Value<0.05) & abs(delta_base_beta_CD4$`2nd_1st`>0.05))
CD4_DEG_3rd_1st <- sum((CD4_MSvsHP_3rd_1st$P.Value<0.05) & abs(delta_base_beta_CD4$`3rd_1st`>0.05))
CD4_DEG_PP_1st <- sum((CD4_MSvsHP_PP_1st$P.Value<0.05) & abs(delta_base_beta_CD4$`PP_1st`>0.05))


CD8_DEG_2nd_1st <- sum((CD8_MSvsHP_2nd_1st$P.Value<0.05) & abs(delta_base_beta_CD8$`2nd_1st`>0.05))
CD8_DEG_3rd_1st <- sum((CD8_MSvsHP_3rd_1st$P.Value<0.05) & abs(delta_base_beta_CD8$`3rd_1st`>0.05))
CD8_DEG_PP_1st <- sum((CD8_MSvsHP_PP_1st$P.Value<0.05) & abs(delta_base_beta_CD8$`PP_1st`>0.05))


df <- rbind( data.frame(n=c(CD4_DEG_2nd_1st, CD4_DEG_3rd_1st, CD4_DEG_PP_1st), time = 1:3, group = rep("CD4",3)),
             data.frame(n=c(CD8_DEG_2nd_1st, CD8_DEG_3rd_1st, CD8_DEG_PP_1st), time =1:3, group = rep("CD8",3)))

dev.off()
par(mar=c(1,1,1,1))

pdf("figures_manus/n_DMPs.pdf", width=5, height=5)
ggplot(df, aes(x=time, y=n, col=group, cex.axis=2)) + geom_line(size=2) + geom_point(size=5) +
  scale_x_continuous(breaks=c(1, 2, 3), labels=c("2nd-1st", "3rd-1st", "PP-1st")) + ylab("n. DMPs") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 
dev.off()








