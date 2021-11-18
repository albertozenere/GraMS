#Here we look at the comparison bewteen MS and HP, at each time point.
#We look at the overlap between (nominally) DEG, DMG, and MS-associated genes.


###################
# Limma analysis
###################

# Created and modified by Alberto Zenere, 2021-04-11. Based on "RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR" by Law et. al.

# 0.1 Setup
# 0.2 Split CD4 and CD8
# 0.3 My functions
# 1.1 Limma on CD4
# 1.1.1 Comparisons between time points in MS
# 1.1.2 Comparisons between time points in HP
# 1.1.3 Comparisons MS and HP
# 1.2 Limma on CD8
# 1.2.1 Comparisons between time points in MS
# 1.2.2 Comparisons between time points in HP
# 1.2.3 Comparisons MS and HP 
# 2.1 Using CD4 and CD8 together
# 2.2 Add activation

rm(list=ls()) # remove all entries in the global environment 

################
## 0.1 Setup  ##
################

#Set the virtual environment folder as R library
pathlist <- .libPaths();
newpath <-  "C:/Users/albze08/.conda/envs/graMS/Lib/R/library"
.libPaths(newpath)
setwd('C:/Users/albze08/Desktop/P4/RNAseq')

#Load packages
pack_R <- c("sva","tidyverse","readxl","assertthat", "gridExtra",
            "edgeR","patchwork","ggfortify","RColorBrewer","gam",
            "mgcv", "tidyr","factoextra")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}
set.seed(541)

#Set important folders
RDS_folder <- "RDS/"
data_folder <- "data/"
figure_folder <- "figures/"
results_folder <- "results/"

par(mar=c(1,1,1,1))

##############
## 0.2 Load ##
##############

#CD4
count_cd4 <- readRDS(file = paste0(RDS_folder, "count_cd4_tmm.RDS"))
metadata_cd4 <- readRDS(file = paste0(RDS_folder, "metadata_cd4.RDS"))
stopifnot(colnames(count_cd4) == metadata_cd4$NGI_ID)

#CD8
count_cd8 <- readRDS(file = paste0(RDS_folder, "count_cd8_tmm.RDS"))
metadata_cd8 <- readRDS(file = paste0(RDS_folder, "metadata_cd8.RDS"))
stopifnot(colnames(count_cd8) == metadata_cd8$NGI_ID)

#MS relevant genes
MS_new_genes <- read.csv(file = "data/complete_MS_annotation_p10e-6.csv")
MS_new_genes <- unique(MS_new_genes$SYMBOL)

#other MS relevant genes
MS_old_genes <- read_excel("data/annotated_MS_SNPs_proc.xlsx")
idx <- MS_old_genes[,2]<0.00001
MS_old_genes <- unique(MS_old_genes$nearest_gene[idx])

#Genes from methylation data
MS1stvsHP1st_CD4_DMGs <- read.csv("data/Methylation_DMGs/MS1stvsHP1st_CD4_DMGs.txt")
MS2ndvsHP2nd_CD4_DMGs <- read.csv("data/Methylation_DMGs/MS2ndvsHP2nd_CD4_DMGs.txt")
MS3rdvsHP3rd_CD4_DMGs <- read.csv("data/Methylation_DMGs/MS3rdvsHP3rd_CD4_DMGs.txt")
MSPPvsHPPP_CD4_DMGs <- read.csv("data/Methylation_DMGs/MSPPvsHPPP_CD4_DMGs.txt")

MS1stvsHP1st_CD8_DMGs <- read.csv("data/Methylation_DMGs/MS1stvsHP1st_CD8_DMGs.txt")
MS2ndvsHP2nd_CD8_DMGs <- read.csv("data/Methylation_DMGs/MS2ndvsHP2nd_CD8_DMGs.txt")
MS3rdvsHP3rd_CD8_DMGs <- read.csv("data/Methylation_DMGs/MS3rdvsHP3rd_CD8_DMGs.txt")
MSPPvsHPPP_CD8_DMGs <- read.csv("data/Methylation_DMGs/MSPPvsHPPP_CD8_DMGs.txt")


######################
## 0.3 My functions 
######################

plot_pvalues <- function(fit){
  
  #plot p-values of each covariate
  covar <- colnames(fit$coefficients)
  Ncovar <- length(covar)
  
  n_col <- 3 #number of rows and columns
  n_row <- ceiling(Ncovar/n_col)
  n_slots <- n_row*n_col
  
  lay <- rep(0,n_slots) #create layout
  lay[1:Ncovar] <- 1:Ncovar
  layout(matrix(lay, ncol=3))
  
  for (n in 1:Ncovar){ #plot distributions of p-values
    if (!any(is.na(fit$p.value[,n]))){
      hist(fit$p.value[,n], xlab="", ylab="", main = covar[n],breaks = seq(from=0, to =1, by=0.05))
    }
  }
}

#This function does an hypergeometric test between the genes in subset_1 and those in the pre-defined lists MS_old_genes and MS_new_genes.
#The second list contains all the genes that have been measured.
hyge_test <- function(subset_1, full_1){
  
  ## Check the input
  
  #convert to character list
  if (!class(subset_1)=="character"){
    stopifnot(ncol(subset_1)==1) #subset_1 needs to have one column only
    stopifnot(class(subset_1[,1])=="character") #full_1 needs contain characters
    subset_1 <- as.character(subset_1)
  }
  if (!class(full_1)=="character"){
    stopifnot(ncol(full_1)==1) #full_1 needs to have one column only
    stopifnot(class(full_1[,1])=="character") #full_1 needs contain characters
    full_1 <- as.character(full_1)
  }
  
  #check subset_1 is a subset of full_1
  stopifnot(all(is.element(subset_1, full_1))) #all the genes in subset_1 need to be in full_1
  stopifnot(length(subset_1)<length(full_1)) #subset_1 needs to be a subset of full_1
  
  stopifnot(length(subset_1)>0) #subset_1 contains 0 zero. Check you code
  stopifnot(length(full_1)>0) #full_1 contains 0 zero. Check you code
  
  #check MS_old_genes and MS_new_genes
  stopifnot(class(MS_old_genes)=="character") #MS_old_genes should be a character array
  stopifnot(class(MS_new_genes)=="character") #MS_new_genes should be a character array
  
  stopifnot(length(MS_old_genes)>0) #MS_old_genes contains 0 zero. Check you code
  stopifnot(length(MS_new_genes)>0) #MS_new_genes contains 0 zero. Check you code
  
  ## Do the hypergeometric test on old list
  q <- sum(is.element(MS_old_genes, subset_1))
  m <- sum(is.element(MS_old_genes, full_1))
  n <- length(full_1) - m
  k <- length(subset_1)
  pval_old <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  odds_old <- (q/k) / ( m/(n+m) ) #ratio in subset_1 compared to the ratio in full_1
  n_genes_old <- q #n. genes
  genes_list_old <- list(intersect(MS_old_genes, subset_1))
  
  if (m==0){
    warning("There is no overlap between all the genes measured and MS_old_genes. Check the data.")
  }
  
  ## Do the hypergeometric test on new list
  q <- sum(is.element(MS_new_genes, subset_1))
  m <- sum(is.element(MS_new_genes, full_1))
  n <- length(full_1) - m
  k <- length(subset_1)
  pval_new <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  odds_new <- (q/k) / ( m/(n+m) ) #ratio in subset_1 compared to the ratio in full_1
  n_genes_new <- q 
  genes_list_new <- list(intersect(MS_new_genes, subset_1))
  
  if (m==0){
    warning("There is no overlap between all the genes measured and MS_new_genes. Check the data.")
  }
  
  p_value <- c(pval_old, pval_new)
  odds_ratio <- c(odds_old, odds_new)
  n_genes <- c(n_genes_old, n_genes_new)
  gene_list <- c(genes_list_old, genes_list_new)
  
  out <- data.frame(p_value, odds_ratio, n_genes)
  out$gene_list <- gene_list 
  
  rownames(out) <- c("old_list", "new_list")
  
  return(out)
}


## Comparison with genes from methylation data
hyge_test_sub <- function(subset_1, full_1, independent_list){
  
  ## Check the input
  #convert to character list
  if (!class(subset_1)=="character"){
    stopifnot(ncol(subset_1)==1) #subset_1 needs to have one column only
    stopifnot(class(subset_1[,1])=="character") #full_1 needs contain characters
    subset_1 <- as.character(subset_1)
  }
  if (!class(full_1)=="character"){
    stopifnot(ncol(full_1)==1) #full_1 needs to have one column only
    stopifnot(class(full_1[,1])=="character") #full_1 needs contain characters
    full_1 <- as.character(full_1)
  }
  
  #check subset_1 is a subset of full_1
  stopifnot(all(is.element(subset_1, full_1))) #all the genes in subset_1 need to be in full_1
  stopifnot(length(subset_1)<length(full_1)) #subset_1 needs to be a subset of full_1
  
  stopifnot(length(subset_1)>0) #subset_1 contains 0 zero. Check you code
  stopifnot(length(full_1)>0) #full_1 contains 0 zero. Check you code
  
  stopifnot(class(independent_list)=="character") #independent_list should be a character array
  stopifnot(length(independent_list)>0) #independent_list contains 0 zero. Check you code
  
  #hypergeometric test
  q <- sum(is.element(independent_list, subset_1))
  m <- sum(is.element(independent_list, full_1))
  n <- length(full_1) - m
  k <- length(subset_1)
  p_value <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  odds <- (q/k) / ( m/(n+m) ) #ratio in subset_1 compared to the ratio in full_1
  n_genes <- q 
  gene_list <- list(intersect(independent_list, subset_1))
  
  out <- data.frame(p_value, odds, n_genes)
  out$gene_list <- gene_list 
  
  return(out)
  
}

######################
## 1.1.a Prepare CD4
######################

#remove HC
is_HC <- which(is.element(metadata_cd4$Disease, "HC"))
count_cd4 <- count_cd4[, -is_HC]
metadata_cd4 <- metadata_cd4[-is_HC, ]

is_resting <- is.element(metadata_cd4$State, "Resting")

count_cd4 <- count_cd4[, is_resting]
metadata_cd4 <- metadata_cd4[is_resting, ]

stopifnot(colnames(count_cd4) == metadata_cd4$NGI_ID)
stopifnot(ncol(count_cd4) == nrow(metadata_cd4))

metadata_cd4$Sample_Group <- gsub(" ", "", metadata_cd4$Sample_Group)

#####################
## 1.1.b Prepare CD8
#####################

#remove HC
is_HC <- which(is.element(metadata_cd8$Disease, "HC"))
count_cd8 <- count_cd8[, -is_HC]
metadata_cd8 <- metadata_cd8[-is_HC, ]

is_resting <- is.element(metadata_cd8$State, "Resting")

count_cd8 <- count_cd8[, is_resting]
metadata_cd8 <- metadata_cd8[is_resting, ]

metadata_cd8$Sample_Group <- gsub(" ", "", metadata_cd8$Sample_Group)

stopifnot(colnames(count_cd8) == metadata_cd8$NGI_ID)
stopifnot(ncol(count_cd8) == nrow(metadata_cd8))

#####################
## 1.2.a DEG in CD4
#####################

#comparison for each time-point
design <- model.matrix(~ 0 + Sample_Group, metadata_cd4 )
corfit <- duplicateCorrelation(count_cd4, design, block=metadata_cd4$Individual) #block the individual effect
fit <- lmFit(count_cd4, design, block=metadata_cd4$Individual, correlation=corfit$consensus)
contr_matrix <- makeContrasts(
  MS1stvsHP1st = Sample_GroupMS_1st_CD4 - Sample_GroupHP_1st_CD4 ,
  MS2ndvsHP2nd = Sample_GroupMS_2nd_CD4 - Sample_GroupHP_2nd_CD4 ,
  MS3rdvsHP3rd = Sample_GroupMS_3rd_CD4 - Sample_GroupHP_3rd_CD4 ,
  MSPPvsHPPP = Sample_GroupMS_PP_CD4 - Sample_GroupHP_PP_CD4 ,
  levels=colnames(design)
)

fit <- contrasts.fit(fit, contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

plot_pvalues(fit)

#extract nominally significant genes at each time-point
DEG_1st <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS1stvsHP1st")
DEG_2nd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS2ndvsHP2nd")
DEG_3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS3rdvsHP3rd")
DEG_PP <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MSPPvsHPPP")

#Save files
saveRDS(DEG_1st, paste0(RDS_folder, "DEG/DEG_1st_CD4.RDS"))
saveRDS(DEG_2nd, paste0(RDS_folder, "DEG/DEG_2nd_CD4.RDS"))
saveRDS(DEG_3rd, paste0(RDS_folder, "DEG/DEG_3rd_CD4.RDS"))
saveRDS(DEG_PP, paste0(RDS_folder, "DEG/DEG_PP_CD4.RDS"))

############################
## 1.2.b DEG and MS in CD4
############################

#extract nominally significant genes at each time-point
DEG_MS_1st <- hyge_test(rownames(DEG_1st), rownames(count_cd4))
DEG_MS_2nd <- hyge_test(rownames(DEG_2nd), rownames(count_cd4))
DEG_MS_3rd <- hyge_test(rownames(DEG_3rd), rownames(count_cd4))
DEG_MS_PP <- hyge_test(rownames(DEG_PP), rownames(count_cd4))

#Save files
saveRDS(DEG_1st, paste0(RDS_folder, "DEG_MS/DEG_MS_1st.RDS"))
saveRDS(DEG_2nd, paste0(RDS_folder, "DEG_MS/DEG_MS_2nd.RDS"))
saveRDS(DEG_3rd, paste0(RDS_folder, "DEG_MS/DEG_MS_3rd.RDS"))
saveRDS(DEG_PP, paste0(RDS_folder, "DEG_MS/DEG_MS_PP.RDS"))

#############################
## 1.2.c DEG and DMG in CD4
#############################

#Overlap with genes from methylation data, time-point by time-point
DEG_DMG_1st_CD4 <- hyge_test_sub(rownames(DEG_1st), rownames(count_cd4), MS1stvsHP1st_CD4_DMGs$Gene )
DEG_DMG_2nd_CD4 <- hyge_test_sub(rownames(DEG_2nd), rownames(count_cd4), MS2ndvsHP2nd_CD4_DMGs$Gene )
DEG_DMG_3rd_CD4 <- hyge_test_sub(rownames(DEG_3rd), rownames(count_cd4), MS3rdvsHP3rd_CD4_DMGs$Gene )
DEG_DMG_PP_CD4 <- hyge_test_sub(rownames(DEG_PP), rownames(count_cd4), MSPPvsHPPP_CD4_DMGs$Gene )

#Save files
saveRDS(DEG_DMG_1st_CD4, paste0(RDS_folder, "DEG_DMG/DEG_DMG_1st_CD4.RDS"))
saveRDS(DEG_DMG_2nd_CD4, paste0(RDS_folder, "DEG_DMG/DEG_DMG_2nd_CD4.RDS"))
saveRDS(DEG_DMG_3rd_CD4, paste0(RDS_folder, "DEG_DMG/DEG_DMG_3rd_CD4.RDS"))
saveRDS(DEG_DMG_PP_CD4, paste0(RDS_folder, "DEG_DMG/DEG_DMG_PP_CD4.RDS"))

#Overlap with genes from methylation data, all time-points merged
DMG_CD4 <- unique(c(MS1stvsHP1st_CD4_DMGs$Gene, MS2ndvsHP2nd_CD4_DMGs$Gene,
                    MS3rdvsHP3rd_CD4_DMGs$Gene, MSPPvsHPPP_CD4_DMGs$Gene))

DEG_CD4 <- unique(c(rownames(DEG_1st), rownames(DEG_2nd),
                    rownames(DEG_3rd), rownames(DEG_PP)))

DEG_DMG_CD4 <- hyge_test_sub(DEG_CD4, rownames(count_cd4), DMG_CD4 )

#Save file
saveRDS(DEG_DMG_CD4, paste0(RDS_folder, "DEG_DMG_CD4.RDS"))

####################################
## 1.2.d DEG and DMG and MS in CD4
####################################

DEG_DMG_1st_MS_CD4 <- hyge_test(DEG_DMG_1st_CD4$gene_list[[1]], rownames(count_cd4) )
DEG_DMG_2nd_MS_CD4 <- hyge_test(DEG_DMG_2nd_CD4$gene_list[[1]], rownames(count_cd4) )
DEG_DMG_3rd_MS_CD4 <- hyge_test(DEG_DMG_3rd_CD4$gene_list[[1]], rownames(count_cd4) )
DEG_DMG_PP_MS_CD4 <- hyge_test(DEG_DMG_PP_CD4$gene_list[[1]], rownames(count_cd4) )

#Save files
saveRDS(DEG_DMG_1st_MS_CD4, paste0(RDS_folder, "DEG_DMG_MS/DEG_DMG_1st_MS_CD4.RDS"))
saveRDS(DEG_DMG_2nd_MS_CD4, paste0(RDS_folder, "DEG_DMG_MS/DEG_DMG_2nd_MS_CD4.RDS"))
saveRDS(DEG_DMG_3rd_MS_CD4, paste0(RDS_folder, "DEG_DMG_MS/DEG_DMG_3rd_MS_CD4.RDS"))
saveRDS(DEG_DMG_PP_MS_CD4, paste0(RDS_folder, "DEG_DMG_MS/DEG_DMG_PP_MS_CD4.RDS"))

#Overlap with genes from methylation data, all time-points merged
DEG_DMG_CD4 <- unique(c(DEG_DMG_1st_MS_CD4$gene_list[[1]], DEG_DMG_2nd_MS_CD4$gene_list[[1]],
                        DEG_DMG_3rd_MS_CD4$gene_list[[1]], DEG_DMG_PP_MS_CD4$gene_list[[1]]))

DEG_DMG_MS_CD4 <- hyge_test(DEG_DMG_CD4, rownames(count_cd4) )

#Save file
saveRDS(DEG_DMG_CD4, paste0(RDS_folder, "DEG_DMG_MS_CD4.RDS"))

#####################
## 1.3.a DEG in CD8
#####################

#comparison for each time-point
design <- model.matrix(~ 0 + Sample_Group, metadata_cd8 )
corfit <- duplicateCorrelation(count_cd8, design, block=metadata_cd8$Individual) #block the individual effect
fit <- lmFit(count_cd8, design, block=metadata_cd8$Individual, correlation=corfit$consensus)
contr_matrix <- makeContrasts(
  MS1stvsHP1st = Sample_GroupMS_1st_CD8 - Sample_GroupHP_1st_CD8 ,
  MS2ndvsHP2nd = Sample_GroupMS_2nd_CD8 - Sample_GroupHP_2nd_CD8 ,
  MS3rdvsHP3rd = Sample_GroupMS_3rd_CD8 - Sample_GroupHP_3rd_CD8 ,
  MSPPvsHPPP = Sample_GroupMS_PP_CD8 - Sample_GroupHP_PP_CD8 ,
  levels=colnames(design)
)

fit <- contrasts.fit(fit, contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

plot_pvalues(fit)

#extract nominally significant genes at each time-point
DEG_1st <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS1stvsHP1st")
DEG_2nd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS2ndvsHP2nd")
DEG_3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS3rdvsHP3rd")
DEG_PP <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MSPPvsHPPP")

#Save files
saveRDS(DEG_1st, paste0(RDS_folder, "DEG/DEG_1st_CD8.RDS"))
saveRDS(DEG_2nd, paste0(RDS_folder, "DEG/DEG_2nd_CD8.RDS"))
saveRDS(DEG_3rd, paste0(RDS_folder, "DEG/DEG_3rd_CD8.RDS"))
saveRDS(DEG_PP, paste0(RDS_folder, "DEG/DEG_PP_CD8.RDS"))

#############################
## 1.3.b DEG and MS in CD8
#############################

#extract nominally significant genes at each time-point
DEG_MS_1st <- hyge_test(rownames(DEG_1st), rownames(count_cd8))
DEG_MS_2nd <- hyge_test(rownames(DEG_2nd), rownames(count_cd8))
DEG_MS_3rd <- hyge_test(rownames(DEG_3rd), rownames(count_cd8))
DEG_MS_PP <- hyge_test(rownames(DEG_PP), rownames(count_cd8))

#Save files
saveRDS(DEG_MS_1st, paste0(RDS_folder, "DEG_MS/DEG_MS_1st_CD8.RDS"))
saveRDS(DEG_MS_2nd, paste0(RDS_folder, "DEG_MS/DEG_MS_2nd_CD8.RDS"))
saveRDS(DEG_MS_3rd, paste0(RDS_folder, "DEG_MS/DEG_MS_3rd_CD8.RDS"))
saveRDS(DEG_MS_PP, paste0(RDS_folder, "DEG_MS/DEG_MS_PP_CD8.RDS"))

##############################
## 1.3.c DEG and DMG in CD8
##############################

#Overlap with genes from methylation data
DEG_DMG_1st_CD8 <- hyge_test_sub(rownames(DEG_1st), rownames(count_cd8), MS1stvsHP1st_CD8_DMGs$Gene )
DEG_DMG_2nd_CD8 <- hyge_test_sub(rownames(DEG_2nd), rownames(count_cd8), MS2ndvsHP2nd_CD8_DMGs$Gene )
DEG_DMG_3rd_CD8 <- hyge_test_sub(rownames(DEG_3rd), rownames(count_cd8), MS3rdvsHP3rd_CD8_DMGs$Gene )
DEG_DMG_PP_CD8 <- hyge_test_sub(rownames(DEG_PP), rownames(count_cd8), MSPPvsHPPP_CD8_DMGs$Gene )

#Save files
saveRDS(DEG_DMG_1st_CD8, paste0(RDS_folder, "DEG_DMG/DEG_DMG_1st_CD8.RDS"))
saveRDS(DEG_DMG_2nd_CD8, paste0(RDS_folder, "DEG_DMG/DEG_DMG_2nd_CD8.RDS"))
saveRDS(DEG_DMG_3rd_CD8, paste0(RDS_folder, "DEG_DMG/DEG_DMG_3rd_CD8.RDS"))
saveRDS(DEG_DMG_PP_CD8, paste0(RDS_folder, "DEG_DMG/DEG_DMG_PP_CD8.RDS"))

#Overlap with genes from methylation data, all time-points merged
DMG_CD8 <- unique(c(MS1stvsHP1st_CD8_DMGs$Gene, MS2ndvsHP2nd_CD8_DMGs$Gene,
                    MS3rdvsHP3rd_CD8_DMGs$Gene, MSPPvsHPPP_CD8_DMGs$Gene))

DEG_CD8 <- unique(c(rownames(DEG_1st), rownames(DEG_2nd),
                    rownames(DEG_3rd), rownames(DEG_PP)))

DEG_DMG_CD8 <- hyge_test_sub(DEG_CD8, rownames(count_cd8), DMG_CD8 )

saveRDS(DEG_DMG_CD8, paste0(RDS_folder, "DEG_DMG/DEG_DMG_CD8.RDS"))

####################################
## 1.3.d DEG and DMG and MS in CD8
####################################

DEG_DMG_1st_MS_CD8 <- hyge_test(DEG_DMG_1st_CD8$gene_list[[1]], rownames(count_cd8) )
DEG_DMG_2nd_MS_CD8 <- hyge_test(DEG_DMG_2nd_CD8$gene_list[[1]], rownames(count_cd8) )
DEG_DMG_3rd_MS_CD8 <- hyge_test(DEG_DMG_3rd_CD8$gene_list[[1]], rownames(count_cd8) )
DEG_DMG_PP_MS_CD8 <- hyge_test(DEG_DMG_PP_CD8$gene_list[[1]], rownames(count_cd8) )

#Save files
saveRDS(DEG_DMG_1st_MS_CD8, paste0(RDS_folder, "DEG_DMG_MS/DEG_DMG_1st_MS_CD8.RDS"))
saveRDS(DEG_DMG_2nd_MS_CD8, paste0(RDS_folder, "DEG_DMG_MS/DEG_DMG_2nd_MS_CD8.RDS"))
saveRDS(DEG_DMG_3rd_MS_CD8, paste0(RDS_folder, "DEG_DMG_MS/DEG_DMG_3rd_MS_CD8.RDS"))
saveRDS(DEG_DMG_PP_MS_CD8, paste0(RDS_folder, "DEG_DMG_MS/DEG_DMG_PP_MS_CD8.RDS"))

#Overlap with genes from methylation data, all time-points merged
DEG_DMG_CD8 <- unique(c(DEG_DMG_1st_MS_CD8$gene_list[[1]], DEG_DMG_2nd_MS_CD8$gene_list[[1]],
                        DEG_DMG_3rd_MS_CD8$gene_list[[1]], DEG_DMG_PP_MS_CD8$gene_list[[1]]))

DEG_DMG_MS_CD8 <- hyge_test(DEG_DMG_CD8, rownames(count_cd8) )

#Save file
saveRDS(DEG_DMG_CD8, paste0(RDS_folder, "DEG_DMG_MS_CD8.RDS"))