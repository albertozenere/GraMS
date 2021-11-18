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
pack_R <- c("sva","tidyverse","readxl","assertthat",
            "edgeR","patchwork","ggfortify","RColorBrewer")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}
set.seed(541)

#Set important folders
RDS_folder <- "RDS/"
data_folder <- "data/"
figure_folder <- "figures/"
results_folder <- "results/"

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
MS_new_genes <- read.csv(file = "data/complete_MS_annotation.csv")
MS_new_genes <- unique(MS_new_genes$SYMBOL)

#other MS relevant genes
MS_old_genes <- read_excel("data/annotated_MS_SNPs_proc.xlsx")
idx <- MS_old_genes[,2]<0.00001
MS_old_genes <- unique(MS_old_genes$nearest_gene[idx])

######################
## 0.3 My functions ##
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

assign_colors <- function(string_array){ #creates a vector the same length as the input and replaces colors in it, useful for MDS
  out <- string_array
  
  u <- unique(string_array)
  if (length(u)==2){
    color_palette <- c("red","blue")
  }
  if (length(u)>2 & length(u)<13){
    color_palette <- brewer.pal(n = length(u), name = "Paired")
  }
  if (length(u)>12){
    color_palette <- c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
  }
  for (i in 1:length(u)){
    idx <- which( string_array == u[i])
    out[idx] <- color_palette[i]
  }
  return(out)
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
  

  ## Do the hypergeometric test

  #hypergeometric test
  q <- sum(is.element(MS_old_genes, subset_1))
  m <- sum(is.element(MS_old_genes, full_1))
  n <- length(full_1) - m
  k <- length(subset_1)
  pval_old <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  odds_old <- (k/q) / ( (k+m)/m ) #ratio in subset_1 compared to the ratio in full_1
  n_genes_old <- q #n. genes
  genes_list_old <- list(intersect(MS_old_genes, subset_1))
  
  if (m==0){
    warning("There is no overlap between all the genes measured and MS_old_genes. Check the data.")
  }
  
  #hypergeometric test
  q <- sum(is.element(MS_new_genes, subset_1))
  m <- sum(is.element(MS_new_genes, full_1))
  n <- length(full_1) - m
  k <- length(subset_1)
  pval_new <- phyper(q, m, n, k, lower.tail = FALSE, log.p = FALSE)
  odds_new <- (k/q) / ( (k+m)/m ) #ratio in subset_1 compared to the ratio in full_1
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


#Calculates how many genes are up or down regulated
kegg_regulation <- function(kegg, count, metadata){
  
  Up_Down <- character(nrow(kegg)) 
  
  for (n in 1:nrow(kegg)){ #one pathway at a time
    genes <- strsplit(kegg$geneID[n], "/")[[1]]
    
    # find which comparison it is
    if (str_detect( kegg$lab[n],"PPvs3rd" ) ){
      idx_2 <- metadata$TimePoint=="Postpartum"
      idx_1 <- metadata$TimePoint=="3rd"
    } else if (str_detect( kegg$lab[n],"3rdvs2nd" ) ){
      idx_2 <- metadata$TimePoint=="3rd"
      idx_1 <- metadata$TimePoint=="2nd"
    } else if (str_detect( kegg$lab[n],"2ndvs1st" ) ){
      idx_2 <- metadata$TimePoint=="2nd"
      idx_1 <- metadata$TimePoint=="1st"
    } else{
      stopifnot(1<0)
    }
    
    stopifnot(any(idx_2))
    stopifnot(any(idx_1))
    
    #extract count
    count_2 <- count[genes,idx_2]
    count_1 <- count[genes,idx_1]
    
    #extract individuals
    Ind_2 <- metadata$Individual[idx_2]
    Ind_1 <- metadata$Individual[idx_1]
    
    colnames(count_1) <- Ind_1
    colnames(count_2) <- Ind_2
    
    #select only common individuals
    common_Ind <- intersect(Ind_1,Ind_2)
    count_1 <- count_1[, is.element(colnames(count_1),common_Ind)]
    count_2 <- count_2[, is.element(colnames(count_2),common_Ind)]
    
    #order in same way
    ord <- match(colnames(count_2), colnames(count_1))
    count_1 <- count_1[,ord]
    stopifnot(all(colnames(count_1)==colnames(count_2)))
    
    str = ""
    #count in how many samples genes go up or down
    for (n_gene in 1:length(genes)){
      
      G <- genes[n_gene]
      
      n_Up <- sum((count_2[G,] - count_1[G,])>0)
      n_Down <- sum((count_2[G,] - count_1[G,])<0)
      
      str <- paste0(str, n_Up, ":", n_Down,"/")
    }
    
    str <- substr(str,1,nchar(str)-1)    
    Up_Down[n] <- str
    
  }
  
  kegg$Up_Down <- Up_Down   
  
  return(kegg)
  
}

######################
## 1.1 Limma on CD4 ##
######################

#For the initial limma analysis we select only complete samples (i.e. all time-points) and resting

#from Paired_Samples.txt
#complete_cd4 <- as.character( c(11,22,26,59,63,901,924,925,928))

#remove HC
is_HC <- which(is.element(metadata_cd4$Disease, "HC"))
count_cd4 <- count_cd4[, -is_HC]
metadata_cd4 <- metadata_cd4[-is_HC, ]

#is_complete <- is.element(metadata_cd4$Individual, complete_cd4) #should be ok to select all
is_complete <- rep(TRUE, nrow(metadata_cd4))
is_resting <- is.element(metadata_cd4$State, "Resting")

count_cd4 <- count_cd4[, is_complete & is_resting]
metadata_cd4 <- metadata_cd4[is_complete & is_resting, ]

stopifnot(colnames(count_cd4) == metadata_cd4$NGI_ID)
stopifnot(ncol(count_cd4) == nrow(metadata_cd4))

metadata_cd4$Sample_Group <- gsub(" ", "", metadata_cd4$Sample_Group)

#split MS and HC
is_MS <- is.element(metadata_cd4$Disease, "MS")
is_HP <- is.element(metadata_cd4$Disease, "HP")

count_MS <- count_cd4[, is_MS]
metadata_MS <- metadata_cd4[is_MS, ]
stopifnot(ncol(count_MS) == nrow(metadata_MS))

count_HP <- count_cd4[, is_HP]
metadata_HP <- metadata_cd4[is_HP, ]
stopifnot(ncol(count_HP) == nrow(metadata_HP))

#################################################
## 1.1.1 Comparisons between time points in MS ##
#################################################

#no contrast with duplicateCorrelation
design <- model.matrix(~ TimePoint, metadata_MS )
corfit <- duplicateCorrelation(count_MS, design, block=metadata_MS$Individual) #block the individual effect
fit <- lmFit(count_MS, design, block=metadata_MS$Individual, correlation=corfit$consensus)
fit <- eBayes(fit)
summary(decideTests(fit))

#Plot distribution of p-values of all covariates in the model
pdf(file=paste0(figure_folder, "pval_cd4_MS.pdf"))
plot_pvalues(fit)
dev.off()

#Plot the first two components
pdf(paste0(figure_folder, "pca12_cd4_MS.pdf"), width = 8, height = 5)
count_pca <- prcomp(t(count_MS))
p1 <- autoplot(count_pca, data = metadata_MS, colour = 'TimePoint', size = 3)
p2 <- autoplot(count_pca, data = metadata_MS, colour = 'Individual', size = 3)
p1 + p2
dev.off()

#with contrast
design <- model.matrix(~ 0 + TimePoint, metadata_MS )
corfit <- duplicateCorrelation(count_MS, design, block=metadata_MS$Individual) #block the individual effect
contr_matrix <- makeContrasts(
  TimePoint_1stvsBP = TimePoint1st-TimePointBefore,
  TimePoint_2ndvs1st = TimePoint2nd-TimePoint1st,
  TimePoint_3rdvs2nd = TimePoint3rd-TimePoint2nd,
  TimePoint_PPvs3rd = TimePointPostpartum-TimePoint3rd,
  levels=colnames(design)
)
fit <- lmFit(count_MS , design, block=metadata_MS$Individual, correlation=corfit$consensus)
fit <- contrasts.fit(fit, contrasts=contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

pdf(file=paste0(figure_folder, "pval_cd4_MS_timepoint.pdf"), width = 9)
plot_pvalues(fit)
dev.off()

#Print nominally expressed genes
tab_1stvsBP <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_1stvsBP",number=Inf)
tab_2ndvs1st <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_2ndvs1st",number=Inf)
tab_3rdvs2nd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_3rdvs2nd",number=Inf)
tab_PPvs3rd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_PPvs3rd",number=Inf)

#print list to file
lapply(rownames(tab_1stvsBP), write, paste0(results_folder, "cd4_MS_1stvsBP.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tab_2ndvs1st), write, paste0(results_folder, "cd4_MS_2ndvs1st.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tab_3rdvs2nd), write, paste0(results_folder, "cd4_MS_3rdvs2nd.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tab_PPvs3rd), write, paste0(results_folder, "cd4_MS_PPvs3rd.txt"), append= TRUE ,ncolumns=1000)

#How many genes go up or down
kegg <- read.csv(file = paste0(figure_folder, "kegg_cd4_MS_timepoint.txt"), sep = " ")
kegg_cd4_MS <- kegg_regulation(kegg,count_cd4,metadata_cd4)
write.table(kegg_cd4_MS, paste0(results_folder, "kegg_cd4_MS_timepoint_.txt"))

#################################################
## 1.1.2 Comparisons between time points in HP ##
#################################################

#no contrast
design <- model.matrix(~ TimePoint, metadata_HP)
corfit <- duplicateCorrelation(count_HP, design, block=metadata_HP$Individual) #block the individual effect
fit <- lmFit(count_HP, design, block=metadata_HP$Individual, correlation=corfit$consensus)
fit <- eBayes(fit)
summary(decideTests(fit))

#Plot distribution of p-values of all covariates in the model
pdf(file=paste0(figure_folder, "pval_cd4_HP.pdf"))
plot_pvalues(fit)
dev.off()

#Plot the first two components
pdf(paste0(figure_folder, "pca12_cd4_HP.pdf"), width = 8, height = 5)
count_pca <- prcomp(t(count_HP))
p1 <- autoplot(count_pca, data = metadata_HP, colour = 'TimePoint', size = 3)
p2 <- autoplot(count_pca, data = metadata_HP, colour = 'Individual', size = 3)
p1 + p2
dev.off()

#with contrast
design <- model.matrix(~ 0 + TimePoint, metadata_HP )
corfit <- duplicateCorrelation(count_HP, design, block=metadata_HP$Individual) #block the individual effect
contr_matrix <- makeContrasts(
  TimePoint_2ndvs1st = TimePoint2nd-TimePoint1st,
  TimePoint_3rdvs2nd = TimePoint3rd-TimePoint2nd,
  TimePoint_PPvs3rd = TimePointPostpartum-TimePoint3rd,
  levels=colnames(design)
)
fit <- lmFit(count_HP , design, block=metadata_HP$Individual, correlation=corfit$consensus)
fit <- contrasts.fit(fit, contrasts=contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

pdf(file=paste0(figure_folder, "pval_cd4_HP_timepoint.pdf"), width = 9, height = 4)
plot_pvalues(fit)
dev.off()

#Print nominally expressed genes
tab_2ndvs1st <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_2ndvs1st",number=Inf)
tab_3rdvs2nd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_3rdvs2nd",number=Inf)
tab_PPvs3rd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_PPvs3rd",number=Inf)

#print list to file
lapply(rownames(tab_2ndvs1st), write, paste0(results_folder, "cd4_HP_2ndvs1st.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tab_3rdvs2nd), write, paste0(results_folder, "cd4_HP_3rdvs2nd.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tab_PPvs3rd), write, paste0(results_folder, "cd4_HP_PPvs3rd.txt"), append= TRUE ,ncolumns=1000)

#How many genes go up or down
kegg <- read.csv(file = paste0(figure_folder, "kegg_cd4_HP_timepoint.txt"), sep = " ")
kegg_cd4_HP <- kegg_regulation(kegg,count_cd4,metadata_cd4)
write.table(kegg_cd4_HP, paste0(results_folder, "kegg_cd4_HP_timepoint_.txt"))

#################################
## 1.1.3 Comparisons MS and HP ##
#################################

pdf(file=paste0(figure_folder, "mds_cd4.pdf"), width = 9, height = 8)
par(mfrow=c(1,3))
mds1 <- plotMDS(count_cd4, top = dim(count_cd4)[1], labels = metadata_cd4$TimePoint, col = assign_colors( metadata_cd4$TimePoint),cex=1.5)
mds2 <- plotMDS(count_cd4, top = dim(count_cd4)[1], labels = metadata_cd4$Individual, col = assign_colors( metadata_cd4$Individual),cex=1.5)
mds3 <- plotMDS(count_cd4, top = dim(count_cd4)[1], labels = metadata_cd4$Disease, col = assign_colors( metadata_cd4$Disease),cex=1.5)
mds4 <- plotMDS(count_cd4, top = dim(count_cd4)[1], labels = metadata_cd4$Sample_Site, col = assign_colors( metadata_cd4$Sample_Site),cex=1.5)
dev.off()

#no constrast with duplicateCorrelation
design <- model.matrix(~ TimePoint + Disease, metadata_cd4 )
corfit <- duplicateCorrelation(count_cd4, design, block=metadata_cd4$Individual) #block the individual effect
fit <- lmFit(count_cd4, design, block=metadata_cd4$Individual, correlation=corfit$consensus)
fit <- eBayes(fit)
summary(decideTests(fit))

#enrichment on Disease term
tabMS <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="DiseaseMS")
hyge_test(rownames(tabMS), rownames(count_cd4))

#enrichment on 3rd trimester term
tab3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="TimePoint3rd")
hyge_test(rownames(tab3rd), rownames(count_cd4))

#Plot distribution of p-values of all covariates in the model
pdf(file=paste0(figure_folder, "pval_cd4.pdf"))
plot_pvalues(fit)
dev.off()

#Print nominally expressed genes
tab3rd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint3rd",number=Inf)
tabMS <- topTable(fit,adjust.method="none", p.value = 0.01,coef="DiseaseMS",number=Inf)

#print list to file
lapply(rownames(tab3rd), write, paste0(results_folder, "cd4_3rd.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tabMS), write, paste0(results_folder, "cd4_MS.txt"), append= TRUE ,ncolumns=1000)

#Plot the first two components
pdf(paste0(figure_folder, "pca12_cd4.pdf"), width = 13, height = 5)
count_pca <- prcomp(t(count_cd4))
p1 <- autoplot(count_pca, data = metadata_cd4, colour = 'TimePoint', size = 3)
p2 <- autoplot(count_pca, data = metadata_cd4, colour = 'Individual', size = 3)
p3 <- autoplot(count_pca, data = metadata_cd4, colour = 'Disease', size = 3)
p4 <- autoplot(count_pca, data = metadata_cd4, colour = 'Sample_Site', size = 3)
p1 + p2 + p3 + p4
dev.off()

#compare each time-point MS vs HP
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

#enrichment on 3rd trimester term
tab2nd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS2ndvsHP2nd")
hyge_test(rownames(tab2nd), rownames(count_cd4))

tab3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS3rdvsHP3rd")
hyge_test(rownames(tab3rd), rownames(count_cd4))

tabPP <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MSPPvsHPPP")
hyge_test(rownames(tabPP), rownames(count_cd4))


#print list to file
tab2nd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="MS2ndvsHP2nd", number=Inf)
tab3rd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="MS3rdvsHP3rd", number=Inf)
tabPP <- topTable(fit,adjust.method="none", p.value = 0.01,coef="MSPPvsHPPP", number=Inf)

#print list to file
lapply(rownames(tab2nd), write, paste0(results_folder, "cd4_MS2nd_vs_HP2nd.txt"), append=TRUE, ncolumns=1000)
lapply(rownames(tab3rd), write, paste0(results_folder, "cd4_MS3rd_vs_HP3rd.txt"), append=TRUE,ncolumns=1000)
lapply(rownames(tabPP), write, paste0(results_folder, "cd4_MSPP_vs_HPPP.txt"), append=TRUE, ncolumns=1000)


#Plot distribution of p-values of all covariates in the model
pdf(file=paste0(figure_folder, "pval_cd4_MSvsHP.pdf"))
plot_pvalues(fit)
dev.off()

#compare each time-point MS vs HP with baseline correction
design <- model.matrix(~ 0 + Sample_Group, metadata_cd4 )
corfit <- duplicateCorrelation(count_cd4, design, block=metadata_cd4$Individual) #block the individual effect
fit <- lmFit(count_cd4, design, block=metadata_cd4$Individual, correlation=corfit$consensus)

contr_matrix <- makeContrasts(
  MS2ndvsHP2nd = (Sample_GroupMS_2nd_CD4 - Sample_GroupMS_1st_CD4) - (Sample_GroupHP_2nd_CD4 - Sample_GroupHP_1st_CD4) ,
  MS3rdvsHP3rd = (Sample_GroupMS_3rd_CD4 - Sample_GroupMS_1st_CD4) - (Sample_GroupHP_3rd_CD4 - Sample_GroupHP_1st_CD4) ,
  MSPPvsHPPP = (Sample_GroupMS_PP_CD4 - Sample_GroupMS_1st_CD4) - (Sample_GroupHP_PP_CD4 - Sample_GroupHP_1st_CD4) ,
  levels=colnames(design)
)
#rownames(contr_matrix)[1] <- "(Intercept)" #if intercept is included in design
fit <- contrasts.fit(fit, contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

pdf(file=paste0(figure_folder, "pval_cd4_MSvsHP_baselinecorrected.pdf"))
plot_pvalues(fit)
dev.off()


#compare sequential time-points
design <- model.matrix(~ 0 + TimePoint + Disease, metadata_cd4)
corfit <- duplicateCorrelation(count_cd4, design, block=metadata_cd4$Individual) #block the individual effect
fit <- lmFit(count_cd4, design, block=metadata_cd4$Individual, correlation=corfit$consensus)

contr_matrix <- makeContrasts(
  TimePoint_2ndvs1st = TimePoint2nd - TimePoint1st ,
  TimePoint_3rdvs2nd = TimePoint3rd - TimePoint2nd ,
  TimePoint_PPvs3rd = TimePointPostpartum - TimePoint3rd ,
  levels=colnames(design)
)

fit <- contrasts.fit(fit, contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

pdf(file=paste0(figure_folder, "pval_cd4_timepoint.pdf"))
plot_pvalues(fit)
dev.off()

#remove BP
is_BP <- which(is.element(metadata_cd4$TimePoint, "Before"))
count_cd4 <- count_cd4[, -is_BP]
metadata_cd4 <- metadata_cd4[-is_BP, ]

#Interaction term
design <- model.matrix(~ Disease:TimePoint, metadata_cd4 )
corfit <- duplicateCorrelation(count_cd4, design, block=metadata_cd4$Individual) #block the individual effect
fit <- lmFit(count_cd4, design, block=metadata_cd4$Individual, correlation=corfit$consensus)
fit <- eBayes(fit)
summary(decideTests(fit))
plot_pvalues(fit)

######################
## 1.2 Limma on CD8 ##
######################

#For the initial limma analysis we select only complete samples (i.e. all time-points) and resting

#from Paired_Samples.txt
#complete_cd8 <- as.character( c(10,11,12,3,62,923,924,925))

#remove HC
is_HC <- which(is.element(metadata_cd8$Disease, "HC"))
count_cd8 <- count_cd8[, -is_HC]
metadata_cd8 <- metadata_cd8[-is_HC, ]

#is_complete <- is.element(metadata_cd8$Individual, complete_cd8) #should be fine to include non complete individuals
is_complete <- rep(TRUE,nrow(metadata_cd8))
is_resting <- is.element(metadata_cd8$State, "Resting")

count_cd8 <- count_cd8[, is_complete & is_resting]
metadata_cd8 <- metadata_cd8[is_complete & is_resting, ]

metadata_cd8$Sample_Group <- gsub(" ", "", metadata_cd8$Sample_Group)

stopifnot(colnames(count_cd8) == metadata_cd8$NGI_ID)
stopifnot(ncol(count_cd8) == nrow(metadata_cd8))

#split MS and HC
is_MS <- is.element(metadata_cd8$Disease, "MS")
is_HP <- is.element(metadata_cd8$Disease, "HP")

count_MS <- count_cd8[, is_MS]
metadata_MS <- metadata_cd8[is_MS, ]
stopifnot(ncol(count_MS) == nrow(metadata_MS))

count_HP <- count_cd8[, is_HP]
metadata_HP <- metadata_cd8[is_HP, ]
stopifnot(ncol(count_HP) == nrow(metadata_HP))

#################################################
## 1.2.1 Comparisons between time points in MS ##
#################################################

#no constrast
design <- model.matrix(~ TimePoint + Individual, metadata_MS )
fit_MS <- lmFit(count_MS , design)
fit_MS <- eBayes(fit_MS)
summary(decideTests(fit_MS))

#Plot distribution of p-values of all covariates in the model
pdf(file=paste0(figure_folder, "pval_CD8_MS.pdf"))
plot_pvalues(fit_MS)
dev.off()

#Plot the first two components of PCA
pdf(paste0(figure_folder, "pca12_cd8_MS.pdf"), width = 8, height = 5)
count_pca <- prcomp(t(count_MS))
p1 <- autoplot(count_pca, data = metadata_MS, colour = 'TimePoint', size = 3)
p2 <- autoplot(count_pca, data = metadata_MS, colour = 'Individual', size = 3)
p1 + p2
dev.off()

#with constrast
design <- model.matrix(~ 0 + TimePoint, metadata_MS )
corfit <- duplicateCorrelation(count_MS, design, block=metadata_MS$Individual) #block the individual effect
contr_matrix <- makeContrasts(
  TimePoint_1stvsBP = TimePoint1st-TimePointBefore,
  TimePoint_2ndvs1st = TimePoint2nd-TimePoint1st,
  TimePoint_3rdvs2nd = TimePoint3rd-TimePoint2nd,
  TimePoint_PPvs3rd = TimePointPostpartum-TimePoint3rd,
  levels=colnames(design)
)
fit <- lmFit(count_MS , design, block=metadata_MS$Individual, correlation=corfit$consensus)
fit <- contrasts.fit(fit, contrasts=contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

pdf(file=paste0(figure_folder, "pval_cd8_MS_timepoint.pdf"), width = 9, height = 4)
plot_pvalues(fit)
dev.off()

#Print nominally expressed genes
tab_1stvsBP <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_1stvsBP",number=Inf)
tab_2ndvs1st <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_2ndvs1st",number=Inf)
tab_3rdvs2nd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_3rdvs2nd",number=Inf)
tab_PPvs3rd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_PPvs3rd",number=Inf)

#print list to file
lapply(rownames(tab_1stvsBP), write, paste0(results_folder, "cd8_MS_1stvsBP.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tab_2ndvs1st), write, paste0(results_folder, "cd8_MS_2ndvs1st.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tab_3rdvs2nd), write, paste0(results_folder, "cd8_MS_3rdvs2nd.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tab_PPvs3rd), write, paste0(results_folder, "cd8_MS_PPvs3rd.txt"), append= TRUE ,ncolumns=1000)

#How many genes go up or down
kegg <- read.csv(file = paste0(figure_folder, "kegg_cd8_MS_timepoint.txt"), sep = " ")
kegg_cd8_MS <- kegg_regulation(kegg,count_cd8,metadata_cd8)
write.table(kegg_cd8_MS, paste0(results_folder, "kegg_cd8_MS_timepoint_.txt"))


#################################################
## 1.2.2 Comparisons between time points in HP ##
#################################################

#no contrast
design <- model.matrix(~ TimePoint + Individual, metadata_HP )
fit_HP <- lmFit(count_HP , design)
fit_HP <- eBayes(fit_HP)
summary(decideTests(fit_HP))

#Plot distribution of p-values of all covariates in the model
pdf(file=paste0(figure_folder, "pval_CD8_HP.pdf"))
plot_pvalues(fit_HP)
dev.off()

#Plot the first two components
pdf(paste0(figure_folder, "pca12_cd8_HP.pdf"), width = 8, height = 5)
count_pca <- prcomp(t(count_HP))
p1 <- autoplot(count_pca, data = metadata_HP, colour = 'TimePoint', size = 3)
p2 <- autoplot(count_pca, data = metadata_HP, colour = 'Individual', size = 3)
p1 + p2
dev.off()

#with contrast
design <- model.matrix(~ 0 + TimePoint, metadata_HP )
corfit <- duplicateCorrelation(count_HP, design, block=metadata_HP$Individual) #block the individual effect
contr_matrix <- makeContrasts(
  TimePoint_2ndvs1st = TimePoint2nd-TimePoint1st,
  TimePoint_3rdvs2nd = TimePoint3rd-TimePoint2nd,
  TimePoint_PPvs3rd = TimePointPostpartum-TimePoint3rd,
  levels=colnames(design)
)
fit <- lmFit(count_HP , design, block=metadata_HP$Individual, correlation=corfit$consensus)
fit <- contrasts.fit(fit, contrasts=contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

pdf(file=paste0(figure_folder, "pval_cd8_HP_timepoint.pdf"), width = 9, height = 4)
plot_pvalues(fit)
dev.off()

#Print nominally expressed genes
tab_2ndvs1st <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_2ndvs1st",number=Inf)
tab_3rdvs2nd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_3rdvs2nd",number=Inf)
tab_PPvs3rd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint_PPvs3rd",number=Inf)

#print list to file
lapply(rownames(tab_2ndvs1st), write, paste0(results_folder, "cd8_HP_2ndvs1st.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tab_3rdvs2nd), write, paste0(results_folder, "cd8_HP_3rdvs2nd.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tab_PPvs3rd), write, paste0(results_folder, "cd8_HP_PPvs3rd.txt"), append= TRUE ,ncolumns=1000)

#How many genes go up or down
kegg <- read.csv(file = paste0(figure_folder, "kegg_cd8_HP_timepoint.txt"), sep = " ")
kegg_cd8_MS <- kegg_regulation(kegg,count_cd8,metadata_cd8)
write.table(kegg_cd8_MS, paste0(results_folder, "kegg_cd8_HP_timepoint_.txt"))

#################################
## 1.2.3 Comparisons MS and HP ##
#################################

#MDS plot
pdf(file=paste0(figure_folder, "mds_cd8.pdf"), width = 9, height = 8)
par(mfrow=c(1,3))
mds1 <- plotMDS(count_cd8, top = dim(count_cd8)[1], labels = metadata_cd8$TimePoint, col = assign_colors( metadata_cd8$TimePoint),cex=1.5)
mds2 <- plotMDS(count_cd8, top = dim(count_cd8)[1], labels = metadata_cd8$Individual, col = assign_colors( metadata_cd8$Individual),cex=1.5)
mds3 <- plotMDS(count_cd8, top = dim(count_cd8)[1], labels = metadata_cd8$Disease, col = assign_colors( metadata_cd8$Disease),cex=1.5)
mds4 <- plotMDS(count_cd8, top = dim(count_cd8)[1], labels = metadata_cd8$Sample_Site, col = assign_colors( metadata_cd8$Sample_Site),cex=1.5)
dev.off()

#no contrast with duplicateCorrelation
design <- model.matrix(~ TimePoint + Disease, metadata_cd8 )
corfit <- duplicateCorrelation(count_cd8, design, block=metadata_cd8$Individual) #block the individual effect
fit <- lmFit(count_cd8, design, block=metadata_cd8$Individual, correlation=corfit$consensus)
fit <- eBayes(fit)
summary(decideTests(fit))

#hypergeometric test
tabMS <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="DiseaseMS")
hyge_test(rownames(tabMS), rownames(count_cd8))

tab3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="TimePoint3rd")
hyge_test(rownames(tab3rd), rownames(count_cd8))

#Plot distribution of p-values of all covariates in the model
pdf(file=paste0(figure_folder, "pval_cd8.pdf"))
plot_pvalues(fit)
dev.off()

#Print nominally expressed genes
tab3rd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint3rd",number=Inf)
tabMS <- topTable(fit,adjust.method="none", p.value = 0.01,coef="DiseaseMS",number=Inf)

#print list to file
lapply(rownames(tab3rd), write, paste0(results_folder, "cd8_3rd.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tabMS), write, paste0(results_folder, "cd8_MS.txt"), append= TRUE ,ncolumns=1000)

#Plot the first two components
pdf(paste0(figure_folder, "pca12_cd8.pdf"), width = 13, height = 5)
count_pca <- prcomp(t(count_cd8))
p1 <- autoplot(count_pca, data = metadata_cd8, colour = 'TimePoint', size = 3)
p2 <- autoplot(count_pca, data = metadata_cd8, colour = 'Individual', size = 3)
p3 <- autoplot(count_pca, data = metadata_cd8, colour = 'Disease', size = 3)
p4 <- autoplot(count_pca, data = metadata_cd8, colour = 'Sample_Site', size = 3)
p1 + p2 + p3 + p4
dev.off()

#compare each time-point MS vs HP
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
#rownames(contr_matrix)[1] <- "(Intercept)" #if intercept is included in design
fit <- contrasts.fit(fit, contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

#hypergeometric test on independent lists of MS-associated genes
tab1st <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS1stvsHP1st")
hyge_test(rownames(tab1st), rownames(count_cd8))

tab2nd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS2ndvsHP2nd")
hyge_test(rownames(tab2nd), rownames(count_cd8))

tab3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS3rdvsHP3rd")
hyge_test(rownames(tab3rd), rownames(count_cd8))

tabPP <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MSPPvsHPPP")
hyge_test(rownames(tabPP), rownames(count_cd8))

#print list to file
tab1st <- topTable(fit,adjust.method="none", p.value = 0.01,coef="MS1stvsHP1st",number=Inf)
tab2nd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="MS2ndvsHP2nd",number=Inf)
tab3rd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="MS3rdvsHP3rd",number=Inf)
tabPP <- topTable(fit,adjust.method="none", p.value = 0.01,coef="MSPPvsHPPP",number=Inf)

#print list to file
lapply(rownames(tab1st), write, paste0(results_folder, "cd8_MS1st_vs_HP1st.txt"), append= TRUE, ncolumns=1000)
lapply(rownames(tab2nd), write, paste0(results_folder, "cd8_MS2nd_vs_HP2nd.txt"), append= TRUE, ncolumns=1000)
lapply(rownames(tab3rd), write, paste0(results_folder, "cd8_MS3rd_vs_HP3rd.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tabPP), write, paste0(results_folder, "cd8_MSPP_vs_HPPP.txt"),  append= TRUE, ncolumns=1000)

pdf(file=paste0(figure_folder, "pval_cd8_MSvsHP.pdf"))
plot_pvalues(fit)
dev.off()

#compare each time-point MS vs HP with baseline correction
design <- model.matrix(~ 0 + Sample_Group, metadata_cd8 )
corfit <- duplicateCorrelation(count_cd8, design, block=metadata_cd8$Individual) #block the individual effect
fit <- lmFit(count_cd8, design, block=metadata_cd8$Individual, correlation=corfit$consensus)

contr_matrix <- makeContrasts(
  MS2ndvsHP2nd = (Sample_GroupMS_2nd_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupHP_2nd_CD8 - Sample_GroupHP_1st_CD8) ,
  MS3rdvsHP3rd = (Sample_GroupMS_3rd_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupHP_3rd_CD8 - Sample_GroupHP_1st_CD8) ,
  MSPPvsHPPP = (Sample_GroupMS_PP_CD8 - Sample_GroupMS_1st_CD8) - (Sample_GroupHP_PP_CD8 - Sample_GroupHP_1st_CD8) ,
  levels=colnames(design)
)
#rownames(contr_matrix)[1] <- "(Intercept)" #if intercept is included in design
fit <- contrasts.fit(fit, contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

pdf(file=paste0(figure_folder, "pval_cd8_MSvsHP_baselinecorrected.pdf"))
plot_pvalues(fit)
dev.off()

#compare sequential time-points
design <- model.matrix(~ 0 + TimePoint + Disease, metadata_cd8)
corfit <- duplicateCorrelation(count_cd8, design, block=metadata_cd8$Individual) #block the individual effect
fit <- lmFit(count_cd8, design, block=metadata_cd8$Individual, correlation=corfit$consensus)

contr_matrix <- makeContrasts(
  TimePoint_2ndvs1st = TimePoint2nd - TimePoint1st ,
  TimePoint_3rdvs2nd = TimePoint3rd - TimePoint2nd ,
  TimePoint_PPvs3rd = TimePointPostpartum - TimePoint3rd ,
  levels=colnames(design)
)

fit <- contrasts.fit(fit, contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

pdf(file=paste0(figure_folder, "pval_cd8_timepoint.pdf"))
plot_pvalues(fit)
dev.off()


##################################
## 2.1 Using CD4 and CD8 together 
##################################

#merge cd4 and cd8
count <- merge(count_cd4,count_cd8, by=0)
rownames(count) <- count$Row.names
count <- count[,-1]

metadata <- rbind(metadata_cd4, metadata_cd8)

stopifnot(all(colnames(count)==metadata$NGI_ID))

#PCA
pdf(paste0(figure_folder, "pca12_merged.pdf"), width = 13, height = 5)
count_pca <- prcomp(t(count))
p1 <- autoplot(count_pca, data = metadata, colour = 'TimePoint', size = 3)
p2 <- autoplot(count_pca, data = metadata, colour = 'Disease', size = 3)
p3 <- autoplot(count_pca, data = metadata, colour = 'Sample_Type', size = 3)
p1 + p2 + p3 
dev.off()

#no contrast with duplicateCorrelation
design <- model.matrix(~ TimePoint + Disease + Sample_Type, metadata )
corfit <- duplicateCorrelation(count, design, block=metadata$Individual) #block the individual effect
fit <- lmFit(count, design, block=metadata$Individual, correlation=corfit$consensus)
fit <- eBayes(fit)
summary(decideTests(fit))

#hypergeometric test
tabCellType <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Sample_TypeCD8")
hyge_test(rownames(tabCellType), rownames(count))

tabMS <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="DiseaseMS")
hyge_test(rownames(tabMS), rownames(count))

tab3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="TimePoint3rd")
hyge_test(rownames(tab3rd), rownames(count))

#Plot distribution of p-values of all covariates in the model
pdf(file=paste0(figure_folder, "pval_merged.pdf"))
plot_pvalues(fit)
dev.off()

#Print nominally expressed genes
tab3rd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint3rd",number=Inf)
tabMS <- topTable(fit,adjust.method="none", p.value = 0.01,coef="DiseaseMS",number=Inf)
tabCellType <- topTable(fit,adjust.method="none", p.value = 0.01,coef="Sample_TypeCD8",number=Inf)

#print list to file
lapply(rownames(tab3rd), write, paste0(results_folder, "merged_3rd.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tabMS), write, paste0(results_folder, "merged_MS.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tabCellType), write, paste0(results_folder, "merged_CellType.txt"), append= TRUE ,ncolumns=1000)


#compare sequential time-points
design <- model.matrix(~ 0 + TimePoint + Disease + Sample_Type, metadata)
corfit <- duplicateCorrelation(count, design, block=metadata$Individual) #block the individual effect
fit <- lmFit(count, design, block=metadata$Individual, correlation=corfit$consensus)

contr_matrix <- makeContrasts(
  TimePoint_2ndvs1st = TimePoint2nd - TimePoint1st ,
  TimePoint_3rdvs2nd = TimePoint3rd - TimePoint2nd ,
  TimePoint_PPvs3rd = TimePointPostpartum - TimePoint3rd ,
  levels=colnames(design)
)

fit <- contrasts.fit(fit, contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

pdf(file=paste0(figure_folder, "pval_merged_timepoint.pdf"))
plot_pvalues(fit)
dev.off()

#replace Before, 1st ad 2nd with Early
metadata$TimePoint[ is.element(metadata$TimePoint, c("Before", "1st","2nd")) ] <- "Early"

#Interaction term
design <- model.matrix(~ Sample_Type + Disease*TimePoint, metadata )
corfit <- duplicateCorrelation(count, design, block=metadata$Individual) #block the individual effect
fit <- lmFit(count, design, block=metadata$Individual, correlation=corfit$consensus)
fit <- eBayes(fit)
summary(decideTests(fit))
plot_pvalues(fit)

###################################################
## 2.2 Using CD4 and CD8 together, add Activation 
###################################################
#Need to load them again, since I have removed the activated samples earlier

#CD8
count_cd4_ <- readRDS(file = paste0(RDS_folder, "count_cd4_tmm.RDS"))
metadata_cd4_ <- readRDS(file = paste0(RDS_folder, "metadata_cd4.RDS"))
stopifnot(colnames(count_cd4_) == metadata_cd4_$NGI_ID)

#CD8
count_cd8_ <- readRDS(file = paste0(RDS_folder, "count_cd8_tmm.RDS"))
metadata_cd8_ <- readRDS(file = paste0(RDS_folder, "metadata_cd8.RDS"))
stopifnot(colnames(count_cd8_) == metadata_cd8_$NGI_ID)

#merge cd4 and cd8
count <- merge(count_cd4_,count_cd8_, by=0)
rownames(count) <- count$Row.names
count <- count[,-1]

metadata <- rbind(metadata_cd4_, metadata_cd8_)
stopifnot(all(colnames(count)==metadata$NGI_ID))

#remove HC
is_HC <- which(is.element(metadata$Disease, "HC"))
count <- count[, -is_HC]
metadata <- metadata[-is_HC, ]

metadata$Sample_Group <- gsub(" ", "", metadata$Sample_Group)

#PCA
pdf(paste0(figure_folder, "pca12_merged_activation.pdf"), width = 13, height = 5)
count_pca <- prcomp(t(count))
p1 <- autoplot(count_pca, data = metadata, colour = 'TimePoint', size = 3)
p2 <- autoplot(count_pca, data = metadata, colour = 'Disease', size = 3)
p3 <- autoplot(count_pca, data = metadata, colour = 'Sample_Type', size = 3)
p4 <- autoplot(count_pca, data = metadata, colour = 'State', size = 3)
p1 + p2 + p3 + p4
dev.off()

#no contrast with duplicateCorrelation
design <- model.matrix(~ TimePoint + Disease + Sample_Type + State + Cell_Viability, metadata )
corfit <- duplicateCorrelation(count, design, block=metadata$Individual) #block the individual effect
fit <- lmFit(count, design, block=metadata$Individual, correlation=corfit$consensus)
fit <- eBayes(fit)
summary(decideTests(fit))

#Plot distribution of p-values of all covariates in the model
pdf(file=paste0(figure_folder, "pval_merged_activation.pdf"))
plot_pvalues(fit)
dev.off()

#hypergeometric test
tabCellType <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="Sample_TypeCD8")
hyge_test(rownames(tabCellType), rownames(count))

tabMS <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="DiseaseMS")
hyge_test(rownames(tabMS), rownames(count))

tab3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="TimePoint3rd")
hyge_test(rownames(tab3rd), rownames(count))

tabState <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="StateResting")
hyge_test(rownames(tabState), rownames(count))

#Print nominally expressed genes
tab3rd <- topTable(fit,adjust.method="none", p.value = 0.01,coef="TimePoint3rd",number=Inf)
tabMS <- topTable(fit,adjust.method="none", p.value = 0.01,coef="DiseaseMS",number=Inf)
tabCellType <- topTable(fit,adjust.method="none", p.value = 0.01,coef="Sample_TypeCD8",number=Inf)
tabState <- topTable(fit,adjust.method="none", p.value = 0.01,coef="StateResting",number=Inf)

#print list to file
lapply(rownames(tab3rd), write, paste0(results_folder, "merged_activation_3rd.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tabMS), write, paste0(results_folder, "merged_activation_MS.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tabCellType), write, paste0(results_folder, "merged_activation_CellType.txt"), append= TRUE ,ncolumns=1000)
lapply(rownames(tabState), write, paste0(results_folder, "merged_activation_State.txt"), append= TRUE ,ncolumns=1000)

#compare sequential time-points
design <- model.matrix(~ 0 + TimePoint + Disease  + State + Cell_Viability, metadata)
corfit <- duplicateCorrelation(count, design, block=metadata$Individual) #block the individual effect
fit <- lmFit(count, design, block=metadata$Individual, correlation=corfit$consensus)

contr_matrix <- makeContrasts(
  TimePoint_2ndvs1st = TimePoint2nd - TimePoint1st ,
  TimePoint_3rdvs2nd = TimePoint3rd - TimePoint2nd ,
  TimePoint_PPvs3rd = TimePointPostpartum - TimePoint3rd ,
  levels=colnames(design)
)

fit <- contrasts.fit(fit, contr_matrix)
fit <- eBayes(fit)
summary(decideTests(fit))

pdf(file=paste0(figure_folder, "pval_merged_activation_timepoint.pdf"))
plot_pvalues(fit)
dev.off()




















