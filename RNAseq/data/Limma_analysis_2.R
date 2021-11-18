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
MS_new_genes <- read.csv(file = "data/complete_MS_annotation_p10-6.csv")
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
  odds_old <- (q/k) / ( m/(n+m) ) #ratio in subset_1 compared to the ratio in full_1
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


do_clustering <- function(count,metadata,n_cl){
  
  #merge sample of the same timepoint by taking the mean
  timepoints <- c("1st","2nd","3rd","Postpartum")
  n_timepoint <- length(timepoints)
  count_mean <- matrix(0, nrow(count),n_timepoint)
  rownames(count_mean) <- rownames(count)
  
  for ( t in 1:length(timepoints)){
    idx <- metadata$TimePoint==timepoints[t]
    stopifnot(any(idx))
    idx <- which(idx)
    
    count_mean[,t] <- rowSums(count[,idx])/length(idx)
    
  }
  
  #Z-Normalize each gene
  for (n in 1:nrow(count_mean)){
    c <- count_mean[n,]
    count_mean[n,] <- (c-mean(c))/sqrt(var(c)) 
  }
  
  #optimal n. clusters
  # f <- fviz_nbclust(count_mean, kmeans, method = "silhouette", verbose = FALSE)
  # measure_cl <- f$data$y
  # opt_n_cl <- which(measure_cl==max(measure_cl))  
  
  #kmeans
  cl <- kmeans(count_mean,centers=n_cl)
  
  #plot
  Ngene <- length(cl$cluster)
  
  myplot <- list()
  p_values <- rep(0, n_cl)
  odds_ratio <- rep(0, n_cl)
  n_genes <- rep(0, n_cl)
  for (n in 1:n_cl){
    idx <- which(cl$cluster==n)
    
    cl_info <- matrix(0,n_timepoint,3) #init
    
    c <- count_mean[idx,] #count matrix
    
    cl_info <- t( apply(c, 2, quantile, probs=c(0.1,0.5,0.9)) )
    
    colnames(cl_info) <- c("Lower","Median","Upper")
    rownames(cl_info) <- colnames(count_mean)
    
    #store plot object
    myplot[[n]] <- ggplot(cl_info) + geom_line(aes(x=1:n_timepoint,y=Median), colour="black") +
      geom_errorbar(aes(x = 1:n_timepoint, ymin=Lower, ymax=Upper), colour="black", width=.1) +
      xlab("Time") + ylab(paste0("Cluster", n, " (n. ", length(idx),")")) +
      scale_x_continuous(breaks = 1:n_timepoint,labels = rownames(cl_info))
    
    
    #hyge on new list
    t <- hyge_test(rownames(count)[idx], rownames(count))
    p_values[n] <- t["new_list", "p_value"]
    odds_ratio[n] <- t["new_list", "odds_ratio"]
    n_genes[n] <- t["new_list", "n_genes"] 
    
  }
  
  do.call("grid.arrange", c(myplot, ncol=4)) 
  
  summary_clusters <- data.frame(p_values, odds_ratio, n_genes)
  rownames(summary_clusters) <- paste0("Cluster", 1:n_cl)
  
  print(summary_clusters)
  
  
}

######################
## 1.1 Prepare CD4
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

######################
## 1.1 Prepare CD8
######################

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

##############################
## 1.1 Development over time 
##############################

key_genes <- c("STAT1", "IFNG","IL2RA","NFKB1","SOCS1","SOCS2","IRF1","IL4R",
               "IL15RA","CCR7","TNFSF8")

stopifnot(all(is.element(key_genes, rownames(count_cd4))))
stopifnot(all(is.element(key_genes, rownames(count_cd8))))


timeseries_genes <- function(key_genes, count, metadata, filename){
  
  Ngene <- length(key_genes)
  time_measured <- c("1st","2nd","3rd","Postpartum")
  
  myplot <- list()
  for (n in 1:Ngene){
    g <- key_genes[n]
    
    #HP
    Ind_HP <- which(metadata$Disease == "HP")
    count_HP <- count[g,Ind_HP]
    
    g_info_HP <- matrix(0,length(time_measured),3)
    colnames(g_info_HP) <- c("Lower", "Median", "Upper")
    rownames(g_info_HP) <- time_measured
    
    
    for (t in 1:length(time_measured)){
      idx <- which(str_detect(metadata$TimePoint[Ind_HP], time_measured[t]))
      
      c <- count_HP[idx]
      
      g_info_HP[t,] <- quantile(c, probs=c(0.1,0.5,0.9))
    }
    g_info_HP <- as.data.frame(g_info_HP)
    
    
    #MS
    Ind_MS <- which(metadata$Disease == "MS")
    count_MS <- count[g,Ind_MS]
    
    g_info_MS <- matrix(0,length(time_measured),3)
    colnames(g_info_MS) <- c("Lower", "Median", "Upper")
    rownames(g_info_MS) <- time_measured
    
    for (t in 1:length(time_measured)){
      idx <- which(str_detect(metadata_cd4$TimePoint[Ind_MS], time_measured[t]))
      
      c <- count_MS[idx]
      
      g_info_MS[t,] <- quantile(c, probs=c(0.1,0.5,0.9))
      
    }
    g_info_MS <- as.data.frame(g_info_MS)
    
    #plot
    myplot[[n]] <- ggplot(g_info_HP) + geom_line(aes(x=1:length(time_measured),y=Median), colour="black") +
      geom_errorbar(aes(x = 1:length(time_measured), ymin=Lower, ymax=Upper), colour="black", width=.1) +
      geom_line(data=g_info_MS, aes(x=1:length(time_measured),y=Median), colour="red") +
      geom_errorbar(data=g_info_MS, aes(x = 1:length(time_measured), ymin=Lower, ymax=Upper), colour="red", width=.1) +
      xlab("Time") + ylab(key_genes[n])
    
    
  }
  
  pdf(file=paste0(figure_folder, filename, ".pdf"), height = 13)
  do.call("grid.arrange", c(myplot, ncol=3))
  dev.off()
  
}

timeseries_genes(key_genes, count_cd4, metadata_cd4, "cd4_key_genes")
timeseries_genes(key_genes, count_cd8, metadata_cd8, "cd8_key_genes")

##############################
## 1.1n DEG over time CD4
##############################

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
DEG_MS_1st <- hyge_test(rownames(DEG_1st), rownames(count_cd4))

DEG_2nd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS2ndvsHP2nd")
DEG_MS_2nd <- hyge_test(rownames(DEG_2nd), rownames(count_cd4))

DEG_3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS3rdvsHP3rd")
DEG_MS_3rd <- hyge_test(rownames(DEG_3rd), rownames(count_cd4))

DEG_PP <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MSPPvsHPPP")
DEG_MS_PP <- hyge_test(rownames(DEG_PP), rownames(count_cd4))

#Save files
saveRDS(DEG_1st, paste0(results_folder, "DEG_1st_CD4.RDS"))
saveRDS(DEG_2nd, paste0(results_folder, "DEG_2nd_CD4.RDS"))
saveRDS(DEG_3rd, paste0(results_folder, "DEG_3rd_CD4.RDS"))
saveRDS(DEG_PP, paste0(results_folder, "DEG_PP_CD4.RDS"))

#N. nominally DEG over time
DEG_over_time_cd4 <- data.frame(ngene = c(nrow(DEG_1st), nrow(DEG_2nd), nrow(DEG_3rd), nrow(DEG_PP)))

#N. genes overlapping with the MS-associated gene lists
MS_over_time_cd4 <- data.frame(ngene = c(MS_1st["new_list","n_genes"], MS_2nd["new_list","n_genes"], 
                                         MS_3rd["new_list","n_genes"], MS_PP["new_list","n_genes"]))

#plot
ggplot(DEG_over_time_cd4) + geom_line(aes(x=1:4,y=ngene)) + xlab("time") + ylab("n. nominally DEG") +
  ggplot(MS_over_time_cd4) + geom_line(aes(x=1:4,y=ngene)) + xlab("time") + ylab("n. MS-associated genes")


#Overlap with genes from methylation data
DEG_DMG_1st_MS_CD4 <- hyge_test_sub(rownames(DEG_1st), rownames(count_cd4), MS1stvsHP1st_CD4_DMGs$Gene )
DEG_DMG_2nd_MS_CD4 <- hyge_test_sub(rownames(DEG_2nd), rownames(count_cd4), MS2ndvsHP2nd_CD4_DMGs$Gene )
DEG_DMG_3rd_MS_CD4 <- hyge_test_sub(rownames(DEG_3rd), rownames(count_cd4), MS3rdvsHP3rd_CD4_DMGs$Gene )
DEG_DMG_PP_MS_CD4 <- hyge_test_sub(rownames(DEG_PP), rownames(count_cd4), MSPPvsHPPP_CD4_DMGs$Gene )

#Save files
saveRDS(DEG_DMG_1st_MS_CD4, paste0(results_folder, "DEG_DMG_1st_MS_CD4.RDS"))
saveRDS(DEG_DMG_2nd_MS_CD4, paste0(results_folder, "DEG_DMG_2nd_MS_CD4.RDS"))
saveRDS(DEG_DMG_3rd_MS_CD4, paste0(results_folder, "DEG_DMG_3rd_MS_CD4.RDS"))
saveRDS(DEG_DMG_PP_MS_CD4, paste0(results_folder, "DEG_DMG_PP_MS_CD4.RDS"))


#Comparison between the union of the DMG and DEG lists
DMG_CD4 <- unique(c(MS1stvsHP1st_CD4_DMGs$Gene, MS2ndvsHP2nd_CD4_DMGs$Gene,
                    MS3rdvsHP3rd_CD4_DMGs$Gene, MSPPvsHPPP_CD4_DMGs$Gene))

DEG_CD4 <- unique(c(rownames(DEG_1st), rownames(DEG_2nd),
                    rownames(DEG_3rd), rownames(DEG_PP)))

DEG_DMG_CD4 <- intersect(DMG_CD4,DEG_CD4)

DEG_DMG_MS_CD4 <- hyge_test(DEG_DMG_CD4, rownames(count_cd4))
saveRDS(DEG_DMG_MS_CD4, paste0(results_folder, "DEG_DMG_MS_CD4.RDS"))
DEG_DMG_MS_CD4 <- DEG_DMG_MS_CD4["new_list","gene_list"][[1]]

#############################
## 1.1n DEG over time CD8
#############################

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
MS_1st <- hyge_test(rownames(DEG_1st), rownames(count_cd8))

DEG_2nd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS2ndvsHP2nd")
MS_2nd <- hyge_test(rownames(DEG_2nd), rownames(count_cd8))

DEG_3rd <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MS3rdvsHP3rd")
MS_3rd <- hyge_test(rownames(DEG_3rd), rownames(count_cd8))

DEG_PP <- topTable(fit,adjust.method="none", p.value = 0.05, number = Inf, coef="MSPPvsHPPP")
MS_PP <- hyge_test(rownames(DEG_PP), rownames(count_cd8))

#N. nominally DEG over time
DEG_over_time_cd8 <- data.frame(ngene = c(nrow(DEG_1st), nrow(DEG_2nd), nrow(DEG_3rd), nrow(DEG_PP)))

#Save files
saveRDS(DEG_1st, paste0(results_folder, "DEG_1st_CD8.RDS"))
saveRDS(DEG_2nd, paste0(results_folder, "DEG_2nd_CD8.RDS"))
saveRDS(DEG_3rd, paste0(results_folder, "DEG_3rd_CD8.RDS"))
saveRDS(DEG_PP, paste0(results_folder, "DEG_PP_CD8.RDS"))

#N. genes overlapping with the MS-associated gene lists
MS_over_time_cd8 <- data.frame(ngene = c(MS_1st["new_list","n_genes"], MS_2nd["new_list","n_genes"], 
                                         MS_3rd["new_list","n_genes"], MS_PP["new_list","n_genes"]))

#plot
ggplot(DEG_over_time_cd8) + geom_line(aes(x=1:4,y=ngene)) + xlab("time") + ylab("n. nominally DEG") +
  ggplot(MS_over_time_cd8) + geom_line(aes(x=1:4,y=ngene)) + xlab("time") + ylab("n. MS-associated genes") 

#Overlap with genes from methylation data
DEG_DMG_1st_MS <- hyge_test_sub(rownames(DEG_1st), rownames(count_cd8), MS1stvsHP1st_CD8_DMGs$Gene )
DEG_DMG_2nd_MS <- hyge_test_sub(rownames(DEG_2nd), rownames(count_cd8), MS2ndvsHP2nd_CD8_DMGs$Gene )
DEG_DMG_3rd_MS <- hyge_test_sub(rownames(DEG_3rd), rownames(count_cd8), MS3rdvsHP3rd_CD8_DMGs$Gene )
DEG_DMG_PP_MS <- hyge_test_sub(rownames(DEG_PP), rownames(count_cd8), MSPPvsHPPP_CD8_DMGs$Gene )

#Save files
saveRDS(DEG_DMG_1st_MS, paste0(results_folder, "DEG_DMG_1st_MS_CD8.RDS"))
saveRDS(DEG_DMG_2nd_MS, paste0(results_folder, "DEG_DMG_2nd_MS_CD8.RDS"))
saveRDS(DEG_DMG_3rd_MS, paste0(results_folder, "DEG_DMG_3rd_MS_CD8.RDS"))
saveRDS(DEG_DMG_PP_MS, paste0(results_folder, "DEG_DMG_PP_MS_CD8.RDS"))

#Comparison between the union of the DMG and DEG lists
DMG_CD8 <- unique(c(MS1stvsHP1st_CD8_DMGs$Gene, MS2ndvsHP2nd_CD8_DMGs$Gene,
                 MS3rdvsHP3rd_CD8_DMGs$Gene, MSPPvsHPPP_CD8_DMGs$Gene))

DEG_CD8 <- unique(c(rownames(DEG_1st), rownames(DEG_2nd),
                 rownames(DEG_3rd), rownames(DEG_PP)))

DEG_DMG_CD8 <- intersect(DMG_CD8,DEG_CD8)

DEG_DMG_MS_CD8 <- hyge_test(DEG_DMG_CD8, rownames(count_cd8))
saveRDS(DEG_DMG_MS_CD8, paste0(results_folder, "DEG_DMG_MS_CD8.RDS"))
DEG_DMG_MS_CD8 <- DEG_DMG_MS_CD8["new_list","gene_list"][[1]]

#######################################
## 1. Time development of DEG DMG MS
#######################################

timeseries_genes(DEG_DMG_MS_CD4, count_cd4, metadata_cd4, "cd4_DEG_DMG_MS")
timeseries_genes(DEG_DMG_MS_CD8, count_cd8, metadata_cd8, "cd8_DEG_DMG_MS")

###################
## 1. GAM on CD4
###################


# n<-1
# 
# g <- count_cd4[n,]
# t <- metadata_cd4$TimePoint
# t <- gsub("Before", 0, t); t <- gsub("1st", 1, t); t <- gsub("2nd", 2, t); t <- gsub("3rd", 3, t); t <- gsub("Postpartum", 4, t); 
# t <- as.numeric
# 
# df <- data.frame(g,t)
# 
# mod_gam <- gam(g ~ s(t, bs="cr"), data=df)
# summary(mod_gam)

##############################
## 1. Clustering
##############################

## clustering on all CD4
do_clustering(count_cd4,metadata_cd4,5)

#HP
is_HP <- which(is.element(metadata_cd4$Disease, "HP"))
count_HP <- count_cd4[, is_HP]
metadata_HP <- metadata_cd4[is_HP,]

do_clustering(count_HP,metadata_HP,4)

#MS
is_MS <- which(is.element(metadata_cd4$Disease, "MS"))
count_MS <- count_cd4[, is_MS]
metadata_MS <- metadata_cd4[is_MS,]

do_clustering(count_MS,metadata_MS,4)

#select only DEG between in MS vs HP
idx <- is.element(rownames(count_cd4),DEG_CD4)

pdf(paste0(figure_folder,"clustering_HP_CD4.pdf"), width = 8, height = 5)
do_clustering(count_HP[idx,],metadata_HP,8) #HP
dev.off()

pdf(paste0(figure_folder,"clustering_MS_CD4.pdf"), width = 8, height = 5)
do_clustering(count_MS[idx,],metadata_MS,8) #MS
dev.off()

## clustering on CD8
do_clustering(count_cd8,metadata_cd8,10)

#HP
is_HP <- which(is.element(metadata_cd8$Disease, "HP"))
count_HP <- count_cd8[, is_HP]
metadata_HP <- metadata_cd8[is_HP,]

#MS
is_MS <- which(is.element(metadata_cd8$Disease, "MS"))
count_MS <- count_cd8[, is_MS]
metadata_MS <- metadata_cd8[is_MS,]

#select only DEG between in MS vs HP
idx <- is.element(rownames(count_cd8),DEG_CD8)

pdf(paste0(figure_folder,"clustering_HP_CD8.pdf"), width = 8, height = 5)
do_clustering(count_HP[idx,],metadata_HP,8) #HP
dev.off()

pdf(paste0(figure_folder,"clustering_MS_CD8.pdf"), width = 8, , height = 5)
do_clustering(count_MS[idx,],metadata_MS,8) #MS
dev.off()
