# This script performes batch correction of RNA-seq data
# Created and modified by Alberto Zenere, 2021-04-11

# 2. Batch correction with Combat_seq ####

# 2.1 Setup
# 2.2 Define functions that will be used later on
# 2.1 Import files and re-format
# 2.2 PCA before batch correction
# 2.3 batch correction with Combat_seq 
# 2.4 PCA after batch correction

rm(list=ls()) # remove all entries in the global environment 


# 2.1 Setup  ####

#Set the virtual environment folder as R library
pathlist <- .libPaths();
newpath <-  "C:/Users/albze08/.conda/envs/graMS/Lib/R/library"
.libPaths(newpath)
setwd('C:/Users/albze08/Desktop/phd/P4/RNAseq')

#Load packages
pack_R <- c("sva","tidyverse","readxl","assertthat","ggfortify","patchwork")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}
set.seed(541)

#Set important folders
RDS_folder <- "RDS/"
data_folder <- "data/"
figure_folder <- "figures/"


# 2.2 Define functions that will be used later on ####

#function to create heatmap taken from https://rdrr.io/bioc/ChAMP/src/R/champ.SVD.R
drawheatmap <- function(svdPV.m)
{
  myPalette <- c("darkred","red","orange","pink","white");
  breaks.v <- c(-10000,-10,-5,-2,log10(0.05),0);
  image(x=1:nrow(svdPV.m), y=1:ncol(svdPV.m), z=log10(svdPV.m), col=myPalette, breaks=breaks.v, xlab="", ylab="", axes=FALSE, main="Principal Component Analysis (PCA)");
  axis(1,at=1:nrow(svdPV.m),labels=paste("PCA",1:nrow(svdPV.m),sep=""),las=2);
  axis(2,cex.axis=0.25);
  suppressWarnings(axis(2,at=1:ncol(svdPV.m),labels=colnames(svdPV.m),las=2));
  legend(x=-5.5,y=6,legend=c(expression("p < 1x"~10^{-10}),expression("p < 1x"~10^{-5}),"p < 0.01", "p < 0.05", "p > 0.05"), fill=c("darkred","red","orange","pink", "white"),par('usr')[2], par('usr')[4], xpd=NA);
}


# 2.1 Import files and re-format ####

#Import
metadata <- read_excel(paste0(data_folder, "Md.xlsx"))
count <- read.csv(paste0(data_folder,"salmon.merged.gene_counts.tsv"), sep = "")
rownames(count) <- count$gene_id
count <- as.matrix(count[,-1])

#remove _R1 suffix
colnames(count) <- gsub("_R1", "", colnames(count) )

#count matrix contains some samples that have low quality
tokeep <- which(is.element(colnames(count), metadata$NGI_ID))
count <- count[,tokeep]

#count is in alphabetical order whereas metadata is not
metadata <- metadata[order(metadata$NGI_ID),]

#check that they have the same order
stopifnot(colnames(count) == metadata$NGI_ID)

#save
saveRDS(metadata, file = paste0(RDS_folder, "metadata.RDS"))
saveRDS(count, file = paste0(RDS_folder, "count.RDS"))


# 2.2 PCA before batch correction ####

#Load
metadata <- readRDS(paste0(RDS_folder, "metadata.RDS"))
count <- readRDS(paste0(RDS_folder, "count.RDS"))

#PCA
count_log <- log(count + 1)
count_pca <- prcomp(t(count_log))

#Plot the first two components
pdf(paste0(figure_folder, "pca12_before_norm.pdf"), width = 12, height = 5)
p1 <- autoplot(count_pca, data = metadata, colour = 'Library_Batch')
p2 <- autoplot(count_pca, data = metadata, colour = 'State')
p3 <- autoplot(count_pca, data = metadata, colour = 'Sample_Group')
p4 <- autoplot(count_pca, data = metadata, colour = 'Sample_Type')
p1 + p2 + p3 + p4
dev.off()

#Each phenotype column is assigned a p-value that correspond to its correlation with each PCA
topPCA <- min(10, ncol(count_pca$x))
pheno <- metadata[,-1]
pca_pval <-  matrix(nrow=topPCA,ncol=ncol(pheno))

for(c in 1:topPCA){
  for(f in 1:ncol(pheno)){
    if(typeof(metadata[[f]])!="numeric"){
      pca_pval[c,f] <- kruskal.test(count_pca$x[,c] ~ as.factor(pheno[[f]]))$p.value 
    }
    else{
      pca_pval[c,f] <- summary(lm(count_pca$x[,c] ~ pheno[[f]]))$coeff[2,4];
    }}}
colnames(pca_pval) <- colnames(pheno)
rownames(pca_pval) <- paste0("PCA",1:topPCA)


# Plot p-values of each phenotype and the variance explained by each PCA ####
#tiff(paste0(figure_folder, "pca_before_norm.tiff"), width=550)
pdf(paste0(figure_folder, "pca_before_norm.pdf"), width=7.3)
layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))

#Pvalue for each phenotype and PCA
par(mar = c(1,16,2,4)) #bottom, left, top, right margins of the plot region
drawheatmap(pca_pval)

#Variance Explained
par(mar = c(4,11,3,2))
plot(count_pca$sdev[1:topPCA]^2/sum(count_pca$sdev^2), type = "b", xlab="", ylab="Variance Explained", xaxt='n')
axis(1,at=1:topPCA,labels=paste("PCA",1:topPCA,sep=""),las=2)

dev.off()


# 1.3 batch correction with Combat_seq ####

covar_mat <- data.frame(cbind(metadata$Sample_Group, metadata$State)) # variables to keep
colnames(covar_mat) <- c("Sample_Group", "State")
count_adjusted <- ComBat_seq(count, batch=metadata$Library_Batch, group=NULL, covar_mod = covar_mat) #run Combat-seq on the count matrix

#save
saveRDS(count_adjusted, file = paste0(RDS_folder, "count_adjusted.RDS"))


# 1.4 PCA after batch correction ####

count_adjusted <- readRDS(paste0(RDS_folder, "count_adjusted.RDS"))

count_log_adjusted <- log(count_adjusted+1)
count_pca <- prcomp(t(count_log_adjusted))

#Plot the first two components
pdf(paste0(figure_folder, "pca12_after_norm.pdf"), width = 12, height = 5)
p1 <- autoplot(count_pca, data = metadata, colour = 'Library_Batch')
p2 <- autoplot(count_pca, data = metadata, colour = 'State')
p3 <- autoplot(count_pca, data = metadata, colour = 'Sample_Group')
p4 <- autoplot(count_pca, data = metadata, colour = 'Sample_Type')
p1 + p2 + p3 + p4
dev.off()

#Each phenotype column is assigned a p-value that correspond to its correlation with each PCA
topPCA <- min(10, ncol(count_pca$x))
pheno <- metadata[,-1]
pca_pval <-  matrix(nrow=topPCA,ncol=ncol(pheno))

for(c in 1:topPCA){
  for(f in 1:ncol(pheno)){
    if(typeof(metadata[[f]])!="numeric"){
      pca_pval[c,f] <- kruskal.test(count_pca$x[,c] ~ as.factor(pheno[[f]]))$p.value 
    }
    else{
      pca_pval[c,f] <- summary(lm(count_pca$x[,c] ~ pheno[[f]]))$coeff[2,4];
    }}}
colnames(pca_pval) <- colnames(pheno)
rownames(pca_pval) <- paste0("PCA",1:topPCA)


#Plot p-values of each phenotype and the variance explained by each PCA

tiff(paste0(figure_folder, "pca_after_norm.tiff"), height = 500, width=495.5)
#pdf(paste0(figure_folder, "pca_after_norm.pdf"), width=6.8)
layout(matrix(c(1,1,2,2), 2, 2, byrow = TRUE))

#Pvalue for each phenotype and PCA
par(mar = c(1,14,2,4)) #bottom, left, top, right margins of the plot region
drawheatmap(pca_pval)

#Variance Explained
par(mar = c(4,11,3,2))
plot(count_pca$sdev[1:topPCA]^2/sum(count_pca$sdev^2), type = "b", xlab="", ylab="Variance Explained", xaxt='n')
axis(1,at=1:topPCA,labels=paste("PCA",1:topPCA,sep=""),las=2)

dev.off()




