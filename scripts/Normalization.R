###############################################################################
# Normalization with TMM
###############################################################################

# Created and modified by Alberto Zenere, 2021-04-11

# 0.1 Setup
# 0.2 Convert id
# 0.3 Split CD4 and CD8
# 0.4 Remove lowly expressed genes in CD4
# 0.5 Remove lowly expressed genes in CD8
# 0.6 Save
# 1.1 Calculate TMM on CD4
# 1.2 Calculate TMM on CD8

rm(list=ls()) # remove all entries in the global environment 

################
## 0.1 Setup  ##
################

#Set the virtual environment folder as R library
pathlist <- .libPaths();
newpath <-  "C:/Users/albze08/.conda/envs/graMS/Lib/R/library"
.libPaths(newpath)
setwd('C:/Users/albze08/Desktop/phd/P4/RNAseq')

#Load packages
pack_R <- c("biomaRt","tidyverse","readxl","patchwork",
            "assertthat","edgeR","sva","ggfortify")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}
set.seed(541)

#Set important folders
RDS_folder <- "RDS/"
data_folder <- "data/"
figure_folder <- "figures/"

####################
## 0.2 Convert id ##
####################

count <- readRDS(paste0(RDS_folder, "count_adjusted.RDS")) 
metadata <- readRDS(paste0(RDS_folder, "metadata.RDS"))

rownames(count) <- gsub("\\..*","",rownames(count)) #the digits after the dot simply indicate the version in Gencode

#convert id using biomaRt
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

ensg_symbol <- getBM(attributes <- c("ensembl_gene_id" ,"hgnc_symbol"), filters = "ensembl_gene_id", values = rownames(count), mart = mart)
colnames(ensg_symbol) <- c("ensg","symbol")
ensg_symbol <- ensg_symbol[!duplicated(ensg_symbol$ensg),]

#remove ensg not in ensg_symbol
count <- count[is.element(rownames(count),ensg_symbol$ensg),]

#make sure that count and ensg_symbol are in the same order
stopifnot( all(rownames(count)==ensg_symbol$ensg) )

#remove genes without a symbol id
idx <- !(ensg_symbol$symbol=="")
count <- count[idx,]
ensg_symbol <- ensg_symbol[idx,]
stopifnot( all(rownames(count)==ensg_symbol$ensg) )

#multiple genes map to the same symbol, we sum them together
rownames(count) <- ensg_symbol$symbol
count <- t(sapply(by(count,rownames(count),colSums),identity))

#PCA
par(mar=c(1,1,1,1))
pdf(paste0(figure_folder, "pca12_after_2.pdf"))
count_pca <- prcomp(t(log2(count+1)))
p1 <- autoplot(count_pca, data = metadata, colour = 'TimePoint', size = 3)
p2 <- autoplot(count_pca, data = metadata, colour = 'Disease', size = 3)
p3 <- autoplot(count_pca, data = metadata, colour = 'Sample_Type', size = 3)
p4 <- autoplot(count_pca, data = metadata, colour = 'State', size = 3)

p1 + p2 + p3 + p4
dev.off()

###########################
## 0.3 Split CD4 and CD8 ##  
###########################

#index
idx_cd4 <- grep("CD4", metadata$Sample_Type)
idx_cd8 <- grep("CD8", metadata$Sample_Type)

#split
count_cd8 <- count[,idx_cd8]
metadata_cd8 <- metadata[idx_cd8,]
metadata_cd8$Individual <- as.character(metadata_cd8$Individual)

count_cd4 <- count[,idx_cd4]
metadata_cd4 <- metadata[idx_cd4,]
metadata_cd4$Individual <- as.character(metadata_cd4$Individual)

#############################################
## 0.4 Remove lowly expressed genes in CD4 ##  
#############################################

par(mar=c(1,1,1,1))
dev.off()

pdf(paste0(figure_folder, "density_beforeafter_filtering_CD4.pdf"))

#density before removing
par(mfrow=c(1,2))
plot(density(log2(count_cd4[,1]+1)), main = "Before Filtering")
for (n in 2:ncol(count_cd4)){
  lines(density(log2(count_cd4[,n]+1)))
}

#remove
tokeep_cd4 <- filterByExpr(count_cd4)
count_cd4 <- count_cd4[tokeep_cd4,]

#density after removing
plot(density(log2(count_cd4[,1]+1)), main = "After Filtering")
for (n in 2:ncol(count_cd4)){
  lines(density(log2(count_cd4[,n]+1)))
}

dev.off()

#############################################
## 0.5 Remove lowly expressed genes in CD8 
#############################################

par(mar=c(1,1,1,1))
dev.off()

pdf(paste0(figure_folder, "density_beforeafter_filtering_CD8.pdf"))

#density before removing
par(mfrow=c(1,2))
plot(density(log2(count_cd8[,1]+1)), main = "Before Filtering")
for (n in 2:ncol(count_cd8)){
  lines(density(log2(count_cd8[,n]+1)))
}

#remove
tokeep_cd8 <- filterByExpr(count_cd8)
count_cd8 <- count_cd8[tokeep_cd8,]

#density after removing
plot(density(log2(count_cd8[,1]+1)), main = "After Filtering")
for (n in 2:ncol(count_cd8)){
  lines(density(log2(count_cd8[,n]+1)))
}

dev.off()

##############
## 0.6 Save ##  
##############

saveRDS(count_cd8, file = paste0(RDS_folder, "count_cd8.RDS"))
saveRDS(metadata_cd8, file = paste0(RDS_folder, "metadata_cd8.RDS"))

saveRDS(count_cd4, file = paste0(RDS_folder, "count_cd4.RDS"))
saveRDS(metadata_cd4, file = paste0(RDS_folder, "metadata_cd4.RDS"))

##############################
## 1.1 Calculate TMM on CD4 ##  
##############################
#TMM 
count_cd4_tmm <- DGEList(count_cd4)
count_cd4_tmm <- calcNormFactors(count_cd4_tmm, method="TMM")
#count_cd4_tmm <- cpm(count_cd4_tmm, log=F) #no log because voom uses non-log counts

#save
saveRDS(count_cd4_tmm, file = paste0(RDS_folder, "count_cd4_tmm.RDS"))

##############################
## 1.2 Calculate TMM on CD8 ##  
##############################

#TMM 
count_cd8_tmm <- DGEList(count_cd8)
count_cd8_tmm <- calcNormFactors(count_cd8_tmm, method="TMM")
#count_cd8_tmm <- cpm(count_cd8_tmm, log=F) #no log because voom uses non-log counts

#save
saveRDS(count_cd8_tmm, file = paste0(RDS_folder, "count_cd8_tmm.RDS"))


