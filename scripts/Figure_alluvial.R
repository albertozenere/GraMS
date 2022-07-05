
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
library("ggalluvial")

set.seed(206)

#Load ProbeFeatures
ProbeFeatures <- readRDS(file=paste0(FOLDER_RDS, "TH2636_ProbeFeatures_211123.RDS"))

# Load methyl DMP ####
CD4_3rd_HP_methyl <- readRDS(paste0(FOLDER_RDS, "DMR/CD4_3rd_HP_methyl.RDS"))
CD4_PP_HP_methyl <- readRDS(paste0(FOLDER_RDS, "DMR/CD4_PP_HP_methyl.RDS"))
CD4_3rd_MS_methyl <- readRDS(paste0(FOLDER_RDS, "DMR/CD4_3rd_MS_methyl.RDS"))
CD4_PP_MS_methyl <- readRDS(paste0(FOLDER_RDS, "DMR/CD4_PP_MS_methyl.RDS"))

CD8_3rd_HP_methyl <- readRDS(paste0(FOLDER_RDS, "DMR/CD8_3rd_HP_methyl.RDS"))
CD8_PP_HP_methyl <- readRDS(paste0(FOLDER_RDS, "DMR/CD8_PP_HP_methyl.RDS"))
CD8_3rd_MS_methyl <- readRDS(paste0(FOLDER_RDS, "DMR/CD8_3rd_MS_methyl.RDS"))
CD8_PP_MS_methyl <- readRDS(paste0(FOLDER_RDS, "DMR/CD8_PP_MS_methyl.RDS"))

# Load all methyl cpgs ####
CD4_3rd_HP_methyl_all <- readRDS(paste0(FOLDER_RDS, "DMR/CD4_3rd_HP_methyl_all.RDS"))
CD4_PP_HP_methyl_all <- readRDS(paste0(FOLDER_RDS, "DMR/CD4_PP_HP_methyl_all.RDS"))
CD4_3rd_MS_methyl_all <- readRDS(paste0(FOLDER_RDS, "DMR/CD4_3rd_MS_methyl_all.RDS"))
CD4_PP_MS_methyl_all <- readRDS(paste0(FOLDER_RDS, "DMR/CD4_PP_MS_methyl_all.RDS"))

CD8_3rd_HP_methyl_all <- readRDS(paste0(FOLDER_RDS, "DMR/CD8_3rd_HP_methyl_all.RDS"))
CD8_PP_HP_methyl_all <- readRDS(paste0(FOLDER_RDS, "DMR/CD8_PP_HP_methyl_all.RDS"))
CD8_3rd_MS_methyl_all <- readRDS(paste0(FOLDER_RDS, "DMR/CD8_3rd_MS_methyl_all.RDS"))
CD8_PP_MS_methyl_all <- readRDS(paste0(FOLDER_RDS, "DMR/CD8_PP_MS_methyl_all.RDS"))


# Load RNA-seq DMG ####
CD4_3rd_HP_rna <- readRDS(paste0(FOLDER_RDS, "DMG/CD4_3rd_HP_rna.RDS"))
CD4_PP_HP_rna <- readRDS(paste0(FOLDER_RDS, "DMG/CD4_PP_HP_rna.RDS"))
CD4_3rd_MS_rna <- readRDS(paste0(FOLDER_RDS, "DMG/CD4_3rd_MS_rna.RDS"))
CD4_PP_MS_rna <- readRDS(paste0(FOLDER_RDS, "DMG/CD4_PP_MS_rna.RDS"))

CD8_3rd_HP_rna <- readRDS(paste0(FOLDER_RDS, "DMG/CD8_3rd_HP_rna.RDS"))
CD8_PP_HP_rna <- readRDS(paste0(FOLDER_RDS, "DMG/CD8_PP_HP_rna.RDS"))
CD8_3rd_MS_rna <- readRDS(paste0(FOLDER_RDS, "DMG/CD8_3rd_MS_rna.RDS"))
CD8_PP_MS_rna <- readRDS(paste0(FOLDER_RDS, "DMG/CD8_PP_MS_rna.RDS"))


# Load all RNA-seq genes ####
CD4_3rd_HP_rna_all <- readRDS(paste0(FOLDER_RDS, "DMG/CD4_3rd_HP_rna_all.RDS"))
CD4_PP_HP_rna_all <- readRDS(paste0(FOLDER_RDS, "DMG/CD4_PP_HP_rna_all.RDS"))
CD4_3rd_MS_rna_all <- readRDS(paste0(FOLDER_RDS, "DMG/CD4_3rd_MS_rna_all.RDS"))
CD4_PP_MS_rna_all <- readRDS(paste0(FOLDER_RDS, "DMG/CD4_PP_MS_rna_all.RDS"))

CD8_3rd_HP_rna_all <- readRDS(paste0(FOLDER_RDS, "DMG/CD8_3rd_HP_rna_all.RDS"))
CD8_PP_HP_rna_all <- readRDS(paste0(FOLDER_RDS, "DMG/CD8_PP_HP_rna_all.RDS"))
CD8_3rd_MS_rna_all <- readRDS(paste0(FOLDER_RDS, "DMG/CD8_3rd_MS_rna_all.RDS"))
CD8_PP_MS_rna_all <- readRDS(paste0(FOLDER_RDS, "DMG/CD8_PP_MS_rna_all.RDS"))



alluvial_plot_methyl <- function(beta_values,cpg){
  stopifnot(length(cpg)==nrow(beta_values))
  
  #make initial data.frame on cpgs properties ####
  
  #3rd-1st
  methyl_3rd_1st <- rep("", nrow(beta_values))
  methyl_3rd_1st[beta_values$Third_First>0] <- "Hyper"
  methyl_3rd_1st[beta_values$Third_First<0] <- "Hypo"
  methyl_3rd_1st <- as.character(methyl_3rd_1st)
  
  #PP-3rd
  methyl_PP_3rd <- rep("", nrow(beta_values))
  methyl_PP_3rd[beta_values$PP_Third>0] <- "Hyper"
  methyl_PP_3rd[beta_values$PP_Third<0] <- "Hypo"
  methyl_PP_3rd <- as.character(methyl_PP_3rd)
  
  #combination of 3rd-1st and PP-3rd
  group <- rep("", nrow(beta_values))
  group[methyl_3rd_1st=="Hyper" & methyl_PP_3rd=="Hypo"] <- "Hyper-Hypo"
  group[methyl_3rd_1st=="Hyper" & methyl_PP_3rd=="Hyper"] <- "Hyper-Hyper"
  group[methyl_3rd_1st=="Hypo" & methyl_PP_3rd=="Hyper"] <- "Hypo-Hyper"
  group[methyl_3rd_1st=="Hypo" & methyl_PP_3rd=="Hypo"] <- "Hypo-Hypo"
  
  #genomic region of cpgs
  ord <- match(cpg, rownames(ProbeFeatures))
  region <- ProbeFeatures$feature[ord] %>% as.character()
  region <- as.character(region)
  
  #create data.frame on cpgs
  dat <- data.frame(Third_First=methyl_3rd_1st, PP_Third=methyl_PP_3rd,
                    region = region, group = group)
  
  
  #Make all combinations of properties ####
  n_var <- ncol(dat)
  
  list_var <- vector(mode = "list", n_var)
  for (n in 1:n_var){
    list_var[[n]] <- unique(dat[,n])
  }
  
  df <- expand.grid(list_var)
  colnames(df) <- colnames(dat)
  df$Freq <- 0
  
  #count frequency of each combination
  for (i in 1:nrow(df)){
    idx_mat <- matrix(FALSE, nrow(dat), ncol(df)-1) #last column is frequency
    
    for (j in 1:ncol(idx_mat)){
      idx_mat[,j] <- dat[,j] == df[i,j]
    }
    df$Freq[i] <- sum(apply(idx_mat,1,all))
  }
  
  
  #Plot
  p <- ggplot(df, aes(y=Freq, axis2=Third_First , axis3=PP_Third)) + 
    geom_alluvium(aes(fill = group), width = 1/12) +
    geom_stratum(width = 1/8, fill = "white", color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              reverse = T) +
    scale_x_discrete(limits = c("3rd-1st", "PP-3rd"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ylab(NULL) + 
    theme_bw() +
    theme(legend.position="none",  
          axis.ticks.y = element_line(),
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() )
  #print(p)
  
  return(p)               
  
}


alluvial_plot_rna <- function(rna_values){

  #make initial data.frame on cpgs properties ####
  
  #3rd-1st
  rna_3rd_1st <- rep("", nrow(rna_values))
  rna_3rd_1st[rna_values$Third_First>0] <- "Upreg"
  rna_3rd_1st[rna_values$Third_First<0] <- "Downreg"
  rna_3rd_1st <- as.character(rna_3rd_1st)
  
  #PP-3rd
  rna_PP_3rd <- rep("", nrow(rna_values))
  rna_PP_3rd[rna_values$PP_Third>0] <- "Upreg"
  rna_PP_3rd[rna_values$PP_Third<0] <- "Downreg"
  rna_PP_3rd <- as.character(rna_PP_3rd)
  
  #combination of 3rd-1st and PP-3rd
  group <- rep("", nrow(rna_values))
  group[rna_3rd_1st=="Upreg" & rna_PP_3rd=="Downreg"] <- "Upreg-Downreg"
  group[rna_3rd_1st=="Upreg" & rna_PP_3rd=="Upreg"] <- "Upreg-Upreg"
  group[rna_3rd_1st=="Downreg" & rna_PP_3rd=="Upreg"] <- "Downreg-Upreg"
  group[rna_3rd_1st=="Downreg" & rna_PP_3rd=="Downreg"] <- "Downreg-Downreg"

  
  #create data.frame on cpgs
  dat <- data.frame(Third_First=rna_3rd_1st, PP_Third=rna_PP_3rd,
                    group = group)
  
  
  #Make all combinations of properties ####
  n_var <- ncol(dat)
  
  list_var <- vector(mode = "list", n_var)
  for (n in 1:n_var){
    list_var[[n]] <- unique(dat[,n])
  }
  
  df <- expand.grid(list_var)
  colnames(df) <- colnames(dat)
  df$Freq <- 0
  
  #count frequency of each combination
  for (i in 1:nrow(df)){
    idx_mat <- matrix(FALSE, nrow(dat), ncol(df)-1) #last column is frequency
    
    for (j in 1:ncol(idx_mat)){
      idx_mat[,j] <- dat[,j] == df[i,j]
    }
    df$Freq[i] <- sum(apply(idx_mat,1,all))
  }
  
  
  #Plot
  p <- ggplot(df, aes(y=Freq, axis2=Third_First , axis3=PP_Third)) + 
    geom_alluvium(aes(fill = group), width = 1/12) +
    geom_stratum(width = 1/8, fill = "white", color = "grey") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)),
              reverse = T) +
    scale_x_discrete(limits = c("3rd-1st", "PP-3rd"), expand = c(.05, .05)) +
    scale_fill_brewer(type = "qual", palette = "Set1") +
    ylab(NULL) + 
    theme_bw() +
    theme(legend.position="none",  
          axis.ticks.y = element_line(),
          axis.text.y=element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank() )

  
  return(p)               
  
}

# RNA-seq ####
#CD4 HP
deg <- union(rownames(CD4_3rd_HP_rna), rownames(CD4_PP_HP_rna))
values_CD4_HP <- cbind(CD4_3rd_HP_rna_all[deg,"logFC"],
                       CD4_PP_HP_rna_all[deg,"logFC"]) %>% data.frame()
colnames(values_CD4_HP) <- c("Third_First", "PP_Third")

p_CD4_HP_rna <- alluvial_plot_rna(values_CD4_HP)

#CD4 MS
deg <- union(rownames(CD4_3rd_MS_rna), rownames(CD4_PP_MS_rna))
values_CD4_MS <- cbind(CD4_3rd_MS_rna_all[deg,"logFC"],
                       CD4_PP_MS_rna_all[deg,"logFC"]) %>% data.frame()
colnames(values_CD4_MS) <- c("Third_First", "PP_Third")

p_CD4_MS_rna <- alluvial_plot_rna(values_CD4_MS)

#CD8 HP
deg <- union(rownames(CD8_3rd_HP_rna), rownames(CD8_PP_HP_rna))
values_CD8_HP <- cbind(CD8_3rd_HP_rna_all[deg,"logFC"],
                       CD8_PP_HP_rna_all[deg,"logFC"]) %>% data.frame()
colnames(values_CD8_HP) <- c("Third_First", "PP_Third")

p_CD8_HP_rna <- alluvial_plot_rna(values_CD8_HP)

#CD8 MS
deg <- union(rownames(CD8_3rd_MS_rna), rownames(CD8_PP_MS_rna))
values_CD8_MS <- cbind(CD8_3rd_MS_rna_all[deg,"logFC"],
                       CD8_PP_MS_rna_all[deg,"logFC"]) %>% data.frame()
colnames(values_CD8_MS) <- c("Third_First", "PP_Third")

p_CD8_MS_rna <- alluvial_plot_rna(values_CD8_MS)




# Methylation ####
#CD4 HP 
cpg <- union(rownames(CD4_3rd_HP_methyl), rownames(CD4_PP_HP_methyl))
beta_values_CD4_HP <- cbind(CD4_3rd_HP_methyl_all[cpg,"logFC"],
                            CD4_PP_HP_methyl_all[cpg,"logFC"]) %>% data.frame()
colnames(beta_values_CD4_HP) <- c("Third_First", "PP_Third")

p_CD4_HP_methyl <- alluvial_plot_methyl(beta_values_CD4_HP,cpg)

#CD4 MS 
cpg <- union(rownames(CD4_3rd_MS_methyl), rownames(CD4_PP_MS_methyl))
beta_values_CD4_MS <- cbind(CD4_3rd_MS_methyl_all[cpg,"logFC"],
                            CD4_PP_MS_methyl_all[cpg,"logFC"]) %>% data.frame()
colnames(beta_values_CD4_MS) <- c("Third_First", "PP_Third")

p_CD4_MS_methyl <- alluvial_plot_methyl(beta_values_CD4_MS,cpg)

#CD8 HP 
cpg <- union(rownames(CD8_3rd_HP_methyl), rownames(CD8_PP_HP_methyl))
beta_values_CD8_HP <- cbind(CD8_3rd_HP_methyl_all[cpg,"logFC"],
                            CD8_PP_HP_methyl_all[cpg,"logFC"]) %>% data.frame()
colnames(beta_values_CD8_HP) <- c("Third_First", "PP_Third")

p_CD8_HP_methyl <- alluvial_plot_methyl(beta_values_CD8_HP,cpg)

#CD8 MS 
cpg <- union(rownames(CD8_3rd_MS_methyl), rownames(CD8_PP_MS_methyl))
beta_values_CD8_MS <- cbind(CD8_3rd_MS_methyl_all[cpg,"logFC"],
                            CD8_PP_MS_methyl_all[cpg,"logFC"]) %>% data.frame()
colnames(beta_values_CD8_MS) <- c("Third_First", "PP_Third")

p_CD8_MS_methyl <- alluvial_plot_methyl(beta_values_CD8_MS,cpg)


# Combine plots ####
pdf("figures_manus/alluvial.pdf", width=19, height=6)
ggarrange(p_CD4_HP_methyl, p_CD4_MS_methyl, p_CD8_HP_methyl, p_CD8_MS_methyl, 
          p_CD4_HP_rna, p_CD4_MS_rna, p_CD8_HP_rna, p_CD8_MS_rna, 
          nrow=2, ncol=4)
dev.off()

