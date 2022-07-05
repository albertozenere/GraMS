# This script contains the pre-processing of the methylation data
# Tutorial: David Martinez-Enguita (2020) and https://www.bioconductor.org/packages/release/bioc/vignettes/ChAMP/inst/doc/ChAMP.html
# Created and modified by Sandra Hellberg, 2021-08-13

# 1_Differential methylation analysis using ChAMP ####

# 1.1 Importation of methylation raw idat files
# 1.2 Quality control and exploratory analysis
# 1.3 Normalization
# 1.4 Singular Value Decomposition Plot
# 1.5 Batch effect correction
# 1.6 Reference-free deconvolution by Houseman et al.
# 1.7 Differential methylation analysis with Limma


rm(list=ls()) # remove all entries in the global environment 

lib_path <- .libPaths("C:/Rpackages/") # specifing the package library
#update.packages(lib.loc = lib_path, ask = FALSE, checkBuilt = TRUE, dependencies = TRUE)
#BiocManager::install(lib.loc = lib_path, update = TRUE, ask = FALSE)
pack_R <- c("ChAMP", "limma", "RefFreeEWAS")
for (i in 1:length(pack_R)) {
  library(pack_R[i], character.only = TRUE)
}
set.seed(541)

# Set directory structure
setwd("P:/LiU/Postdok/Projekt/GraMS/Methylation/TH-2636")
main_dir <- "P:/LiU/Postdok/Projekt/GraMS/Methylation/TH-2636"

# 1.1 Importation of methylation files ####

# Data set ID, data tag 
array_type <- "EPIC"
dset_tag <- "all"
dset_name <- "TH2636"

# Importation of idat methylation files using champ.load
methyl_raw <- champ.load(directory=(paste0(main_dir,"/RawData/IDAT")), 
                         arraytype      = array_type,
                         methValue      = "B", # Beta-values
                         autoimpute     = TRUE, # Impute values if NA after filtering
                         filterDetP     = TRUE, # All probes above the detPcut will be filtered out
                         ProbeCutoff    = 0.2, # The NA ratio threshhold for probes
                         SampleCutoff   = 0.1, # The failed p value (or NA) threshhold for samples 
                         detPcut        = 0.01, # The detection p-value threshhold
                         filterBeads    = TRUE, # Probes with a beadcount less than 3 will be removed depending on the beadCutoff value
                         beadCutoff     = 0.05, # The fraction of samples that must have a beadcount less than 3 before the probe is remove
                         filterNoCG     = TRUE, # Remove non-cg probes
                         filterSNPs     = TRUE, # Probes in which the probed CpG falls near a SNP as defined in Nordlund et al are removed
                         filterMultiHit = TRUE, # Probes in which the probe aligns to multiple locations wawith bwa as defined in Nordlund et al are removed
                         filterXY       = TRUE, # Probes from X and Y chromosomes are removed
                         force          = FALSE) 

saveRDS(methyl_raw, file = paste0(main_dir, "/Results/RDS_files/Methylation_files/", 
                                  dset_name, "_methyl_raw.RDS"))
methyl_raw <- readRDS(file = paste0(main_dir, "/Results/RDS_files/Methylation_files/", 
                                    dset_name, "_methyl_raw.RDS"))


# 1.2 Quality control and exploratory analysis ####

# QC plot, mds plot, density plot, dendrogram
champ.QC(beta        = methyl_raw$beta, 
         pheno       = methyl_raw$pd$Sample_Group, 
         mdsPlot     = TRUE,
         densityPlot = TRUE,
         dendrogram  = TRUE, 
         PDFplot     = TRUE, 
         Rplot       = FALSE, 
         Feature.sel = "None", 
         resultsDir  = paste0(main_dir, "/Results/Figures"))


mdsPlot(methyl_raw$beta, sampNames = methyl_raw$pd$Sample_Type)


# 1.3 Normalization ####

# Normalization (BMIQ: Beta-Mixture Quantile method from Teschendorff et al. Bioinformatics, 2013
methyl_norm <- champ.norm(beta       = methyl_raw$beta, 
                          rgSet      = NULL,
                          mset       = NULL,
                          method     = "BMIQ", 
                          plotBMIQ   = TRUE,
                          arraytype  = array_type, 
                          cores      = 1,
                          resultsDir = paste0(main_dir, "/Figures"))

# Storage of normalized methylation inputs as RDS files
saveRDS(methyl_norm, file = paste0(main_dir, "/Results/RDS_files/Methylation_files/", 
                                   dset_name, "_methyl_norm.RDS"))

methyl_norm <- readRDS(file = paste0(main_dir, "/Results/RDS_files/Methylation_files/", 
                                     dset_name, "_methyl_norm.RDS"))


# 1.4 Singular value decomposition ####

# Singular Value Decomposition (SVD) analysis for batch effects prediction
# Heatmaps: darker colors indicate stronger correlation of SVD component
# with a factor of interest. If technical factors account for a substantial
# part of variation, batch effect correction (ComBat) is needed.

# Set all variables to character
methyl_raw$pd[] <- lapply(methyl_raw$pd, as.character)

# Set class for the different variables
methyl_raw$pd$Age <- as.numeric(methyl_raw$pd$Age)
methyl_raw$pd$Cell.viability <- as.numeric(methyl_raw$pd$Cell.viability)
methyl_raw$pd$Delivery_Week <- as.numeric(methyl_raw$pd$Delivery_Week)
methyl_raw$pd$Pregnancy_Week <- as.numeric(methyl_raw$pd$Pregnancy_Week)


champ.SVD(beta       = methyl_norm,
          rgSet      = TRUE, 
          pd         = methyl_raw$pd, 
          RGEffect   = FALSE, # Include 18 internal probe controls (only for idat files)
          PDFplot    = TRUE,
          Rplot      = TRUE,
          resultsDir = paste0(main_dir, "/Figures"))



# 1.5 Batch effect correction  ####

# If Batch is detected, run champ.runCombat()
methyl_combat <- champ.runCombat(beta      = methyl_norm,
                                 pd        = methyl_raw$pd,
                                 variablename = "Sample_Group", 
                                 logitTrans = TRUE,
                                 batchname = c("Slide"))

mdsPlot(methyl_combat, sampNames = methyl_raw$pd$ID)
mdsPlot(methyl_combat, sampNames = methyl_raw$pd$Sample_Type)


saveRDS(methyl_combat, file = paste0(main_dir, "/Results/RDS_files/Methylation_files/",
                                     dset_name, "_methyl_combat.RDS"))

methyl_combat <- readRDS(file = paste0(main_dir, "/Results/RDS_files/Methylation_files/", 
                                       dset_name, "_methyl_combat.RDS"))

