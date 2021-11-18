# GraMS

Description of the GraMS (graviditet in MS).


## Introduction 

This document is meant to be a journal about the RNAseq processing. 

In the first section we describe the preprocessing of the data (remove batch effects and normalize). In the second we present our differential analysis. 

 

## Preprocessing 

Remove batch effects 

Many factors could influence the data, such as disease, cell type, etc. 

To quantify the effect of each we calculate p-values that express if a given factor is very similar to a given PCA component. This plot can be found in ‘figures/pca_before_norm.pdf’. 

From this plot we notice that ‘Library_Batch’ has a significant p-value wrt the second PCA, thus it is worth to adjust its effect. 

We use Combat_seq (https://rdrr.io/bioc/sva/man/ComBat_seq.html) to do batch correction, which was built precisely for RNAseq data. 

 
During the normalization, we make sure that the core of the signal (i.e., variables ‘Disease’, ‘State’ and ‘Sample_type’ in the metadata) is conserved by using the option ‘covar_mod’ in Combat_seq.  

![](RNAseq/figures/pca_before_norm.png?raw=true)
