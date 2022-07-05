probe_to_gene <- function(cpg){

  stopifnot(class(cpg)=="character")
  
  # Load probes annotation ####
  ProbeFeatures <- readRDS(file = "C:/Users/albze08/Desktop/phd/P4/methylation/RDS_files/TH2636_ProbeFeatures.RDS")

  cpg_gene <- data.frame(cpg = rownames(ProbeFeatures), gene = ProbeFeatures$gene)

  rownames(cpg_gene) <- cpg_gene$cpg
  cpg_gene$cpg <- NULL

  # Convert cpg to gene ####
  out <- data.frame( cpg = cpg, gene = rep("", length(cpg)))
  out$gene <- cpg_gene[cpg,"gene"] %>% as.character()
  
  return(out)
  
  
}
  