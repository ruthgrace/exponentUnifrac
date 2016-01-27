#!/usr/bin/env Rscript 
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata
# run: nohup Rscript ./tongue_rarefaction_script.r > tongue_rarefaction_script_nohup.out 2>&1&

source("UniFrac_memory_efficient.r")
library(ape)
library(phangorn)
library(vegan)

# read and root tree (rooted tree is required)
tree <- read.tree("data/tongue_rarefaction_data/fasttree_all_seed_OTUs.tre")
tree <- midpoint(tree)
tree$tip.label <- gsub("'","",tree$tip.label)

for (i in 51:101) {
  data <- read.table(paste("data/tongue_rarefaction_data/otus_rarefied_",i,".txt",sep=""),sep="\t",header=TRUE,check.names=FALSE,comment.char="",quote="",skip=1,row.names=1)
  
  data <- t(data)
  
  all_distance_matrices <- getDistanceMatrix(data,tree,method="all",verbose=TRUE)
  
  information <- all_distance_matrices[["information"]]
  write.table(information,file=paste("tongue_rarefaction_output/otus_rarefied_",i,"_information_distance_matrix.txt",sep=""),sep="\t",quote=FALSE)
  ratio_no_log <- all_distance_matrices[["ratio_no_log"]]
  write.table(ratio_no_log,file=paste("tongue_rarefaction_output/otus_rarefied_",i,"_ratio_no_log_normalize_distance_matrix.txt",sep=""),sep="\t",quote=FALSE)
  
}
