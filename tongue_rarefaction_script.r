#!/usr/bin/env Rscript 
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata
# run: nohup Rscript ./tongue_rarefaction_script.r > tongue_rarefaction_script_nohup.out 2>&1&

source("UniFrac.r")
library(ape)
library(phangorn)
library(vegan)

# read and root tree (rooted tree is required)
tree <- read.tree("data/tongue_rarefaction_data/fasttree_all_seed_OTUs.tre")
tree <- midpoint(tree)
tree$tip.label <- gsub("'","",tree$tip.label)

for (i in c(1:101)) {
  data <- read.table(paste("data/tongue_rarefaction_data/otus_rarefied_",i,".txt",sep=""),sep="\t",header=TRUE,check.names=FALSE,comment.char="",quote="",skip=1,row.names=1)
  data <- t(data)
    
  tips <- tree$tip.label
  tips <- tips[which(!(tree$tip.label %in% colnames(data)))]
  new.tree <- tree
  new.tree <- drop.tip(new.tree, tips)
  
  all_distance_matrices <- getDistanceMatrix(data,new.tree,method="all",verbose=TRUE)
  
  ratio <- all_distance_matrices[["ratio"]]
  write.table(ratio,file=paste("tongue_rarefaction_output/otus_rarefied_",i,"_ratio_distance_matrix.txt",sep=""),sep="\t",quote=FALSE)
  
}
