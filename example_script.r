#!/usr/bin/env Rscript
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata

source("UniFrac.r")
library(ape)
library(phangorn)

otu.tab <- read.table("data/td_OTU_tag_mapped_lineage.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- otu.tab$taxonomy
otu.tab <- otu.tab[-length(colnames(otu.tab))]
otu.tab <- t(as.matrix(otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
otu.tab <- otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
tree <- read.tree("data/fasttree_all_seed_OTUs.tre")
tree <- midpoint(tree)

# read metadata
MyMeta<- read.table("data/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

# filter OTU table and metadata so that only samples which appear in both are retained
otu_indicies <- match(rownames(MyMeta),rownames(otu.tab))
otu_indicies <- otu_indicies[!is.na(otu_indicies)]
otu.tab <- otu.tab[otu_indicies,]
MyMetaOrdered <- MyMeta[match(rownames(otu.tab),rownames(MyMeta)),]

#calculate distance matrix
unweighted <- getDistanceMatrix(otu.tab,tree,method="unweighted",verbose=TRUE)
weighted <- getDistanceMatrix(otu.tab,tree,method="weighted",verbose=TRUE)
information <- getDistanceMatrix(otu.tab,tree,method="information",verbose=TRUE)
exponent <- getDistanceMatrix(otu.tab,tree,method="exponent",verbose=TRUE,normalize=FALSE)

#output distance matrices
write.table(unweighted,file="output/unweighted_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(weighted,file="output/weighted_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(information,file="output/information_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(exponent,file="output/exponent_distance_matrix.txt",sep="\t",quote=FALSE)


groups <- MyMetaOrdered$diagnosis

unweighted.pcoa <- pcoa(unweighted)
weighted.pcoa <- pcoa(weighted)
information.pcoa <- pcoa(information)
exponent.pcoa <- pcoa(exponent)


#function to get variance explained for the PCOA component labels
getVarExplained <- function(vector) {
	rawVarEx <- apply(vector,2,function(x) sd(x)*sd(x))
	totalVarExplained <- sum(rawVarEx)
	varEx <- rawVarEx/totalVarExplained
	return(varEx)
}


unweighted.varEx <- getVarExplained(unweighted.pcoa$vectors)
weighted.varEx <- getVarExplained(weighted.pcoa$vectors)
information.varEx <- getVarExplained(information.pcoa$vectors)
exponent.varEx <- getVarExplained(exponent.pcoa$vectors)

#get vector version of distance matrices for correlation plots below
unweighted.vector <- unlist(unweighted[lower.tri(unweighted,diag=TRUE)])
weighted.vector <- unlist(weighted[lower.tri(weighted,diag=TRUE)])
information.vector <- unlist(information[lower.tri(information,diag=TRUE)])
exponent.vector <- unlist(exponent[lower.tri(exponent,diag=TRUE)])



pdf("output/pcoa_plots.pdf")

#plot pcoa plots
plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(exponent.pcoa$vectors[,1],exponent.pcoa$vectors[,2], col=groups,main="Exponent UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(exponent.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(exponent.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

#plot correlation between different UniFrac modes

plot(unweighted.vector,information.vector,main="unweighted vs. information UniFrac")
plot(weighted.vector,information.vector,main="weighted vs. information UniFrac")
plot(unweighted.vector,weighted.vector,main="unweighted vs. weighted UniFrac")
plot(unweighted.vector,exponent.vector,main="unweighted vs. exponent UniFrac")
plot(weighted.vector,exponent.vector,main="weighted vs. exponent UniFrac")
plot(information.vector,exponent.vector,main="information vs. exponent UniFrac")


dev.off()

