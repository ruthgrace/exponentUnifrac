#!/usr/bin/env Rscript
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata

source("UniFrac.r")
library(ape)
library(phangorn)
library(vegan)

otu.tab <- read.table("data/human_mouth_data/hmp_mouth_data.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- rownames(otu.tab)
#make row names samples, and col names OTUs
otu.tab <- t(as.matrix(otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
otu.tab <- otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
tree <- read.tree("data/human_mouth_data/hmp_mouth_subtree.tre")
tree <- midpoint(tree)

# read metadata
MyMetaOrdered <- rownames(otu.tab)
MyMetaOrdered <- gsub("_.*$","",MyMetaOrdered)

#rarefy data for unweighted unifrac
otu.tab.rarefy <- rrarefy(otu.tab, min(apply(otu.tab,1,sum)))

#calculate distance matrix
unweighted <- getDistanceMatrix(otu.tab.rarefy,tree,method="unweighted",verbose=TRUE)
#output distance matrices
write.table(unweighted,file="output/unweighted_distance_matrix.txt",sep="\t",quote=FALSE)
weighted <- getDistanceMatrix(otu.tab,tree,method="weighted",verbose=TRUE)
write.table(weighted,file="output/weighted_distance_matrix.txt",sep="\t",quote=FALSE)
information <- getDistanceMatrix(otu.tab,tree,method="information",verbose=TRUE)
write.table(information,file="output/information_distance_matrix.txt",sep="\t",quote=FALSE)
exponent <- getDistanceMatrix(otu.tab,tree,method="exponent",verbose=TRUE)
write.table(exponent,file="output/exponent_normalize_distance_matrix.txt",sep="\t",quote=FALSE)

#conditions (bv - bacterial vaginosis as scored by nugent/amsel, i - intermediate, n - normal/healthy)
groups <- MyMetaOrdered
originalgroups <- groups
# change conditions so that samples which are more than 50% one taxa are colored by that taxa
otuSum <- apply(otu.tab,1,sum)
otuMax <- apply(otu.tab,1,max)
otuWhichMax <- apply(otu.tab,1,which.max)
otuDominated <- which(otuMax > otuSum/2)

otuMaxTax <- taxonomy[otuWhichMax]
#otuDominated <- c(otuDominated[which(as.numeric(otuMaxTax[otuDominated])==32)],otuDominated[which(as.numeric(otuMaxTax[otuDominated])==33)])

taxonomyGroups <- as.character(groups)
taxonomyGroups[otuDominated] <- as.character(otuMaxTax[otuDominated])

taxonomyGroups <- as.factor(taxonomyGroups)

groups <- taxonomyGroups

# assign appropriate names to single taxa dominated groups
newLevels <- levels(taxonomyGroups)
splittaxa <- strsplit(levels(taxonomyGroups),split=";")

for (i in 1:length(splittaxa)) {
	if (length(splittaxa[[i]])>1) {
		newLevels[i] <- paste(splittaxa[[i]][length(splittaxa[[i]])-1],splittaxa[[i]][length(splittaxa[[i]])])
	}
	else {
		newLevels[i] <- splittaxa[[i]][1]
	}
}

levels(taxonomyGroups) <- newLevels

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

pdf("output/human_mouth_output/pcoa_plots.pdf")

#plot pcoa plots
plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.2,0.32,levels(taxonomyGroups),col=palette(),pch=19)
plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(exponent.pcoa$vectors[,1],exponent.pcoa$vectors[,2], col=groups,main="Exponent Normalized UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(exponent.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(exponent.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

#plot correlation between different UniFrac modes
plot(unweighted.vector,information.vector,main="unweighted vs. information UniFrac")
plot(weighted.vector,information.vector,main="weighted vs. information UniFrac")
plot(unweighted.vector,weighted.vector,main="unweighted vs. weighted UniFrac")
plot(exponent.vector,information.vector,main="normalized exponent vs. information UniFrac")
plot(exponent.vector,weighted.vector,main="normalized exponent vs. weighted UniFrac")

dev.off()
