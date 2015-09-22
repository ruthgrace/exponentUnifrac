#!/usr/bin/env Rscript
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata

source("UniFrac.r")
library(ape)
library(phangorn)
library(vegan)


# read OTU table and format appropriately for input into UniFrac methods
breastmilk.otu.tab <- read.table("data/camilla_data/td_OTU_tag_mapped_lineage.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- breastmilk.otu.tab$taxonomy
breastmilk.otu.tab <- breastmilk.otu.tab[-length(colnames(breastmilk.otu.tab))]
breastmilk.otu.tab <- t(as.matrix(breastmilk.otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(breastmilk.otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
breastmilk.otu.tab <- breastmilk.otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
breastmilk.tree <- read.tree("data/camilla_data/fasttree_all_seed_OTUs.tre")
breastmilk.tree <- midpoint(breastmilk.tree)

# read metadata
MyMeta<- read.table("data/camilla_data/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove infected sample S38I
#MyMeta <- MyMeta[(which(rownames(MyMeta)!="S38I")),]

# filter OTU table and metadata so that only samples which appear in both are retained
otu_indicies <- match(rownames(MyMeta),rownames(breastmilk.otu.tab))
otu_indicies <- otu_indicies[!is.na(otu_indicies)]
breastmilk.otu.tab <- breastmilk.otu.tab[otu_indicies,]
MyMetaOrdered <- MyMeta[match(rownames(breastmilk.otu.tab),rownames(MyMeta)),]

#rarefy data for unweighted unifrac
breastmilk.otu.tab.rarefy <- rrarefy(breastmilk.otu.tab, min(apply(otu.tab,1,sum)))

#calculate distance matrix
unweighted <- getDistanceMatrix(breastmilk.otu.tab.rarefy,breastmilk.tree,method="unweighted",verbose=TRUE)
weighted <- getDistanceMatrix(breastmilk.otu.tab,breastmilk.tree,method="weighted",verbose=TRUE)
information <- getDistanceMatrix(breastmilk.otu.tab,breastmilk.tree,method="information",verbose=TRUE)
exponent.normalize <- getDistanceMatrix(breastmilk.otu.tab,breastmilk.tree,method="exponent",verbose=TRUE,normalize=TRUE)

#output distance matrices
write.table(unweighted,file="experimental_output/unweighted_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(weighted,file="experimental_output/weighted_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(information,file="experimental_output/information_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(exponent.normalize,file="experimental_output/exponent_normalize_distance_matrix.txt",sep="\t",quote=FALSE)

#conditions (bv - bacterial vaginosis as scored by nugent/amsel, i - intermediate, n - normal/healthy)
groups <- MyMetaOrdered$Gestation
originalgroups <- groups
groups <- c("No OTU > 50%")
# change conditions so that samples which are more than 50% one taxa are colored by that taxa
otuSum <- apply(breastmilk.otu.tab,1,sum)
otuMax <- apply(breastmilk.otu.tab,1,max)
otuWhichMax <- apply(breastmilk.otu.tab,1,which.max)
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
exponent.normalize.pcoa <- pcoa(exponent.normalize)


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
exponent.normalize.varEx <- getVarExplained(exponent.normalize.pcoa$vectors)

#get vector version of distance matrices for correlation plots below
unweighted.vector <- unlist(unweighted[lower.tri(unweighted,diag=TRUE)])
weighted.vector <- unlist(weighted[lower.tri(weighted,diag=TRUE)])
information.vector <- unlist(information[lower.tri(information,diag=TRUE)])
exponent.normalize.vector <- unlist(exponent.normalize[lower.tri(exponent.normalize,diag=TRUE)])



pdf("experimental_output/pcoa_plots.pdf")

#plot pcoa plots
plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend(0.1,0.17,levels(taxonomyGroups),col=palette(),pch=19)
plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
plot(exponent.normalize.pcoa$vectors[,1],exponent.normalize.pcoa$vectors[,2], col=groups,main="Exponent Normalized UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(exponent.normalize.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(exponent.normalize.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

#plot correlation between different UniFrac modes
plot(unweighted.vector,information.vector,main="unweighted vs. information UniFrac")
plot(weighted.vector,information.vector,main="weighted vs. information UniFrac")
plot(unweighted.vector,weighted.vector,main="unweighted vs. weighted UniFrac")
plot(unweighted.vector,exponent.normalize.vector,main="unweighted vs. exponent UniFrac")
plot(weighted.vector,exponent.normalize.vector,main="weighted vs. exponent UniFrac")
plot(information.vector,exponent.normalize.vector,main="information vs. exponent UniFrac")

dev.off()
