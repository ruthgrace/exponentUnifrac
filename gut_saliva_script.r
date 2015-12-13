#!/usr/bin/env Rscript
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata
# run in bash like: nohup Rscript gut_saliva_script.r > gut_saliva_script_nohup.out 2>&1&

source("UniFrac.r")
library(ape)
library(phangorn)
library(vegan)
#function to get variance explained for the PCOA component labels
getVarExplained <- function(vector) {
	rawVarEx <- apply(vector,2,function(x) sd(x)*sd(x))
	totalVarExplained <- sum(rawVarEx)
	varEx <- rawVarEx/totalVarExplained
	return(varEx)
}

plot_all_gut_saliva_unifrac <- function(count_file, tree_file, output_file) {

	otu.tab <- read.table(count_file, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
	otu.tab <- data.matrix(otu.tab)
	# read and root tree (rooted tree is required)
	tree <- read.tree(tree_file)
	tree <- midpoint(tree)

	# read metadata
	MyMeta <- colnames(otu.tab)
	MyMeta[grepl("Saliva_",colnames(otu.tab))] <- "Saliva"
	MyMeta[grepl("Stool_",colnames(otu.tab))] <- "Stool"
	
	#rarefy data for unweighted unifrac
	otu.tab.rarefy <- rrarefy(otu.tab, min(apply(otu.tab,1,sum)))
	otu.tab.rarefy <- data.matrix(otu.tab.rarefy)
	#calculate distance matrix
	unweighted <- getDistanceMatrix(otu.tab.rarefy,tree,method="unweighted",verbose=TRUE)
	#output distance matrices
	write.table(unweighted,file="unweighted_distance_matrix.txt",sep="\t",quote=FALSE)
	weighted <- getDistanceMatrix(otu.tab,tree,method="weighted",verbose=TRUE)
	write.table(weighted,file="weighted_distance_matrix.txt",sep="\t",quote=FALSE)
	information <- getDistanceMatrix(otu.tab,tree,method="information",verbose=TRUE)
	write.table(information,file="information_distance_matrix.txt",sep="\t",quote=FALSE)
	exponent <- getDistanceMatrix(otu.tab,tree,method="exponent",verbose=TRUE)
	write.table(exponent,file="exponent_normalize_distance_matrix.txt",sep="\t",quote=FALSE)

	#conditions (bv - bacterial vaginosis as scored by nugent/amsel, i - intermediate, n - normal/healthy)
	groups <- MyMetaOrdered$SSvsNASH
	originalgroups <- groups

	# healthy is represented by 0, SS/NASH is represented by 1
	groups <- groups + 1;
	groups[which(is.na(MyMetaOrdered$SSvsNASH))] <- 0

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

	unweighted.varEx <- getVarExplained(unweighted.pcoa$vectors)
	weighted.varEx <- getVarExplained(weighted.pcoa$vectors)
	information.varEx <- getVarExplained(information.pcoa$vectors)
	exponent.varEx <- getVarExplained(exponent.pcoa$vectors)

	#get vector version of distance matrices for correlation plots below
	unweighted.vector <- unlist(unweighted[lower.tri(unweighted,diag=TRUE)])
	weighted.vector <- unlist(weighted[lower.tri(weighted,diag=TRUE)])
	information.vector <- unlist(information[lower.tri(information,diag=TRUE)])
	exponent.vector <- unlist(exponent[lower.tri(exponent,diag=TRUE)])



	pdf("nash_output/pcoa_plots.pdf")

	#plot pcoa plots
	plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend(-0.37,-0.05,c("Healthy","SS","NASH"),col=palette(),pch=19)
	plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	plot(exponent.pcoa$vectors[,1],exponent.pcoa$vectors[,2], col=groups,main="Exponent Normalized UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(exponent.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(exponent.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)

	#plot correlation between different UniFrac modes
	plot(unweighted.vector,information.vector,main="unweighted vs. information UniFrac")
	plot(weighted.vector,information.vector,main="weighted vs. information UniFrac")
	plot(unweighted.vector,weighted.vector,main="unweighted vs. weighted UniFrac")
	plot(exponent.vector,information.vector,main="normalized exponent vs. information UniFrac")
	plot(exponent.vector,weighted.vector,main="normalized exponent vs. weighted UniFrac")

	dev.off()
}

plot_all_gut_saliva_unifrac("data/gut_saliva_data/low_sequencing_depth_hmp_data.txt", "data/gut_saliva_data/low_sequencing_depth_subtree.tre","gut_vs_saliva_sequencing_depth_less_than_3000.pdf")
plot_all_gut_saliva_unifrac("data/gut_saliva_data/med_sequencing_depth_hmp_data.txt", "data/gut_saliva_data/med_sequencing_depth_subtree.tre","gut_vs_saliva_sequencing_depth_3000_to_6000.pdf")
plot_all_gut_saliva_unifrac("data/gut_saliva_data/high_sequencing_depth_hmp_data.txt", "data/gut_saliva_data/high_sequencing_depth_subtree.tre","gut_vs_saliva_sequencing_depth_6000_plus.pdf")

