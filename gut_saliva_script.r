#!/usr/bin/env Rscript
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata
# run in bash like: nohup Rscript gut_saliva_script.r > gut_saliva_script_nohup.out 2>&1&


#function to get variance explained for the PCOA component labels
getVarExplained <- function(vector) {
	rawVarEx <- apply(vector,2,function(x) sd(x)*sd(x))
	totalVarExplained <- sum(rawVarEx)
	varEx <- rawVarEx/totalVarExplained
	return(varEx)
}

plot_all_gut_saliva_unifrac <- function(count_file, tree_file, output_folder) {
	source("UniFrac.r")
	library(ape)
	library(phangorn)
	library(vegan)
	otu.tab <- read.table(count_file, header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
	otu.tab <- data.matrix(otu.tab)
	#put samples in rows and OTUs in columns
	otu.tab <- t(otu.tab)
	# read and root tree (rooted tree is required)
	tree <- read.tree(tree_file)
	tree <- midpoint(tree)

	# read metadata
	MyMeta <- rownames(otu.tab)
	MyMeta[grepl("Saliva_",rownames(otu.tab))] <- "Saliva"
	MyMeta[grepl("Stool_",rownames(otu.tab))] <- "Stool"
	
	#rarefy data for unweighted unifrac
	otu.tab.rarefy <- rrarefy(otu.tab, min(apply(otu.tab,1,sum)))
	otu.tab.rarefy <- data.matrix(otu.tab.rarefy)
	
	#calculate distance matrix
	
	## for debug / testing
	# otuTable <- otu.tab
	# tree.original <- tree
	# method = "ratio"
	# verbose = TRUE
	# pruneTree=FALSE
	# normalize=TRUE
	
	all_distance_matrices <- getDistanceMatrix(otu.tab,tree,method="all",verbose=TRUE)
	
	unweighted <- all_distance_matrices[["unweighted"]]
	#output distance matrices
	write.table(unweighted,file=paste(output_folder,"unweighted_distance_matrix.txt",sep=""),sep="\t",quote=FALSE)
	weighted <- all_distance_matrices[["weighted"]]
	write.table(weighted,file=paste(output_folder,"weighted_distance_matrix.txt",sep=""),sep="\t",quote=FALSE)
	information <- all_distance_matrices[["information"]]
	write.table(information,file=paste(output_folder,"information_distance_matrix.txt",sep=""),sep="\t",quote=FALSE)
	ratio <- all_distance_matrices[["ratio"]]
	write.table(ratio,file=paste(output_folder,"ratio_normalize_distance_matrix.txt",sep=""),sep="\t",quote=FALSE)
	ratio_no_log <- all_distance_matrices[["ratio_no_log"]]
	write.table(ratio_no_log,file=paste(output_folder,"ratio_no_log_normalize_distance_matrix.txt",sep=""),sep="\t",quote=FALSE)

	# unweighted <- read.table("gut_saliva_output/low_sequencing_depth/unweighted_distance_matrix.txt",sep="\t",quote="",check.names=FALSE)
	# weighted <- read.table("gut_saliva_output/low_sequencing_depth/weighted_distance_matrix.txt",sep="\t",quote="",check.names=FALSE)
	# information <- read.table("gut_saliva_output/low_sequencing_depth/information_distance_matrix.txt",sep="\t",quote="",check.names=FALSE)
	# ratio <- read.table("gut_saliva_output/low_sequencing_depth/ratio_normalize_distance_matrix.txt",sep="\t",quote="",check.names=FALSE)
	# ratio_no_log <- read.table("gut_saliva_output/low_sequencing_depth/ratio_no_log_normalize_distance_matrix.txt",sep="\t",quote="",check.names=FALSE)
	
	#conditions (bv - bacterial vaginosis as scored by nugent/amsel, i - intermediate, n - normal/healthy)
	groups <- MyMeta
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
	ratio.pcoa <- pcoa(ratio)
	ratio_no_log.pcoa <- pcoa(ratio_no_log)

	unweighted.varEx <- getVarExplained(unweighted.pcoa$vectors)
	weighted.varEx <- getVarExplained(weighted.pcoa$vectors)
	information.varEx <- getVarExplained(information.pcoa$vectors)
	ratio.varEx <- getVarExplained(ratio.pcoa$vectors)
	ratio_no_log.varEx <- getVarExplained(ratio_no_log.pcoa$vectors)

	#get vector version of distance matrices for correlation plots below
	unweighted.vector <- unlist(unweighted[lower.tri(unweighted,diag=TRUE)])
	weighted.vector <- unlist(weighted[lower.tri(weighted,diag=TRUE)])
	information.vector <- unlist(information[lower.tri(information,diag=TRUE)])
	ratio.vector <- unlist(ratio[lower.tri(ratio,diag=TRUE)])
	ratio_no_log.vector <- unlist(ratio_no_log[lower.tri(ratio_no_log,diag=TRUE)])

	pdf(paste(output_folder,"pcoa_plots.pdf",sep=""))
	# par(mar=c(3,3,3,2.1), oma=c(1,0,0,5))
	par(oma=c(1,1,1,5))
	#plot pcoa plots
	plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("right", levels(groups), pch=c(19,19), col=palette()[c(1:2)], xpd=NA, inset=c(-0.25,0))
	plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("right", levels(groups), pch=c(19,19), col=palette()[c(1:2)], xpd=NA, inset=c(-0.25,0))
	plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("right", levels(groups), pch=c(19,19), col=palette()[c(1:2)], xpd=NA, inset=c(-0.25,0))
	plot(ratio.pcoa$vectors[,1],ratio.pcoa$vectors[,2], col=groups,main="ratio Normalized UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(ratio.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(ratio.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("right", levels(groups), pch=c(19,19), col=palette()[c(1:2)], xpd=NA, inset=c(-0.25,0))
	plot(ratio_no_log.pcoa$vectors[,1],ratio_no_log.pcoa$vectors[,2], col=groups,main="ratio no log Normalized UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(ratio_no_log.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(ratio_no_log.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
	legend("right", levels(groups), pch=c(19,19), col=palette()[c(1:2)], xpd=NA, inset=c(-0.25,0))

	#plot correlation between different UniFrac modes
	plot(unweighted.vector,information.vector,main="unweighted vs. information UniFrac")
	plot(weighted.vector,information.vector,main="weighted vs. information UniFrac")
	plot(unweighted.vector,weighted.vector,main="unweighted vs. weighted UniFrac")
	plot(ratio.vector,information.vector,main="normalized ratio vs. information UniFrac")
	plot(ratio.vector,weighted.vector,main="normalized ratio vs. weighted UniFrac")
	plot(ratio_no_log.vector,weighted.vector,main="normalized no log ratio vs. weighted UniFrac")

	dev.off()
}

## for testing/debugging purposes
# count_file <- "data/gut_saliva_data/low_sequencing_depth_hmp_data.txt"
# tree_file <- "data/gut_saliva_data/low_sequencing_depth_subtree.tre"
# output_folder <- "gut_saliva_output/low_sequencing_depth/less_than_3000_readcounts_"

plot_all_gut_saliva_unifrac("data/gut_saliva_data/low_sequencing_depth_hmp_data.txt", "data/gut_saliva_data/low_sequencing_depth_subtree.tre","gut_saliva_output/low_sequencing_depth/less_than_3000_readcounts_")
plot_all_gut_saliva_unifrac("data/gut_saliva_data/med_sequencing_depth_hmp_data.txt", "data/gut_saliva_data/med_sequencing_depth_subtree.tre","gut_saliva_output/low_sequencing_depth/3000_to_6000_readcounts_")
plot_all_gut_saliva_unifrac("data/gut_saliva_data/high_sequencing_depth_hmp_data.txt", "data/gut_saliva_data/high_sequencing_depth_subtree.tre","gut_saliva_output/low_sequencing_depth/more_than_6000_readcounts_")

