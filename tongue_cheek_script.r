#!/usr/bin/env Rscript 
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata
# run: nohup Rscript ./tongue_cheek_script.r > tongue_cheek_script_nohup.out 2>&1&

source("UniFrac.r")
library(ape)
library(phangorn)
library(vegan)

otu.tab <- read.table("data/tongue_cheek_data/hmp_tongue_cheek_data.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- rownames(otu.tab)
#make row names samples, and col names OTUs
otu.tab <- t(as.matrix(otu.tab))

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(otu.tab,2,sum)))
taxonomy <- taxonomy[taxaOrder]
otu.tab <- otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
tree <- read.tree("data/tongue_cheek_data/hmp_tongue_cheek_subtree.tre")
tree <- midpoint(tree)

# read metadata
MyMetaOrdered <- rownames(otu.tab)
MyMetaOrdered <- gsub("_.*$","",MyMetaOrdered)

#rarefy data for unweighted unifrac
otu.tab.rarefy <- rrarefy(otu.tab, min(apply(otu.tab,1,sum)))

# do Unweighted UniFrac separately with rarefied data
unweighted <- getDistanceMatrix(otu.tab.rarefy,tree,method="unweighted",verbose=TRUE)
write.table(unweighted,file="tongue_cheek_output/unweighted_distance_matrix.txt",sep="\t",quote=FALSE)

#calculate distance matrix
all_distance_matrices <- getDistanceMatrix(otu.tab,tree,method="all",verbose=TRUE)

#output distance matrices
weighted <- all_distance_matrices[["weighted"]]
write.table(weighted,file="tongue_cheek_output/weighted_distance_matrix.txt",sep="\t",quote=FALSE)
information <- all_distance_matrices[["information"]]
write.table(information,file="tongue_cheek_output/information_distance_matrix.txt",sep="\t",quote=FALSE)
ratio <- all_distance_matrices[["ratio"]]
write.table(ratio,file="tongue_cheek_output/ratio_normalize_distance_matrix.txt",sep="\t",quote=FALSE)
ratio_no_log <- all_distance_matrices[["ratio_no_log"]]
write.table(ratio_no_log,file="tongue_cheek_output/ratio_no_log_normalize_distance_matrix.txt",sep="\t",quote=FALSE)

# unweighted <- read.table("tongue_cheek_output/unweighted_distance_matrix.txt", sep = "\t", header = TRUE, row.names=1, quote = "")
# weighted <- read.table("tongue_cheek_output/weighted_distance_matrix.txt", sep = "\t", header = TRUE, row.names=1, quote = "")
# information <- read.table("tongue_cheek_output/information_distance_matrix.txt", sep = "\t", header = TRUE, row.names=1, quote = "")
# ratio <- read.table("tongue_cheek_output/ratio_normalize_distance_matrix.txt", sep = "\t", header = TRUE, row.names=1, quote = "")
# ratio_no_log <- read.table("tongue_cheek_output/ratio_no_log_normalize_distance_matrix.txt", sep = "\t", header = TRUE, row.names=1, quote = "")


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
ratio.pcoa <- pcoa(ratio)
ratio_no_log.pcoa <- pcoa(ratio_no_log)


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
ratio.varEx <- getVarExplained(ratio.pcoa$vectors)
ratio_no_log.varEx <- getVarExplained(ratio_no_log.pcoa$vectors)

#get vector version of distance matrices for correlation plots below
unweighted.vector <- unlist(unweighted[lower.tri(unweighted,diag=TRUE)])
weighted.vector <- unlist(weighted[lower.tri(weighted,diag=TRUE)])
information.vector <- unlist(information[lower.tri(information,diag=TRUE)])
ratio.vector <- unlist(ratio[lower.tri(ratio,diag=TRUE)])
ratio_no_log.vector <- unlist(ratio_no_log[lower.tri(ratio_no_log,diag=TRUE)])

# replace abbreviations with full body site names (there aren't actually any dominant taxa in this data set)
taxonomyGroups <- as.factor(c("Cheek", "Tongue"))

palette(c("black", "blue", "red"))
dev.off()

pdf("tongue_cheek_output/pcoa_plots.pdf")
par(oma=c(1,1,1,5))
#plot pcoa plots
plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Coordinate", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Coordinate", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("right", levels(taxonomyGroups), pch=rep(19,length(taxonomyGroups)), col=palette()[1:length(taxonomyGroups)], xpd=NA, inset=c(-0.25,0))
plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Coordinate", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Coordinate", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("right", levels(taxonomyGroups), pch=rep(19,length(taxonomyGroups)), col=palette()[1:length(taxonomyGroups)], xpd=NA, inset=c(-0.25,0))
plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Coordinate", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Coordinate", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("right", levels(taxonomyGroups), pch=rep(19,length(taxonomyGroups)), col=palette()[1:length(taxonomyGroups)], xpd=NA, inset=c(-0.25,0))
plot(ratio.pcoa$vectors[,1],ratio.pcoa$vectors[,2], col=groups,main="Centered Log Ratio UniFrac\nprincipal coordinates analysis",xlab=paste("First Coordinate", round(ratio.varEx[1],digits=3),"variance explained"),ylab=paste("Second Coordinate", round(ratio.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("right", levels(taxonomyGroups), pch=rep(19,length(taxonomyGroups)), col=palette()[1:length(taxonomyGroups)], xpd=NA, inset=c(-0.25,0))
plot(ratio_no_log.pcoa$vectors[,1],ratio_no_log.pcoa$vectors[,2], col=groups,main="Centered Ratio UniFrac\nprincipal coordinates analysis",xlab=paste("First Coordinate", round(ratio_no_log.varEx[1],digits=3),"variance explained"),ylab=paste("Second Coordinate", round(ratio_no_log.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("right", levels(taxonomyGroups), pch=rep(19,length(taxonomyGroups)), col=palette()[1:length(taxonomyGroups)], xpd=NA, inset=c(-0.25,0))

#plot correlation between different UniFrac modes
plot(unweighted.vector,information.vector,main="unweighted vs. information UniFrac")
plot(weighted.vector,information.vector,main="weighted vs. information UniFrac")
plot(unweighted.vector,weighted.vector,main="unweighted vs. weighted UniFrac")
plot(ratio.vector,information.vector,main="normalized ratio vs. information UniFrac")
plot(ratio.vector,weighted.vector,main="normalized ratio vs. weighted UniFrac")
plot(ratio_no_log.vector,weighted.vector,main="normalized no log ratio vs. weighted UniFrac")

dev.off()

### BIPLOT
library(ALDEx2)
library(compositions)

otu.tab <- read.table("data/tongue_cheek_data/hmp_tongue_cheek_data.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

groups <- colnames(otu.tab)
groups <- gsub("_.*$","",groups)

x <- aldex.clr(otu.tab)
x.e <- aldex.effect(x, groups)
x.t <- aldex.ttest(x, groups)
x.all <- data.frame(x.t, x.e)

# duplicate the data for the biplot and clustering
otu.tab.prior <- otu.tab
otu.tab.prior[otu.tab.prior == 0] <- 0.5

bi <- acomp(t(otu.tab.prior))
clr <- t(apply(otu.tab.prior, 2, function(x){log(x) - mean(log(x))}))
pcx <- prcomp(clr)
sum(pcx$sdev[1]^2)/mvar(clr)
# [1] 0.1926549
sum(pcx$sdev[2]^2)/mvar(clr)
# [1] 0.06741566
sum(pcx$sdev[3]^2)/mvar(clr)
# [1] 0.0477463

save(pcx,file="tongue_cheek_output/biplot_data.rdat")
## load with
# load(file="tongue_cheek_output/biplot_data.rdat")

pdf("tongue_cheek_output/biplot.pdf")
biplot(pcx, cex=0.6, col=c("black", "red"), scale=0, arrow.len=0, var.axes=F)
dev.off()

