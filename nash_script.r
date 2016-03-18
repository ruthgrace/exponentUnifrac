#!/usr/bin/env Rscript
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata

source("UniFrac.r")
library(ape)
library(phangorn)
library(vegan)

otu.tab <- read.table("data/nash_data/summed_data_gg.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(otu.tab,2,sum)))
otu.tab <- otu.tab[,taxaOrder]

# read and root tree (rooted tree is required)
tree <- read.tree("data/nash_data/fasttree_all_seed_OTUs.tre")
tree <- midpoint(tree)

# read metadata
MyMeta<- read.table("data/nash_data/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
metadata <- MyMeta[grepl("a$",rownames(MyMeta)),]
rownames(metadata) <- gsub("a$","",rownames(metadata))

metagenomic_samples <- c("CL-119", "CL-139-6mo-2", "CL-141-BL-R2", "CL-144-2", "CL-160", "CL-165", "CL-166-BL", "CL-169-BL", "CL-173-2", "CL-177", "HLD-100", "HLD-102", "HLD-111-2", "HLD-112", "HLD-23", "HLD-28", "HLD-47", "HLD-72-2", "HLD-80", "HLD-85")

## sanity check to make sure all your counts have metadata
# which(!(colnames(otu.tab) %in% rownames(metadata)))

otu.tab <- t(otu.tab)

#rarefy data for unweighted unifrac
otu.tab.rarefy <- rrarefy(otu.tab, min(apply(otu.tab,1,sum)))

#calculate distance matrix

# do unweighted separately with rarefied data set
unweighted <- getDistanceMatrix(otu.tab.rarefy,tree,method="unweighted",verbose=TRUE)

all_distance_matrices <- getDistanceMatrix(otu.tab,tree,method="all",verbose=TRUE)

weighted <- all_distance_matrices[["weighted"]]
information <- all_distance_matrices[["information"]]
ratio <- all_distance_matrices[["ratio"]]
ratio_no_log <- all_distance_matrices[["ratio_no_log"]]

#output distance matrices
write.table(unweighted,file="nash_output/unweighted_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(weighted,file="nash_output/weighted_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(information,file="nash_output/information_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(ratio,file="nash_output/ratio_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(ratio_no_log,file="nash_output/ratio_no_log_distance_matrix.txt",sep="\t",quote=FALSE)

# unweighted <- read.table("nash_output/unweighted_distance_matrix.txt", sep = "\t", quote = "", row.names = 1, check.names = FALSE)
# weighted <- read.table("nash_output/weighted_distance_matrix.txt", sep = "\t", quote = "", row.names = 1, check.names = FALSE)
# information <- read.table("nash_output/information_distance_matrix.txt", sep = "\t", quote = "", row.names = 1, check.names = FALSE)
# ratio <- read.table("nash_output/ratio_distance_matrix.txt", sep = "\t", quote = "", row.names = 1, check.names = FALSE)
# ratio_no_log <- read.table("nash_output/ratio_no_log_distance_matrix.txt", sep = "\t", quote = "", row.names = 1, check.names = FALSE)

# conditions: Originally 0 meant steatohepatosis, and 1 meant NASH
groups <- metadata$SSvsNASH[match(rownames(otu.tab),rownames(metadata))]
originalgroups <- groups

# Make healthy represented by 0, SS by 1, NASH by 2
groups <- groups + 1;
groups[which(is.na(groups))] <- 0

# make healthy 1, ss 2, nash 3 (healthy metagenomic will be 0 and nash metagenomic will be 4)
groups <- groups + 1

# mark healthy samples selected for metagenomic study
groups[which(rownames(otu.tab) %in% metagenomic_samples & groups == 1)] <- 0

# mark nash samples selected for metagenomic study
groups[which(rownames(otu.tab) %in% metagenomic_samples & groups == 3)] <- 4

# # remove SS (intermediate group)
# otu.tab.original <- otu.tab
# otu.tab <- otu.tab[which(groups != 1),]
# groups <- groups[which(groups != 1)]

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

groups[which(groups == 0)] <- "Healthy Metagenomic"
groups[which(groups == 1)] <- "Healthy"
groups[which(groups == 2)] <- "SS"
groups[which(groups == 3)] <- "NASH"
groups[which(groups == 4)] <- "NASH Metagenomic"

groups <- as.factor(groups)

originalpalette <- palette()

palette(c("blue4", "blue", "purple", "red", "red4"))

pdf("nash_output/pcoa_plots.pdf",height=7,width=9)

originalpar <- par()
par(mar=c(5.1, 5.1, 5.1, 14.1),xpd=TRUE)

#plot pcoa plots
plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(ratio.pcoa$vectors[,1],ratio.pcoa$vectors[,2], col=groups,main="Centered Log Ratio UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(ratio.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(ratio.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(ratio_no_log.pcoa$vectors[,1],ratio_no_log.pcoa$vectors[,2], col=groups,main="Centered Ratio UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(ratio_no_log.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(ratio_no_log.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))

#plot correlation between different UniFrac modes
plot(unweighted.vector,information.vector,main="unweighted vs. information UniFrac")
plot(weighted.vector,information.vector,main="weighted vs. information UniFrac")
plot(unweighted.vector,weighted.vector,main="unweighted vs. weighted UniFrac")
plot(ratio.vector,information.vector,main="normalized ratio vs. information UniFrac")
plot(ratio.vector,weighted.vector,main="normalized ratio vs. weighted UniFrac")
plot(ratio_no_log.vector,information.vector,main="normalized ratio (no log) vs. information UniFrac")
plot(ratio_no_log.vector,weighted.vector,main="normalized ratio (no log) vs. weighted UniFrac")

dev.off()

palette(originalpalette)
par(originalpar)

# ALDEx effect size plots

library(ALDEx2)

pdf("nash_output/ALDEx_effect_size_plots.pdf")

h.ss <- t(otu.tab)
h.ss.cond <- groups
h.ss <- h.ss[,which(h.ss.cond == "Healthy Metagenomic" | h.ss.cond == "SS")]
h.ss.cond <- h.ss.cond[which(h.ss.cond == "Healthy Metagenomic" | h.ss.cond == "SS")]

h.ss.aldex <- aldex(data.frame(h.ss),as.character(h.ss.cond))
aldex.plot(h.ss.aldex,type="MA",main="Bland-Altman style plot for healthy vs. SS",ylab="ratio",xlab="average")
aldex.plot(h.ss.aldex,type="MW",main="Difference within vs. difference between for healthy vs. SS",xlab="Difference within",ylab="Difference between")

h.nash <- t(otu.tab)
h.nash.cond <- groups
h.nash.cond[which(h.nash.cond == "NASH Metagenomic")] <- "NASH"
h.nash <- h.nash[,which(h.nash.cond == "Healthy Metagenomic" | h.nash.cond == "NASH")]
h.nash.cond <- h.nash.cond[which(h.nash.cond == "Healthy Metagenomic" | h.nash.cond == "NASH")]

h.nash.aldex <- aldex(data.frame(h.nash),as.character(h.nash.cond))
aldex.plot(h.nash.aldex,type="MA",main="Bland-Altman style plot for healthy vs. NASH",ylab="ratio",xlab="average")
aldex.plot(h.nash.aldex,type="MW",main="Difference within vs difference between for healthy vs. NASH",xlab="Difference within",ylab="Difference between")

h.metnash <- t(otu.tab)
h.metnash.cond <- groups
h.metnash <- h.metnash[,which(h.metnash.cond == "Healthy Metagenomic" | h.metnash.cond == "NASH Metagenomic")]
h.metnash.cond <- h.metnash.cond[which(h.metnash.cond == "Healthy Metagenomic" | h.metnash.cond == "NASH Metagenomic")]

h.metnash.aldex <- aldex(data.frame(h.metnash),as.character(h.metnash.cond))
aldex.plot(h.metnash.aldex,type="MA",main="Bland-Altman style plot for healthy vs. extreme NASH",ylab="ratio",xlab="average")
aldex.plot(h.metnash.aldex,type="MW",main="Difference within vs difference between for healthy vs. extreme NASH",xlab="Difference within",ylab="Difference between")

mycolor <- c(col2rgb("turquoise4"))
red <- mycolor[1]
green <- mycolor[2]
blue <- mycolor[3]
mycolor <- rgb(red/255, green/255, blue/255, 0.3)

plot(h.ss.aldex$effect, h.nash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs SS compared to healthy vs NASH",xlab="Healthy vs. SS",ylab="Healthy vs. NASH")
plot(h.ss.aldex$effect, h.metnash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs SS compared to healthy vs extreme NASH",xlab="Healthy vs. SS",ylab="Healthy vs. extreme NASH")
plot(h.nash.aldex$effect, h.metnash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs NASH compared to healthy vs extreme NASH",xlab="Healthy vs. NASH",ylab="Healthy vs. extreme NASH")

dev.off()

print("Average difference in distance from zero of healthy vs. nash compared to healthy vs. ss:")
summary(abs(h.nash.aldex$effect) - abs(h.ss.aldex$effect))
print("Average difference in distance from zero of healthy vs. extreme nash compared to healthy vs. ss:")
summary(abs(h.metnash.aldex$effect) - abs(h.ss.aldex$effect))
print("Average difference in distance from zero of healthy vs. extreme nash compared to healthy vs. nash:")
summary(abs(h.metnash.aldex$effect) - abs(h.nash.aldex$effect))
