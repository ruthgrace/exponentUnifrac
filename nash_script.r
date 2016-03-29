#!/usr/bin/env Rscript
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata

source("UniFrac.r")
library(ape)
library(phangorn)
library(vegan)

otu.tab <- read.table("data/nash_data/summed_data_baseline_only.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

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

palette(c("blue", "blue4", "red", "red4", "purple"))

pdf("nash_output/pcoa_plots.pdf",height=7,width=9)

originalpar <- par()
par(mar=c(5.1, 5.1, 5.1, 14.1),xpd=TRUE)

#plot pcoa plots
plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac",xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac",xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac",xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(ratio.pcoa$vectors[,1],ratio.pcoa$vectors[,2], col=groups,main="Centered Log Ratio UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(ratio.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(ratio.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(ratio_no_log.pcoa$vectors[,1],ratio_no_log.pcoa$vectors[,2], col=groups,main="Ratio UniFrac",xlab=paste("First Component", round(ratio_no_log.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(ratio_no_log.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
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
cor(h.ss.aldex$effect, y = h.nash.aldex$effect, use = "everything", method = "spearman")
# [1] 0.7692734

plot(h.ss.aldex$effect, h.metnash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs SS compared to healthy vs extreme NASH",xlab="Healthy vs. SS",ylab="Healthy vs. extreme NASH")
cor(h.ss.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman")
# [1] 0.7193963

plot(h.nash.aldex$effect, h.metnash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs NASH compared to healthy vs extreme NASH",xlab="Healthy vs. NASH",ylab="Healthy vs. extreme NASH")
cor(h.nash.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman")
# [1] 0.8704643

# make plots with top and bottom decile effect sizes from metnash comparison colored
orchid <- c(col2rgb("darkorchid3"))
red <- orchid[1]
green <- orchid[2]
blue <- orchid[3]
orchid <- rgb(red/255, green/255, blue/255, 0.3)

myblue <- c(col2rgb("deepskyblue3"))
red <- myblue[1]
green <- myblue[2]
blue <- myblue[3]
myblue <- rgb(red/255, green/255, blue/255, 0.3)

mygray <- c(col2rgb("gray18"))
red <- mygray[1]
green <- mygray[2]
blue <- mygray[3]
mygray <- rgb(red/255, green/255, blue/255, 0.3)

effectorder <- order(h.metnash.aldex$effect, decreasing=TRUE)
h.ss.aldex <- h.ss.aldex[effectorder,]
h.nash.aldex <- h.nash.aldex[effectorder,]
h.metnash.aldex <- h.metnash.aldex[effectorder,]

palette(c(myblue,mygray,orchid))

effectgroup <- rep("other",nrow(h.metnash.aldex))
decile <- round(length(effectgroup)/10)
effectgroup[1:decile] <- "top"
effectgroup[(length(effectgroup) - decile + 1):length(effectgroup)] <- "bottom"
effectgroup <- as.factor(effectgroup)

plot(h.ss.aldex$effect, h.nash.aldex$effect, pch=19,col=effectgroup, main="Effect sizes of healthy vs SS compared to healthy vs NASH",xlab="Healthy vs. SS",ylab="Healthy vs. NASH")
cor(h.ss.aldex$effect, y = h.nash.aldex$effect, use = "everything", method = "spearman")
# [1] 0.7692734

plot(h.ss.aldex$effect, h.metnash.aldex$effect, pch=19,col=effectgroup, main="Effect sizes of healthy vs SS compared to healthy vs extreme NASH",xlab="Healthy vs. SS",ylab="Healthy vs. extreme NASH")
cor(h.ss.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman")
# [1] 0.7193963

plot(h.nash.aldex$effect, h.metnash.aldex$effect, pch=19,col=effectgroup, main="Effect sizes of healthy vs NASH compared to healthy vs extreme NASH",xlab="Healthy vs. NASH",ylab="Healthy vs. extreme NASH")
cor(h.nash.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman")
# [1] 0.8704643

dev.off()

print("Average difference in distance from zero of healthy vs. nash compared to healthy vs. ss:")
summary(abs(h.nash.aldex$effect) - abs(h.ss.aldex$effect))

# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.489900 -0.120600 -0.005547 -0.014710  0.072460  0.440900 

print("Average difference in distance from zero of healthy vs. extreme nash compared to healthy vs. ss:")
summary(abs(h.metnash.aldex$effect) - abs(h.ss.aldex$effect))

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.61550 -0.09591  0.01094  0.01128  0.11820  0.60140 

print("Average difference in distance from zero of healthy vs. extreme nash compared to healthy vs. nash:")
summary(abs(h.metnash.aldex$effect) - abs(h.nash.aldex$effect))

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.34440 -0.04699  0.03105  0.02599  0.08771  0.47590 

rownames(h.metnash.aldex)[1:decile]
#  [1] "39"   "89"   "194"  "226"  "141"  "49"   "501"  "936"  "68"   "984" 
# [11] "1170" "215"  "181"  "154"  "178"  "51"   "47"   "139"  "1477" "93"  
# [21] "13"   "1243" "92"  

rownames(h.metnash.aldex)[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)]
#  [1] "129"  "888"  "150"  "88"   "1127" "149"  "25"   "99"   "1854" "11"  
# [11] "1"    "161"  "23"   "48"   "31"   "41"   "30"   "19"   "982"  "24"  
# [21] "8"    "61"   "17"  

h.metnash.aldex$effect[1:decile]
# [1] 0.8223052 0.7456387 0.5331967 0.5144386 0.5044868 0.4834193 0.4819016
# [8] 0.4606236 0.4480371 0.4441347 0.4365944 0.4358672 0.4221069 0.4188573
# [15] 0.3743713 0.3717831 0.3589893 0.3443864 0.3406586 0.3405844 0.3300911
# [22] 0.3297504 0.3282504
paste(round(h.metnash.aldex$effect[1:decile], digits=3),collapse="\n")

h.metnash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)]
# [1] -0.4853962 -0.4879183 -0.4880809 -0.5063250 -0.5111801 -0.5122526
# [7] -0.5131539 -0.5241678 -0.5247356 -0.5281960 -0.5341553 -0.5468414
# [13] -0.5547340 -0.5590358 -0.5604316 -0.5655144 -0.5744744 -0.5772855
# [19] -0.5949985 -0.6344935 -0.7214833 -0.7758869 -0.8734976
paste(round(h.metnash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n")

h.ss.aldex$effect[1:decile]
# [1]  0.22091446  0.60189782  0.30495031  0.41812395  0.22621953  0.21242478
# [7]  0.70310602  0.19935897  0.33122063  0.11923063  0.50000251  0.50989240
# [13]  0.10249190  0.31735013  0.55114706  0.08650403 -0.07052043  0.19028603
# [19]  0.24666039  0.75471668  0.26309971  0.34314307  0.18927969
paste(round(h.ss.aldex$effect[1:decile], digits=3),collapse="\n")

h.ss.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)]
# [1] -0.18260123 -0.30078646 -0.14862096 -0.49557550  0.04783647 -0.16860430
# [7] -0.28465360 -0.50404282  0.05277693 -0.15168854 -0.18848113 -0.34969275
# [13] -0.40642904 -0.64746624 -0.57357971 -0.36201428 -0.27432791 -0.55385864
# [19] -0.41021008 -0.33896505 -0.61918668 -0.50220860 -0.30620984
paste(round(h.ss.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n")

h.nash.aldex$effect[1:decile]
# [1]  0.34643719  0.60623335  0.37745840  0.26498522  0.51937513  0.16391584
# [7]  0.30746304  0.23857476  0.44600004  0.37740534  0.33005539  0.28094810
# [13]  0.18885896  0.36630199  0.21151861  0.29063807 -0.05445759  0.34947815
# [19]  0.31118779  0.54750210  0.59361189  0.25211330  0.26096296
paste(round(h.nash.aldex$effect[1:decile], digits=3),collapse="\n")

h.nash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)]
# [1] -0.1628749 -0.2789568 -0.1134037 -0.5477269 -0.4887603 -0.2206858
# [7] -0.2630476 -0.6045677 -0.3342644 -0.2935091 -0.3686036 -0.3816436
# [13] -0.5527248 -0.3336079 -0.6785748 -0.5144209 -0.6289683 -0.6420113
# [19] -0.5909743 -0.5954985 -0.5888082 -0.7700357 -0.7201200
paste(round(h.nash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n")

# BARPLOT

# ALDEX BY GENUS TO COMPARE WITH QPCR
