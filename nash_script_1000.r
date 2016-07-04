#!/usr/bin/env Rscript
options(error=recover)

#this script prints out PDF pcoa plots and distance matrices, given an OTU table, phylogenetic tree, and metadata

source("UniFrac.r")
library(ape)
library(phangorn)
library(vegan)
library(stringr)

otu.tab <- read.table("data/nash_data/td_OTU_tag_mapped_lineage_ab.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
# take out taxonomy
otu.tab <- otu.tab[,c(1:(ncol(otu.tab)-1))]

taxonomy.table <- read.table("data/nash_data/td_OTU_tag_mapped_lineage_ab.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)

taxonomy <- taxonomy.table$taxonomy
names(taxonomy) <- rownames(taxonomy.table)

#sort taxa from most to least abundant
taxaOrder <- rev(order(apply(otu.tab,1,sum)))
otu.tab <- otu.tab[taxaOrder,]
colnames(otu.tab) <- gsub("\\.","-",colnames(otu.tab))

# read and root tree (rooted tree is required)
tree <- read.tree("data/nash_data/OTU_seeds.tre")
tree <- midpoint(tree)

# fix tree labels to match OTUs
newlabels <- tree$tip.label
for (i in c(1:length(newlabels))) {
	newlabels[i] <- strsplit(newlabels[i],"\\|")[[1]][6]
}
tree$tip.label <- newlabels

# read metadata
MyMeta<- read.table("data/nash_data/metadata.txt", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
metadata <- MyMeta[grepl("a$",rownames(MyMeta)),]
rownames(metadata) <- gsub("a$","",rownames(metadata))
metadata <- metadata[!(grepl("m",rownames(metadata))),]
rownames(metadata) <- str_extract(rownames(metadata), "^[A-Z]+-[0-9]+")

# metagenomic_samples <- c("CL-119", "CL-139-6mo-2", "CL-141-BL-R2", "CL-144-2", "CL-160", "CL-165", "CL-166-BL", "CL-169-BL", "CL-173-2", "CL-177", "HLD-100", "HLD-102", "HLD-111-2", "HLD-112", "HLD-23", "HLD-28", "HLD-47", "HLD-72-2", "HLD-80", "HLD-85")
metagenomic_samples <- c("CL-119", "CL-139", "CL-141", "CL-144", "CL-160", "CL-165", "CL-166", "CL-169", "CL-173", "CL-177", "HLD-100", "HLD-102", "HLD-111", "HLD-112", "HLD-23", "HLD-28", "HLD-47", "HLD-72", "HLD-80", "HLD-85")

## sanity check to make sure all your counts have metadata
# which(!(colnames(otu.tab) %in% rownames(metadata)))

# remove OTUs that don't have at least 1% abundance in at least 1 sample
sample.sum <- apply(otu.tab,2,sum)
otu.filter <- apply(otu.tab,1,function(x) { return(length(which(x >= (sample.sum*0.01)))) } )
otu.tab <- otu.tab[which(otu.filter > 0),]

# take OTUs that were removed out of tree
absent <- which(!(tree$tip.label %in% rownames(otu.tab)))
if (length(absent > 0)) {
	tree <- drop.tip(tree, absent)
}

taxonomy <- taxonomy[which(names(taxonomy) %in% rownames(otu.tab))]

otu.tab <- t(otu.tab)

#rarefy data for unweighted unifrac
otu.tab.rarefy <- rrarefy(otu.tab, min(apply(otu.tab,1,sum)))

#calculate distance matrix

# do unweighted separately with rarefied data set
unweighted <- getDistanceMatrix(otu.tab.rarefy,tree,method="unweighted",verbose=TRUE)

all_distance_matrices <- getDistanceMatrix(otu.tab,tree,method="all",verbose=TRUE)

weighted <- all_distance_matrices[["weighted"]]
information <- all_distance_matrices[["information"]]
ratio <- all_distance_matrices[["ratio_no_log"]]

#output distance matrices
write.table(unweighted,file="nash_output/1000/unweighted_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(weighted,file="nash_output/1000/weighted_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(information,file="nash_output/1000/information_distance_matrix.txt",sep="\t",quote=FALSE)
write.table(ratio,file="nash_output/1000/ratio_distance_matrix.txt",sep="\t",quote=FALSE)

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

unweighted.pcoa <- pcoa(unweighted)
weighted.pcoa <- pcoa(weighted)
information.pcoa <- pcoa(information)
ratio.pcoa <- pcoa(ratio)

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

#get vector version of distance matrices for correlation plots below
unweighted.vector <- unlist(unweighted[lower.tri(unweighted,diag=TRUE)])
weighted.vector <- unlist(weighted[lower.tri(weighted,diag=TRUE)])
information.vector <- unlist(information[lower.tri(information,diag=TRUE)])
ratio.vector <- unlist(ratio[lower.tri(ratio,diag=TRUE)])

groups[which(groups == 0)] <- "Healthy Metagenomic"
groups[which(groups == 1)] <- "Healthy"
groups[which(groups == 2)] <- "SS"
groups[which(groups == 3)] <- "NASH"
groups[which(groups == 4)] <- "NASH Metagenomic"

groups <- as.factor(groups)

originalpalette <- palette()

palette(c("blue", "blue4", "red", "red4", "purple"))

pdf("nash_output/1000/pcoa_plots.pdf",height=7,width=9)

originalpar <- par()
par(mar=c(5.1, 5.1, 5.1, 14.1),xpd=TRUE)

#plot pcoa plots
plot(unweighted.pcoa$vectors[,1],unweighted.pcoa$vectors[,2], col=groups,main="Unweighted UniFrac",xlab=paste("First Component", round(unweighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(unweighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(weighted.pcoa$vectors[,1],weighted.pcoa$vectors[,2], col=groups,main="Weighted UniFrac",xlab=paste("First Component", round(weighted.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(weighted.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(information.pcoa$vectors[,1],information.pcoa$vectors[,2], col=groups,main="Information UniFrac",xlab=paste("First Component", round(information.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(information.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))
plot(ratio.pcoa$vectors[,1],ratio.pcoa$vectors[,2], col=groups,main="Ratio UniFrac\nprincipal coordinates analysis",xlab=paste("First Component", round(ratio.varEx[1],digits=3),"variance explained"),ylab=paste("Second Component", round(ratio.varEx[2],digits=3),"variance explained"),pch=19,cex.lab=1.4,cex.main=2)
legend("topright", levels(groups), pch=c(19,19), col=palette(), xpd=NA, inset=c(-0.45,0))

# #plot correlation between different UniFrac modes
# plot(unweighted.vector,information.vector,main="unweighted vs. information UniFrac")
# plot(weighted.vector,information.vector,main="weighted vs. information UniFrac")
# plot(unweighted.vector,weighted.vector,main="unweighted vs. weighted UniFrac")
# plot(ratio.vector,information.vector,main="normalized ratio vs. information UniFrac")
# plot(ratio.vector,weighted.vector,main="normalized ratio vs. weighted UniFrac")

dev.off()

palette(originalpalette)
par(originalpar)

# ALDEx effect size plots

library(ALDEx2)

pdf("nash_output/1000/ALDEx_effect_size_plots.pdf")

# randomly pick which healthy sample will be controls for the SS or the NASH comparison
# healthy.index <- which(groups == "Healthy")
# control.size <- floor(length(healthy.index)/2)
# healthy.ss <- sample(healthy.index,control.size)
# healthy.nash <- healthy.index[which(!(healthy.index %in% healthy.ss))]
# 
# groups <- as.character(groups)
# 
# groups[healthy.ss] <- "Healthy SS"
# groups[healthy.nash] <- "Healthy NASH"

h.ss <- t(otu.tab)
h.ss.cond <- as.character(groups)
# h.ss <- h.ss[,which(h.ss.cond == "Healthy SS" | h.ss.cond == "SS")]
# h.ss.cond <- h.ss.cond[which(h.ss.cond == "Healthy SS" | h.ss.cond == "SS")]
h.ss <- h.ss[,which(h.ss.cond == "Healthy" | h.ss.cond == "Healthy Metagenomic" | h.ss.cond == "SS")]
h.ss.cond <- h.ss.cond[which(h.ss.cond == "Healthy" | h.ss.cond == "Healthy Metagenomic" | h.ss.cond == "SS")]
h.ss.cond <- str_extract(h.ss.cond, "^[^ ]+")

h.ss.aldex <- aldex(data.frame(h.ss),as.character(h.ss.cond),mc.samples=1000)
aldex.plot(h.ss.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. SS")
aldex.plot(h.ss.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs. difference between for healthy vs. SS")
write.table(h.ss.aldex,file="nash_output/1000/H_vs_SS_aldex_output_OTU_level_1000_MC_samples.txt",sep="\t",quote=FALSE)

h.nash <- t(otu.tab)
h.nash.cond <- as.character(groups)
# h.nash <- h.nash[,which(h.nash.cond == "Healthy NASH" | h.nash.cond == "NASH")]
# h.nash.cond <- h.nash.cond[which(h.nash.cond == "Healthy NASH" | h.nash.cond == "NASH")]
h.nash <- h.nash[,which(h.nash.cond == "Healthy" | h.nash.cond == "Healthy Metagenomic" | h.nash.cond == "NASH" | h.nash.cond == "NASH Metagenomic")]
h.nash.cond <- h.nash.cond[which(h.nash.cond == "Healthy" | h.nash.cond == "Healthy Metagenomic" | h.nash.cond == "NASH" | h.nash.cond == "NASH Metagenomic")]
h.nash.cond <- str_extract(h.nash.cond, "^[^ ]+")

h.nash.aldex <- aldex(data.frame(h.nash),as.character(h.nash.cond),mc.samples=1000)
aldex.plot(h.nash.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. NASH")
aldex.plot(h.nash.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs difference between for healthy vs. NASH")
write.table(h.nash.aldex,file="nash_output/1000/H_vs_NASH_aldex_output_OTU_level_1000_MC_samples.txt",sep="\t",quote=FALSE)

h.metnash <- t(otu.tab)
h.metnash.cond <- as.character(groups)
h.metnash <- h.metnash[,which(h.metnash.cond == "Healthy Metagenomic" | h.metnash.cond == "NASH Metagenomic")]
h.metnash.cond <- h.metnash.cond[which(h.metnash.cond == "Healthy Metagenomic" | h.metnash.cond == "NASH Metagenomic")]
h.metnash.cond <- str_extract(h.metnash.cond, "^[^ ]+")

h.metnash.aldex <- aldex(data.frame(h.metnash),as.character(h.metnash.cond),mc.samples=1000)
aldex.plot(h.metnash.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. extreme NASH")
aldex.plot(h.metnash.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs difference between for healthy vs. extreme NASH")
write.table(h.metnash.aldex,file="nash_output/1000/H_vs_extreme_NASH_aldex_output_OTU_level_1000_MC_samples.txt",sep="\t",quote=FALSE)

h.nafld <- t(otu.tab)
h.nafld.cond <- as.character(groups)
h.nafld.cond <- str_extract(h.nafld.cond, "^[^ ]+")
h.nafld.cond[which(h.nafld.cond == "NASH")] <- "NAFLD"
h.nafld.cond[which(h.nafld.cond == "SS")] <- "NAFLD"

h.nafld.aldex <- aldex(data.frame(h.nafld),as.character(h.nafld.cond),mc.samples=1000)
aldex.plot(h.nafld.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. NAFLD")
aldex.plot(h.nafld.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs difference between for healthy vs. NAFLD")
write.table(h.nafld.aldex,file="nash_output/1000/H_vs_NAFLD_aldex_output_OTU_level_1000_MC_samples.txt",sep="\t",quote=FALSE)


write.table(taxonomy,file="nash_output/1000/taxonomy_map.txt",sep="\t",quote=FALSE,col.names=FALSE)

mycolor <- c(col2rgb("turquoise4"))
red <- mycolor[1]
green <- mycolor[2]
blue <- mycolor[3]
mycolor <- rgb(red/255, green/255, blue/255, 0.3)

# remove and print out OTUs that are not in metnash comparison
print(h.ss.aldex[which(! rownames(h.ss.aldex) %in% rownames(h.metnash.aldex)),])
# (none)
print(h.nash.aldex[which(! rownames(h.nash.aldex) %in% rownames(h.metnash.aldex)),])
# (none)

shared.otus <- intersect(intersect(rownames(h.ss.aldex),rownames(h.nash.aldex)),rownames(h.metnash.aldex))

h.ss.aldex <- h.ss.aldex[which(rownames(h.ss.aldex) %in% shared.otus),]
h.nash.aldex <- h.nash.aldex[which(rownames(h.nash.aldex) %in% shared.otus),]
h.metnash.aldex <- h.metnash.aldex[which(rownames(h.metnash.aldex) %in% shared.otus),]

# sanity check that things are int he same order
print(paste("OTUs in the same order in all ALDEX comparisons:",all.equal(rownames(h.ss.aldex),rownames(h.nash.aldex),rownames(h.metnash.aldex))))

plot(h.ss.aldex$effect, h.nash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs SS compared to healthy vs NASH",xlab="Healthy vs. SS",ylab="Healthy vs. NASH")
cor(h.ss.aldex$effect, y = h.nash.aldex$effect, use = "everything", method = "spearman")
# [1] 0.7192259

plot(h.ss.aldex$effect, h.metnash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs SS compared to healthy vs extreme NASH",xlab="Healthy vs. SS",ylab="Healthy vs. extreme NASH")
cor(h.ss.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman")
# [1] 0.4517627

plot(h.nash.aldex$effect, h.metnash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs NASH compared to healthy vs extreme NASH",xlab="Healthy vs. NASH",ylab="Healthy vs. extreme NASH")
cor(h.nash.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman")
# [1] 0.6735478

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

effectgroup <- rep("other",nrow(h.metnash.aldex))
decile <- round(length(effectgroup)/10)
effectgroup[1:decile] <- "top"
effectgroup[(length(effectgroup) - decile + 1):length(effectgroup)] <- "bottom"
effectgroup <- as.factor(effectgroup)

palette(c(myblue, mygray, orchid))

plot(h.ss.aldex$effect, h.nash.aldex$effect, pch=19,col=effectgroup, main="Effect sizes of healthy vs SS compared to healthy vs NASH",xlab="Healthy vs. SS",ylab="Healthy vs. NASH")
print("Spearman correlation for ss vs nash")
print(cor(h.ss.aldex$effect, y = h.nash.aldex$effect, use = "everything", method = "spearman"))
# [1] 0.7192259

plot(h.ss.aldex$effect, h.metnash.aldex$effect, pch=19,col=effectgroup, main="Effect sizes of healthy vs SS compared to healthy vs extreme NASH",xlab="Healthy vs. SS",ylab="Healthy vs. extreme NASH")
print("Spearman correlation for ss vs metnash")
print(cor(h.ss.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman"))
# [1] 0.4517627

plot(h.nash.aldex$effect, h.metnash.aldex$effect, pch=19,col=effectgroup, main="Effect sizes of healthy vs NASH compared to healthy vs extreme NASH",xlab="Healthy vs. NASH",ylab="Healthy vs. extreme NASH")
print("Spearman correlation for nash vs metnash")
print(cor(h.nash.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman"))
# [1] 0.6735478

dev.off()

print("Average difference in distance from zero of healthy vs. nash compared to healthy vs. ss:")
print(summary(abs(h.nash.aldex$effect) - abs(h.ss.aldex$effect)))

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.53880 -0.12360 -0.02656 -0.04624  0.04314  0.43230

print("Average difference in distance from zero of healthy vs. extreme nash compared to healthy vs. ss:")
print(summary(abs(h.metnash.aldex$effect) - abs(h.ss.aldex$effect)))

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.53930 -0.08681  0.01607  0.03540  0.16870  0.75490 

print("Average difference in distance from zero of healthy vs. extreme nash compared to healthy vs. nash:")
print(summary(abs(h.metnash.aldex$effect) - abs(h.nash.aldex$effect)))

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.30810 -0.03790  0.06341  0.08165  0.18070  0.98370 

print("aldex metnash top decile rownames:")
print(paste(rownames(h.metnash.aldex)[1:decile],collapse=", "))
# [1] "36, 84, 65, 15, 318, 538, 178, 30, 1031, 64, 61, 72, 170, 80, 100, 160"
print("aldex metnash bottom decile rownames:")
print(paste(rownames(h.metnash.aldex)[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)],collapse=", "))
# [1] "40, 16, 41, 28, 2, 17, 383, 0, 12, 1445, 18, 54, 23, 14, 8, 31"

print("aldex metnash top decile effect sizes")
print(h.metnash.aldex$effect[1:decile])
# [1] 1.0188783 0.7799609 0.6575043 0.5881891 0.5739999 0.5218084 0.4393049
# [8] 0.4345936 0.4340999 0.4340738 0.3666616 0.3566753 0.3543986 0.3353891
# [15] 0.3317879 0.3093140
print(paste(round(h.metnash.aldex$effect[1:decile], digits=3),collapse="\n"))
# [1] "1.019\n0.78\n0.658\n0.588\n0.574\n0.522\n0.439\n0.435\n0.434\n0.434\n0.367\n0.357\n0.354\n0.335\n0.332\n0.309"

print("aldex metnash bottom decile effect sizes:")
print(h.metnash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)])
# [1] -0.4457842 -0.4541306 -0.4784393 -0.4925164 -0.5089918 -0.5536379
# [7] -0.5541652 -0.5703956 -0.6225181 -0.6825579 -0.7144518 -0.7197835
# [13] -0.7255898 -0.7930968 -0.8011380 -0.8621549
print(paste(round(h.metnash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n"))
# [1] "-0.446\n-0.454\n-0.478\n-0.493\n-0.509\n-0.554\n-0.554\n-0.57\n-0.623\n-0.683\n-0.714\n-0.72\n-0.726\n-0.793\n-0.801\n-0.862"

print("aldex ss top decile effect sizes")
print(h.ss.aldex$effect[1:decile])
# [1] -0.26395091  0.47156929  0.33675393  0.09773811  0.20240223  0.24970541
# [7]  0.07239912 -0.23653322 -0.43659622  0.20427478  0.14179955  0.19110348
# [13]  0.22378771  0.25877459  0.03292758  0.29565657
print(paste(round(h.ss.aldex$effect[1:decile], digits=3),collapse="\n"))
# [1] "-0.264\n0.472\n0.337\n0.098\n0.202\n0.25\n0.072\n-0.237\n-0.437\n0.204\n0.142\n0.191\n0.224\n0.259\n0.033\n0.296"

print("aldex ss bottom decile effect sizes:")
print(h.ss.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)])
# [1] -0.36064479 -0.29070205 -0.48702868  0.02983030 -0.48101876  0.23698117
# [7] -0.40115437 -0.17173978 -0.30421444 -0.07988267 -0.17284943 -0.30323925
# [13] -0.55768030 -0.22311316 -0.06577780 -0.56203728
print(paste(round(h.ss.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n"))
# [1] "-0.361\n-0.291\n-0.487\n0.03\n-0.481\n0.237\n-0.401\n-0.172\n-0.304\n-0.08\n-0.173\n-0.303\n-0.558\n-0.223\n-0.066\n-0.562"

print("aldex nash top decile effect sizes")
print(h.nash.aldex$effect[1:decile])
# [1]  0.035193852  0.604990880  0.250966678  0.530065918  0.140834648
# [6]  0.138317063  0.180932733  0.229116084 -0.028673620  0.197925882
# [11]  0.245302189 -0.004954869  0.266927646  0.271976324  0.024941330
# [16]  0.324412304
print(paste(round(h.nash.aldex$effect[1:decile], digits=3),collapse="\n"))
# [1] "0.035\n0.605\n0.251\n0.53\n0.141\n0.138\n0.181\n0.229\n-0.029\n0.198\n0.245\n-0.005\n0.267\n0.272\n0.025\n0.324"
print("aldex nash bottom decile effect sizes:")
print(h.nash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)])
# [1] -0.3066982 -0.3994963 -0.2917866 -0.3244448 -0.3970767 -0.1493649
# [7] -0.3082160 -0.2208430 -0.3500466 -0.2950210 -0.3250654 -0.3559420
# [13] -0.5347665 -0.3458652 -0.1368042 -0.4541884
print(paste(round(h.nash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n"))
# [1] "-0.307\n-0.399\n-0.292\n-0.324\n-0.397\n-0.149\n-0.308\n-0.221\n-0.35\n-0.295\n-0.325\n-0.356\n-0.535\n-0.346\n-0.137\n-0.454"

top <- c(1:decile)
bottom <- c((length(effectgroup) - decile + 1):length(effectgroup))

print(summary(h.nash.aldex$effect[top] - h.ss.aldex$effect[top]))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.19610 -0.02138  0.03595  0.09790  0.17490  0.46560 
print(summary(h.metnash.aldex$effect[top] - h.ss.aldex$effect[top]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.01366 0.21000 0.30360 0.38090 0.40130 1.28300 
print(summary(h.metnash.aldex$effect[top] - h.nash.aldex$effect[top]))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.0151  0.1129  0.2473  0.2830  0.3893  0.9837 

print(summary(h.nash.aldex$effect[bottom] -  h.ss.aldex$effect[bottom]))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.38630 -0.13010 -0.05090 -0.06258  0.06145  0.19520 
print(summary(h.metnash.aldex$effect[bottom] -  h.ss.aldex$effect[bottom]))
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.790600 -0.548700 -0.358500 -0.361600 -0.160800  0.008589 
print(summary(h.metnash.aldex$effect[bottom] -  h.nash.aldex$effect[bottom]))
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.66430 -0.39310 -0.31100 -0.29900 -0.18200 -0.05463 

otu.genus <- c(as.character(taxonomy))
otu.family <- c(as.character(taxonomy))
bootstrap <- c(as.character(taxonomy))

for (i in c(1:length(taxonomy))) {
  otu.genus[i] <- strsplit(otu.genus[i],c(";"))[[1]][6]
	otu.family[i] <- strsplit(otu.family[i],c(";"))[[1]][5]
	bootstrap[i] <- strsplit(bootstrap[i],c(";"))[[1]][7]
}

names(otu.genus) <- names(taxonomy)
names(otu.family) <- names(taxonomy)
names(bootstrap) <- names(taxonomy)

print("top decile bootstrap")
print(paste(bootstrap[match(rownames(h.metnash.aldex)[top],names(bootstrap))], collapse="\n"))
print("top decile OTU genus")
print(paste(otu.genus[match(rownames(h.metnash.aldex)[top],names(otu.genus))], collapse="\n"))
print("top decile OTU family")
print(paste(otu.family[match(rownames(h.metnash.aldex)[top],names(otu.family))], collapse="\n"))

print("bottom decile bootstrap")
print(paste(bootstrap[match(rownames(h.metnash.aldex)[bottom],names(bootstrap))], collapse="\n"))
print("bottom decile OTU genus")
print(paste(otu.genus[match(rownames(h.metnash.aldex)[bottom],names(otu.genus))], collapse="\n"))
print("bottom decile OTU family")
print(paste(otu.family[match(rownames(h.metnash.aldex)[bottom],names(otu.family))], collapse="\n"))


# GENUS LEVEL

otu.genus <- c(as.character(taxonomy))

for (i in c(1:length(taxonomy))) {
  otu.genus[i] <- strsplit(otu.genus[i],c(";"))[[1]][6]
}

otu.tab.genus <- aggregate(t(otu.tab),list(otu.genus),sum)
rownames(otu.tab.genus) <- otu.tab.genus$Group.1
otu.tab.genus <- otu.tab.genus[,c(2:ncol(otu.tab.genus))]

pdf("nash_output/1000/ALDEx_genus_level_plots.pdf")

h.ss.genus <- otu.tab.genus
h.ss.genus.cond <- groups
# h.ss.genus <- h.ss.genus[,which(h.ss.genus.cond == "Healthy SS" | h.ss.genus.cond == "SS")]
# h.ss.genus.cond <- h.ss.genus.cond[which(h.ss.genus.cond == "Healthy SS" | h.ss.genus.cond == "SS")]
h.ss.genus <- h.ss.genus[,which(h.ss.genus.cond == "Healthy" | h.ss.genus.cond == "Healthy Metagenomic" | h.ss.genus.cond == "SS")]
h.ss.genus.cond <- h.ss.genus.cond[which(h.ss.genus.cond == "Healthy" | h.ss.genus.cond == "Healthy Metagenomic" | h.ss.genus.cond == "SS")]
h.ss.genus.cond <- str_extract(h.ss.genus.cond, "^[^ ]+")

h.ss.genus.aldex <- aldex(data.frame(h.ss.genus),as.character(h.ss.genus.cond),mc.samples=1000)
aldex.plot(h.ss.genus.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. SS")
aldex.plot(h.ss.genus.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs. difference between for healthy vs. SS")
write.table(h.ss.genus.aldex,file="nash_output/1000/H_vs_SS_aldex_output_genus_level_1000_MC_samples.txt",sep="\t",quote=FALSE)

h.nash.genus <- otu.tab.genus
h.nash.genus.cond <- groups
# h.nash.genus <- h.nash.genus[,which(h.nash.genus.cond == "Healthy NASH" | h.nash.genus.cond == "NASH")]
# h.nash.genus.cond <- h.nash.genus.cond[which(h.nash.genus.cond == "Healthy NASH" | h.nash.genus.cond == "NASH")]
h.nash.genus <- h.nash.genus[,which(h.nash.genus.cond == "Healthy" | h.nash.genus.cond == "Healthy Metagenomic" | h.nash.genus.cond == "NASH" | h.nash.genus.cond == "NASH Metagenomic")]
h.nash.genus.cond <- h.nash.genus.cond[which(h.nash.genus.cond == "Healthy" | h.nash.genus.cond == "Healthy Metagenomic" | h.nash.genus.cond == "NASH" | h.nash.genus.cond == "NASH Metagenomic")]
h.nash.genus.cond <- str_extract(h.nash.genus.cond, "^[^ ]+")

h.nash.genus.aldex <- aldex(data.frame(h.nash.genus),as.character(h.nash.genus.cond),mc.samples=1000)
aldex.plot(h.nash.genus.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. NASH")
aldex.plot(h.nash.genus.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs difference between for healthy vs. NASH")
write.table(h.nash.genus.aldex,file="nash_output/1000/H_vs_NASH_aldex_output_genus_level_1000_MC_samples.txt",sep="\t",quote=FALSE)

h.metnash.genus <- otu.tab.genus
h.metnash.genus.cond <- groups
h.metnash.genus <- h.metnash.genus[,which(h.metnash.genus.cond == "Healthy Metagenomic" | h.metnash.genus.cond == "NASH Metagenomic")]
h.metnash.genus.cond <- h.metnash.genus.cond[which(h.metnash.genus.cond == "Healthy Metagenomic" | h.metnash.genus.cond == "NASH Metagenomic")]

h.metnash.genus.aldex <- aldex(data.frame(h.metnash.genus),as.character(h.metnash.genus.cond),mc.samples=1000)
aldex.plot(h.metnash.genus.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. extreme NASH")
aldex.plot(h.metnash.genus.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs difference between for healthy vs. extreme NASH")
write.table(h.metnash.genus.aldex,file="nash_output/1000/H_vs_extreme_NASH_aldex_output_genus_level_1000_MC_samples.txt",sep="\t",quote=FALSE)

h.nafld.genus <- otu.tab.genus
h.nafld.genus.cond <- as.character(groups)
h.nafld.genus.cond <- str_extract(h.nafld.genus.cond, "^[^ ]+")
h.nafld.genus.cond[which(h.nafld.genus.cond == "NASH")] <- "NAFLD"
h.nafld.genus.cond[which(h.nafld.genus.cond == "SS")] <- "NAFLD"

h.nafld.genus.aldex <- aldex(data.frame(h.nafld.genus),as.character(h.nafld.genus.cond),mc.samples=1000)
aldex.plot(h.nafld.genus.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. NAFLD")
aldex.plot(h.nafld.genus.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs difference between for healthy vs. NAFLD")
write.table(h.nafld.genus.aldex,file="nash_output/1000/H_vs_NAFLD_aldex_output_genus_level_1000_MC_samples.txt",sep="\t",quote=FALSE)


dev.off()

# FAMILY LEVEL

otu.family <- c(as.character(taxonomy))

for (i in c(1:length(taxonomy))) {
  otu.family[i] <- strsplit(otu.family[i],c(";"))[[1]][5]
}

otu.tab.family <- aggregate(t(otu.tab),list(otu.family),sum)
rownames(otu.tab.family) <- otu.tab.family$Group.1
otu.tab.family <- otu.tab.family[,c(2:ncol(otu.tab.family))]

pdf("nash_output/1000/ALDEx_family_level_plots.pdf")

h.ss.family <- otu.tab.family
h.ss.family.cond <- groups
# h.ss.family <- h.ss.family[,which(h.ss.family.cond == "Healthy SS" | h.ss.family.cond == "SS")]
# h.ss.family.cond <- h.ss.family.cond[which(h.ss.family.cond == "Healthy SS" | h.ss.family.cond == "SS")]
h.ss.family <- h.ss.family[,which(h.ss.family.cond == "Healthy" | h.ss.family.cond == "Healthy Metagenomic" | h.ss.family.cond == "SS")]
h.ss.family.cond <- h.ss.family.cond[which(h.ss.family.cond == "Healthy" | h.ss.family.cond == "Healthy Metagenomic" | h.ss.family.cond == "SS")]
h.ss.family.cond <- str_extract(h.ss.family.cond, "^[^ ]+")

h.ss.family.aldex <- aldex(data.frame(h.ss.family),as.character(h.ss.family.cond),mc.samples=1000)
aldex.plot(h.ss.family.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. SS")
aldex.plot(h.ss.family.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs. difference between for healthy vs. SS")
write.table(h.ss.family.aldex,file="nash_output/1000/H_vs_SS_aldex_output_family_level_1000_MC_samples.txt",sep="\t",quote=FALSE)

h.nash.family <- otu.tab.family
h.nash.family.cond <- groups
# h.nash.family <- h.nash.family[,which(h.nash.family.cond == "Healthy NASH" | h.nash.family.cond == "NASH")]
# h.nash.family.cond <- h.nash.family.cond[which(h.nash.family.cond == "Healthy NASH" | h.nash.family.cond == "NASH")]
h.nash.family <- h.nash.family[,which(h.nash.family.cond == "Healthy" | h.nash.family.cond == "Healthy Metagenomic" | h.nash.family.cond == "NASH" | h.nash.family.cond == "NASH Metagenomic")]
h.nash.family.cond <- h.nash.family.cond[which(h.nash.family.cond == "Healthy" | h.nash.family.cond == "Healthy Metagenomic" | h.nash.family.cond == "NASH" | h.nash.family.cond == "NASH Metagenomic")]
h.nash.family.cond <- str_extract(h.nash.family.cond, "^[^ ]+")

h.nash.family.aldex <- aldex(data.frame(h.nash.family),as.character(h.nash.family.cond),mc.samples=1000)
aldex.plot(h.nash.family.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. NASH")
aldex.plot(h.nash.family.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs difference between for healthy vs. NASH")
write.table(h.nash.family.aldex,file="nash_output/1000/H_vs_NASH_aldex_output_family_level_1000_MC_samples.txt",sep="\t",quote=FALSE)

h.metnash.family <- otu.tab.family
h.metnash.family.cond <- groups
h.metnash.family <- h.metnash.family[,which(h.metnash.family.cond == "Healthy Metagenomic" | h.metnash.family.cond == "NASH Metagenomic")]
h.metnash.family.cond <- h.metnash.family.cond[which(h.metnash.family.cond == "Healthy Metagenomic" | h.metnash.family.cond == "NASH Metagenomic")]

h.metnash.family.aldex <- aldex(data.frame(h.metnash.family),as.character(h.metnash.family.cond),mc.samples=1000)
aldex.plot(h.metnash.family.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. extreme NASH")
aldex.plot(h.metnash.family.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs difference between for healthy vs. extreme NASH")
write.table(h.metnash.family.aldex,file="nash_output/1000/H_vs_extreme_NASH_aldex_output_family_level_1000_MC_samples.txt",sep="\t",quote=FALSE)

h.nafld.family <- otu.tab.family
h.nafld.family.cond <- as.character(groups)
h.nafld.family.cond <- str_extract(h.nafld.family.cond, "^[^ ]+")
h.nafld.family.cond[which(h.nafld.family.cond == "NASH")] <- "NAFLD"
h.nafld.family.cond[which(h.nafld.family.cond == "SS")] <- "NAFLD"

h.nafld.family.aldex <- aldex(data.frame(h.nafld.family),as.character(h.nafld.family.cond),mc.samples=1000)
aldex.plot(h.nafld.family.aldex,type="MA",ylab="ratio",xlab="average")
title(main="Bland-Altman style plot for healthy vs. NAFLD")
aldex.plot(h.nafld.family.aldex,type="MW",xlab="Difference within",ylab="Difference between")
title(main="Difference within vs difference between for healthy vs. NAFLD")
write.table(h.nafld.family.aldex,file="nash_output/1000/H_vs_NAFLD_aldex_output_family_level_1000_MC_samples.txt",sep="\t",quote=FALSE)


dev.off()

#output sample groupings

names(groups) <- rownames(otu.tab)
write.table(groups,file="nash_output/1000/sample_grouping.txt",sep="\t",quote=FALSE,col.names=FALSE)
