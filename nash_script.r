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

# randomly pick which healthy sample will be controls for the SS or the NASH comparison
healthy.index <- which(groups == "Healthy")
control.size <- floor(length(healthy.index)/2)
healthy.ss <- sample(healthy.index,control.size)
healthy.nash <- healthy.index[which(!(healthy.index %in% healthy.ss))]

groups <- as.character(groups)

groups[healthy.ss] <- "Healthy SS"
groups[healthy.nash] <- "Healthy NASH"

h.ss <- t(otu.tab)
h.ss.cond <- groups
h.ss <- h.ss[,which(h.ss.cond == "Healthy SS" | h.ss.cond == "SS")]
h.ss.cond <- h.ss.cond[which(h.ss.cond == "Healthy SS" | h.ss.cond == "SS")]

h.ss.aldex <- aldex(data.frame(h.ss),as.character(h.ss.cond))
aldex.plot(h.ss.aldex,type="MA",main="Bland-Altman style plot for healthy vs. SS",ylab="ratio",xlab="average")
aldex.plot(h.ss.aldex,type="MW",main="Difference within vs. difference between for healthy vs. SS",xlab="Difference within",ylab="Difference between")

h.nash <- t(otu.tab)
h.nash.cond <- groups
h.nash <- h.nash[,which(h.nash.cond == "Healthy NASH" | h.nash.cond == "NASH")]
h.nash.cond <- h.nash.cond[which(h.nash.cond == "Healthy NASH" | h.nash.cond == "NASH")]

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

# remove and print out OTUs that are not in metnash comparison
print(h.ss.aldex[which(! rownames(h.ss.aldex) %in% rownames(h.metnash.aldex)),])
# rab.all rab.win.Healthy.SS rab.win.SS  diff.btw diff.win diff.btw.025
# 370 -4.748640          -5.663842  -4.110681 1.5973947 4.797527    -7.496665
# 492 -5.019791          -5.288625  -4.877867 0.3968992 4.089898    -7.819287
# diff.btw.975 diff.win.025 diff.win.975     effect effect.025 effect.975
# 370     12.24924    0.8289231     13.91703 0.30732376  -2.391286   5.309102
# 492     10.06937    0.7443974     11.93355 0.08233881  -3.626621   4.674094
# overlap     we.ep    we.eBH     wi.ep    wi.eBH
# 370 0.3310226 0.2971612 0.7034279 0.3034283 0.6633413
# 492 0.4523397 0.5161322 0.8289618 0.5295436 0.8245989
print(h.nash.aldex[which(! rownames(h.nash.aldex) %in% rownames(h.metnash.aldex)),])
# rab.all rab.win.Healthy.NASH rab.win.NASH   diff.btw diff.win
# 370 -4.986099            -4.771754    -5.083839 -0.3905306 4.223909
# diff.btw.025 diff.btw.975 diff.win.025 diff.win.975      effect effect.025
# 370    -9.059532     9.192058     0.807064     13.39104 -0.08618866  -3.717694
# effect.975   overlap     we.ep    we.eBH     wi.ep   wi.eBH
# 370   3.995423 0.4513889 0.4862919 0.9473435 0.5417632 0.962178

h.ss.aldex <- h.ss.aldex[which(rownames(h.ss.aldex) %in% rownames(h.metnash.aldex)),]
h.nash.aldex <- h.nash.aldex[which(rownames(h.nash.aldex) %in% rownames(h.metnash.aldex)),]

# sanity check that things are int he same order
print(paste("OTUs in the same order in all ALDEX comparisons:",all.equal(rownames(h.ss.aldex),rownames(h.nash.aldex),rownames(h.metnash.aldex))))

plot(h.ss.aldex$effect, h.nash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs SS compared to healthy vs NASH",xlab="Healthy vs. SS",ylab="Healthy vs. NASH")
cor(h.ss.aldex$effect, y = h.nash.aldex$effect, use = "everything", method = "spearman")

plot(h.ss.aldex$effect, h.metnash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs SS compared to healthy vs extreme NASH",xlab="Healthy vs. SS",ylab="Healthy vs. extreme NASH")
cor(h.ss.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman")

plot(h.nash.aldex$effect, h.metnash.aldex$effect, pch=19,col=mycolor, main="Effect sizes of healthy vs NASH compared to healthy vs extreme NASH",xlab="Healthy vs. NASH",ylab="Healthy vs. extreme NASH")
cor(h.nash.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman")

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
# [1] 0.1760099

plot(h.ss.aldex$effect, h.metnash.aldex$effect, pch=19,col=effectgroup, main="Effect sizes of healthy vs SS compared to healthy vs extreme NASH",xlab="Healthy vs. SS",ylab="Healthy vs. extreme NASH")
print("Spearman correlation for ss vs metnash")
print(cor(h.ss.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman"))
# [1] 0.2034683

plot(h.nash.aldex$effect, h.metnash.aldex$effect, pch=19,col=effectgroup, main="Effect sizes of healthy vs NASH compared to healthy vs extreme NASH",xlab="Healthy vs. NASH",ylab="Healthy vs. extreme NASH")
print("Spearman correlation for nash vs metnash")
print(cor(h.nash.aldex$effect, y = h.metnash.aldex$effect, use = "everything", method = "spearman"))
# [1] 0.1633695

dev.off()

print("Average difference in distance from zero of healthy vs. nash compared to healthy vs. ss:")
print(summary(abs(h.nash.aldex$effect) - abs(h.ss.aldex$effect)))

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -1.00000 -0.19000 -0.03317 -0.04747  0.10070  0.94550 

print("Average difference in distance from zero of healthy vs. extreme nash compared to healthy vs. ss:")
print(summary(abs(h.metnash.aldex$effect) - abs(h.ss.aldex$effect)))

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.92850 -0.21580 -0.03809 -0.05042  0.11870  0.74830 

print("Average difference in distance from zero of healthy vs. extreme nash compared to healthy vs. nash:")
print(summary(abs(h.metnash.aldex$effect) - abs(h.nash.aldex$effect)))

# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -1.013000 -0.126900 -0.009027 -0.002950  0.141300  0.764500 

print("aldex metnash top decile rownames:")
print(paste(rownames(h.metnash.aldex)[1:decile],collapse=", "))
# [1] "36, 84, 120, 422, 538, 635, 15, 318, 64, 937, 142, 174, 185, 88, 65, 1447, 1000, 178, 170, 742, 211, 1331, 332, 155, 1245, 197, 1279, 100, 280, 878, 286, 116, 350, 66, 112, 293, 297, 261, 225, 758, 96, 160, 631, 359, 1097, 80, 30, 191, 205, 990, 177, 236, 387"
print("aldex metnash bottom decile rownames:")
print(paste(rownames(h.metnash.aldex)[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)],collapse=", "))
# [1] "257, 308, 7, 1289, 127, 553, 202, 91, 1303, 50, 254, 287, 103, 846, 83, 714, 259, 10, 329, 1118, 5, 70, 110, 195, 122, 157, 1005, 17, 1101, 40, 0, 180, 16, 145, 815, 558, 23, 383, 79, 2, 996, 224, 320, 8, 54, 14, 18, 12, 1202, 111, 718, 31, 176"

print("aldex metnash top decile effect sizes")
print(h.metnash.aldex$effect[1:decile])
# [1] 0.9985202 0.8754650 0.8369680 0.6492463 0.6217012 0.6125587 0.5952462
# [8] 0.5851389 0.5403535 0.5363887 0.5330805 0.5285141 0.5134040 0.4997956
# [15] 0.4785943 0.4690567 0.4664723 0.4659647 0.4555221 0.4549954 0.4520423
# [22] 0.4311886 0.4253795 0.4154577 0.4153501 0.4114170 0.4110055 0.4056718
# [29] 0.4043549 0.4024815 0.3949997 0.3917501 0.3713569 0.3670031 0.3659648
# [36] 0.3652170 0.3467735 0.3419681 0.3399944 0.3371219 0.3370802 0.3358125
# [43] 0.3330868 0.3327869 0.3304090 0.3193641 0.3149830 0.3145493 0.3129694
# [50] 0.3093094 0.3036032 0.3018303 0.2982272
print(paste(round(h.metnash.aldex$effect[1:decile], digits=3),collapse="\n"))
# [1] "0.999\n0.875\n0.837\n0.649\n0.622\n0.613\n0.595\n0.585\n0.54\n0.536\n0.533\n0.529\n0.513\n0.5\n0.479\n0.469\n0.466\n0.466\n0.456\n0.455\n0.452\n0.431\n0.425\n0.415\n0.415\n0.411\n0.411\n0.406\n0.404\n0.402\n0.395\n0.392\n0.371\n0.367\n0.366\n0.365\n0.347\n0.342\n0.34\n0.337\n0.337\n0.336\n0.333\n0.333\n0.33\n0.319\n0.315\n0.315\n0.313\n0.309\n0.304\n0.302\n0.298"

print("aldex metnash bottom decile effect sizes:")
print(h.metnash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)])
# [1] -0.3422843 -0.3464672 -0.3497061 -0.3503735 -0.3547473 -0.3618035
# [7] -0.3730018 -0.3780739 -0.3787961 -0.3788051 -0.3814072 -0.3858562
# [13] -0.3928465 -0.3941044 -0.3942228 -0.3993499 -0.4122165 -0.4137053
# [19] -0.4162458 -0.4216362 -0.4336658 -0.4397961 -0.4416509 -0.4493063
# [25] -0.4599566 -0.4729102 -0.4746872 -0.4870829 -0.4892602 -0.5026108
# [31] -0.5070114 -0.5139160 -0.5188113 -0.5285926 -0.5430143 -0.5459578
# [37] -0.5497792 -0.5510761 -0.5514221 -0.5659363 -0.5842574 -0.6213035
# [43] -0.6213849 -0.6543683 -0.6548169 -0.6694540 -0.6858669 -0.6913024
# [49] -0.6959176 -0.6988220 -0.7451949 -0.7501492 -0.9543258
print(paste(round(h.metnash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n"))
# [1] "-0.342\n-0.346\n-0.35\n-0.35\n-0.355\n-0.362\n-0.373\n-0.378\n-0.379\n-0.379\n-0.381\n-0.386\n-0.393\n-0.394\n-0.394\n-0.399\n-0.412\n-0.414\n-0.416\n-0.422\n-0.434\n-0.44\n-0.442\n-0.449\n-0.46\n-0.473\n-0.475\n-0.487\n-0.489\n-0.503\n-0.507\n-0.514\n-0.519\n-0.529\n-0.543\n-0.546\n-0.55\n-0.551\n-0.551\n-0.566\n-0.584\n-0.621\n-0.621\n-0.654\n-0.655\n-0.669\n-0.686\n-0.691\n-0.696\n-0.699\n-0.745\n-0.75\n-0.954"

print("aldex ss top decile effect sizes")
print(h.ss.aldex$effect[1:decile])
# [1] -0.70086179  0.32963148 -0.08868995 -0.20361595  0.02229984 -0.14598706
# [7] -0.26973240  0.02387802  0.13599104 -0.04740110 -0.07151937  0.13416052
# [13] -0.11287012  0.05268304  0.27767967  0.42584547 -0.50362351 -0.17939955
# [19]  0.15502683  0.85864333 -0.08690670  0.31419080  0.12091298 -0.06030971
# [25]  0.76749361 -0.03706545 -0.15388961 -0.26290312  0.58915251  0.49025269
# [31] -0.33754824  0.02128101  0.10833983  0.03534174  0.02267682  0.01158103
# [37]  0.30037382  0.20429782 -0.06144975 -0.04297356  0.20769389  0.48728114
# [43] -0.13231014 -0.19702417 -0.06206644  0.03135385 -0.70783514 -0.35980225
# [49]  0.07325836  0.61953163 -0.09703600 -1.16856011  0.11726201
print(paste(round(h.ss.aldex$effect[1:decile], digits=3),collapse="\n"))
# [1] "-0.701\n0.33\n-0.089\n-0.204\n0.022\n-0.146\n-0.27\n0.024\n0.136\n-0.047\n-0.072\n0.134\n-0.113\n0.053\n0.278\n0.426\n-0.504\n-0.179\n0.155\n0.859\n-0.087\n0.314\n0.121\n-0.06\n0.767\n-0.037\n-0.154\n-0.263\n0.589\n0.49\n-0.338\n0.021\n0.108\n0.035\n0.023\n0.012\n0.3\n0.204\n-0.061\n-0.043\n0.208\n0.487\n-0.132\n-0.197\n-0.062\n0.031\n-0.708\n-0.36\n0.073\n0.62\n-0.097\n-1.169\n0.117"

print("aldex ss bottom decile effect sizes:")
print(h.ss.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)])
# [1] -0.04528287 -0.17270944 -0.38052010  0.04915930 -0.23166717  0.18786675
# [7] -0.41259164  0.23194486 -0.04597739 -0.38202939 -0.22205060 -0.06387055
# [13] -0.47557425 -0.19777167  0.18453737  0.21187640 -0.67108178  0.47256428
# [19]  0.23127634 -0.25402727 -0.83350300  0.03836738 -0.20930442 -0.35200615
# [25] -0.33401859 -0.12159441 -0.11276543  0.25688158  0.08566400 -0.18530112
# [31] -0.15128794 -0.27668869 -0.15969842 -0.09396876 -0.26705460  0.16606799
# [37] -0.65018577 -0.92772805  0.27635236 -0.22715163 -0.47864137 -0.58019043
# [43] -0.09122315  0.25573645 -0.12805737 -0.29077792 -0.04722916 -0.24982571
# [49] -0.30816195  0.48395380  0.06170587 -1.01008779 -0.31173879
print(paste(round(h.ss.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n"))
# [1] "-0.045\n-0.173\n-0.381\n0.049\n-0.232\n0.188\n-0.413\n0.232\n-0.046\n-0.382\n-0.222\n-0.064\n-0.476\n-0.198\n0.185\n0.212\n-0.671\n0.473\n0.231\n-0.254\n-0.834\n0.038\n-0.209\n-0.352\n-0.334\n-0.122\n-0.113\n0.257\n0.086\n-0.185\n-0.151\n-0.277\n-0.16\n-0.094\n-0.267\n0.166\n-0.65\n-0.928\n0.276\n-0.227\n-0.479\n-0.58\n-0.091\n0.256\n-0.128\n-0.291\n-0.047\n-0.25\n-0.308\n0.484\n0.062\n-1.01\n-0.312"

print("aldex nash top decile effect sizes")
print(h.nash.aldex$effect[1:decile])
# [1] -0.3290460243  0.7265017902  0.1922670224 -0.0035312656  0.0376863921
# [6] -0.0011982329  0.9575208941 -0.0938086906 -0.0926458871  0.3125021301
# [11]  0.1104757669  0.1403730822  0.0107851188  0.1583540484  0.0032257969
# [16] -0.0431842087 -0.2836717785  0.0877664347  0.1708481717  0.1149355946
# [21] -0.0883735586  0.1010410156  0.0009973402  0.2142651449  0.2361030131
# [26] -0.1252414778  0.1570956588 -0.0253328039 -0.0615638995  0.1880988709
# [31] -0.2352463124  0.2410225711  0.4436129551 -0.1751130617  0.2061204967
# [36] -0.3633166575 -0.0330833512 -0.4621218503  0.0862620196  0.1111759835
# [41] -0.3471449044  0.1292514491  0.0621746792 -0.0007338209  0.0755392035
# [46]  0.5122518515  0.3754535652 -0.1263016459 -0.3814788055  0.0879933369
# [51]  0.0390634534 -0.2572648034  0.0580825292
print(paste(round(h.nash.aldex$effect[1:decile], digits=3),collapse="\n"))
# [1] "-0.329\n0.727\n0.192\n-0.004\n0.038\n-0.001\n0.958\n-0.094\n-0.093\n0.313\n0.11\n0.14\n0.011\n0.158\n0.003\n-0.043\n-0.284\n0.088\n0.171\n0.115\n-0.088\n0.101\n0.001\n0.214\n0.236\n-0.125\n0.157\n-0.025\n-0.062\n0.188\n-0.235\n0.241\n0.444\n-0.175\n0.206\n-0.363\n-0.033\n-0.462\n0.086\n0.111\n-0.347\n0.129\n0.062\n-0.001\n0.076\n0.512\n0.375\n-0.126\n-0.381\n0.088\n0.039\n-0.257\n0.058"
print("aldex nash bottom decile effect sizes:")
print(h.nash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)])
# [1]  0.0319145510 -0.0002590834  0.1346245649  0.2874443089 -0.0325989631
# [6]  0.1788515663 -0.0990311867 -0.1292862613 -0.3556286887  0.0143885466
# [11] -0.1206216732 -0.1888588974 -0.5002189808 -0.0409717576 -0.6073349531
# [16] -0.1455875242 -0.1148739663 -0.0882708626  0.2093166365 -0.3823249336
# [21]  0.2244937292 -0.0816762237 -0.3096203650  0.2338474224 -0.1037390394
# [26]  0.0075137161 -0.2175400308 -0.0541820502  0.1603304627 -0.3734777213
# [31]  0.0465939446 -0.2504298021 -0.6161754850  0.2151494223 -0.3079874137
# [36]  0.0404993577 -0.3904718763  0.0851001273  0.1277322213 -0.6214785841
# [41]  0.2038593841  0.0800074619 -0.5416815878 -0.0659487594 -0.4172416805
# [46] -0.2339427659 -0.1681916123 -0.1645333467 -0.2795182356 -0.2454305623
# [51] -0.1674095434  0.0183521300  0.1898556764
print(paste(round(h.nash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n"))
# [1] "0.032\n0\n0.135\n0.287\n-0.033\n0.179\n-0.099\n-0.129\n-0.356\n0.014\n-0.121\n-0.189\n-0.5\n-0.041\n-0.607\n-0.146\n-0.115\n-0.088\n0.209\n-0.382\n0.224\n-0.082\n-0.31\n0.234\n-0.104\n0.008\n-0.218\n-0.054\n0.16\n-0.373\n0.047\n-0.25\n-0.616\n0.215\n-0.308\n0.04\n-0.39\n0.085\n0.128\n-0.621\n0.204\n0.08\n-0.542\n-0.066\n-0.417\n-0.234\n-0.168\n-0.165\n-0.28\n-0.245\n-0.167\n0.018\n0.19"

print(summary(abs(h.nash.aldex$effect[deciles]) - abs(h.ss.aldex$effect[deciles])))
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.99170 -0.22770 -0.05858 -0.06974  0.10540  0.68780 
print(summary(abs(h.metnash.aldex$effect[deciles]) - abs(h.ss.aldex$effect[deciles])))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.8667  0.1185  0.2657  0.2091  0.3691  0.7483 
print(summary(abs(h.metnash.aldex$effect[deciles]) - abs(h.nash.aldex$effect[deciles])))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.3623  0.1843  0.2724  0.2788  0.4211  0.7645 


# BARPLOT

# ALDEX BY GENUS TO COMPARE WITH QPCR
