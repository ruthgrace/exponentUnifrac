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

# remove and print out OTUs that are not in metnash comparison
print(h.ss.aldex[which(! rownames(h.ss.aldex) %in% rownames(h.metnash.aldex)),])
print(h.nash.aldex[which(! rownames(h.nash.aldex) %in% rownames(h.metnash.aldex)),])

h.ss.aldex <- h.ss.aldex[which(rownames(h.ss.aldex) %in% rownames(h.metnash.aldex)),]
h.nash.aldex <- h.nash.aldex[which(rownames(h.nash.aldex) %in% rownames(h.metnash.aldex)),]

# sanity check that things are int he same order
print(paste("OTUs in the same order in all ALDEX comparisons:",all.equal(rownames(h.ss.aldex),rownames(h.nash.aldex),rownames(h.metnash.aldex))))

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

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.72130 -0.14550 -0.02135 -0.04785  0.05492  0.62890 

print("Average difference in distance from zero of healthy vs. extreme nash compared to healthy vs. ss:")
summary(abs(h.metnash.aldex$effect) - abs(h.ss.aldex$effect))

# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -0.96430 -0.13720 -0.02211 -0.03677  0.06846  0.87010 

print("Average difference in distance from zero of healthy vs. extreme nash compared to healthy vs. nash:")
summary(abs(h.metnash.aldex$effect) - abs(h.nash.aldex$effect))

# Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# -0.386000 -0.064030  0.007811  0.011070  0.084120  0.583700 

paste(rownames(h.metnash.aldex)[1:decile],collapse=", ")
# [1] "36, 84, 120, 538, 635, 422, 15, 88, 174, 318, 178, 142, 1000, 937, 64, 65, 211, 1447, 1279, 742, 280, 1245, 332, 1331, 286, 155, 170, 185, 878, 116, 100, 197, 30, 338, 160, 66, 297, 350, 225, 293, 96, 387, 631, 112, 1097, 261, 72, 236, 359, 61, 758, 990, 80"

paste(rownames(h.metnash.aldex)[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)],collapse=", ")
# [1] "28, 10, 316, 57, 254, 344, 7, 127, 5, 1289, 1118, 259, 257, 83, 103, 50, 91, 846, 329, 202, 122, 287, 70, 195, 110, 558, 157, 17, 383, 1101, 815, 1005, 79, 996, 180, 40, 0, 23, 54, 16, 2, 145, 224, 14, 320, 8, 1202, 18, 111, 12, 718, 31, 176"

h.metnash.aldex$effect[1:decile]
# [1] 0.9808669 0.9594850 0.7962744 0.6447331 0.6374616 0.6225576 0.5719270
# [8] 0.5665163 0.5540923 0.5252060 0.5180265 0.5130678 0.5097834 0.5058127
# [15] 0.5025663 0.4849291 0.4843138 0.4551219 0.4452383 0.4451407 0.4224283
# [22] 0.4220231 0.4146682 0.4078054 0.4049878 0.3967917 0.3927877 0.3926212
# [29] 0.3830912 0.3807892 0.3782888 0.3750775 0.3744567 0.3735914 0.3669307
# [36] 0.3659202 0.3638281 0.3588629 0.3567970 0.3521021 0.3500685 0.3371275
# [43] 0.3362127 0.3357944 0.3305064 0.3296209 0.3226165 0.3194716 0.3162951
# [50] 0.3082643 0.3080926 0.3073452 0.3071650
paste(round(h.metnash.aldex$effect[1:decile], digits=3),collapse="\n")
# [1] "0.981\n0.959\n0.796\n0.645\n0.637\n0.623\n0.572\n0.567\n0.554\n0.525\n0.518\n0.513\n0.51\n0.506\n0.503\n0.485\n0.484\n0.455\n0.445\n0.445\n0.422\n0.422\n0.415\n0.408\n0.405\n0.397\n0.393\n0.393\n0.383\n0.381\n0.378\n0.375\n0.374\n0.374\n0.367\n0.366\n0.364\n0.359\n0.357\n0.352\n0.35\n0.337\n0.336\n0.336\n0.331\n0.33\n0.323\n0.319\n0.316\n0.308\n0.308\n0.307\n0.307"

h.metnash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)]
# [1] -0.3498652 -0.3533812 -0.3575882 -0.3578354 -0.3594253 -0.3605675
# [7] -0.3608379 -0.3740034 -0.3753918 -0.3803591 -0.3827707 -0.3836090
# [13] -0.3897053 -0.3952533 -0.4023195 -0.4041430 -0.4078187 -0.4115223
# [19] -0.4145232 -0.4214440 -0.4221561 -0.4244282 -0.4287520 -0.4357271
# [25] -0.4436833 -0.4718773 -0.4748695 -0.4924988 -0.4933874 -0.4935576
# [31] -0.5005012 -0.5043554 -0.5176838 -0.5243671 -0.5308861 -0.5309035
# [37] -0.5316331 -0.5316476 -0.5510731 -0.5523870 -0.5682283 -0.5768735
# [43] -0.6109385 -0.6178225 -0.6214374 -0.6269958 -0.6533644 -0.6660126
# [49] -0.6940952 -0.7325584 -0.7326152 -0.7328583 -0.9619965
paste(round(h.metnash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n")
# [1] "-0.35\n-0.353\n-0.358\n-0.358\n-0.359\n-0.361\n-0.361\n-0.374\n-0.375\n-0.38\n-0.383\n-0.384\n-0.39\n-0.395\n-0.402\n-0.404\n-0.408\n-0.412\n-0.415\n-0.421\n-0.422\n-0.424\n-0.429\n-0.436\n-0.444\n-0.472\n-0.475\n-0.492\n-0.493\n-0.494\n-0.501\n-0.504\n-0.518\n-0.524\n-0.531\n-0.531\n-0.532\n-0.532\n-0.551\n-0.552\n-0.568\n-0.577\n-0.611\n-0.618\n-0.621\n-0.627\n-0.653\n-0.666\n-0.694\n-0.733\n-0.733\n-0.733\n-0.962"

h.ss.aldex$effect[1:decile]
# [1]  0.110726245  0.495565275  0.193620053  0.406084164  0.236220107
# [6]  0.615120077  0.194071647  0.734783310  0.169281601  0.345123624
# [11]  0.209752541  0.425601285 -0.434988408  0.348590432  0.438079297
# [16]  0.589146581  0.441429366  0.436995379  0.211683902  0.922213342
# [21]  0.468890784  0.600150283  0.292888788  0.281431090  0.380489062
# [26]  0.138884184  0.227708983  0.313050674  0.396890197  0.133311047
# [31]  0.137016254  0.162182142 -0.160722709  0.178335182  0.358051943
# [36]  0.167190302  0.279849135 -0.011763482  0.271448972  0.166680871
# [41] -0.018004158  0.188882475  0.137704102  0.131628494  0.066220707
# [46]  0.317270941  0.027808416 -0.262842386  0.119506191  0.045841738
# [51] -0.006272643  0.383008730  0.375762786
paste(round(h.ss.aldex$effect[1:decile], digits=3),collapse="\n")
# [1] "0.111\n0.496\n0.194\n0.406\n0.236\n0.615\n0.194\n0.735\n0.169\n0.345\n0.21\n0.426\n-0.435\n0.349\n0.438\n0.589\n0.441\n0.437\n0.212\n0.922\n0.469\n0.6\n0.293\n0.281\n0.38\n0.139\n0.228\n0.313\n0.397\n0.133\n0.137\n0.162\n-0.161\n0.178\n0.358\n0.167\n0.28\n-0.012\n0.271\n0.167\n-0.018\n0.189\n0.138\n0.132\n0.066\n0.317\n0.028\n-0.263\n0.12\n0.046\n-0.006\n0.383\n0.376"

h.ss.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)]
# [1] -0.17687181 -0.11538454 -0.50496981 -0.52064715 -0.26654065 -0.27270409
# [7] -0.47562424 -0.46370905 -0.49592026  0.04375069 -0.03919277 -0.21113176
# [13] -0.17641491 -0.57831395 -0.56710033 -0.69337220  0.03814349 -0.10907167
# [19] -0.03689019 -0.61970393 -0.34832434 -0.29235083 -0.28737849 -0.43225579
# [25] -1.02362334 -0.62066910 -0.34211771  0.08265796 -0.70755499 -0.53477480
# [31] -0.59689793 -0.59294705 -0.30720780 -0.47776984 -0.38709346 -0.55585706
# [37] -0.39506415 -0.77728673 -0.68099018 -0.59832250 -0.87759943 -0.06819383
# [43] -0.53262630 -0.48407214 -0.53586024 -0.54577679 -0.72637908 -0.55450177
# [49] -0.55245968 -0.59128075 -0.60170345 -0.95828023 -0.51920813
paste(round(h.ss.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n")
# [1] "-0.177\n-0.115\n-0.505\n-0.521\n-0.267\n-0.273\n-0.476\n-0.464\n-0.496\n0.044\n-0.039\n-0.211\n-0.176\n-0.578\n-0.567\n-0.693\n0.038\n-0.109\n-0.037\n-0.62\n-0.348\n-0.292\n-0.287\n-0.432\n-1.024\n-0.621\n-0.342\n0.083\n-0.708\n-0.535\n-0.597\n-0.593\n-0.307\n-0.478\n-0.387\n-0.556\n-0.395\n-0.777\n-0.681\n-0.598\n-0.878\n-0.068\n-0.533\n-0.484\n-0.536\n-0.546\n-0.726\n-0.555\n-0.552\n-0.591\n-0.602\n-0.958\n-0.519"

h.nash.aldex$effect[1:decile]
# [1]  0.39715289  0.69005252  0.54514417  0.21596893  0.34353115  0.34035126
# [7]  0.82299962  0.77456805  0.30976597  0.30331341  0.33639054  0.53548576
# [13]  0.22923001  0.39639806  0.32276564  0.45968732  0.34067365  0.20973358
# [19]  0.26658249  0.41722621  0.19824482  0.22107496  0.33844001  0.30738894
# [25]  0.17000819  0.17780029  0.39237596  0.24529739  0.27126359  0.32166406
# [31]  0.14139397  0.25992817  0.35355448  0.11344520  0.42809363  0.16390138
# [37]  0.15534075  0.22453853  0.28213731  0.25907459  0.15608735  0.07416651
# [43]  0.31248701  0.27414291  0.30331671 -0.04140424 -0.06004984  0.04453913
# [49]  0.18934530  0.17186440  0.22404918  0.13673769  0.44886034
paste(round(h.nash.aldex$effect[1:decile], digits=3),collapse="\n")
# [1] "0.397\n0.69\n0.545\n0.216\n0.344\n0.34\n0.823\n0.775\n0.31\n0.303\n0.336\n0.535\n0.229\n0.396\n0.323\n0.46\n0.341\n0.21\n0.267\n0.417\n0.198\n0.221\n0.338\n0.307\n0.17\n0.178\n0.392\n0.245\n0.271\n0.322\n0.141\n0.26\n0.354\n0.113\n0.428\n0.164\n0.155\n0.225\n0.282\n0.259\n0.156\n0.074\n0.312\n0.274\n0.303\n-0.041\n-0.06\n0.045\n0.189\n0.172\n0.224\n0.137\n0.449"

h.nash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)]
# [1] -0.4569313 -0.1980081 -0.4449075 -0.3261535 -0.3911496 -0.3003279
# [7] -0.2564907 -0.2263294 -0.2383183 -0.1483125 -0.4035740 -0.2761080
# [13] -0.4030956 -0.4492730 -0.5955755 -0.4808998 -0.4690371 -0.2692206
# [19] -0.1144343 -0.5290339 -0.3562260 -0.4004923 -0.3031141 -0.2332924
# [25] -0.6510059 -0.1717429 -0.5700789 -0.4448545 -0.4996060 -0.4107511
# [31] -0.5355268 -0.4999002 -0.4183439 -0.2592910 -0.3425145 -0.4956334
# [37] -0.3920621 -0.6284506 -0.7418515 -0.6391082 -0.6100064 -0.2127948
# [43] -0.4217157 -0.6081564 -0.6131261 -0.6400437 -0.5638736 -0.6366474
# [49] -0.6671395 -0.5959037 -0.7744150 -0.5819000 -0.5832378
paste(round(h.nash.aldex$effect[(nrow(h.metnash.aldex)-decile+1):nrow(h.metnash.aldex)], digits=3),collapse="\n")
# [1] "-0.457\n-0.198\n-0.445\n-0.326\n-0.391\n-0.3\n-0.256\n-0.226\n-0.238\n-0.148\n-0.404\n-0.276\n-0.403\n-0.449\n-0.596\n-0.481\n-0.469\n-0.269\n-0.114\n-0.529\n-0.356\n-0.4\n-0.303\n-0.233\n-0.651\n-0.172\n-0.57\n-0.445\n-0.5\n-0.411\n-0.536\n-0.5\n-0.418\n-0.259\n-0.343\n-0.496\n-0.392\n-0.628\n-0.742\n-0.639\n-0.61\n-0.213\n-0.422\n-0.608\n-0.613\n-0.64\n-0.564\n-0.637\n-0.667\n-0.596\n-0.774\n-0.582\n-0.583"

# BARPLOT

# ALDEX BY GENUS TO COMPARE WITH QPCR
