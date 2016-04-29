taxonomy <- read.table("taxonomy_map.txt",header=FALSE,quote="",row.names=1,comment.char="",stringsAsFactors=FALSE)

ss <- read.table("H_vs_SS_aldex_output_OTU_level_128_MC_samples.txt",header=TRUE,row.names=1,quote="",comment.char="")
nash <- read.table("H_vs_NASH_aldex_output_OTU_level_128_MC_samples.txt",header=TRUE,row.names=1,quote="",comment.char="")
metnash <- read.table("H_vs_extreme_NASH_aldex_output_OTU_level_128_MC_samples.txt",header=TRUE,row.names=1,quote="",comment.char="")

metnash <- metnash[order(metnash$effect,decreasing=TRUE),]

decile <- round(nrow(metnash)/10)

top <- data.frame(matrix(nrow=decile,ncol=6))
bottom <- data.frame(matrix(nrow=decile,ncol=6))

colnames(top) <- c("OTU family","OTU genus","SILVA bootstrap value","H vs. SS effect","H vs. NASH effect","H vs. extreme NASH effect")
colnames(bottom) <- colnames(top)

rownames(top) <- rownames(metnash)[c(1:decile)]
rownames(bottom) <- rownames(metnash)[c((nrow(metnash)-decile+1):nrow(metnash))]

top.taxa <- taxonomy[match(rownames(top),rownames(taxonomy)),]
bottom.taxa <- taxonomy[match(rownames(bottom),rownames(taxonomy)),]

top.genus <- c(as.character(top.taxa))
top.family <- c(as.character(top.taxa))
top.bootstrap <- c(as.character(top.taxa))
bottom.genus <- c(as.character(bottom.taxa))
bottom.family <- c(as.character(bottom.taxa))
bottom.bootstrap <- c(as.character(bottom.taxa))

for (i in c(1:decile)) {
  top.genus[i] <- strsplit(top.genus[i],c(";"))[[1]][6]
	top.family[i] <- strsplit(top.family[i],c(";"))[[1]][5]
	top.bootstrap[i] <- strsplit(top.bootstrap[i],c(";"))[[1]][7]
  top.bootstrap[i] <- gsub("^.","",top.bootstrap[i])
  bottom.genus[i] <- strsplit(bottom.genus[i],c(";"))[[1]][6]
  bottom.family[i] <- strsplit(bottom.family[i],c(";"))[[1]][5]
  bottom.bootstrap[i] <- strsplit(bottom.bootstrap[i],c(";"))[[1]][7]
  bottom.bootstrap[i] <- gsub("^.","",bottom.bootstrap[i])
}

top[,1] <- top.family
top[,2] <- top.genus
top[,3] <- top.bootstrap

bottom[,1] <- bottom.family
bottom[,2] <- bottom.genus
bottom[,3] <- bottom.bootstrap

top[,4] <- ss$effect[match(rownames(top),rownames(ss))]
top[,5] <- nash$effect[match(rownames(top),rownames(nash))]
top[,6] <- metnash$effect[match(rownames(top),rownames(metnash))]

bottom[,4] <- ss$effect[match(rownames(bottom),rownames(ss))]
bottom[,5] <- nash$effect[match(rownames(bottom),rownames(nash))]
bottom[,6] <- metnash$effect[match(rownames(bottom),rownames(metnash))]

top.round <- top
top.round[,c(4:6)] <- round(top.round[,c(4:6)],digits=3)
bottom.round <- bottom
bottom.round[,c(4:6)] <- round(bottom.round[,c(4:6)],digits=3)

top.round[,4] <- paste("$",as.character(top.round[,4]),"$",sep="")
top.round[,5] <- paste("$",as.character(top.round[,5]),"$",sep="")
top.round[,6] <- paste("$",as.character(top.round[,6]),"$",sep="")
bottom.round[,4] <- paste("$",as.character(bottom.round[,4]),"$",sep="")
bottom.round[,5] <- paste("$",as.character(bottom.round[,5]),"$",sep="")
bottom.round[,6] <- paste("$",as.character(bottom.round[,6]),"$",sep="")

write.table(top.round,file="top_decile.txt",row.names=FALSE,quote=FALSE,col.names=TRUE,sep="\t")
write.table(bottom.round,file="bottom_decile.txt",row.names=FALSE,quote=FALSE,col.names=TRUE,sep="\t")

cutoff <- 0.4

top.ss <- which(top[,4] >= cutoff)
top.nash <- which(top[,5] >= cutoff)
top.metnash <- which(top[,6] >= cutoff)
top.consistent <- top[intersect(intersect(top.ss,top.nash),top.metnash),]

cutoff <- -1*cutoff

bottom.ss <- which(bottom[,4] <= cutoff)
bottom.nash <- which(bottom[,5] <= cutoff)
bottom.metnash <- which(bottom[,6] <= cutoff)
bottom.consistent <- bottom[intersect(intersect(bottom.ss,bottom.nash),bottom.metnash),]

top.consistent.round <- top.round[intersect(intersect(top.ss,top.nash),top.metnash),]
bottom.consistent.round <- bottom.round[intersect(intersect(bottom.ss,bottom.nash),bottom.metnash),]

write.table(top.consistent.round,file="top_decile_consistent.txt",row.names=FALSE,quote=FALSE,col.names=TRUE,sep="\t")
write.table(bottom.consistent.round,file="bottom_decile_consistent.txt",row.names=FALSE,quote=FALSE,col.names=TRUE,sep="\t")

library(ALDEx2)

pdf("ALDEx_OTU_level_plots_top_consistent.pdf")

aldex.plot(ss,type="MW",xlab="Dispersion within groups", ylab="Difference between groups")
title("Healthy vs. SS")
points(ss$diff.win[match(rownames(top.consistent),rownames(ss))],ss$diff.btw[match(rownames(top.consistent),rownames(ss))],pch=1,col="red")
points(ss$diff.win[match(rownames(bottom.consistent),rownames(ss))],ss$diff.btw[match(rownames(bottom.consistent),rownames(ss))],pch=1,col="black")

aldex.plot(nash,type="MW",xlab="Dispersion within groups", ylab="Difference between groups")
title("Healthy vs. NASH")
points(nash$diff.win[match(rownames(top.consistent),rownames(nash))],nash$diff.btw[match(rownames(top.consistent),rownames(nash))],pch=1,col="red")
points(nash$diff.win[match(rownames(bottom.consistent),rownames(nash))],nash$diff.btw[match(rownames(bottom.consistent),rownames(nash))],pch=1,col="black")

aldex.plot(metnash,type="MW",xlab="Dispersion within groups", ylab="Difference between groups")
title("Healthy vs. extreme NASH")
points(metnash$diff.win[match(rownames(top.consistent),rownames(metnash))],metnash$diff.btw[match(rownames(top.consistent),rownames(metnash))],pch=1,col="red")
points(metnash$diff.win[match(rownames(bottom.consistent),rownames(metnash))],metnash$diff.btw[match(rownames(bottom.consistent),rownames(metnash))],pch=1,col="black")

dev.off()
