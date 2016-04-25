d <- read.table("td_OTU_tag_mapped_lineage.txt",sep="\t", header=TRUE,quote="",row.names=1)

#remove taxonomy column to make otu count matrix numeric
taxonomy <- d$taxonomy
d.mat <- d[-length(colnames(d))]
d.mat <- t(as.matrix(d.mat))

rownames(d.mat)[which(rownames(d.mat) ==  "HLD.102bba")] <- "HLD.102b"
rownames(d.mat)[which(rownames(d.mat) ==  "Cl.141.6mo.NASH_FO.C")] <- "CL.141.6mo.NASH_FO.C"

rownames(d.mat) <- gsub("NASH_FO.","",rownames(d.mat))
rownames(d.mat) <- gsub("NASH.","",rownames(d.mat))
rownames(d.mat) <- gsub("Healthy.","",rownames(d.mat))
rownames(d.mat) <- gsub(".C","c",rownames(d.mat))
rownames(d.mat) <- gsub("SS[^.]*.","",rownames(d.mat))

d.mat <- d.mat[grepl(".*[ab]$",rownames(d.mat)),]

samples <- substr(rownames(d.mat),1,nchar(rownames(d.mat))-1)
samples <- list(samples)
d.df <- as.data.frame(d.mat)
d.agg <- aggregate(d.df,samples,sum)
rownames(d.agg) <- d.agg$Group.1
d.agg <- d.agg[,c(2:ncol(d.agg))]

MyMeta<- read.table("../../../metadoot.tsv", header=T, sep="\t", row.names=1, comment.char="", check.names=FALSE)
metadata.baseline <- MyMeta[which(MyMeta$stage == "baseline"),]
metadata.baseline <- metadata.baseline[which(metadata.baseline$exclude.for.IM == 1),]

rownames(metadata.baseline) <- gsub("-",".",rownames(metadata.baseline))
rownames(metadata.baseline) <- gsub("\\(",".",rownames(metadata.baseline))
rownames(metadata.baseline) <- gsub("\\)","",rownames(metadata.baseline))

rownames(metadata.baseline)[which(rownames(metadata.baseline) == "CL.151.BL.2")] <- "CL.151.BL.R2"
rownames(metadata.baseline)[which(rownames(metadata.baseline) == "CL.141.BL.2")] <- "CL.141.BL.R2"
rownames(metadata.baseline)[which(rownames(metadata.baseline) == "CL.166.BL.2")] <- "CL.166.BL"

d.agg <- d.agg[which(rownames(d.agg) %in% rownames(metadata.baseline)),]

# sanity check
print(paste("taxonomy still in order:", all.equal(colnames(d.agg),rownames(d))))

d.baseline <- data.frame(t(d.agg))
colnames(d.baseline) <- gsub("\\.","-",colnames(d.baseline))

write.table(d.baseline,file="summed_data_baseline_only.txt",sep="\t",quote=FALSE)

d.baseline$taxonomy <- taxonomy

write.table(d.baseline,file="summed_data_baseline_only_with_taxonomy.txt",sep="\t",quote=FALSE)
