options(error=recover)

#GETTING THE DATA FROM Human Microbiome Project OTU table
library(phangorn)

hmpData <- "../fodor/otu_table_psn_v35.txt"
hmpMetadata <- "../fodor/v35_map_uniquebyPSN.txt"
treeFile <- "../fodor/rep_set_v35.tre"

id <- read.table(hmpMetadata, header=TRUE, sep="\t", row.names=1)
otu <- t( read.table(hmpData, header=T, sep="\t", row.names=1, check.names=FALSE) )
tree <- read.tree(treeFile)

# BODY SITES
#  [1] "Anterior_nares"               "Attached_Keratinized_gingiva"
#  [3] "Buccal_mucosa"                "Hard_palate"                 
#  [5] "Left_Antecubital_fossa"       "Left_Retroauricular_crease"  
#  [7] "Mid_vagina"                   "Palatine_Tonsils"            
#  [9] "Posterior_fornix"             "Right_Antecubital_fossa"     
# [11] "Right_Retroauricular_crease"  "Saliva"                      
# [13] "Stool"                        "Subgingival_plaque"          
# [15] "Supragingival_plaque"         "Throat"                      
# [17] "Tongue_dorsum"                "Vaginal_introitus"  

bm <- "Buccal_mucosa" # cheek
td <- "Tongue_dorsum" # tongue

bm.id <- rownames(id)[which(id$HMPbodysubsite==bm)]
td.id <- rownames(id)[which(id$HMPbodysubsite==td)]

bm.otu <- site <- otu[rownames(otu) %in% bm.id,]
td.otu <- site <- otu[rownames(otu) %in% td.id,]

otuIDs <- colnames(bm.otu)

bm.otu <- apply(bm.otu, 1, function(x){as.numeric(x)})
td.otu <- apply(td.otu, 1, function(x){as.numeric(x)})

rownames(bm.otu) <- otuIDs
rownames(td.otu) <- otuIDs

# remove all samples with read count lower than 10,000
bm.sum <- apply(bm.otu,2,sum)
td.sum <- apply(td.otu,2,sum)

bm.otu <- bm.otu[,which(bm.sum>5000 & bm.sum<10000)]
td.otu <- td.otu[,which(td.sum>5000 & td.sum<10000)]

#pick random samples from each category
nSamples <- 30

bm.rand <- bm.otu[,as.integer(sample(seq(1,length(colnames(bm.otu)),1),nSamples,replace=FALSE))]
td.rand <- td.otu[,as.integer(sample(seq(1,length(colnames(td.otu)),1),nSamples,replace=FALSE))]

#concatenate

data <- data.frame(bm.rand,td.rand)
colnames(data) <- sub("^X", "", colnames(data))

#make condition vector

groups <- as.factor(c(rep("bm",nSamples),rep("td",nSamples)))

#get rid of zero sum rows

data.sum <- apply(data,1,sum)
data.0 <- data[data.sum > 0,]
data <- data.0

# get rid of extra OTUs in tree
tree$tip.label <- gsub("'","",tree$tip.label)
absent <- tree$tip.label[!(tree$tip.label %in% rownames(data))]
if (length(absent) != 0) {
		tree <- drop.tip(tree, absent)
}
write.tree(tree,file="data/tongue_cheek_data/hmp_tongue_cheek_subtree.tre")

#make sample names that contain condition

colnames(data) <- paste(groups,colnames(data),sep="_")

#write otu counts into table
write.table(data,file="data/tongue_cheek_data/hmp_tongue_cheek_data.txt",sep="\t",quote=FALSE)
# read in with read.table("hmp_mouth_data.txt",sep="\t",header=TRUE,row.names=1)