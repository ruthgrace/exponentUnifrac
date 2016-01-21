# ruth's unifrac method


library(phangorn)
library(ape)
library(zCompositions)

gm_mean = function(x, na.rm=TRUE){
  exp(mean(log(x), na.rm=na.rm) )
}

gm_vector = function(node_count, other_counts, geometric_mean) {
  gm_vec <- array(0,length(node_count))
  # a count of -1 means that the OTU has been amalgamated to other OTUs in the node_count
  zero_means <- list()
  counter = 1
  for (i in 1:length(node_count)) {
  # CHECK IF THIS IS CORRECT - indexing for other_counts
    gm_vec[i] <- gm_mean(c(node_count[i],other_counts[i,][which(other_counts[i,]!=-1)]))
    if (gm_vec[i] == 0) {
      gm_vec[i] = geometric_mean[i]
    }
  }
  return(gm_vec)
}

#valid methods are unweighted, weighted, information. Any other method will result in a warning and the unweighted analysis
#pruneTree option prunes the tree for each comparison to exclude branch lengths not present in both samples
#normalize divides the value at each node by sum of weights to guarantee output between 0 and 1 (breaks the triangle inequality)

#otuTable must have samples as rows, OTUs as columns
#tree must be phylo tree object from ape package (can use read.tree method to read in a newick tree to this format)

getDistanceMatrix <- function(otuTable,tree,method="weighted",verbose=FALSE,pruneTree=FALSE,normalize=TRUE)  {

	if (length(which(is.na(otuTable))) > 0) {
		stop("OTU count table has NA")
	}

	if (!is.rooted(tree)) {
		tree <- midpoint(tree)
		if(verbose) { print("Rooting tree by midpoint") }
	}

	if (attributes(tree)$order!="postorder") {
		tree <- reorder(tree,order="postorder")
		if (verbose) { print("Reordering tree as postorder for distance calculation algorithm") }
	}

	# get proportions
	readsPerSample <- apply(otuTable,1,sum)
	otu.prop <- otuTable/readsPerSample
	otu.prop <- as.matrix(otu.prop)
	rownames(otu.prop) <- rownames(otuTable)
	colnames(otu.prop) <- colnames(otuTable)
	if(verbose) {	print("calculated proportional abundance")	}

	# add priors to zeros based on bayesian approach
	otuTable.adjustedZeros <- cmultRepl(otuTable, method="CZM", output="counts")
	# calculate geometric mean & geometric sum for exponent weighted UniFrac
	geometric_mean <- apply(otuTable.adjustedZeros, 1, gm_mean)
	geometric_sum <- apply(otuTable.adjustedZeros, 1, sum) / geometric_mean
	if(verbose) {	print("calculated geometric mean per sample")	}

	##get cumulative proportional abundance for the nodes (nodes are ordered same as in the phylo tree representation)

	#cumulative proportional abundance stored in weights
	weights <- matrix(NA,ncol=(length(tree$edge.length) + 1),nrow=nrow(otuTable))
	#cumulative abundance
	absolute_weights <- matrix(NA,ncol=(length(tree$edge.length) + 1),nrow=nrow(otuTable))
  geometric_means <- matrix(NA,ncol=(length(tree$edge.length) + 1),nrow=nrow(otuTable))
  otuCountsPerNode <- array(0, c(length(rownames(otuTable)), length(tree$edge.length)+1, length(colnames(otuTable))));
	#each row is a sample
	rownames(weights) <- rownames(otuTable)
	rownames(absolute_weights) <- rownames(otuTable)
	rownames(geometric_means) <- rownames(otuTable)
  dimnames(otuCountsPerNode)[[1]] <- rownames(otuTable)
	#each column is the abundance weighting for a node in the phylogenetic tree
	colnames(weights) <- c(1:(length(tree$edge.length)+1))
	colnames(absolute_weights) <- c(1:(length(tree$edge.length)+1))
	colnames(geometric_means) <- c(1:(length(tree$edge.length)+1))
  dimnames(otuCountsPerNode)[[2]] <- c(1:(length(tree$edge.length)+1))
  dimnames(otuCountsPerNode)[[3]] <- colnames(otuTable)

	treeLeaves <- c(1:length(tree$tip.label))

	#loop through edges
	#if child node of edge (of which there is only one) has not been seen before, it is a leaf
	#	(this is a property of a postorder tree -- children are listed before parents)
	if(verbose) {	print("calculating weights...") }
	for (i in c(1:nrow(tree$edge))) {
		parentNode <- tree$edge[i,1]
		childNode <- tree$edge[i,2]

		#if node is all NA, node has not been seen before, and node is a leaf (ie. an OTU)
		if (is.na(weights[1,childNode])) {
			#put OTU abundance in weights
			otuName <- tree$tip.label[childNode]
			otuIndex <- which(colnames(otu.prop) == otuName)[1]
			weights[,childNode] <- otu.prop[,otuIndex]
			absolute_weights[,childNode] <- otuTable[,otuIndex]
      # MAKE CHILD ZERO IN OTU COUNTS
      otuCountsPerNode[,childNode,] <- as.matrix(otuTable.adjustedZeros)
      otuCountsPerNode[,childNode,otuIndex] <- -1
      geometric_means[,childNode] <- geometric_mean
		}

		if (is.na(weights[1,parentNode])) {
			# initialize parentNode with counts of zero
			weights[,parentNode] <- 0
			absolute_weights[,parentNode] <- 0
      otuCountsPerNode[,parentNode,] <- otuCountsPerNode[,childNode,]
		} else {
        otuCountsPerNode[,parentNode,][which(otuCountsPerNode[,childNode,]==-1)] <- -1
    }

		#add child node abundance to parent node abundance
		weights[,parentNode] <- weights[,parentNode] + weights[,childNode]
		absolute_weights[,parentNode] <- absolute_weights[,parentNode] + absolute_weights[,childNode]
    geometric_means[,parentNode] <- gm_vector(absolute_weights[,parentNode],otuCountsPerNode[,parentNode,],geometric_mean)
	}

	if(verbose) {	print("done calculating weights")	}

	if (method=="information") {
		if(verbose) {	print("information entropy transform")	}
		#information entropy transform
		weights[] <- (-1) * weights[] * log2(weights[])
		weights <- as.matrix(weights)
		weights[which(is.na(weights))] <- 0
	}

	if (method=="exponent") {
		if(verbose) {	print("CLR exponent transform")	}
		weights[] <- (absolute_weights[] / geometric_means[])
		weights <- as.matrix(weights)
		weights[which(is.na(weights))] <- 0
	}

	nSamples <- length(rownames(otuTable))
	distanceMatrix <- data.frame(matrix(ncol=nSamples,nrow=nSamples))
	rownames(distanceMatrix) <- rownames(otuTable)
	colnames(distanceMatrix) <- rownames(otuTable)
	branchLengths <- tree$edge.length

	branchLengths <- branchLengths[order(tree$edge[,2])]
	weights <- weights[,which(!is.na(match(colnames(weights),tree$edge[,2])))]


	if(verbose) {	print("calculating pairwise distances...")	}

	for (i in 1:nSamples) {
		for (j in i:nSamples) {

				if (method=="weighted" || method=="information" || method=="exponent") {
					# the formula is sum of (proportional branch lengths * | proportional abundance for sample A - proportional abundance for sample B| )
					if (pruneTree==TRUE){
						includeBranchLengths <- which( (weights[i,] > 0) | (weights[j,] > 0) )
						if (normalize==TRUE && method!="exponent") {
							distance <- sum( branchLengths[includeBranchLengths] * abs(weights[i,includeBranchLengths] - weights[j,includeBranchLengths]) )/sum( branchLengths[includeBranchLengths]* (weights[i,includeBranchLengths] + weights[j,includeBranchLengths]) )
						}
						else {
							distance <- sum( branchLengths[includeBranchLengths] * abs(weights[i,includeBranchLengths] - weights[j,includeBranchLengths]) )/sum( branchLengths[includeBranchLengths])
						}
					}
					else {
						distance <- sum( branchLengths * abs(weights[i,] - weights[j,]) )/sum(branchLengths)
						if (normalize==TRUE) {
							distance <- sum( branchLengths * abs(weights[i,] - weights[j,]) )/sum(branchLengths * (weights[i,] + weights[j,]))
						}
					}
				}
			else {
				if (method!="unweighted") {
					warning(paste("Invalid method",method,", using unweighted Unifrac"))
				}
				# the formula is sum of (branch lengths * (1 if one sample has counts and not the other, 0 otherwise) )
				#	i call the (1 if one sample has counts and not the other, 0 otherwise) xorBranchLength
				xorBranchLength <- as.numeric(xor( weights[i,] > 0, weights[j,] > 0))
				if (pruneTree==TRUE) {
					includeBranchLengths <- which( (weights[i,] > 0) | (weights[j,] > 0) )
					distance <- sum( branchLengths[includeBranchLengths] *  xorBranchLength[includeBranchLengths])/sum(branchLengths[includeBranchLengths])
				}
				else {
					distance <- sum( branchLengths *  xorBranchLength)/sum(branchLengths)
				}

			}
			distanceMatrix[i,j] <- distance
			distanceMatrix[j,i] <- distance

		}
	}

	if(verbose) {	print("done")	}

	return(distanceMatrix)

}
