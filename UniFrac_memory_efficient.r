# ruth's unifrac method


library(phangorn)
library(ape)
library(zCompositions)

# global variables - used to save memory when being passed between functions
# otuPropsPerNode <- "global"
# otuPropsPerNode.adjustedZeros <- "global"
# weightsPerNode <- "global"
# unifrac.tree <- "global"
# unifrac.method <- "global"
# unifrac.treeLeaves <- "global"

gm_mean = function(x, na.rm=TRUE){
  exp(mean(log(x), na.rm=na.rm) )
}

# gm_vector = function(node_count, other_counts, geometric_mean) {
#   gm_vec <- array(0,length(node_count))
#   # a count of -1 means that the OTU has been amalgamated to other OTUs in the node_count
#   zero_means <- list()
#   counter = 1
#   for (i in 1:length(node_count)) {
#   # CHECK IF THIS IS CORRECT - indexing for other_counts
#     gm_vec[i] <- gm_mean(c(node_count[i],other_counts[i,][which(other_counts[i,]!=-1)]))
#     if (gm_vec[i] == 0) {
#       gm_vec[i] = geometric_mean[i]
#     }
#   }
#   return(gm_vec)
# }

gm_from_props = function(rownum, root) {
  return(log2(otuPropsPerNode.adjustedZeros[rownum,root]) - mean(log2(otuPropsPerNode.adjustedZeros[rownum,which(otuPropsPerNode.adjustedZeros[rownum,]>0)])))
}

build_weights = function(root) {
  print(paste("building weights for node ",root))
  # make all weight values positive
  otuPropsPerNode <- abs(otuPropsPerNode)
  otuPropsPerNode.adjustedZeros <- abs(otuPropsPerNode.adjustedZeros)
  # if root is leaf, then weight has already been calculated and we can skip this
  if (root > length(unifrac.tree$tip.label)) {
    children <- unifrac.tree$edge[which(unifrac.tree$edge[,1]==root),2]
    
    # construct values for left child
    build_weights(children[1])
    leftChildProportions <- otuPropsPerNode
        
    # construct values for right child
    build_weights(children[2])
    rightChildProportions <- otuPropsPerNode
    
    # the negative values in leftChildProportions are for all the nodes in the subtree of the left child
    # the negative values in rightChildProportions are for all the nodes in the subtree of the right child
    
    leafColumns <- which(colnames(otuPropsPerNode) %in% unifrac.treeLeaves)
    
    leftChildSubtree <- colnames(leftChildProportions)[which(leftChildProportions[1,] < 0)]
    rightChildSubtree <- colnames(rightChildProportions)[which(rightChildProportions[1,] < 0)]
    
    leafColumnsForWeightCalculation <- leafColumns
    leafColumnsForWeightCalculation <- which(!(leafColumnsForWeightCalculation %in% leftChildSubtree))
    leafColumnsForWeightCalculation <- which(!(leafColumnsForWeightCalculation %in% rightChildSubtree))
    
    #calculate proportion and weight of current node
    otuPropsPerNode[,root] <- otuPropsPerNode[,children[1]] + otuPropsPerNode[,children[2]]
    otuPropsPerNode[,children[1]] <- (-1)*abs(otuPropsPerNode[,children[1]])
    otuPropsPerNode[,children[2]] <- (-1)*abs(otuPropsPerNode[,children[2]])
    
    otuPropsPerNode.adjustedZeros[,root] <- otuPropsPerNode.adjustedZeros[,children[1]] + otuPropsPerNode.adjustedZeros[,children[2]]
    otuPropsPerNode.adjustedZeros[,children[1]] <- (-1)*abs(otuPropsPerNode.adjustedZeros[,children[1]])
    otuPropsPerNode.adjustedZeros[,children[2]] <- (-1)*abs(otuPropsPerNode.adjustedZeros[,children[2]])

    # TODO fix - the new weightsPerNode is all negative, and the CLR values aren't zeros
    
    weightsPerNode[,root] <- sapply(c(1:nrow(otuPropsPerNode.adjustedZeros)),function(x) { gm_from_props(x, root) })
  }
  # make root node values negative
  otuPropsPerNode[,root] <- (-1)*abs(otuPropsPerNode[,root])
  otuPropsPerNode.adjustedZeros[,root] <- (-1)*abs(otuPropsPerNode.adjustedZeros[,root]
}


#valid methods are unweighted, weighted, information. Any other method will result in a warning and the unweighted analysis
#pruneTree option prunes the tree for each comparison to exclude branch lenxgths not present in both samples
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
  
  #make globally available copy of tree
  unifrac.tree <<- tree
  
  #make globally available copy of method
  unifrac.method <<- method
  
  #make globally available copy of leaves
  unifrac.treeLeaves <<- c(1:length(tree$tip.label))

	# get proportions
	readsPerSample <- apply(otuTable,1,sum)
	otu.prop <- otuTable/readsPerSample
	otu.prop <- as.matrix(otu.prop)
	rownames(otu.prop) <- rownames(otuTable)
	colnames(otu.prop) <- colnames(otuTable)
	if(verbose) {	print("calculated proportional abundance")	}

	# add priors to zeros based on bayesian approach
	otuTable.adjustedZeros <- cmultRepl(otuTable, method="CZM", output="counts")
  readsPerSample <- apply(otuTable.adjustedZeros,1,sum)
  otu.prop.adjustedZeros <- otuTable.adjustedZeros/readsPerSample
  otu.prop.adjustedZeros <- as.matrix(otu.prop.adjustedZeros)
	rownames(otu.prop.adjustedZeros) <- rownames(otuTable)
	colnames(otu.prop.adjustedZeros) <- colnames(otuTable)
  
	# calculate geometric mean & geometric sum for exponent weighted UniFrac
	geometric_mean <- apply(otuTable.adjustedZeros, 1, gm_mean)
	geometric_sum <- apply(otuTable.adjustedZeros, 1, sum) / geometric_mean
	if(verbose) {	print("calculated geometric mean per sample")	}

	##get cumulative proportional abundance for the nodes (nodes are ordered same as in the phylo tree representation)

  otuPropsPerNode <<- matrix(NA, ncol=length(tree$edge.length)+1, nrow=length(colnames(otuTable)))
  otuPropsPerNode.adjustedZeros <<- matrix(NA, ncol=length(tree$edge.length)+1, nrow=length(colnames(otuTable)))
  weightsPerNode <<- matrix(NA, ncol=length(tree$edge.length)+1, nrow=length(colnames(otuTable)))
  
	#each row is a sample
  rownames(otuPropsPerNode) <- colnames(otuTable)
  rownames(otuPropsPerNode.adjustedZeros) <- colnames(otuTable)
  rownames(weightsPerNode) <- colnames(otuTable)
  #each column is a node in the tree
  colnames(otuPropsPerNode) <- c(1:(length(tree$edge.length)+1))
  colnames(otuPropsPerNode.adjustedZeros) <- c(1:(length(tree$edge.length)+1))
  colnames(weightsPerNode) <- c(1:(length(tree$edge.length)+1))
  
  leafNodes <- c(1:length(tree$tip.label))
  leafOrder <- match(tree$tip.label,rownames(otu.prop))
  
  otuPropsPerNode[,leafNodes] <- t(otu.prop[leafOrder,])
  otuPropsPerNode.adjustedZeros[,leafNodes] <- t(otu.prop.adjustedZeros[leafOrder,])
  weightsPerNode[,leafNodes] <- log2(t(otu.prop.adjustedZeros[leafOrder,]))
  weightsPerNode[,leafNodes] <- t(apply(weightsPerNode[,leafNodes],1,function(x) { return(x - mean(x))}))
  
  # the tree is in postorder, so the last edges belong to the root
  root <- tree$edge[nrow(tree$edge),1]

  if(verbose) {	print("calculating weights...")	}
  
  #TODO function doesn't do upper nodes :(
  
  build_weights(root)

  if(verbose) {	print("calculating pairwise distances...")	}
  
  nSamples <- nrow(otuTable)
  distanceMatrix <- matrix(NA,nrow=nSamples,ncol=nSamples)
  rownames(distanceMatrix) <- rownames(otuTable)
  colnames(distanceMatrix) <- rownames(otuTable)
  
  branchLengths <- tree$edge.length
  
  #convert table according to weight
  if (method=="information") {
    weights <- otuPropsPerNode*log2(otuPropsPerNode)
  } else if (method=="exponent") {
    weights <- weightsPerNode
  } else {
    weights <- otuPropsPerNode
  }
  
	for (i in 1:nSamples) {
		for (j in i:nSamples) {

				if (method=="weighted" || method=="information" || method=="exponent") {
					# the formula is sum of (proportional branch lengths * | proportional abundance for sample A - proportional abundance for sample B| )
					if (pruneTree==TRUE){
						includeBranchLengths <- which( (otuPropsPerNode[i,] > 0) | (otuPropsPerNode[j,] > 0) )
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
					warning(paste("Invalid method",method,", using unweighted Unifrac instead"))
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
