#' minHGT permutation test
#'
#' This function takes a list of indicator values (0s and 1s) ranked  by some orthogonal measure,
#' for example if you sorted people by height and then went down the row asking if each individual
#' was male (1) or female (0)
#'
#' This permutation test determines if the successes in that ranked list are predisposed beyond
#' random chance (as in the case with height and sex, which are associated) to be skewed to the front of the list 
#'
#' More specifically, this test identifies the minimum hypergeometric test (minHGT) p-value for a given vector as described here:
#' https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-10-48
#' then performs a specified number of permutations of the list to determine how frequently an equal or lesser minHGT p-value can be obtained randomly.
#' the reported p-value from this script, is then simply just n/N, where n is the number of randomly obtained minHGT scores which were equal to or less
#' than the value obtained by the given list, then divided by the total number of permutations attempted.
#'
#' @param thisVector A vector of individual Bernoulli trials (https://en.wikipedia.org/wiki/Bernoulli_trial). Each value must be a 0 or 1. 
#' @param numPerms The number of permutations to run (defaults to 1000)
#' @param verbose Determines whether to print the number of permutation the run is on (default is "False", prints every 100 permutations)
#' @export
#' @examples
#' v1 <- round(runif(100,0,.7),0)
#' v2 <- round(runif(100,0.3,1),0)
#' testVector <- c(v2,v1)
#' testOut <-  minHGT(testVector,1000,verbose=1)
#' testOut$pval  # the probability of observing a given minimum mypergeometric value given random shuffling of the vector
#' testOut$crit  # the set size which yeilds the greatest minHGT for the given vector

minHGT <- function(thisVector, numPerms=1000, verbose=1) {
	thisVector[which(thisVector>0)] <- 1
	
	getMinHGT <- function(thisVector) {
		thisM <- sum(thisVector==1)
		thisN <- sum(thisVector==0)
		
		endLimit <- 100
		currentMin <- endLimit
		currentMax <- length(thisVector)-endLimit
		oMag <- log10(length(thisVector[currentMin:currentMax]))		
		
		numLoops = 0
		thisMinHGT = 1
		criticalPoint = 1
		
		while (oMag>0.6 & numLoops < oMag) {
	
			
			cutoffs <- currentMin + unique(round(10^ (c(1:100) * oMag/100),0))
			cutoffs <- cutoffs[which(cutoffs > endLimit & cutoffs < currentMax)]
			
			for (i in 1:length(cutoffs)) {
				thisQ = sum(thisVector[1:cutoffs[i]]==1)
				thisK = cutoffs[i]
				tmpPval = phyper(thisQ, thisM, thisN, thisK, lower.tail=F)
				if (tmpPval <= thisMinHGT) {
					thisMinHGT    <-  tmpPval
					criticalPoint <- cutoffs[i]
					
					if (i > 1) { 
						currentMin <- cutoffs[i-1]
					} else { 
						currentMin <- endLimit 
					}
					if (i < length(cutoffs)) { 
						currentMax <- cutoffs[i+1] 
					} else { 
						currentMax <- length(thisVector)-endLimit
					}
					oMag <- log10(length(thisVector[currentMin:currentMax]))
				}
			}
			
			numLoops <- numLoops + 1
		}
		op <- c()
		op$minHGT <- thisMinHGT
		op$crit   <- criticalPoint
		return(op)
	}
	
	thisResult <- getMinHGT(thisVector)
	op <- c()
	op$crit <- thisResult$crit
	
	numLesser = 0
	for (i in 1:numPerms) {
		if (verbose & i %% 100 == 0) { print (paste("perm",i)) }
		tmpMinHGT <- getMinHGT(sample(thisVector))
		if (tmpMinHGT$minHGT <= thisResult$minHGT) {numLesser <- numLesser + 1}
	}
	op$pval <- numLesser/numPerms
	return(op)
}
