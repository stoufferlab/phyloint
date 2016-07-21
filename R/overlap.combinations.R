#' @title Calculate the overlap based on combinations of traits
#' @param adj The adjacency matrix
#' @param omat The overlap matrix for the different traits
#' @references A Eklöf, U Jacob, J Kopp, J Bosch, R Castro-Urgal, NP Chacoff,
#'    B Dalsgaard, C de Sassi, M Galetti, PR Guimarães, SB Lomáscolo, AM Martín
#'    González, MA Pizo, R Rader, A Rodrigo, JM Tylianakis, DP Vázquez, and S
#'    Allesina (2013) The dimensionality of ecological networks. Ecology
#'    Letters 16:577-583 (doi: 10.1111/ele.12081)
#' @return The function returns a table with the overlap produced for all possible combinations
#' of the traits analyzed.
#' @export
#' 
#' @examples
#' combinations <- overlap.combinations(eklof$network, overlap.matrices)

overlap.combinations <- function(adj, omat){
# number of traits
numTraits <- dim(omat)[3]
traitNames <- dimnames(omat)[3][[1]]
# number of combinations of traits
numberComb <- numTraits*(94+numTraits*(5+numTraits*(25+numTraits*(-5+numTraits))))/120

# build an empty result table
resultTable <- matrix(NA,numberComb,2)
labels <- rep("aa",numberComb)
current <- 1
  for (k in 1:numTraits){
    # compute all combinations of numTraits traits taken k at a time
    combinations <- combn(numTraits,k)
    # for each combination
    for(j in 1:dim(combinations)[2]){
      # start with the matrix corresponding to the first trait
      finalMat <- omat[,,combinations[1,j]]
      # and multiply for the matrices corresponding to the other traits
      if (k>1){
        for (z in 2:k){
          finalMat <- finalMat*omat[,,combinations[z,j]]
        }
      }
      # save the results in the table
      # first, save number of traits used in this combination
      resultTable[current,1] <- k
      # second, compute overlap
      resultTable[current,2] <- sum(adj)/sum(finalMat)
      
      labels[current] <- paste(traitNames[combinations[,j]],collapse="+")
      current<-current+1
    }
  }

rownames(resultTable) <- labels
colnames(resultTable) <- c("NumTraits", "Overlap")
print(resultTable)
return(resultTable)
}