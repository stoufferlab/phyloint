#' @title Build a predicted adjacency matrix based on trait values
#' @param adj The adjacency matrix
#' @param corr.traits A list with species as rows, the values of the corrected traits as columns, and 
#' the different traits as components
#' @references A Eklöf, U Jacob, J Kopp, J Bosch, R Castro-Urgal, NP Chacoff,
#'    B Dalsgaard, C de Sassi, M Galetti, PR Guimarães, SB Lomáscolo, AM Martín
#'    González, MA Pizo, R Rader, A Rodrigo, JM Tylianakis, DP Vázquez, and S
#'    Allesina (2013) The dimensionality of ecological networks. Ecology
#'    Letters 16:577-583 (doi: 10.1111/ele.12081)
#' @return The function returns a three-dimensional array where the first and second dimension describes 
#' the species, i.e., the adjacency matrix based on the corrected trait values, and the third dimension 
#' corresponds to the trait analyzed.
#' 
#' @export
#' 
#' @examples
#' data(eklof)
#' overlap.matrices <- build.matrix.from.traits(eklof$network, corr.traits)
build.matrix.from.traits <- function(adj, corr.traits){

# number of traits being analyzed
  numTraits <- length(corr.traits)
  
  # build a container for predictions based on an intervals defined by trait values
  B_all <- array(dim=c(nrow(adj), ncol(adj), length(corr.traits)))

  # loop over the different traits
  for(c in 1:numTraits){
    print(c)
    Trait_temp<-as.matrix(corr.traits[[c]])

    # build a matrix that repeats the row trait value across all columns
    B_temp <- array(dim=c(nrow(adj), ncol(adj), ncol(Trait_temp)))
    # loop over the different columns of the corrected traits (columns are produced for categorical traits)
    for(t in 1:dim(Trait_temp)[2]){
    MatrixTraits <- as.matrix(Trait_temp[,t])[,rep(1,nrow(adj))]

    # build a dummy adjacency matrix but that has trait values instead of 1's and NA's instead of 0's
    Adummy <- MatrixTraits * adj
    Adummy[which(adj == 0)] <- NA

    # determine the min and max trait values when only considering matches to non-NA elements
    MinTrait <- apply(Adummy, 2, function(x){if(sum(is.na(x))==length(x)){NA}else{min(x,na.rm=TRUE)}})
    MaxTrait <- apply(Adummy, 2, function(x){if(sum(is.na(x))==length(x)){NA}else{max(x,na.rm=TRUE)}})

    # get the final predicted adjacency matrix
    B <- (t(t(MatrixTraits)<=MaxTrait) & t(t(MatrixTraits)>=MinTrait))*1

    # turn all NA's back into zeros (since they just mean there are no zeros [and hence min/max values] in that column)
    B[which(is.na(B))] <- 0
    B_temp[,,t] <- B
    }
    # if catecorical trait "summarize" all columns
    FinalMat <- B_temp[,,1]
    
    k <- dim(Trait_temp)[2]
    if(k>1){
      for(z in 2:k){
        FinalMat <- FinalMat + B_temp[,,z]
      }
    }
    # make sure it is still a binary matrix
    FinalMat[FinalMat>0]<-1
    # add this B to the collection of B's
    B_all[,,c] <- FinalMat
  }
   
  rownames(B_all)<-rownames(corr.traits$trait.1)
  colnames(B_all)<-rownames(corr.traits$trait.1)
  dimnames(B_all)[[3]] <- names(corr.traits)

  return(B_all)
}
