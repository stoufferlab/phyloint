#' @title Build a predicted adjacency matrix based on trait values
#' @param adj The adjacency matrix
#' @param trait A vector or matrix of continuous trait values
#' @param direction Whether to build based on rows (1) or columns (2)
#' @references A Eklöf, U Jacob, J Kopp, J Bosch, R Castro-Urgal, NP Chacoff,
#'    B Dalsgaard, C de Sassi, M Galetti, PR Guimarães, SB Lomáscolo, AM Martín
#'    González, MA Pizo, R Rader, A Rodrigo, JM Tylianakis, DP Vázquez, and S
#'    Allesina (2013) The dimensionality of ecological networks. Ecology Letters (doi: 10.1111/ele.12081)
#' @export
build.matrix.from.traits <- function(adj, traits, direction=1){
  # Transpose the matrix
  if(direction==2){
      adj <- t(adj)
  }

  # build a container for predictions based on an intervals defined by trait values
  B_all <- array(dim=c(nrow(adj), ncol(adj), ncol(traits)))

  # loop over the different traits
  for(c in 1:ncol(traits)){
    Trait_temp<-traits[colnames(traits)[c]]

    # build a matrix that repeats the row trait value across all columns
    MatrixTraits <- as.matrix(Trait_temp)[,rep(1,nrow(adj))]

    # build a dummy adjacency matrix but that has trait values instead of 1's and NA's instead of 0's
    Adummy <- MatrixTraits * adj
    Adummy[which(adj == 0)] <- NA

    # determine the min and max trait values when only considering matches to non-NA elements
    MinTrait <- apply(Adummy, 2, min, na.rm = TRUE)
    MaxTrait <- apply(Adummy, 2, max, na.rm = TRUE)

    # do some linear algebra and R magic to get the final predicted adjacency matrix
    B <- (t(t(MatrixTraits)<=MaxTrait) & t(t(MatrixTraits)>=MinTrait))*1

    # turn all NA's back into zeros (since they just mean there are no zeros [and hence min/max values] in that column)
    B[which(is.na(B))] <- 0

    # add this B to the collection of B's
    B_all[,,c] <- B
  }
  
  # start with the first trait as the basis for the predicted matrix
  FinalMat <- B_all[,,1]
  
  # add the matrices corresponding to the other traits
  # (if the trait analyzed is categorical B_all will consist of one matrix for each category (k matrices) and these have to be combined with each other)
  k <- ncol(traits)
  if(k>1){
    for(z in 2:k){
      FinalMat <- FinalMat + B_all[,,z]
    }
  }
      
  # make sure it is still a binary matrix
  FinalMat[FinalMat>0]<-1
  
  # prepare the final return value    
  B<-FinalMat
  rownames(B)<-rownames(traits)
  
  # retranspose in case we started off that way
  if(direction==2){
      B <- t(B)
  }
  
  return(B)
}
