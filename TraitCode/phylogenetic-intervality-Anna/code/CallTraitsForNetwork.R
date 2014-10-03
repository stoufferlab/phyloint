require(combinat)
source("BuildMatrixFromTrait.R")


SaveResultsCorretedTraits<- function(Network, A, Direction, CorrectedTraits, ListTraits){
    ## Get the adjacency matrix
    print("Doing Traits")
   
    ## Number of traits
    NumTraits = dim(ListTraits)[1]
 
    ## for each trait we build the matrix and store it in an array
    AllMatrices <- array(0,c(dim(A),NumTraits))
    for (i in 1:NumTraits){

        # read the corected trait value
        TraitValue<-CorrectedTraits[[i]]

        # if we want to analyze both from the prey perspective and predator perspective
        B <- BuildMatrixTrait(A=A,MatrixType=as.character(ListTraits[i,4]),Direction,Trait=TraitValue)
        
        # B is then the matrix bulit on the corrected values of one of the traits
        
        print(paste(as.character(ListTraits[i,1]),sum(B*A)/sum(B)))
        AllMatrices[,,i]<-B
        
    }

    return(AllMatrices)
}