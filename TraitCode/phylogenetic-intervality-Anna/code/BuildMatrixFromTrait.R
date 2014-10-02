## Builds a matrix using a trait:
## Basically, each col (row) interacts with all the rows (cols)
## whose trait value is between a min and a max
## Trait must be a vector of real numbers

BuildMatrixTrait <- function(A,MatrixType,Direction,Trait){
  ## prey choose predators
  print(paste("Direction=",Direction,"now"))
  if (Direction==2){
      A <- t(A)
  }
  ## Dimensions of the matrix
  NR <- dim(A)[1]
  NC <- dim(A)[2]

  B_all <- array(dim=c(NR, NC, dim(Trait)[2]))

# here we need to loop over the different categories
      for(c in 1:dim(Trait)[2]){
      
        Trait_temp<-Trait[[c]]

      ## ADummy is a matrix in which each 1 is multiplied by the
      ## corresponding Trait value
        ADummy <- A*Trait_temp
      ## for each column, find the min and the max
      ## finding the maximum is easy
        MaxTrait <- apply(ADummy,2,max)
      ## the minimum could potentially return 0 for all species
      ## to prevent this, assign a large value
      ## to all the elements that are zero
      ## and then take the min
        ADummy[ADummy==0] <- max(MaxTrait)+1.0
        MinTrait <- apply(ADummy,2,min)
      ## MatrixTraits is a matrix repeating the
      ## trait value by column
        MatrixTraits <- matrix(Trait_temp,NR,NR) ## we can do this because the matrix is squared
      ## now do some linear algebra magic to get the final matrix
        B_all[,,c] <- (t(t(MatrixTraits)<=MaxTrait)& t(t(MatrixTraits)>=MinTrait))*1
  }
    
      FinalMat <- B_all[,,1]
      k<-dim(Trait)[2]
      ## and multiply for the matrices corresponding to the other traits
      # if the trait analyzed is categorical B_all will consist of one matrix for each category (k matrices) and these has to be combined with each other.
      if (k>1){
          for (z in 2:k){
              FinalMat <- FinalMat+B_all[,,z]
          }
      }
      
      FinalMat[FinalMat>0]<-1
      
  B<-FinalMat
  rownames(B)<-rownames(Trait)
  
  if (Direction==2){
      B <- t(B)
  }
  
  return(B)

}


