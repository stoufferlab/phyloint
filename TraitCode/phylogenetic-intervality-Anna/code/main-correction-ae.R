# load the libraries and source files
require(ape)
require(picante)
source("taxonomy.phylo.tree.R")
source("phylo.correction.R")
source("CallTraitsForNetwork.R")



SaveResultsTableAllCombinationTraits <- function(Network,Direction){

    # call the function creating the tree and the modified adj.matrix
    TaxTreeRes <- TaxonomyTree(Network)
    
     tree <- TaxTreeRes$tree
     adj <- TaxTreeRes$adj
     removeSp <- TaxTreeRes$removeSp

    # create a dataframe with the traits
    
    # load the list of traits for the specific network
    ListTraits <- as.matrix(read.csv(paste("../data/",Network,"/ListOfTraits.txt",sep=""),header=FALSE))
    NumTraits <- dim(ListTraits)[1]
   
    
    for(i in 1:dim(ListTraits)[1]){
        
         # just to make sure the categorical trait is counted as catecorical even if it is listed as a number
        if (ListTraits[i,3]=="Categ"){
            X[,1]<-as.character(X[,1])
        }
        
            X<-as.data.frame(read.csv(paste("../data/",Network,"/traits/",ListTraits[i,2],sep=""),header=FALSE, sep=" "))
            
            #rownames(X)<- as.character(paste("sp_",(1:dim(X)[1]), sep=""))
            X<-X[-removeSp,]
            # we put the rownames to 1:numberOfSpecies as that is what is done in the TaxonomyTable
            X<-data.frame(sp=rownames(adj),t=X)
            X<-X[match(tree$tip.label, X$sp),]
            X<-X$t
            X<-data.frame(row.names=tree$tip.label, t=X)
           

            # calculate the phylogeneically-corrected traits
            U_cont <- phylo.correction(X,tree)
            assign(paste0("corr_", ListTraits[i,1]), U_cont)
        }
    }

    # construct a list consisting of all phylogenetic corrected traits
    CorrectedTraits <- list(corr_BM=corr_BM, corr_HB=corr_HB, corr_FT=corr_FT, corr_FM=corr_FM,
    corr_MC=corr_MC,corr_MB=corr_MB)
    
    # reorder adj so the species ordering is the same as tree
    A<-adj[match(tree$tip.label, rownames(adj)), match(tree$tip.label, rownames(adj))]
    
    # calculate the A' matrices (interaction matrices basd on the corrected traits)
        AllMatrices <- SaveResultsCorretedTraits(Network,A, Direction, CorrectedTraits, ListTraits)


# combine traits
TraitNames <- ListTraits[,1]
## MaxNumComb
MaxNumComb <- min(c(5,NumTraits))
NumberComb <- NumTraits*(94+NumTraits*(5+NumTraits*(25+NumTraits*(-5+NumTraits))))/120+1
if (NumTraits<6){
    NumberComb <- 2^MaxNumComb-1
    MaxNumComb <- NumTraits-1
}
ResultTable <- matrix(NA,NumberComb,2)
Labels <- rep("aa",NumberComb)
print(paste("MaxNumComb",MaxNumComb,"NumTraits",NumTraits,"NumberComb",NumberComb))
Current <- 1
if (NumberComb>1){
    for (k in 1:MaxNumComb){
        ## compute all combinations of NumTraits traits taken k at a time
        Combinations <- combn(NumTraits,k)#(7,1)
        ## for each combination
        for(j in 1:dim(Combinations)[2]){
            ## start with the matrix corresponding to the first trait
            FinalMat <- AllMatrices[,,Combinations[1,j]]
            ## and multiply for the matrices corresponding to the other traits
            if (k>1){
                for (z in 2:k){
                    FinalMat <- FinalMat*AllMatrices[,,Combinations[z,j]]
                }
            }
            ## save the results in the table
            ## first, save number of traits used in this combination
            ResultTable[Current,1] <- k
            ## second, compute overlap
            ResultTable[Current,2] <- sum(A)/sum(FinalMat)
           
            Labels[Current] <- paste(TraitNames[Combinations[,j]],collapse="+")
            Current<-Current+1
        }
	}
}

## the last one is all traits taken together
FinalMat<-AllMatrices[,,1]
if (NumTraits>1){
    for (z in 2:NumTraits){
        FinalMat <- FinalMat*AllMatrices[,,z]
    }
}

ResultTable[Current,1] <- NumTraits
## second, compute overlap
ResultTable[Current,2] <- sum(A)/sum(FinalMat)
rownames(ResultTable) <- Labels
colnames(ResultTable) <- c("NumTraits", "Overlap")

## write and save all results
write.table(ResultTable, paste("../results/AllCombTraits_",Network,"_dir",Direction,".txt", sep=""), sep=" ")
}


#Networks <- as.matrix(read.csv("../Data/TraitWebs/ListNetworks.csv", header=TRUE))[,1]
#Networks<-c("kongsfjorden", "loughhyne", "stmarks", "weddell", "ythanjacob")
Networks<-c("kongsfjorden")
Directions <- c(1,2)

for(i in 1:length(Networks)){
    print(Networks[i])
    for(d in 1:length(Directions)){
        print(Directions[d])
        SaveResultsTableAllCombinationTraits(Network=Networks[i], Direction=Directions[d])
    }
}











