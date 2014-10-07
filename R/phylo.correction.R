#' @title Remove potential phylogenetic signal from a set of continuous or discrete traits
#' @param traits a data frame with species as rows and traits as columns
#' @param phy a phylogenetic tree of class 'phylo'
#' @references MA Butler, TW Schoener, and JB Losos (2000) The relationship between
#'    sexual size dimorphism and habitat use in Greater Antillean Anolis
#'    lizards. Evolution 54:259-272
#' @export
`phylo.correction` <- function(traits, phy)
{
    # calculate the Gmatrix based on the phylogeny   
    Gmatrix <- ape::vcv(phy);

    # determine the phylogenetic 'correction factor' as defined in Butler et al. (2000)
    correction.factor <- chol(solve(Gmatrix));

    # determine the original order of the traits
    rnames.orig <- rownames(traits);

    # make sure that the order of the traits corresponds to the correction.factor matrix
    traits <- traits[as.character(rownames(correction.factor)), ,drop=FALSE];

    # build a container for the corrected trait values
    U = data.frame(row.names=rownames(correction.factor));

    # correct everything trait by trait (useful if we have a mixture of continuous and discrete traits)
    for(t in colnames(traits)){
        if(is.factor(traits[,t])){
          traits[,t] <- droplevels(traits[,t])

          # just in case dropping levels leaves us with only one
          if(nlevels(traits[,t]) == 1)
            traits[,1] <- 1
        }

        # build a model matrix based on the trait values (this helps automate the separation of categorical traits)
        M <- model.matrix(as.formula(paste0("~0+",t)),traits)

        # use the correction factor to determine the corrected trait values
        corrected.traits <- data.frame(correction.factor %*% M);

        # populate the U container
        for(ct in colnames(corrected.traits)){
            U[,ct] <- corrected.traits[,ct]
        }
    }

    # reorder the corrected traits as they originally came
    # TODO: should make sure this works when a species was dropped because it isn't in the tree
    U <- U[rnames.orig, ,drop=FALSE];
 
    return(U);
}
