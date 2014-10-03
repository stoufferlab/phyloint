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
    # vcv computes the expected varainces and covariances of a continuous trait assuming it evolved under a given model (default=Brownian)
    
    Gmatrix <- ape::vcv(phy);

    # determine the phylogenetic 'correction factor' as defined in Butler et al. (2000)
    correction.factor <- chol(solve(Gmatrix));

    U = data.frame(row.names=rownames(traits));
    for(t in colnames(traits)){
        # so this is if the trait is categorical?
        if(is.factor(traits[,t])){
          traits[,t] <- droplevels(traits[,t])
        }
        
        # a matrix with number of cols equal the number of levels of the categorical trait is produced with ones only in the column for the trait-value the species has
        M <- model.matrix(as.formula(paste0("~0+",t)),traits)
        U[,t] <- correction.factor %*% M;
    }

    td <- geiger::treedata(phy,U);
    return(data.frame(td$data));
}
