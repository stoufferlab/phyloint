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

    U = data.frame(row.names=rownames(traits));
    for(t in colnames(traits)){
        if(is.factor(traits[,t])){
          traits[,t] <- droplevels(traits[,t])
        }

        M <- model.matrix(as.formula(paste0("~0+",t)),traits)
        U[,t] <- correction.factor %*% M;
    }

    td <- geiger::treedata(phy,U);
    return(data.frame(td$data));
}
