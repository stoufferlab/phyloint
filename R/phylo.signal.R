#' @title Function to compute phylogenetic signal of a discrete trait based on comparisons of Pagel's lambda
#' @param trait vector with same order as phy$tip.label
#' @param phy a phylogenetic tree of class 'phylo'
#' @references M Pagel (1999) Inferring the historical patterns of biological
#'   evolution. Nature 401:877-884
#' @export
`phylo.signal` <- function(trait, phy, rep = 999) {
    if (length(attributes(factor(trait))$levels) == length(trait)) 
        stop("Are you sure this variable is categorical?")

    # calculate likelihood corresponding to maximum likelihood value of lambda
    obs <- fitDiscrete(phylo, trait, transform="lambda")

    # calculate likelihood of model with no phylogenetic signal
    null <- fitDiscrete(transform(phylo, "lambda", 0), trait)

    # calculate the likelihood ratio between the two models
    LLR <- -2*(null$opt$lnL - obs$opt$lnL)

    # what is the p value of this likelihood ratio?
    p <- pchisq(LLR, df=1, lower.tail=FALSE)

    return(data.frame(row.names=NULL, lambda=obs$opt$lambda, obs=obs$opt$lnL, null=null$opt$lnL, LLR=LLR, p=p))
}
