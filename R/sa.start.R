#' @title Estimate a starting temperature for simulated annealing
#' @param adj The adjacency matrix
#' @param phy The phylo object
#' @param prob The target initial probability of accepting a move; defaults to 0.50
#' @param reps The number of reps to estimate across; defaults to 100
#' @export
`sa.start` <- function(adj, phy, prob = 0.50, reps = 100)
{
    gaps <- sum.phylo.gaps(adj,phy)

    de <- rep(NA,reps)
    for(i in 1:reps){
      	phy.new <- swap.tree(phy,1)
        gaps.new <- sum.phylo.gaps(adj,phy)

  	    de[i] <- abs(gaps-gaps.new)

        phy <- phy.new
        gaps <- gaps.new
    }

    t.start <- mean(de)/(-log(prob))
    return(t.start)
}
