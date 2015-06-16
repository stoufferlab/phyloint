#' @title Estimate a starting temperature for simulated annealing
#' @param adj The adjacency matrix
#' @param phy The phylo object
#' @param prob The target initial probability of accepting a move; defaults to 0.50
#' @param reps The number of reps to estimate across; defaults to 100
#' @param start.type Whether the target probability refers to accepting the min, mean, or max change; defaults to mean
#' 
`sa.start` <- function(adj, phy, prob = 0.50, reps = 100, start.type = 'mean')
{
    gaps <- sum.phylo.gaps(adj,phy)

    de <- rep(NA,reps)
    for(i in 1:reps){
      	phy.new <- swap.tree(phy)
        gaps.new <- sum.phylo.gaps(adj,phy)

  	    de[i] <- abs(gaps.new-gaps)

        phy <- phy.new
        gaps <- gaps.new
    }

    t.start <- switch(start.type,
                      mean = mean(de)/(-log(prob)),
                      max = max(de)/(-log(prob)),
                      min = min(de)/(-log(prob)),
                      stop("Invalid 'type' given to phyloint::sa.start function."))
    
    return(t.start)
}
