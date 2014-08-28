#' @title Minimize diet gaps in an adjacency matrix based on prey phylogeny
#' @param adj an adjacency matrix of class 'matrix'
#' @param phy a phylogenetic tree of class 'phylo'
#' @export
`minimize.phylo.gaps` <- function(adj, phy, t.min=1E-5, t.cooling=0.99, swap.factor=1.0, verbose=FALSE){
    # find a good minimum temperature
    t <- sa.start(adj, phy)

    # how many swaps at each temp?
    nswaps <- round(swap.factor * dim(adj)[1]**2)

    # get the initial state and energy
    gaps <- sum.phylo.gaps(adj, phy)

    # set the optima at the intial conditions
    phy.optim <- phy
    gaps.optim <- gaps

    if(verbose){
        cat(paste("Temp","Gaps","Optim",sep="\t"))
        cat("\n")
        cat(paste(t,gaps,gaps.optim,sep="\t"))
        cat("\n")
    }

    # keep going so long as we are above the minimum temp
    while(t > t.min){
        s.count <- 0
        while(s.count < nswaps){
            # suggest a swap of the tree
            phy.new <- swap.tree(phy,1)

            # calculate the new energy
            gaps.new <- sum.phylo.gaps(adj, phy)

            # check whether we accept
            if(gaps.new < gaps || runif(1) < exp(-(gaps.new-gaps)/t)){
                phy <- phy.new
                gaps <- gaps.new

                # is the new state the best ever?
                if(gaps < gaps.optim){
                    phy.optim <- phy
                    gaps.optim <- gaps
                }
            }

            s.count <- s.count + 1
        }
        t <- t * t.cooling
        
        if(verbose){
            cat(paste(t,gaps,gaps.optim,sep="\t"))
            cat("\n")
        }
        
    }

    obj <- list(adj=adj,phy=phy.optim,gaps=gaps.optim)
    class(obj) <- 'phylo.gaps'
    obj
}
