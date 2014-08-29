#' @title Minimize diet gaps in an adjacency matrix based on prey phylogeny
#' @param adj an adjacency matrix of class 'matrix'
#' @param phy a phylogenetic tree of class 'phylo'
#' @export
`minimize.phylo.gaps` <- function(adj, phy, p.start=0.5, t.start=NA, t.min=1E-5, t.cooling=0.99, swap.factor=1.0, verbose=FALSE){
    obj <- list(adj=adj,phy.init=phy)
    class(obj) <- 'phylo.gaps'

    # start at a specified initial temperature
    if(!is.na(t.start)){
        t <- t.start
    }else{
        # find a good initial temperature
        t <- sa.start(adj, phy, prob=p.start)
    }

    # how many swaps at each temp?
    nswaps <- round(swap.factor * nrow(adj)^2)

    # get the initial state and energy
    obj$gaps.init <- gaps <- sum.phylo.gaps(adj, phy)

    # set the optima at the intial conditions
    phy.optim <- phy
    gaps.optim <- gaps

    if(verbose){
        cat(paste("Temp","Swaps","Better","Worse","Gaps","Optim",sep="\t"))
        cat("\n")
        cat(paste(t,0,0,0,gaps,gaps.optim,sep="\t"))
        cat("\n")
    }

    # keep going so long as we are above the minimum temp
    while(t > t.min){
        b.count <- 0
        w.count <- 0
        s.count <- 0
        while(s.count < nswaps){
            # suggest a swap of the tree
            phy.new <- swap.tree(phy,1)

            # calculate the new energy
            gaps.new <- sum.phylo.gaps(adj, phy)

            # check whether we accept
            if(gaps.new < gaps || runif(1) < exp(-(gaps.new-gaps)/t)){
                if(gaps.new <= gaps) b.count <- b.count + 1
                if(gaps.new > gaps) w.count <- w.count + 1

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
            cat(paste(t,b.count+w.count,b.count,w.count,gaps,gaps.optim,sep="\t"))
            cat("\n")
        }
        
    }

    obj$phy.optim <- phy.optim
    obj$gaps.optim <- gaps.optim
    
    obj
}
