#' @title Minimize diet gaps in an adjacency matrix based on prey phylogeny using steepest descent
#' @param adj an adjacency matrix of class 'matrix'
#' @param phy a phylogenetic tree of class 'phylo'
#' @export
`minimize.phylo.gaps.sd` <- function(adj, phy, verbose=FALSE){
    # make sure we shuffle the tree at the start to give us "random" initial conditions
    for(i in seq_len(Nnode(phy))) phy <- swap.tree(phy)

    # create a container for everything
    obj <- list(adj=adj,phy.init=phy)
    class(obj) <- 'phylo.gaps'

    # get the initial state and energy
    obj$gaps.init <- gaps <- sum.phylo.gaps(adj, phy)

    # set the optima at the intial conditions
    phy.optim <- phy
    gaps.optim <- gaps

    # a bit of verbiage for the people
    if(verbose){
        cat(paste("Step","Worse","Same","Better","Gaps",sep="\t"))
        cat("\n")
        cat(paste(0,NA,NA,NA,gaps,sep="\t"))
        cat("\n")
    }

    # get the labels of internal nodes
    ns <- unique(phy$edge[,1])

    # keep trying to swap internal nodes as long as we can improve things
    m.count <- 0
    moved <- TRUE
    while(moved){
        moved <- FALSE
        
        # perform all swaps
        phy.new <- lapply(ns, function(x)(swap.tree(phy,x)))

        # calculate the gaps for all swaps
        gaps.new <- unlist(lapply(phy.new, function(x)(sum.phylo.gaps(adj,x))))

        # how many are better?
        b.count <- sum(gaps.new < gaps)

        # how many are equivalent?
        e.count <- sum(gaps.new == gaps)

        # how many are worse?
        w.count <- sum(gaps.new > gaps)


        if(b.count > 0){
            # give us another shot to try again
            m.count <- m.count + 1
            moved <- TRUE

            # find the best move
            gaps <- gaps.optim <- min(gaps.new)

            # select the best at random (in case there are ties)
            phy <- phy.optim <- phy.new[[sample(which(gaps.new == gaps.optim),1)]]
        }

        # keep the people placid
        if(verbose){
            cat(paste(m.count,w.count,e.count,b.count,gaps,sep="\t"))
            cat("\n")
        }
        
    }

    # return an object that corresponds to the optimum
    obj$phy.optim <- phy.optim
    obj$gaps.optim <- gaps.optim
    obj
}
