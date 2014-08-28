#' @title Calculate the number of diet gaps in an adjacency matrix given an ordering compatible with the phylogeny
#' @param adj The adjacency matrix
#' @param tc The tree coordinate object
#' @export
`sum.phylo.gaps` <- function(adj, phy){
    tc <- tree.coord(phy)
    nodes <- tc[rownames(adj),]
    nodes$node <- rownames(nodes)
    comm <- Community(nodes=nodes,
	    			  properties=list(title='sa.start'),
		    		  trophic.links=PredationMatrixToLinks(adj))
    comm <- OrderCommunity(comm,'y')
  	gaps <- SumDietGaps(comm)
    return(gaps)
}
