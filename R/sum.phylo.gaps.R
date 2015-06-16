#' @title Calculate the number of diet gaps in an adjacency matrix given an ordering compatible with the phylogeny
#' @param adj The adjacency matrix
#' @param phy A phylogenetic tree of class 'phylo'
#' @export
`sum.phylo.gaps` <- function(adj, phy){
    tc <- tree.coord(phy)
    nodes <- tc[rownames(adj),]
    nodes$node <- rownames(nodes)
    comm <- Cheddar::Community(nodes=nodes,
	    				      properties=list(title='varza'),
	    				      trophic.links=PredationMatrixToLinks(adj))
    comm <- Cheddar::OrderCommunity(comm,'y')
  	gaps <- Cheddar::SumDietGaps(comm)
    return(gaps)
}
