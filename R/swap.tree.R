#' @title Permute a random pair of sister clades in a phylogenetic tree
#' @param phy The phylo object
#' @export
`swap.tree` <- function(phy)
{
  # find the existing xy-location of all tips and nodes
  tc <- tree.coord(phy)

  # pick a random node to swap around
  n <- sample(phy$edge[,1],1)

  # find the children of this node
  c <- phangorn::Children(phy, n)

  # find all tips that are descendants of the various children
  ts <- sapply(c,function(x)(unlist(phangorn::Descendants(phy,x,type="tips"))))
  
  # reorder the children so that they have increasing y-values
  c <- c[order(tc[c,"y"])]

  # order the tips within each group
  ts <- lapply(ts,function(x)(x[order(tc[x,"y"])]))

  # what will the new order be?
  ts.new <- sample(ts)

  # make absolutely sure things actually change
  while(sum(unlist(ts.new) == unlist(ts)) != 0){
    ts.new <- sample(ts)
  }

  # set the new location of the tips
  tc[unlist(ts.new),"y"] <- tc[unlist(ts),"y"]

  # figure out the new top to bottom order
  new.order <- rownames(tc)[order(tc$y[1:Ntip(phy)])];

  # reorder the phylo object
  phy.swapped <- ape::rotateConstr(phy, new.order);

  return(phy.swapped);
}
