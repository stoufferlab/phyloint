#' @title Permute a random pair of sister clades in a phylogenetic tree
#' @param phy The phylo object
#' @param rep The number of swaps to perform
#' @export
`swap.tree` <- function(phy, rep = 1)
{
  tc <- tree.coord(phy)
  
  for(k in 1:rep){
      # pick a random node to swap around
      n <- sample(phy$edge[,1],1);

      # find all of the subsets of tips that are below this node
      c <- phangorn::Children(phy, n);
      ts <- sapply(c,function(x)(unlist(phangorn::Descendants(phy,x,type="tips"))))

      # find the existing location of the tips
      tu <- unlist(ts)
      y.old <- tc[tu,"y"]

      # find a new location for them
      tu.new <- unlist(sample(ts))
      while(sum(tu.new != tu) == 0) tu.new <- unlist(sample(ts))

      # set the new location of the tips
      tc[tu.new,"y"] <- y.old
    
      # and fix all parents to line up
      anc <- unique(phangorn::Ancestors(phy, c, 'parent'));
      repeat{
        for(a in anc){
          tc[a, "y"] <- mean(tc[phangorn::Children(phy, a), "y"]);
        }
        anc <- unique(phangorn::Ancestors(phy, anc, 'parent'));
        if(sum(anc == 0) == 1) break;
      }
  }
  
  # figure out the new top to bottom order
  new.order <- rownames(tc)[order(tc$y[1:Ntip(phy)])];

  # reorder the phylo object
  phy.swapped <- ape::rotateConstr(phy, new.order);

  return(phy.swapped);
}
