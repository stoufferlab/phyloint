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

      # swapping the edges is equivalent to swapping descendants of all children
      c <- phangorn::Children(phy, n);

      # and now we need to rearrange their descendants
      t1 <- unlist(phangorn::Descendants(phy, c[1], type='tips'));
      t2 <- unlist(phangorn::Descendants(phy, c[2], type='tips'));

      # since everything is arbitrary, set d1 to be the "lowest" set
      if(min(tc[t2,"y"]) < min(tc[t1,"y"])){
        tmp <- t2;
        t2 <- t1;
        t1 <- tmp;
      }

      # everything relates to a shift based on the size of the corresponding groups
      tc[t1, "y"] <- tc[t1, "y"] + length(t2);
      tc[t2, "y"] <- tc[t2, "y"] - length(t1);
    
      # and fix all parents to line up
      anc <- unique(phangorn::Ancestors(phy, c(t1,t2), 'parent'));
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
