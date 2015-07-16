#' @title Plot a tree coordinate object in the xy plane and an ordered adjacency matrix by its side
#'
#' @param x an object of class 'phylo.gaps'
#' @param edge.width a numeric vector giving the width of the branches of the plotted phylogeny. These are taken to be in the same order as the component 'edge' of the phylogeny. If fewer widths are given than the length of 'edge', the values are recycled.
#' @param link.size a numeric value giving the size of the squares that indicate trophic links in the visual representation of the adjacency matrix
#'
#' @export
#'
`plot.phylo.gaps` <- function(x,type='optim',tree.lwd=1,link.size=1)
{
  adj <- x$adj
  phy <- x$phy.optim

  par(mfrow=c(1,2))

  # plot the tree itself (from left to right)
  par(pty='s')
  plot(phy,edge.width=tree.lwd,show.tip.label=F)
    
  # build a community to plot the matrix using cheddar
  tc <- tree.coord(phy)
  nodes <- tc[rownames(adj),]
  nodes$node <- rownames(nodes)
  comm <- cheddar::Community(nodes=nodes,
                             properties=list(title='barsom'),
                             trophic.links=cheddar::PredationMatrixToLinks(adj))

  # order the community as defined by the phylogeny
  comm <- cheddar::OrderCommunity(comm,'y')
  
  # plot the matrix so that it's lined up with the tree
  par(pty='s')
  cheddar::PlotPredationMatrix(comm,cex=link.size,pch=15,main='',xlab='',ylab='',asp=1,consumer.order=order(cheddar::NumberOfResources(comm)))
  box(lwd=tree.lwd)
}
