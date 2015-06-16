#' @title Calculate the xy-coordinates of nodes in a phylogenetic tree
#' @param phy A phylogenetic tree of class 'phylo'
#' @param type Specify whether to return coordinates for all nodes or just tips
#' @export
`tree.coord` <- function(phy, type=c('all','tips'))
{
    # make sure that a valid type has been passed to this function
    type <- match.arg(type)

    # make a faux plot (just so that we can grab the coordinates from it)
    pdf(file=NULL)
    ape::plot.phylo(phy, type='phylogram', plot=FALSE)
    dev.off()
    last.plot <- get("last_plot.phylo", envir = .PlotPhyloEnv)

    # save the x,y coordinates of everything
    x <- last.plot$xx
    y <- last.plot$yy

    # make sure things map to the default labelling scheme from plot.phylo
    labels <- c(phy$tip.label,(length(phy$tip.label)+1):(length(phy$tip.label)+phy$Nnode))

    tc <- data.frame(row.names=labels, x=x, y=y)

    if(type == 'tips') tc <- subset(tc, rownames(tc) %in% phy$tip.label)

    return(tc);
}
