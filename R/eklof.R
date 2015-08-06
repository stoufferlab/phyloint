#' Sample food-web dataset with which to study phylogenetic intervality
#'
#' The object 'eklof' has three items: a 'network' object, a 'tree' object, 
#' and a 'traits' object. The 'network' object is a binary SxS matrix where
#' a value of 1 indicates that the row species consumes the column species.
#' The 'tree' object is a phylogenetic tree of class 'phylo' (see 
#' 'ape::read.tree'). The 'traits' object is a data frame with species names
#' as row names and columns for different species traits (which can be both
#' continuous or categorical). These data constitute a synthetic dataset 
#' created to demonstrate how to use the 'phyloint' package and hence do
#' not have an associated reference.
#'
"eklof"
