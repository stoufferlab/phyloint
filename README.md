## phyloint

This `R` package provides functions to calculate properties related to food-web intervality (Stouffer et al. 2006; Eklöf et al. 2013) while mixing species phylogenies with interaction matrices as described in:

Anna Eklöf and Daniel B. Stouffer (2016) "The phylogenetic component of food-web structure and intervality." *Theoretical Ecology* 9(1) 107-115 doi:[10.1007/s12080-015-0273-9][doi].

#### How to install directly from github

    require(devtools)
    install_github("stoufferlab/phyloint")

#### How to conduct key analyses in Eklöf & Stouffer

    require(ape)
    require(phyloint)
        
    # load in some sample data
    data(eklof)
        
    # use steepest descent to minimize the number of phylogenetic gaps in the data
    gaps.phylo <- minimize.phylo.gaps.sd(adj=eklof$network, phy=eklof$tree)
        
    # calculate trait values when accounting for the expected covariance between species
    corrected.traits <- phylo.correction(traits=eklof$traits, phy=eklof$tree)

    # produce the overlap matrices based on the phylogentically corrected trait values 
    overlap.matrices <- build.matrix.from.traits(eklof$network, corrected.traits)

    # calculate the overlap values for all different combinations of traits and produce a result table 
    combinations <- overlap.combinations(eklof$network, overlap.matrices)

[doi]: http://dx.doi.org/10.1007/s12080-015-0273-9
