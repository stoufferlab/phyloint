## This code construct a "phylogenetic tree" based on the taxonomy of the species
## Since there are species that are of the same genus individual identities are given as a species number

TaxonomyTree<-function(Network){
# load the data
adj <-as.matrix(read.table(paste("../data/",Network,"/",Network,".txt",sep=""),header=FALSE))

nspecies <- dim(adj)[1] #the rows are resources and the columns are predators

# the matrix MUST have non-numeric row and column names
# name the rows and cols based on phylogeny
TaxonomyTable<-as.data.frame(read.csv(paste("../data/",Network,"/taxonomy/taxonomyall.txt",sep=""),header=FALSE, sep=" "))
TaxonomyTable[,9]<- read.csv(paste("../data/",Network,"/taxonomy/species.txt",sep=""), header=FALSE)
TaxonomyTable[,10]<-as.numeric(rownames(TaxonomyTable))
colnames(TaxonomyTable)<-c("name","nr","kingdom","phyla","class", "order", "family", "genus","species","sp_number")

paste0("sp_",as.character(1:nspecies))
rownames(adj) <- as.character(paste("sp_",TaxonomyTable$sp_number, sep=""))
colnames(adj) <- as.character(paste("sp_", TaxonomyTable$sp_number, sep=""))

# find the rows that has na's
row.has.na <- apply(TaxonomyTable, 1, function(x){any(is.na(x))})
removeSp <- which(row.has.na)

# remove all species that not has a full taxonomy described
adj <- adj[-removeSp,-removeSp]
TaxonomyTable<- TaxonomyTable[-removeSp,]

# and change also the sp_number in TaxonomyTable
TaxonomyTable$sp_number <-rownames(adj)

# convert to factors
TaxonomyTable$sp_number<-as.factor(TaxonomyTable$sp_number)
TaxonomyTable$species<-as.factor(TaxonomyTable$species)
TaxonomyTable$genus<-as.factor(TaxonomyTable$genus)
TaxonomyTable$family<-as.factor(TaxonomyTable$family)
TaxonomyTable$order<-as.factor(TaxonomyTable$order)
TaxonomyTable$class<-as.factor(TaxonomyTable$class)
TaxonomyTable$phyla<-as.factor(TaxonomyTable$phyla)
TaxonomyTable$kingdom<-as.factor(TaxonomyTable$kingdom)
TaxonomyTable$nr<-as.factor(TaxonomyTable$nr)
TaxonomyTable$name<-as.factor(TaxonomyTable$name)

# generate a "phylogenetic" tree for these species based on taxonomy
tree<-(as.phylo(~kingdom/phyla/class/order/family/genus/species/sp_number, data=TaxonomyTable))

# root the tree
tree$root.edge <- 0

# compute branch-lengths
tree <- compute.brlen(tree)

TaxTreeRes <- list(tree=tree, adj=adj, removeSp=removeSp)

}

## Intervality code 
# calculate a decent SA starting temperature
#gaps.min <- minimize.phylo.gaps(adj, tree, t.min=0.01, t.cooling=0.9, swap.factor=1.0, verbose=TRUE)

# plot things all nice and ugly (i.e., beware that I haven't tweak things much at all!)
#plot.phylo.gaps(adj,gaps.min$tc,tree.lwd=2,link.size=1)

