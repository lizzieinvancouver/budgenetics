############################
#
# Function : duplicate.tips
#
# Author: Simon Joly, March 2016
#
# Description:
#
#   Given a tree 'tree', the function duplicates the tips 'n' times,
#   adding a number from 2:n after the tip name. The number is 
#   separated from the tip names by the 'sep.char' character.
#

duplicate.tips <- function(tree,n=2,sep.char="_") {
  require(phytools)
  if (n<2) stop("\"n\" needs to be at lest 2")
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  species = tree$tip.label
  for (i in 1:length(species)) {
    aspecies = species[i]
    for (j in 2:n){
      atip = which(tree$tip.label==aspecies)
      suff <- formatC(j, width = (floor(log(n,10))+1), format = "d", flag = "0")
      tree <- bind.tip(tree,paste(aspecies,sep.char,suff,sep=""),edge.length=0,
                       where=atip, position=0, interactive=FALSE)
    }
    #replace original tip name
    suff <- formatC(1, width = (floor(log(n,10))+1), format = "d", flag = "0")
    tree$tip.label <- gsub(paste("^",aspecies,"$",sep=""),
                           paste(aspecies,sep.char,suff,sep=""),tree$tip.label)
  }
  return(tree)
}
