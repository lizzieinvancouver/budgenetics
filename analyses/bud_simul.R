
library(TreeSim)

###
# Parameters for the simulations

nspecies=10
nindividuals=10
speciationrate=0.5
extinctionrate=0.2
ratefactor=10 # This is approx. the rate increase for the tree within species
              # compared to the rate of the species tree

###
# Simulate species tree

spetree <- sim.bd.taxa(n=nspecies,numbsim=1,lambda=speciationrate,mu=extinctionrate,complete=FALSE)
#rename tips to be able to find them below...
tips<-c("sp1","sp2","sp3","sp4","sp5","sp6","sp7","sp8","sp9","sp10")
spetree[[1]]$tip.label <- tips
plot(spetree[[1]]);add.scale.bar()

###
# Simulate individual trees, one per species

# here, extinction = 0
poptrees <- sim.bd.taxa(n=nindividuals,numbsim=nspecies,
                        lambda=speciationrate*ratefactor,mu=0,complete=FALSE)

###
# Stich trees together

combtree<-spetree[[1]]
for(i in 1:length(spetree[[1]]$tip.label)){
  newtips <- paste(letters[i],seq(1,length(poptrees[[i]]$tip.label)),sep="")
  poptrees[[i]]$tip.label <- newtips
  tipname <- tips[i]
  tipnumber <- seq(1,length(combtree$tip.label))[combtree$tip.label==tipname]
  poptrees[[i]]$root.edge<-0
  combtree<- bind.tree(combtree,poptrees[[i]],where=tipnumber,position=0)
}
plot(combtree,cex=0.5);add.scale.bar()
