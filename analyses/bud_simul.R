###
# Scripts for simulating a variable accorging to a linear model formula
# with the predictor and/or the residual error being phylogenetically 
# correlated and with infraspecific sampling

# Simon Joly, sept. 2016

###
# Load packages

require(nlme)
require(ape)
require(phytools)
require(geiger)


# Phylogenetic simulations ------------------------------------------------

###
# Parameters for the phylogenetic simulations

nspecies=10
nindividuals=10
# Give a relative scale to the species tree vs. the infraspecific tree.
# The total tree length is fixed to 1.
species.scale=0.9 #Total
infrasp.scale=1-species.scale
# Pagel's lambda parameter, to give more or less importance to the
# phylogenetic relationships
lambda = 1.0

###
# Simulate species tree

spetree <- pbtree(n=nspecies,nsim=1,b=1,complete=FALSE,scale=species.scale)
#rename tips to be able to find them below...
tips<-c("sp1","sp2","sp3","sp4","sp5","sp6","sp7","sp8","sp9","sp10")
spetree$tip.label <- tips
plot(spetree);add.scale.bar()

###
# Simulate individual trees, one per species

# Here, extinction = 0
poptrees <- pbtree(n=nindividuals,nsim=nspecies,b=1,scale=infrasp.scale)

###
# Stich trees together

combtree<-spetree
for(i in 1:length(spetree$tip.label)){
  newtips <- paste(letters[i],seq(1,length(poptrees[[i]]$tip.label)),sep="")
  poptrees[[i]]$tip.label <- newtips
  tipname <- tips[i]
  tipnumber <- seq(1,length(combtree$tip.label))[combtree$tip.label==tipname]
  poptrees[[i]]$root.edge<-0
  combtree<- bind.tree(combtree,poptrees[[i]],where=tipnumber,position=0)
}
plot(combtree,cex=0.5);add.scale.bar()

###
# Create tree with no infraspecific info
combtree0<-spetree
for(i in 1:length(spetree$tip.label)){
  newtips <- paste(letters[i],seq(1,length(poptrees[[i]]$tip.label)),sep="")
  poptrees[[i]]$tip.label <- newtips
  #collapse branch lengths
  poptrees[[i]]$edge.length <- rep(0.0001,length.out=length(poptrees[[i]]$edge.length))
  tipname <- tips[i]
  tipnumber <- seq(1,length(combtree0$tip.label))[combtree0$tip.label==tipname]
  poptrees[[i]]$root.edge<-0
  combtree0<- bind.tree(combtree0,poptrees[[i]],where=tipnumber,position=0)
}
# Scale the whole tree to have length 1
combtree0$edge.length <- combtree0$edge.length/max(nodeHeights(combtree0)[,2])*1
plot(combtree0,cex=0.5);add.scale.bar()

###
# Rescale trees

combtree <- rescale(combtree,"lambda",lambda)
combtree0 <- rescale(combtree0,"lambda",lambda)

# Character simulations -------------------------------------------------

# These follow Revell et al. 2010. MEE

###
# Parameters

n.sp <- length(combtree$tip.label) #number of species
sigma.sq.x = 1 # Rate of evolution in X
sigma.sq.e = 0.7 # Rate of residual variation
B = 0.75 # Regression slope
C <- vcv.phylo(combtree) # VCV matrix
u <- rnorm(n=n.sp,mean=10) # Random deviates from the normal distribution 
v <- rnorm(n=n.sp,mean=0) # Random deviates from the normal distribution 

###
# Simulate variables

# Get values for the X variable that are phylogeneticaly correlated
x <- t(chol(C*sigma.sq.x)) %*% u

# Generate residual error that is phylogenetically correlated
e <- t(chol(C*sigma.sq.e)) %*% v

# Get response variable
y <- x %*% B + e


###
# Other options for character simulations

# Get values for the X variable that are NOT phylogeneticaly correlated
# x <- sqrt(sigma.sq.x) * u



# Model fitting -----------------------------------------------------------

require(nlme)

# data frame
thedata <- data.frame(y=y,x=x,sp=gsub("[0-9]*","",combtree$tip.label))

# formula
f1 <- formula(y~x)

# Phylo correlation
phylo.corr <- corBrownian(phy=combtree)
phylo.corr0 <- corBrownian(phy=combtree0) # Tree without infrapsecific branch lengths

# Fit phylogenetic regression with species as random effect
gls.mod0 <- lme(f1,random = ~1 | sp, data=thedata)
gls.mod1 <- lme(f1,random = ~1 | sp, correlation = phylo.corr, data=thedata)
gls.mod2 <- lme(f1,random = ~1 | sp, correlation = phylo.corr0, data=thedata)

# Fit phylogenetic regression WITHOUT species as random effect
gls.mod3 <- gls(f1,correlation = phylo.corr, data=thedata)
gls.mod4 <- gls(f1,correlation = phylo.corr0, data=thedata)

# Compare all models
AIC(gls.mod0,gls.mod1,gls.mod2,gls.mod3,gls.mod4)

