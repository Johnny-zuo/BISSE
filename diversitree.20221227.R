# Load library
library(ape)
library(diversitree)
library(picante)
library(dplyr)

states1<- read.table("group/zotu/diversitree/states.zotu.16S.txt",header=TRUE,row.names=1,sep="\t")
states <- unlist((states1))

# Read tree
#phy <- ape::read.tree(file="phylogenetic_tree.newick")
tree<- read.tree("group/zotu/diversitree/zotu.16S.tree") #读取微生物类群对应的系统发育树
tree <- prune.sample(states1, tree)
new.tree<-multi2di(tree)
#define function
force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

phy <- force.ultrametric(new.tree,method="extend")
is.ultrametric(phy)

# Run BiSSE
lik <- make.bisse(phy,states,strict=TRUE)			## full BiSSE model
lik.ll <- constrain(lik, lambda0 ~ lambda1, mu0 ~ mu1) 	## null model: same lambda and mu
p1 <- starting.point.bisse(phy)				## estimate initial parameters
p.ll <- c(p1[2], p1[3], p1[5], p1[6])			## extract parameters for null model
fit.ll <- find.mle(lik.ll,p.ll)				## fit null model
p2 <- p1
p2[1] <- coef(fit.ll)[1]
p2[2] <- coef(fit.ll)[1]
p2[3] <- coef(fit.ll)[2]
p2[4] <- coef(fit.ll)[2]
p2[5] <- coef(fit.ll)[3]
p2[6] <- coef(fit.ll)[4]
fit2 <- find.mle(lik,p2)					## fit full BiSSE model with null model as starting point

# Display MLE estimates
print("alpha = 0.5")
round(coef(fit2), 3)
round(coef(fit.ll), 3)
anova.result2 <- anova(fit2, equal.lambda_mu=fit.ll) 	## compare BiSSE model vs null hypothesis
head(anova.result2)

# MCMC
tmp <- mcmc(lik, fit2$par, nsteps=1000, prior=NULL, lower=0, upper=100, w=rep(1, 6), print.every=0)		## perform MCMC to sample theparameter space
lapply(tmp, write, "group/zotu/diversitree/mcmc_output.20221227.txt", append=TRUE, ncolumns=1000)

apply(tmp[,-1],2,mean)

write.csv(tmp,"group/zotu/diversitree/diversitree.zotu.20221227.csv")

save.image("group/zotu/diversitree/diversitree.zotu.20221227.RData")

