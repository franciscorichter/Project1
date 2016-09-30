
library(P1)

#set.seed(11)
s <- phyl2(tt=15,seed=runif(1,1,100000))
p <- mle_dd(s,draw=F)
c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)
#lm(p$setoflambda~p$setofbeta)

plot(s$newick)
dropex <- drop.fossil(s$newick)
plot(dropex)
s2 <- phylo2p(tree=dropex,ct=ct)
p2 <- mle_dd(s2,draw=F)
c(p2$lambda,p2$beta,p2$mu)

rec_tree = reconst_tree(bt=s2$t,pars=c(0.8,0.0175,0.1),tt=15) #simulate a random recostructed tree with given parameters
p3 <- mle_dd(s)
c(p3$lambda,(p3$lambda-p3$mu)/p3$beta,p3$mu) # estimation of parameters of the reconstructed tree

library(DDD)
s3 <- dd_sim(c(0.8,0.1,40),15)
s1 <-L2p(L=s3$L)
p2 <- mle_dd(s1)
c(p2$lambda,(p$lambda-p$mu)/p$beta,p2$mu)
#dd_ML(brts = s3$L[,1])


### parallel

library(parallel)
library(foreach)
library(doParallel)

ptm <- proc.time()
u01 <- dd_par(tt=15,mu=0.1,it=100)
proc.time() - ptm

int=NULL
slo=NULL
lambda=NULL
beta=NULL
mu=NULL

for (i in 1:100){
  u=u01[[i]]
  slo[i] = lm(u$setofbeta ~ u$setofbeta)$coef[2]
  int[i] = lm(u$setofbeta ~ u$setofbeta)$coef[1]
  lambda[i] = u$lambda
  beta[i] = u$beta
  mu[i] = u$mu
}

hist(lambda)
hist(mu)
hist((lambda-mu)/beta)


######  EM

s <- phyl2()  # simulate phylogenetic tree
p <- mle_dd(n=s$n,E=s$E,t=s$t) # estimate parameters of complete tree
c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)

plot(s$newick) #plot complete tree
dropex <- drop.fossil(s$newick) # drop extinct species
plot(dropex) #plot observed tree
s2 <- phylo2p(tree=dropex,ct=ct)

p2 <- mle_dd(n=s2$n,E=s2$E,t=s2$t) # estimate parameters of incomplete tree
c(p2$lambda,(p2$lambda-p2$mu)/p2$beta,p2$mu) #we can see no nice estimation and mu=0 always

rec_tree = reconst_tree(bt=s2$t,pars=c(0.8,0.0175,0.1),tt=15) #simulate a random recostructed tree with given parameters

p3 <- mle_dd(n=rec_tree$n,E=rec_tree$E,t=rec_tree$t)
c(p3$lambda,(p3$lambda-p3$mu)/p3$beta,p3$mu) # estimation of parameters of the reconstructed tree


## exp
ptm <- proc.time()
init_pars = c(0.8,0.0175,0.1)
pars=init_pars
n_it = 10
for (i in 1:n_it){
  print(c(pars[1],pars[2],pars[3]))
  n_it = 10
  M = matrix(nrow=n_it,ncol=3)
  for (j in 1:n_it){
    rec_tree = reconst_tree(bt=s2$t,pars=pars,tt=15)
    p3 <- mle_dd(n=rec_tree$n,E=rec_tree$E,t=rec_tree$t)
    M[j,]=c(p3$lambda,p3$beta,p3$mu)
  }
  pars = colMeans(M)
}
proc.time() - ptm

setwd("~/Documents/Code/P1/data")
dendroica =  read.csv(file = 'Dendroica_branchtimes.csv',header = F,sep = ',')
dendroica = rev(as.numeric(matrix(dendroica)))
bt = c(bt[1],diff(dendroica))
init_pars = c(3,0.03,0.1)
pars=init_pars
n_it = 40
N = matrix(nrow=n_it,ncol=4)
for (i in 1:n_it){
  print(c(pars[1],pars[2],pars[3]))
  n_it = 100
  M = matrix(nrow=n_it,ncol=3)
  for (j in 1:n_it){
    rec_tree = reconst_tree_op1(bt=bt,pars=pars,tt=5)
    p3 <- mle_dd(rec_tree,draw=F)
    M[j,]=c(p3$lambda,p3$beta,p3$mu)
  }
  pars = colMeans(M)
  print(i)
  N[i,]=c(pars[1],pars[2],pars[3],(pars[1]-pars[3])/pars[2])
}






######

s = phyl2()
p <- mle_dd(s) # estimate parameters of complete tree
c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)

plot(s$newick) #plot complete tree
dropex <- drop.fossil(s$newick) # drop extinct species
plot(dropex) #plot observed tree
s2 <- phylo2p(tree=dropex,ct=ct)

ptm <- proc.time()
init_pars = c(0.8,0.0175,0.1)
#init_pars = c(8,0.175,0.3)
pars=init_pars
n_it = 50
M = matrix(nrow=n_it,ncol=3)
for (i in 1:n_it){
  print(c(pars[1],pars[2],pars[3]))
  n_it = 10
  S=list()
  for (j in 1:n_it){
    rec_tree = reconst_tree(bt=s2$t,pars=pars,tt=15)
    S[[j]] = rec_tree
    #p3 <- mle_dd(n=rec_tree$n,E=rec_tree$E,t=rec_tree$t)
    #M[j,]=c(p3$lambda,p3$beta,p3$mu)
  }
  pars = mle_dd_setoftrees(S)
  pars = c(pars$lambda,pars$beta,pars$mu)
  M[i,]=c(pars[1],pars[2],pars[3])
  print(proc.time() - ptm)
  ptm <- proc.time()
}
proc.time() - ptm

## another exercise

s = phyl2()
p <- mle_dd(s) # estimate parameters of complete tree
c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)

plot(s$newick) #plot complete tree
dropex <- drop.fossil(s$newick) # drop extinct species
plot(dropex) #plot observed tree
s2 <- phylo2p(tree=dropex,ct=ct)

setwd("~/Documents/Code/P1/data")
dendroica =  read.csv(file = 'Dendroica_branchtimes.csv',header = F,sep = ',')
dendroica = rev(as.numeric(matrix(dendroica)))
#bt = c(dendroica[1],diff(dendroica))
bt = diff(dendroica)
init_pars = c(6,0.3,0.1)  #
init_pars = c(0.8,0.0175,0.1) # converge a
init_pars = c(4, 0.1, 1)
ptm <- proc.time()
#init_pars = c(0.8,0.0175,0.1)
#init_pars = c(8,0.175,0.3)
pars=init_pars
n_it = 50
M = matrix(nrow=n_it,ncol=3)
for (i in 1:n_it){
  print(c(pars[1],pars[2],pars[3],(pars[1]-pars[3])/pars[2]))
  n_it = 10
  S = vector('list',length=n_it)
  for (j in 1:n_it){
    rec_tree = reconst_tree_op1(bt=bt,pars=pars,tt=5)
    S[[j]] = rec_tree
    #p3 <- mle_dd(n=rec_tree$n,E=rec_tree$E,t=rec_tree$t)
    #M[j,]=c(p3$lambda,p3$beta,p3$mu)
  }
  pars = mle_dd_setoftrees(S)
  pars = c(pars$lambda,pars$beta,pars$mu)
  M[i,]=c(pars[1],pars[2],pars[3])
  print(proc.time() - ptm)
  ptm <- proc.time()
}
proc.time() - ptm

foraminifera = read.csv(file='foraminifera_branchtimes.csv',header=F,sep=',')
foraminifera = rev(as.numeric(matrix(foraminifera)))

tt=64.95
#foraminifera = rev(as.numeric(matrix(foraminifera)))
#bt = foraminifera
#bt = c(bt[1],diff(foraminifera))
bt = diff(foraminifera)
init_pars = c(4, 0.1, 1)
ptm <- proc.time()
#init_pars = c(0.8,0.0175,0.1)
#init_pars = c(8,0.175,0.3)
pars=init_pars
n_it = 10
M = matrix(nrow=n_it,ncol=3)
for (i in 1:n_it){
  print(c(pars[1],pars[2],pars[3],(pars[1]-pars[3])/pars[2]))
  n_it = 200
  S = vector('list',length=n_it)
  for (j in 1:n_it){
    rec_tree = reconst_tree_op1(bt=bt,pars=pars,tt=tt)
    S[[j]] = rec_tree
    #p3 <- mle_dd(n=rec_tree$n,E=rec_tree$E,t=rec_tree$t)
    #M[j,]=c(p3$lambda,p3$beta,p3$mu)
  }
  pars = mle_dd_setoftrees(S)
  pars = c(pars$lambda,pars$beta,pars$mu)
  M[i,]=c(pars[1],pars[2],pars[3])
  print(proc.time() - ptm)
  ptm <- proc.time()
}
proc.time() - ptm
