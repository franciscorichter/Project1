# Exercises
#1. get back parameters of a tree
n_it = 100
MM = matrix(nrow=n_it,ncol=3)
for (i in 1:n_it){
  s <- phyl2(seed=runif(1,1,10000000))  # simulate phylogenetic tree
  p <- mle_dd(s) # estimate parameters of complete tree
  MM[i,]=c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)
}

#2. given a reconstructed tree, simulate a complete tree and get back parameters
sM = matrix(nrow=100,ncol=3)
for(j in 1:100){
  print(j)
s = phyl2(seed=j)
plot(s$newick) #plot complete tree
dropex <- drop.fossil(s$newick) # drop extinct species
plot(dropex) #plot observed tree
s2 <- phylo2p(tree=dropex,ct=ct)
n_it = 10
M = matrix(nrow=n_it,ncol=3)
pars=c(0.8,0.0175,0.1)
#pars = c(mle_dd(s)$lambda,mle_dd(s)$beta,mle_dd(s)$mu)
M2 = vector(mode='numeric',length=n_it)
for (i in 1:n_it){
  rec_tree = reconst_tree_op1(bt=s2$t,pars=pars,tt=15) # simulate phylogenetic tree
  p <- mle_dd(rec_tree,draw=F) # estimate parameters of complete tree
  M2[i]=length(rec_tree$t)
  M[i,]=c(p$lambda,p$beta,p$mu)
}
sM[j,] = colMeans(M)
}

plot1 = ggplot(NULL, aes(x=sM[,1])) + geom_histogram()+xlab(paste(expression(lambda),"0"))
plot2 = ggplot(NULL, aes(x=M[,2])) + geom_histogram()+xlab(paste(expression(beta),"0"))
plot3 = ggplot(NULL, aes(x=M[,3])) + geom_histogram()+xlab(paste(expression(mu),"0"))
library(grid)
grid.newpage()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pushViewport(viewport(layout = grid.layout(1, 3)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))
print(plot3, vp = vplayout(1, 3))
colMeans(M)
apply(M,2,summary)

#3. Given a reconst tree, simulate a sets of complete trees, and get back the parameters from the whole llik of the set.
s = phyl2(seed=10)
p <- mle_dd(s) # estimate parameters of complete tree
c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)

plot(s$newick) #plot complete tree
dropex <- drop.fossil(s$newick) # drop extinct species
plot(dropex) #plot observed tree
s2 <- phylo2p(tree=dropex,ct=ct)
ds
plot1 = ggplot(NULL, aes(x=MMM[,1])) + geom_histogram()+xlab(paste(expression(lambda),"0"))
plot2 = ggplot(NULL, aes(x=MMM[,2])) + geom_histogram()+xlab(paste(expression(beta),"0"))
plot3 = ggplot(NULL, aes(x=MMM[,3])) + geom_histogram()+xlab(paste(expression(mu),"0"))
library(grid)
grid.newpage()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pushViewport(viewport(layout = grid.layout(1, 3)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))
print(plot3, vp = vplayout(1, 3))
colMeans(M)
a = apply(MMM,2,summary)
xtable(t(a),digits=4)


#### EM
s = phyl2(seed=8)
p <- mle_dd(s) # estimate parameters of complete tree
c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)

plot(s$newick) #plot complete tree
dropex <- drop.fossil(s$newick) # drop extinct species
plot(dropex) #plot observed tree
s2 <- phylo2p(tree=dropex,ct=ct)

ptm <- proc.time()
init_pars = c(8,0.175,0.05)
pars=init_pars
n_it = 50
MMM = matrix(nrow=n_it,ncol=3)
for (i in 1:n_it){
  n_it = 100
  S = vector('list',length=n_it)
  for (j in 1:n_it){
    rec_tree = reconst_tree_op1(bt=bt,pars=pars,tt=15)
    print(length(rec_tree$t))
    S[[j]] = rec_tree
  }
  pars = mle_dd_setoftrees(S)
  pars = c(pars$lambda,pars$beta,pars$mu)
  MMM[i,]=c(pars[1],pars[2],pars[3])
  print(proc.time() - ptm)
  ptm <- proc.time()
  print(pars)
}
proc.time() - ptm
