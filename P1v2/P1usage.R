library(reshape2)
library(ggplot2)
library(subplex)
library(gridExtra)


rm(list=ls())
#load(file = 'AA.RData')
library(P1)

#set.seed(11)
estimations = T
seed = round(runif(1,1,100000))
#seed = 7
s <- phyl2(tt=15,seed=seed)
plot(s$newick)
p <- subplex(par = c(8,0.175,0.9),fn = llik,n = s$n, E = s$E, t = s$t)
c(p$par[1],(p$par[1]-p$par[3])/p$par[2],p$par[3])
pa=c(p$par[1],p$par[2],p$par[3])
#lm(p$setoflambda~p$setofbeta)

dropex <- drop.fossil(s$newick)
plot(dropex)
s2 <- phylo2p(tree=dropex,ct=ct)

rt = reconst_tree3(bt=s2$t,pars=pa)
p3 <- subplex(par = c(8,0.175,0.9),fn = llik,n = rt$n, E = rt$E, t = rt$t)$par
c(p3[1],(p3[1]-p3[3])/p3[2],p3[3]) # estimation of parameters of the reconstructed tree
#P = matrix(nrow=n_it,ncol=4)
a = data.frame(t=seq(0,15,by=0.1))
a[[2]] = approx(cumsum(s$t),s$n,xou=seq(0,15,by=0.1))$y
n_it = 100
if(estimations) P = matrix(nrow=n_it,ncol=4)
S = vector('list',length=n_it)
for (i in 1:n_it){
  rt = reconst_tree2(bt=s2$t,pars=pa)
  a[[i+2]] = approx(cumsum(rt$t),rt$n,xou=seq(0,15,by=0.1))$y
  if(estimations){
  pars = subplex(par = c(8,0.175,0.9),fn = llik,n = rt$n, E = rt$E, t = rt$t)$par
  pars2 = c(p$par[1],pars[2],p$par[3],(p$par[1]-p$par[3])/p$par[2])
  P[i,] = pars2}
  S[[i]] = rt
}

melted = melt(a, id.vars="t")
p <- ggplot(data=melted, aes(x=t, y=value, group=variable)) + geom_line() + ylab('# Lineages') + xlab('Time') + geom_line(data=melted[melted$variable=='V2',],aes(x=t, y=value),color='blue') +ggtitle(paste('seed',seed,'(mle)'))#  + geom_line(data=data.frame(t=seq(0,15,by=0.1),variable="mean",value=rowMeans(AA, na.rm = TRUE, dims = 1)),aes(x=t, y=value),color='green')
print(p)
pars = mle_dd_setoftrees(S,draw=F)
c(pars$lambda,pars$beta,pars$mu,(pars$lambda-pars$mu)/pars$beta)

#par(mfrow=c(1,3))
if(estimations){
plot1 = ggplot(NULL, aes(x=P[,1])) + geom_histogram()+xlab(paste(expression(lambda),"0"))+geom_vline(xintercept=pa[1])
plot2 = ggplot(NULL, aes(x=P[,2])) + geom_histogram()+xlab(paste(expression(beta),"0"))+geom_vline(xintercept=pa[2])
plot3 = ggplot(NULL, aes(x=P[,3])) + geom_histogram()+xlab(paste(expression(mu),"0"))+geom_vline(xintercept=pa[3])
plot4 = ggplot(NULL, aes(x=P[,4])) + geom_histogram()+xlab(paste(expression(K)))+geom_vline(xintercept=(pa[1]-pa[3])/pa[2])
library(grid)
grid.newpage()
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
pushViewport(viewport(layout = grid.layout(1, 4)))
print(plot1, vp = vplayout(1, 1))
print(plot2, vp = vplayout(1, 2))
print(plot3, vp = vplayout(1, 3))
print(plot4, vp = vplayout(1, 4))
}


## EM
s <- phyl2(seed=8)
plot(s$newick)
p <- subplex(par = c(8,0.175,0.9),fn = llik,n = s$n, E = s$E, t = s$t)
c(p$par[1],(p$par[1]-p$par[3])/p$par[2],p$par[3])
pa=c(p$par[1],p$par[2],p$par[3])


dropex <- drop.fossil(s$newick)
plot(dropex)
s2 <- phylo2p(tree=dropex)

ptm <- proc.time()
init_pars = pa
pars=init_pars
npar=init_pars
n_it = 100
P = matrix(nrow=n_it,ncol=4)
bt = s2$t
for (i in 1:n_it){
  print(pars)
  n_it = 20
  S = vector('list',length=n_it)
  P = matrix(nrow=n_it,ncol=4)
  for (j in 1:n_it){
    rec_tree = reconst_tree2(bt=bt,pars=npar)
    pars = subplex(par = c(8,0.175,0.9),fn = llik,n = rec_tree$n, E = rec_tree$E, t = rec_tree$t)$par
    #pars2 = c(p$par[1],pars[2],p$par[3],(p$par[1]-p$par[3])/p$par[2])
    P[j,] = c(pars,(pars[1]-pars[3])/pars[2])
    #print(length(rec_tree$t))
    S[[j]] = rec_tree
  }
  npar = colMeans(P)[1:3]
  print(paste('npar es',npar))
  #print(proc.time() - ptm)
  #ptm <- proc.time()
  #print('optimization')
  pars = subplex(par = c(8,0.175,0.9),fn = llik_st , setoftrees = S)$par
  #mle_dd_setoftrees(S)
  print(paste('-loglikelihood equal to',llik_st(pars,setoftrees = S)))
  pars1 = c(pars[1],pars[2],pars[3],(pars[1]-pars[3])/pars[2])
  P[i,] = pars1
}


  plot1 = ggplot(NULL, aes(x=P[,1])) + geom_histogram()+xlab(paste(expression(lambda),"0"))+geom_vline(xintercept=pa[1])
  plot2 = ggplot(NULL, aes(x=P[,2])) + geom_histogram()+xlab(paste(expression(beta),"0"))+geom_vline(xintercept=pa[2])
  plot3 = ggplot(NULL, aes(x=P[,3])) + geom_histogram()+xlab(paste(expression(mu),"0"))+geom_vline(xintercept=pa[3])
  plot4 = ggplot(NULL, aes(x=P[,4])) + geom_histogram()+xlab(paste(expression(K)))+geom_vline(xintercept=(pa[1]-pa[3])/pa[2])
  library(grid)
  grid.newpage()
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
  pushViewport(viewport(layout = grid.layout(1, 4)))
  print(plot1, vp = vplayout(1, 1))
  print(plot2, vp = vplayout(1, 2))
  print(plot3, vp = vplayout(1, 3))
  print(plot4, vp = vplayout(1, 4))



ptm = proc.time()
mle_dd_setoftrees(S)
proc.time() - ptm

ptm = proc.time()
subplex(par = c(8,0.175,0.9),fn = llik_st , setoftrees = S)$par
proc.time() - ptm


#exp2

n_it = 1000
HH = matrix(nrow=n_it,ncol=4)
P = matrix(nrow=n_it,ncol=4)
PR = matrix(nrow=n_it,ncol=4)
for (i in 1:n_it){
  s = phyl2(seed=runif(1,1,1000000))
  p <- subplex(par = c(8,0.175,0.9),fn = llik,n = s$n, E = s$E, t = s$t)
  pa=c(p$par[1],p$par[2],p$par[3])
  dropex <- drop.fossil(s$newick)
  PR[i,] = c(pa[1],pa[2],pa[3],(pa[1]-pa[3])/pa[2])
  s2 <- phylo2p(tree=dropex)
  n_it = 100
  bt = s2$t
  S = vector('list',length=n_it)
  H = matrix(nrow=n_it,ncol=3)
    for (j in 1:n_it){
      rec_tree = reconst_tree2(bt=bt,pars=pa)
      pars = subplex(par = c(8,0.175,0.9),fn = llik,n = rec_tree$n, E = rec_tree$E, t = rec_tree$t)$par
      H[j,] = pars
      S[[j]] = rec_tree
    }
    npar = colMeans(H)
    pars = subplex(par = c(8,0.175,0.9),fn = llik_st , setoftrees = S)$par
    HH[i,] = c(npar[1],npar[2],npar[3],(npar[1]-npar[3])/npar[2])
    pars1 = c(pars[1],pars[2],pars[3],(pars[1]-pars[3])/pars[2])
    P[i,] = pars1
}


a = data.frame(t=1:n_it)
a[[2]] = P[,1]
a[[3]] = HH[,1]
melted = melt(a, id.vars="t")
p <- ggplot(data=melted, aes(x=t, color=variable)) + geom_histogram() + ylab('# Lineages') + xlab('Time')# + geom_line(data=melted[melted$variable=='V2',],aes(x=t, y=value),color='blue') +ggtitle(paste('seed',seed,'(mle)'))#  + geom_line(data=data.frame(t=seq(0,15,by=0.1),variable="mean",value=rowMeans(AA, na.rm = TRUE, dims = 1)),aes(x=t, y=value),color='green')
print(p)

plot(PR[,1],P[,1])
abline(0,1)


par_est_vis <- function(P,par){
  if (par == 1){
    int = 0.8 #change for general case
    parname = 'lambda'
  }
  if (par==2){
    int= 0.0175
    parname = 'beta'
  }
  if (par == 3){
    int = 0.1
    parname = 'mu'
  }
  if (par ==4){
    int = 40
    parname = 'K'
    P = P[P[,4]<100,]
  }
hist_top <- ggplot()+geom_histogram(aes(P[,par])) + geom_vline(xintercept=int) + xlab('MLE from incomplete tree')
empty <- ggplot()+geom_point(aes(1,1), colour="white")+
  theme(axis.ticks=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank())

scatter <- ggplot()+geom_point(aes(P[,par], PR[,par]))+ geom_abline(intercept = 0, slope = 1) + ylab(TeX(paste('$\\hat{\\',parname,'}_{C}$',sep='')))+xlab(TeX(paste('$\\hat{\\',parname,'}_{I}$',sep='')))
hist_right <- ggplot()+geom_histogram(aes(PR[,par]))+coord_flip()+ geom_vline(xintercept=int) +xlab('MLE of complete tree')

grid.arrange(hist_top, empty, scatter, hist_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
}


## exp0.1
    # calculating expected value
    # itt=100000
    # Vals = vector(mode='list',length = itt)
    # for (k in 1:itt){
    #   s1 = phyl2(seed=runif(1,1,10000))
    #   Vals[[k]]=data.frame(time=cumsum(s1$t),n=s1$n)
    # }
    # AA = matrix(nrow=length(seq(0,15,by=0.1)),ncol=itt)
    # for (i in 1:itt){
    #   AA[,i]=approx(Vals[[i]]$time,Vals[[i]]$n,xout = seq(0,15,by=0.1))$y
    # }
#    save(AA, file='AA.RData')
    #


# for (j in 1:30){
#     s = phyl2(seed = j)
#     plot(s$newick)
#     a = data.frame(t=seq(0,15,by=0.1))
#     a[[2]] = approx(cumsum(s$t),s$n,xou=seq(0,15,by=0.1))$y
#     dropex <- drop.fossil(s$newick)
#     plot(dropex)
#     s2 <- phylo2p(tree=dropex,ct=ct)
#     n_it = 100
#     for (i in 1:n_it){
#       rt = reconst_tree2(bt=s2$t,pars=c(0.8,0.0175,0.1))
#       a[[i+2]] = approx(cumsum(rt$t),rt$n,xou=seq(0,15,by=0.1))$y
#     }
#     #library('reshape2')
#     #library('ggplot2')
#     melted = melt(a, id.vars="t")
#     p <- ggplot(data=melted, aes(x=t, y=value, group=variable)) + geom_line() + ylab('# Lineages') + xlab('Time') + geom_line(data=melted[melted$variable=='V2',],aes(x=t, y=value),color='blue')  + geom_line(data=data.frame(t=seq(0,15,by=0.1),variable="mean",value=rowMeans(AA, na.rm = TRUE, dims = 1)),aes(x=t, y=value),color='green')+ggtitle(paste('seed',j))
#     print(p)
# }


## exp1
# ptm <- proc.time()
# init_pars = c(0.8,0.0175,0.1)
# pars=init_pars
# n_it = 10
# P = matrix(nrow=n_it,ncol=3)
# PP = matrix(nrow=n_it,ncol=3)
# VP = matrix(nrow=n_it,ncol=3)
# for (i in 1:n_it){
#   print(c(pars[1],pars[2],pars[3]))
#   s <- phyl2(seed=runif(1,1,100000))
#   p2 <- mle_dd(s,draw=F)
#   pars_mle <- c(p2$lambda,p2$beta,p2$mu)
#   dropex <- drop.fossil(s$newick)
#   plot(dropex)
#   s2 <- phylo2p(tree=dropex,ct=ct)
#   n_it = 100
#   M = matrix(nrow=n_it,ncol=3)
#   MM = matrix(nrow=n_it,ncol=3)
#   for (j in 1:n_it){
#     rec_tree = reconst_tree2(bt=s2$t,pars=init_pars)
#     rec_tree2 = reconst_tree2(bt=s2$t,pars=pars_mle)
#     #ltt = approx(cumsum(rt$t),rt$n,xou=seq(0,15,by=0.1))$y
#     p3 <- mle_dd(rec_tree,draw=F)
#     M[j,]=c(p3$lambda,p3$beta,p3$mu)
#     p4 <- mle_dd(rec_tree2,draw=F)
#     MM[j,]=c(p4$lambda,p4$beta,p4$mu)
#   }
#   pars = colMeans(M)
#   P[i,] = pars
#   parsM = colMeans(MM)
#   PP[i,] = parsM
#   VP[i,] = pars_mle
# }
# proc.time() - ptm



#library(DDD)
#s3 <- dd_sim(c(0.8,0.1,40),15)
#s1 <-L2p(L=s3$L)
#p2 <- mle_dd(s1)
#c(p2$lambda,(p$lambda-p$mu)/p$beta,p2$mu)
##dd_ML(brts = s3$L[,1])


### parallel

# library(parallel)
# library(foreach)
# library(doParallel)
#
# ptm <- proc.time()
# u01 <- dd_par(tt=15,mu=0.1,it=100)
# proc.time() - ptm
#
# int=NULL
# slo=NULL
# lambda=NULL
# beta=NULL
# mu=NULL
#
# for (i in 1:100){
#   u=u01[[i]]
#   slo[i] = lm(u$setofbeta ~ u$setofbeta)$coef[2]
#   int[i] = lm(u$setofbeta ~ u$setofbeta)$coef[1]
#   lambda[i] = u$lambda
#   beta[i] = u$beta
#   mu[i] = u$mu
# }
#
# hist(lambda)
# hist(mu)
# hist((lambda-mu)/beta)


######  EM
#
# s <- phyl2()  # simulate phylogenetic tree
# p <- mle_dd(s) # estimate parameters of complete tree
# c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)
#
# plot(s$newick) #plot complete tree
# dropex <- drop.fossil(s$newick) # drop extinct species
# plot(dropex) #plot observed tree
# s2 <- phylo2p(tree=dropex,ct=ct)
# p2 <- mle_dd(s2) # estimate parameters of incomplete tree
# c(p2$lambda,(p2$lambda-p2$mu)/p2$beta,p2$mu) #we can see no nice estimation and mu=0 always
# rec_tree = reconst_tree(bt=s2$t,pars=c(0.8,0.0175,0.1),tt=15) #simulate a random recostructed tree with given parameters
# p3 <- mle_dd(rec_tree)
# c(p3$lambda,(p3$lambda-p3$mu)/p3$beta,p3$mu) # estimation of parameters of the reconstructed tree
#
#
# ## exp
# ptm <- proc.time()
# init_pars = c(0.8,0.0175,0.1)
# pars=init_pars
# n_it = 10
# for (i in 1:n_it){
#   print(c(pars[1],pars[2],pars[3]))
#   n_it = 10
#   M = matrix(nrow=n_it,ncol=3)
#   for (j in 1:n_it){
#     rec_tree = reconst_tree(bt=s2$t,pars=pars,tt=15)
#     p3 <- mle_dd(n=rec_tree$n,E=rec_tree$E,t=rec_tree$t)
#     M[j,]=c(p3$lambda,p3$beta,p3$mu)
#   }
#   pars = colMeans(M)
# }
# proc.time() - ptm
#
# setwd("~/Documents/Code/P1/data")
# dendroica =  read.csv(file = 'Dendroica_branchtimes.csv',header = F,sep = ',')
# dendroica = rev(as.numeric(matrix(dendroica)))
# bt = c(bt[1],diff(dendroica))
# init_pars = c(3,0.03,0.1)
# pars=init_pars
# n_it = 40
# N = matrix(nrow=n_it,ncol=4)
# for (i in 1:n_it){
#   print(c(pars[1],pars[2],pars[3]))
#   n_it = 100
#   M = matrix(nrow=n_it,ncol=3)
#   for (j in 1:n_it){
#     rec_tree = reconst_tree_op1(bt=bt,pars=pars,tt=5)
#     p3 <- mle_dd(rec_tree,draw=F)
#     M[j,]=c(p3$lambda,p3$beta,p3$mu)
#   }
#   pars = colMeans(M)
#   print(i)
#   N[i,]=c(pars[1],pars[2],pars[3],(pars[1]-pars[3])/pars[2])
# }
#
#
# ######
#
# s = phyl2()
# p <- mle_dd(s) # estimate parameters of complete tree
# c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)
#
# plot(s$newick) #plot complete tree
# dropex <- drop.fossil(s$newick) # drop extinct species
# plot(dropex) #plot observed tree
# s2 <- phylo2p(tree=dropex,ct=ct)
#
# ptm <- proc.time()
# init_pars = c(0.8,0.0175,0.1)
# #init_pars = c(8,0.175,0.3)
# pars=init_pars
# n_it = 50
# M = matrix(nrow=n_it,ncol=3)
# for (i in 1:n_it){
#   print(c(pars[1],pars[2],pars[3]))
#   n_it = 10
#   S=list()
#   for (j in 1:n_it){
#     rec_tree = reconst_tree(bt=s2$t,pars=pars,tt=15)
#     S[[j]] = rec_tree
#     #p3 <- mle_dd(n=rec_tree$n,E=rec_tree$E,t=rec_tree$t)
#     #M[j,]=c(p3$lambda,p3$beta,p3$mu)
#   }
#   pars = mle_dd_setoftrees(S)
#   pars = c(pars$lambda,pars$beta,pars$mu)
#   M[i,]=c(pars[1],pars[2],pars[3])
#   print(proc.time() - ptm)
#   ptm <- proc.time()
# }
# proc.time() - ptm
#
# ## another exercise
#
# s = phyl2()
# p <- mle_dd(s) # estimate parameters of complete tree
# c(p$lambda,(p$lambda-p$mu)/p$beta,p$mu)
#
# plot(s$newick) #plot complete tree
# dropex <- drop.fossil(s$newick) # drop extinct species
# plot(dropex) #plot observed tree
# s2 <- phylo2p(tree=dropex,ct=ct)
#
# setwd("~/Documents/Code/P1/data")
# dendroica =  read.csv(file = 'Dendroica_branchtimes.csv',header = F,sep = ',')
# dendroica = rev(as.numeric(matrix(dendroica)))
# #bt = c(dendroica[1],diff(dendroica))
# bt = diff(dendroica)
# init_pars = c(6,0.3,0.1)  #
# init_pars = c(0.8,0.0175,0.1) # converge a
# init_pars = c(4, 0.1, 1)
# ptm <- proc.time()
# #init_pars = c(0.8,0.0175,0.1)
# #init_pars = c(8,0.175,0.3)
# pars=init_pars
# n_it = 50
# M = matrix(nrow=n_it,ncol=3)
# for (i in 1:n_it){
#   print(c(pars[1],pars[2],pars[3],(pars[1]-pars[3])/pars[2]))
#   n_it = 10
#   S = vector('list',length=n_it)
#   for (j in 1:n_it){
#     rec_tree = reconst_tree_op1(bt=bt,pars=pars,tt=5)
#     S[[j]] = rec_tree
#     #p3 <- mle_dd(n=rec_tree$n,E=rec_tree$E,t=rec_tree$t)
#     #M[j,]=c(p3$lambda,p3$beta,p3$mu)
#   }
#   pars = mle_dd_setoftrees(S)
#   pars = c(pars$lambda,pars$beta,pars$mu)
#   M[i,]=c(pars[1],pars[2],pars[3])
#   print(proc.time() - ptm)
#   ptm <- proc.time()
# }
# proc.time() - ptm
#
# foraminifera = read.csv(file='foraminifera_branchtimes.csv',header=F,sep=',')
# foraminifera = rev(as.numeric(matrix(foraminifera)))
#
# tt=64.95
# #foraminifera = rev(as.numeric(matrix(foraminifera)))
# #bt = foraminifera
# #bt = c(bt[1],diff(foraminifera))
# bt = diff(foraminifera)
# init_pars = c(4, 0.1, 1)
# ptm <- proc.time()
# #init_pars = c(0.8,0.0175,0.1)
# #init_pars = c(8,0.175,0.3)
# pars=init_pars
# n_it = 10
# M = matrix(nrow=n_it,ncol=3)
# for (i in 1:n_it){
#   print(c(pars[1],pars[2],pars[3],(pars[1]-pars[3])/pars[2]))
#   n_it = 200
#   S = vector('list',length=n_it)
#   for (j in 1:n_it){
#     rec_tree = reconst_tree_op1(bt=bt,pars=pars,tt=tt)
#     S[[j]] = rec_tree
#     #p3 <- mle_dd(n=rec_tree$n,E=rec_tree$E,t=rec_tree$t)
#     #M[j,]=c(p3$lambda,p3$beta,p3$mu)
#   }
#   pars = mle_dd_setoftrees(S)
#   pars = c(pars$lambda,pars$beta,pars$mu)
#   M[i,]=c(pars[1],pars[2],pars[3])
#   print(proc.time() - ptm)
#   ptm <- proc.time()
# }
# proc.time() - ptm


a = matrix(nrow=100,ncol=8)
for (j in 1:100){

s = phyl2(seed = round(runif(1,1,10000)))
pars = subplex(par = c(8,0.175,0.9),fn = llik,n = s$n, E = s$E, t = s$t)$par
#a1 = matrix(nrow = 100, ncol = 4)
a2 = matrix(nrow = 100, ncol = 4)
#a3 = matrix(nrow = 100, ncol = 4)
#a4 = matrix(nrow = 100, ncol = 4)
s2 = drop.fossil(s$newick)
s2 = phylo2p(s2)
for (i in 1:100){
#  rt1 = reconst_tree(bt = s2$t,pars = pars)
#  a1[i,] = c(length(rt1$t),subplex(par = c(8,0.175,0.9),fn = llik,n = rt1$n, E = rt1$E, t = rt1$t)$par)
  rt2 = reconst_tree2(bt = s2$t,pars = pars)
  a2[i,] = c(length(rt2$t),subplex(par = c(8,0.175,0.9),fn = llik,n = rt2$n, E = rt2$E, t = rt2$t)$par)
#  rt3 = reconst_tree3(bt = s2$t,pars = pars)
#  a3[i,] = c(length(rt3$t),subplex(par = c(8,0.175,0.9),fn = llik,n = rt3$n, E = rt3$E, t = rt3$t)$par)
}

#ggplot(data = data.frame(a1), aes(x=X1,y=X2)) + geom_point() + geom_vline(xintercept = length(s$t)) + geom_hline(yintercept = pars[1])
#ggplot(data = data.frame(a2), aes(x=X1,y=X2)) + geom_point()+ geom_vline(xintercept = length(s$t)) + geom_hline(yintercept = pars[1])+ geom_hline(yintercept = 0.8,color='blue')+geom_hline(yintercept=mean(a2[,2]),color='green')
#ggplot(data = data.frame(a3), aes(x=X1,y=X2)) + geom_point()+ geom_vline(xintercept = length(s$t)) + geom_hline(yintercept = pa[1])
a[j,] = c(colMeans(a2),pars,length(s$t))

}

ggplot(data = data.frame(a),aes(x=X2,X5)) + geom_point() + geom_abline(intercept = 0, slope = 1)
ggplot(data = data.frame(a),aes(x=X4,X7)) + geom_point() + geom_abline(intercept = 0, slope = 1)
ggplot(data = data.frame(a),aes(x=X3,X6)) + geom_point() + geom_abline(intercept = 0, slope = 1)

ggplot(data = data.frame(a),aes(x=X1,X8)) + geom_point() + geom_abline(intercept = 0, slope = 1)







