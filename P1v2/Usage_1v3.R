st = sim_phyl(seed=7)
plot(st$newick)
st2 = drop.fossil(st$newick)
plot(st2)
st2 = phylo2p(st2)
# Recontruction and estimations
rec = rec_tree(wt = st2$t)
p <- subplex(par = c(8,0.175,0.9),fn = llik,n = st$n, E = st$E, t = st$t)
c(p$par[1],(p$par[1]-p$par[3])/p$par[2],p$par[3])
p2 <- subplex(par = c(8,0.175,0.9),fn = llik,n = rec$n, E = rec$E, t = rec$wt)
c(p2$par[1],(p2$par[1]-p2$par[3])/p2$par[2],p2$par[3])
rec = rec_tree(wt = st2$t,pars = p$par)
p3 <- subplex(par = p$par,fn = llik,n = rec$n, E = rec$E, t = rec$wt)
c(p3$par[1],(p3$par[1]-p3$par[3])/p3$par[2],p3$par[3])
###

# Simulations
n_sim = 1000
n_trees = 10
MP = matrix(nrow=n_sim,ncol=3)
MPf = matrix(nrow=n_sim,ncol=3)
RP = matrix(nrow=n_sim,ncol=3)
setofweith = NULL
for(j in 1:n_sim){
  st = sim_phyl(seed=j)
  p <- subplex(par = c(8,0.175,0.9),fn = llik,n = st$n, E = st$E, t = st$t)$par
  RP[j,] = p
  sit = drop.fossil(st$newick)
  sit = phylo2p(sit)
  P = matrix(nrow=n_trees,ncol=3)
  trees = list()
  for (i in 1:n_trees){
    rec = rec_tree(w=sit$t,pars=p)
    setofweith[i] = rec$prob
    p2 <- subplex(par = c(8,0.175,0.9),fn = llik,n = rec$n, E = rec$E, t = rec$wt)$par
    trees[[i]] = rec
    P[i,] = p2
  }

  #####
  MPf[j,] = colMeans(P)
  pars = subplex(par = c(8,0.175,0.9),fn = llik_st , setoftrees = trees, impsam = TRUE, weith=setofweith)$par
  MP[j,] = pars
}


par_est_vis(P=MP,par=1,PR=RP)
