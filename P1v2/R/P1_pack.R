#generate phylogenetic tree under dd model
phyl2 <- function(tt=15, lambda0=0.8, mu0=0.1, K=40, draw=TRUE, model="dd",printEv=FALSE,seed=1){
  set.seed(seed)
  reboot = 0
  N = 1 # Number of species
  i = 1
  Tm = NULL
  sumt = 0
  sigma = 0
  E = NULL # vector with 0 if extinction and 1 if speciation
  n = NULL # vector with number of species at time t_i
  newick = paste(sl[1],";",sep="")  # Newick tree
  identf = data.frame(Spec="aa",Time=0)
  while (sumt<tt){
    if (model == "dd"){  #diversity-dependence model
      lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
      mu = mu0
      lambda = rep(lambda,N)
      mu = rep(mu,N)
    }
    s = sum(lambda)+sum(mu)
    if (s == 0){break}
    tm = rexp(1,s)  # waiting time of iteration i
    if(tm+sumt>tt){break}
    sumt = tm + sumt
    prob = c(lambda,mu)/s  # Probability of extinctions and speciations
    BD = sample(2*N,1,prob=prob)  # speciation/extinction & identification of the species.
    n[i] = N
    if(BD > N){   # Extinction
      E[i] = 0
      ## for newick output
      species = identf[BD-N,1]
      ind = regexpr(species,newick)[1] + 2
      atm=sumt-identf[which(identf[,1]==species),2]
      identf = identf[-(BD-N),]
      newick = paste(substr(newick,1,ind),as.character(atm),substring(newick,ind+2),sep="")
      #
      N = N-1
      if(printEv){print(paste("extinction in time",sumt, sep=" "))}
    }else{  # Speciation
      E[i] = 1
      ## for newick output
      species = as.character(identf[BD,1])
      ind = regexpr(species,newick)[1]-1
      atm=sumt-identf[which(identf[,1]==species),2]
      newick = paste(substr(newick,1,ind),"(",substr(newick,ind+1,ind+4),",",sl[i+1],"):",as.character(atm),substring(newick,ind+5),sep="")
      identf = rbind(identf,data.frame(Spec=substr(sl[i+1],1,2),Time=sumt))
      identf[identf$Spec == species,2] = sumt
      #
      N = N+1
      if(printEv){print(paste("speciation in time",sumt,sep=" "))}
    }
    if (N==0){ # In case all species got extinct: restart
      reboot = reboot + 1
      N = 1 # Number of species
      i = 1
      Tm = NULL
      sumt = 0
      sigma = 0
      E = NULL # vector with 0 if extinction and 1 if speciation
      n = NULL # vector with number of species at time t_i
      newick = paste(sl[1],";",sep="")  # Newick tree
      identf = data.frame(Spec="aa",Time=0)
    }else { # Otherwise, update values and go to next iteration

      Tm[i] = tm
      sigma[i] = s  # sirve?
      i<-i+1
    }
  }
  vals = data.frame(time=cumsum(Tm),n=n)
  newick = compphyl(newi=newick,identf=identf,sumt=sumt)
  newick = read.tree(text=newick)
  treeD = list(t=Tm, E=E, r=reboot, i=i, n=n, vals=vals, newick=newick)
}



## truncated exponential distribution
itexp <- function(x, lambda, R) { -log(1-x*(1-exp(-R*lambda)))/lambda }
rtexp <- function(n, lambda, R) { itexp(runif(n), lambda, R) }

# reconstruct tree

reconst_tree6 <- function(bt,pars,model="dd",seed = F){
  #bt = c(bt,tt-sum(bt))
  bt = c(0,bt)
  tp = sum(bt)
  n_bt = length(bt)
  E = rep(1,n_bt)
  Nu = 2:(n_bt+1)
  lambda0 = pars[1]
  mu0 = pars[3]
  K = (lambda0-mu0)/pars[2]
  i=1
  t_fake = 0
  tm_ext = 0
  while (i < length(bt)){
    if (seed) set.seed(i)
    N = Nu[i] #esta no se necesita si se escribe asi abajo
    if (model == "dd"){  #diversity-dependence model
      lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
      lambda = rep(lambda,N)
      mu = rep(mu0,N)
    }
    s = sum(lambda)+sum(mu)
    if (s == 0){
      #print('S=0!')
      break}
    ts = rexp(1,s)
    tm=ts+t_fake
    if (tm < bt[i+1]){
      tm_ext = rexp(1,mu0)
      if ((sum(bt[1:i])+tm+tm_ext) >= tp){
        #print(paste('fake species at', (sum(bt[1:i])+tm)))
        t_fake = tm
      }else{t_fake = 0}
    }else{t_fake=0}
    while(tm < bt[i+1] & (sum(bt[1:i])+tm+tm_ext) < tp){
      #print(paste('new spec-ext at times',sum(bt[1:i])+tm,sum(bt[1:i])+tm+tm_ext))
      up = update_tree(bt=bt, t_spe=tm, t_ext=tm_ext, pointer=i, E=E, Nu=Nu)
      i = i+1
      E = up$E
      #print(Nu)
      Nu = up$Nu
      bt = up$bt
      N = Nu[i]
      if (model == "dd"){  #diversity-dependence model
        lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
        lambda = rep(lambda,N)
        mu = rep(mu0,N)
      }
      s = sum(lambda)+sum(mu)
      if (s == 0){
        #print('S=0!')
        break}
      ts = rexp(1,s)
      tm=ts+t_fake
      if (tm < bt[i+1]){
        tm_ext = rexp(1,mu0)
        if ((sum(bt[1:i])+tm_ext) >= tp){
          # print(paste('fake species at', sum(bt[1:i])+tm))
          t_fake = tm
        }else {t_fake=0}
      }
    }
    if(t_fake == 0) i=i+1
    #print(paste('nbranch=',length(bt),'iter',i))
  }
  return(list(t=bt,n=Nu,E=E))
}


#####
reconst_tree2 <- function(bt,pars,model="dd",seed = F){
  #bt = c(bt,tt-sum(bt))
  bt = c(0,bt)
  tp = sum(bt)
  n_bt = length(bt)
  E = rep(1,n_bt)
  #Nu = c(2:n_bt,n_bt)
  Nu = 2:(n_bt+1)
  lambda0 = pars[1]
  mu0 = pars[3]
  K = (lambda0-mu0)/pars[2]
  i=1 # this is related to all species
  t_fake = 0
  tm_ext = 0
  while (i < length(bt)){
    if (seed) set.seed(i)
    N = Nu[i] #esta no se necesita si se escribe asi abajo
    if (model == "dd"){  #diversity-dependence model
      lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
      lambda = rep(lambda,N)
      mu = rep(mu0,N)
    }
    s = sum(lambda)+sum(mu)
    if (s == 0){
      #print('S=0!')
      break}
    ts = rexp(1,s)
    tm=ts+t_fake
    if (tm < bt[i+1]){
      tm_ext = rexp(1,mu0)
      if ((sum(bt[1:i])+tm+tm_ext) >= tp){
        #print(paste('fake species at', (sum(bt[1:i])+tm)))
        t_fake = tm
      }else{t_fake = 0}
    }else{t_fake=0}
    while(tm < bt[i+1] & (sum(bt[1:i])+tm+tm_ext) < tp){
      #print(paste('new spec-ext at times',sum(bt[1:i])+tm,sum(bt[1:i])+tm+tm_ext))
      up = update_tree(bt=bt, t_spe=tm, t_ext=tm_ext, pointer=i, E=E, Nu=Nu)
      i = i+1
      E = up$E
      #print(Nu)
      Nu = up$Nu
      bt = up$bt
      N = Nu[i]
      if (model == "dd"){  #diversity-dependence model
        lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
        lambda = rep(lambda,N)
        mu = rep(mu0,N)
      }
      s = sum(lambda)+sum(mu)
      if (s == 0){
        #print('S=0!')
        break}
      ts = rexp(1,s)
      tm=ts+t_fake
      if (tm < bt[i+1]){
        tm_ext = rexp(1,mu0)
        if ((sum(bt[1:i])+tm_ext) >= tp){
         # print(paste('fake species at', sum(bt[1:i])+tm))
          t_fake = tm
        }else {t_fake=0}
      }
    }
    if(t_fake == 0) i=i+1
    #print(paste('nbranch=',length(bt),'iter',i))
  }
  return(list(t=bt,n=Nu,E=E))
}

#mle for a set of trees
mle_dd_setoftrees <- function(setoftrees,draw=T){
  num = NULL
  dem = NULL
  for (i in 1:length(setoftrees)){
    s = setoftrees[[i]]
    num[i] = sum(1-s$E)
    dem[i] = sum(s$t*s$n)
  }
  mu = sum(num)/sum(dem)
  lambda=seq(0.01,2,by=0.02)
  grid = length(lambda)
  y = NULL
  Beta = NULL
  g = NULL
  for (i in 1:grid){
    beta = seq(0,0.025 ,by=0.0005) #this 40 need to be changed
    #beta = beta[2:length(beta)]
    for (j in 1:length(beta)){
      y[j] = llik_st(setoftrees,pars=c(lambda[i],beta[j],mu))
    }
    g[i] = min(y[!is.na(y)])
    #g[i] = optim(0,llik,)
    Beta[i] = beta[y==min(y[!is.na(y)])]

  }
  plot(lambda,Beta)
  plot(lambda,g)
  l = lambda
  lambda = lambda[g==min(g)]
  beta = Beta[g==min(g)]
  return(list(lambda=lambda,beta=beta,mu=mu,setoflambda=l,setofbeta=Beta))
}

reconst_tree <-function(bt,pars,model="dd",seed = F){
  #DUPLICATED, EQUIVALENT TO rt5
  bt = c(0,bt)
  tp = sum(bt)
  n_bt = length(bt)
  E = rep(1,n_bt)
  Nu = 2:(n_bt+1)
  lambda0 = pars[1]
  mu0 = pars[3]
  K = (lambda0-mu0)/pars[2]
  i=1
  # t_fake = 0
  tm_ext = 0
  while (i < length(bt)){
    if (seed) set.seed(i)
    N = Nu[i] #esta no se necesita si se escribe asi abajo
    if (model == "dd"){  #diversity-dependence model
      lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
      lambda = rep(lambda,N)
    }
    s = sum(lambda)
    if (s == 0){
      break}
    tm = rexp(1,s)
    if (tm < bt[i+1]){tm_ext = rexp(1,mu0)}
    while(tm < bt[i+1] & (sum(bt[1:i])+tm+tm_ext) < tp){
      #print(paste('new spec-ext at times',sum(bt[1:i])+tm,sum(bt[1:i])+tm+tm_ext))
      up = update_tree(bt=bt, t_spe=tm, t_ext=tm_ext, pointer=i, E=E, Nu=Nu)
      i = i+1
      E = up$E
      Nu = up$Nu
      bt = up$bt
      N = Nu[i]
      if (model == "dd"){  #diversity-dependence model
        lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
        lambda = rep(lambda,N)
      }
      s = sum(lambda)
      if (s == 0){
        break}
      tm = rexp(1,s)
      if (tm < bt[i+1]){tm_ext = rexp(1,mu0)}
    }
    i=i+1
    #print(paste('nbranch=',length(bt),'iter',i))
  }
  return(list(t=bt,n=Nu,E=E))
}

reconst_tree3 <- function(bt,pars,model="dd",seed = F){
  ##DUPLICATED, EQUIVALENT TO rt4.
  #bt = c(bt,tt-sum(bt))
  bt = c(0,bt)
  tp = sum(bt)
  n_bt = length(bt)
  E = rep(1,n_bt)
  #Nu = c(2:n_bt,n_bt)
  Nu = 2:(n_bt+1)
  lambda0 = pars[1]
  mu0 = pars[3]
  K = (lambda0-mu0)/pars[2]
  i=1 # this is related to all species
  tm_ext = 0
  while (i < length(bt)){
    if (seed) set.seed(i)
    N = Nu[i] #esta no se necesita si se escribe asi abajo
    if (model == "dd"){  #diversity-dependence model
      lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
      lambda = rep(lambda,N)
      mu = rep(pars[3],N)
    }
    s = sum(lambda)+sum(mu)
    if (s == 0){
      #print('S=0!')
      break}
    ts = rexp(1,s)
    tm=ts
    prob = c(lambda,mu)/s  # Probability of extinctions and speciations
    BD = sample(2*N,1,prob=prob)
    if (tm < bt[i+1] & BD < (N+1)){
      tm_ext = rexp(1,mu0)
    }
    while(tm < bt[i+1] & (sum(bt[1:i])+tm+tm_ext) < tp & BD < (N+1)){
      #print(paste('new spec-ext at times',sum(bt[1:i])+tm,sum(bt[1:i])+tm+tm_ext))
      up = update_tree(bt=bt, t_spe=tm, t_ext=tm_ext, pointer=i, E=E, Nu=Nu)
      i = i+1
      E = up$E
      #print(Nu)
      Nu = up$Nu
      bt = up$bt
      N = Nu[i]
      if (model == "dd"){  #diversity-dependence model
        lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
        lambda = rep(lambda,N)
        mu = rep(pars[3],N)
      }
      s = sum(lambda)+sum(mu)
      if (s == 0){
        #print('S=0!')
        break}
      ts = rexp(1,s)
      tm=ts
      prob = c(lambda,mu)/s  # Probability of extinctions and speciations
      BD = sample(2*N,1,prob=prob)
      if (tm < bt[i+1] & BD < (N+1)){
        tm_ext = rexp(1,mu0)
      }
    }
    i=i+1
    #print(paste('nbranch=',length(bt),'iter',i))
  }
  return(list(t=bt,n=Nu,E=E))
}


update_tree <- function(bt, t_spe, t_ext, pointer, E, Nu, B){
  #TODO add newick format to output
  i = pointer
  #adding speciation
  if(length(bt)>(i+1)){
    last_bit = bt[(i+2):length(bt)]
  }
  else{
    last_bit = NULL
  }
  bt = c(bt[1:i],t_spe,bt[(i+1)]-t_spe,last_bit)
  ## ACA VOY
  B = c(B[1:i],1,0,B[(i+2):length(B)])
  E = c(E[1:i],1,E[(i+1):length(E)])
  Nu = c(Nu[1:i],Nu[i:length(Nu)]+1)
  #adding extinction
  n_ext_cum = sum(bt[1:i+1]) + t_ext
  i = i + 1
  E = c(E[cumsum(bt)<n_ext_cum],0,E[cumsum(bt)>n_ext_cum])
  B = c(B[cumsum(bt)<n_ext_cum],1,0,B[(length(B[cumsum(bt)<n_ext_cum])+2):length(B)])
  Nu = c(Nu[cumsum(bt)<n_ext_cum],Nu[length(Nu[cumsum(bt)<n_ext_cum])]-1,Nu[cumsum(bt)>n_ext_cum] - 1)
  if(length(bt[cumsum(bt)>n_ext_cum])>1){
    last_bit = bt[cumsum(bt)>n_ext_cum][2:length(bt[cumsum(bt)>n_ext_cum])]
  }
  else{
    last_bit = NULL
  }
  bt = c(bt[cumsum(bt)<n_ext_cum],n_ext_cum-sum(bt[cumsum(bt)<n_ext_cum]),bt[cumsum(bt)>n_ext_cum][1]-(n_ext_cum-sum(bt[cumsum(bt)<n_ext_cum])),last_bit) #including the new extinction time
  return(list(bt=bt,E=E,Nu=Nu,B=B))
}

llik_st = function(pars,setoftrees){
  m = length(setoftrees)
  l = NULL
  for(i in 1:m){
    s = setoftrees[[i]]
    l[i] = llik(b=pars,n=s$n,E=s$E,t=s$t)
  }
  L = sum(l)
  return(L)
}



llik = function(b,n,E,t){
  sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
  rho = pmax(b[1]*E-b[2]*n*E+b[3]*(1-E),0)
  #if(min(rho)<0){print(paste(min(rho),'help'))}
  l = -sum(-sigma*t+log(rho))
  if(min(b)<0){l = Inf}
  #print(l)
  #print(b)
  return(l)
}

cond1 <- function(b){
  return(b[3]>0)
}

cond2 <- function(b,n){
  if (b[2]>0){return(b[1]/b[2]>max(n))}
  if (b[2]<0){return(b[1]/b[2]<min(n))}
  if (b[2]==0) {return(FALSE)}
}


mle_dd <- function(s,draw=T){
  n=s$n
  E=s$E
  t=s$t
  u=sum(1-E)/(sum(n*t))
  lambda=seq(0.01,2,by=0.01)
  grid = length(lambda)
  y = NULL
  Beta = NULL
  g = NULL
  for (i in 1:grid){
    beta = seq(0,lambda[i]/max(n),by=0.0001)
    beta = beta[2:length(beta)]
    for (j in 1:length(beta)){
      y[j] = llik(c(lambda[i],beta[j],u),n,E,t)
    }
    g[i] = min(y[!is.na(y)])
    #g[i] = optim(0,llik,)
    Beta[i] = beta[y==min(y[!is.na(y)])]
  }
  if(draw) plot(lambda,Beta)
  #plot(lambda,g)
  l = lambda
  lambda = lambda[g==min(g)]
  beta = Beta[g==min(g)]
  return(list(lambda=lambda,beta=beta,mu=u,setoflambda=l,setofbeta=Beta))
}

L2p <- function(L,ct=15){
  rr = L
  if (length(L[,1])<5){
    # message of error
  }
  spec = data.frame(t=(ct-L[3:length(L[,4]),1]),E=1)
  if(length((ct-L[L[,4]!=(-1),4]))>0){
    ext = data.frame(t=(ct-L[L[,4]!=(-1),4]),E=0)
    tree = rbind(spec,ext)
  } else{
    tree = spec}
  tree = tree[order(tree$t),]
  dt = diff(tree$t)
  tree$t = c(tree$t[1],dt)
  n=NULL
  n=1
  for (i in 1:dim(tree)[1]){
    n[i+1] = ifelse(tree$E[i]==1,n[i]+1,n[i]-1)
  }
  n = n[1:(length(n)-1)]
  tr = cbind(tree,n)
  n = tr$n
  E = tr$E
  t=tr$t
  return(list(n=n,E=E,t=t))
}

dd_par <- function(tt=15,mu=0.1,it=10){
  foreach(h = 1:it,
          .combine = list,
          .multicombine = TRUE) %dopar% {
            s1 = phyl(tt=tt, mu0=mu,printEv=FALSE)
            n = s1$n
            E = s1$E
            t=s1$t
            mle = mle_dd(n,E,t)
            mle
          }
}

###


# get_b <- function(b_est,b,n){
#   g=0.5
#   if (cond1(b_est) & cond2(b_est,n)){
#     g = 1
#     bb = b_est
#     #print("a")
#   }
#   if (!cond1(b_est) & cond2(b_est,n)){
#     gamma = b_est[3]/(b_est[3]-b[3])
#     bb = b*gamma + b_est*(1-gamma)
#     #print( "b")
#   }
#   if (cond1(b_est) & !cond2(b_est,n)){
#     gamma = get_gammac2(b_est = b_est, b=b, n=n)
#     bb =  b*gamma + b_est*(1-gamma)
#     #print("c")
#   }
#   if (!cond1(b_est) & !cond2(b_est,n)){
#     gamma = max(b_est[3]/(b_est[3]-b[3]),get_gammac2(b_est = b_est, b=b, n=n))
#     bb = b*gamma + b_est*(1-gamma)
#     #print("d")
#   }
#   bp = bb*g + b*(1-g)  # potential b
#   g = 0.5
#   i=1
#   while (llik(b,n,E,t)<llik(bp,n,E,t) | i>100){
#     bp = bp*g + b*(1-g)
#     i=i+1}
#   ##
#   b = bp
#   return(b)
# }
#
#
# get_gammac2 <- function(b_est,b,n){
#   if (b_est[2]>0){
#     gamma = (max(n)*b_est[2]-b_est[1])/(b[1]-b_est[1]-max(n)*b[2]+max(n)*b_est[2])
#     bb =  b*gamma + b_est*(1-gamma)
#     if (bb[2]<0){
#       gamma = (min(n)*b_est[2]-b_est[1])/(b[1]-b_est[1]-min(n)*b[2]+min(n)*b_est[2])
#       #bb =  b*gamma + b_est*(1-gamma)
#     }
#   }
#   if (b_est[2]<0){
#     gamma = (min(n)*b_est[2]-b_est[1])/(b[1]-b_est[1]-min(n)*b[2]+min(n)*b_est[2])
#     bb =  b*gamma + b_est*(1-gamma)
#     if (bb[2]>0){
#       gamma = (max(n)*b_est[2]-b_est[1])/(b[1]-b_est[1]-max(n)*b[2]+max(n)*b_est[2])
#       #bb =  b*gamma + b_est*(1-gamma)
#     }
#   }
#   return(gamma)
# }
#
#
# ### Newton rapson estimation
# # incomplete
# # jacobian <- function(E,n,b){
# #   J=matrix(nrow=3,ncol=3)
# #   J[1,1] = sum(-E/((E*b[1]-n*E*b[2]+(1-E)*b[3])^2))
# #   J[1,2] = sum((E*n)/((E*b[1]-n*E*b[2]+(1-E)*b[3])^2))
# #   J[1,3] = 0
# #   J[2,1] = sum((E*n)/((E*b[1]-n*E*b[2]+(1-E)*b[3])^2))
# #   J[2,2] = sum((-E*n*n)/((E*b[1]-n*E*b[2]+(1-E)*b[3])^2))
# #   J[2,3] = 0
# #   J[3,1] = 0
# #   J[3,2] = 0
# #   J[3,3] = sum((E-1)/((E*b[1]-n*E*b[2]+(1-E)*b[3])^2))
# #   return(J)
# # }
# #
# # gradient <- function(b,E,n,t){
# #   G = NULL
# #   coe = (b[1]*E-b[2]*n*E+b[3]*(1-E))^2
# #   G[1] = sum(-n*t-E/coe)
# #   G[2] = sum(n*t+(n*E)/coe)
# #   G[3] = sum(-n*t-(1-E)/coe)
# #   return(G)
# # }
#
#
#
# ### GLM estimations
#
# W_bt = function(b,n){
#   sigmas = b[1]*n-b[2]*n*n+n*b[3]
#   invs = 1/(sigmas*sigmas)
#   W = diag(invs)
#   return(W)
# }
#
# W_to = function(b,n){
#   sigmas = b[1]*n-b[2]*n*n+n*b[3]
#   rho = E*(b[1]-b[2]*n)+(1-E)*b[3]
#   W = diag((sigmas-rho)/(rho*sigmas*sigmas))
#   return(W)
# }
#
# W = function(b,n){
#   W = bdiag(W_bt(b,n),W_to(b,n))
#   return(W)
# }
#
# z_bt = function(b,n,t){
#   sigmas = b[1]*n-b[2]*n*n+n*b[3]
#   z = sigmas*(2-sigmas*t)
#   return(z)
# }
#
# z_to = function(b,n){
#   sigmas = b[1]*n-b[2]*n*n+n*b[3]
#   rho = E*(b[1]-b[2]*n)+(1-E)*b[3]
#   z = rho + sigmas
#   return(z)
# }
#
# z = function(b,n,t){
#   z = c(z_bt(b,n,t),z_to(b,n))
#   return(z)
# }

### Parallel
# incompleta, testearla primero

#stopCluster(cl)



### Newick functions
#
sl = paste(letters[1],letters,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],letters,":0",sep="")
  sl = c(sl,ll)
}
#
compphyl <- function(newi,identf,sumt){
  #test if the newick phylo is consistent with the identif matrix
  identf[,1]=as.character(identf[,1])
  identf[,2]=sumt-identf[,2]
  for(i in 1:length(identf[,1])){
    ind = regexpr(identf[i,1],newi)[1] + 2
    newi = paste(substr(newi,1,ind),as.character(identf[i,2]),substring(newi,ind+2),sep="")
  }
  return(newi)
}

phylo2p <- function(tree,ct){
  ltt = ltt.plot.coords(tree)
  n = ltt[,2]
  E = diff(n)
  n = n[2:(length(n)-1)]
  E=E[2:length(E)]
  t = diff(ltt[,1])
  # last t = t[length(t)]
  t = t[1:(length(t)-1)]
  return(list(t=t,E=E,n=n))
}


############

reconst_tree4 <- function(bt, pars, model="dd", seed = F){
  # lambda+mu, no fake species indicator
  if (seed) set.seed(i)
  bt = c(0,bt)
  tp = sum(bt)
  n_bt = length(bt)
  E = rep(1,n_bt)
  N = 2:(n_bt+1)
  lambda0 = pars[1]
  mu0 = pars[3]
  K = (lambda0-mu0)/pars[2]
  i=1
  tm_ext = 0
  while (i < length(bt)){
    if (model == "dd"){  #diversity-dependence model
      lambda = max(0,lambda0 - (lambda0-mu0)*N[i]/K)
      lambda = rep(lambda,N[i])
      mu = rep(pars[3],N[i])
    }
    s = sum(lambda)+sum(mu)
    if (s == 0){
      #print('S=0!')
      break}
    tm = rexp(1,s)
    if (tm < bt[i+1]){
      prob = c(lambda,mu)/s  # Probability of extinctions and speciations
      BD = sample(2*N[i],1,prob=prob)
      if (BD < (N[i]+1)){
        tm_ext = rexp(1,mu0)
        #print('new species')
        if ((sum(bt[1:i])+tm+tm_ext) < tp){
          up = update_tree(bt=bt, t_spe=tm, t_ext=tm_ext, pointer=i, E=E, Nu=N)
          E = up$E
          N = up$Nu
          bt = up$bt
        }
      }#else{print('fue extincion no permitida')}
    }
    i = i+1
  }
  return(list(t=bt,n=N,E=E))
}


reconst_tree5 <- function(bt, pars, model="dd", seed = F){
  # s=sum of lambda only, no fake species indicator
  if (seed) set.seed(i)
  bt = c(0,bt)
  tp = sum(bt)
  n_bt = length(bt)
  E = rep(1,n_bt)
  B = rep(1,n_bt)
  N = 2:(n_bt+1)
  lambda0 = pars[1]
  mu0 = pars[3]
  K = (lambda0-mu0)/pars[2]
  i=1
  tm_ext = 0
  while (i < length(bt)){
    if (model == "dd"){  #diversity-dependence model
      lambda = max(0,lambda0 - (lambda0-mu0)*N[i]/K)
      lambda = rep(lambda,N[i])
      mu = rep(pars[3],N[i])
    }
    s = sum(lambda)#+sum(mu)
    if (s == 0){
      #print('S=0!')
      break}
    tm = rexp(1,s)
    if (tm < bt[i+1]){
        tm_ext = rexp(1,mu0)
        if ((sum(bt[1:i])+tm+tm_ext) < tp){
          up = update_tree(bt=bt, t_spe=tm, t_ext=tm_ext, pointer=i, E=E, Nu=N,B = B)
          B = up$B
          E = up$E
          N = up$Nu
          bt = up$bt
        }
    }
    i = i+1
  }
  return(list(t=bt,n=N,E=E,B=B))
}



    while(tm < bt[i+1] & (sum(bt[1:i])+tm+tm_ext) < tp & BD < (N+1)){
      #print(paste('new spec-ext at times',sum(bt[1:i])+tm,sum(bt[1:i])+tm+tm_ext))
      up = update_tree(bt=bt, t_spe=tm, t_ext=tm_ext, pointer=i, E=E, Nu=Nu)
      i = i+1
      E = up$E
      #print(Nu)
      Nu = up$Nu
      bt = up$bt
      N = Nu[i]
      if (model == "dd"){  #diversity-dependence model
        lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
        lambda = rep(lambda,N)
        mu = rep(pars[3],N)
      }
      s = sum(lambda)+sum(mu)
      if (s == 0){
        #print('S=0!')
        break}
      ts = rexp(1,s)
      tm=ts
      prob = c(lambda,mu)/s  # Probability of extinctions and speciations
      BD = sample(2*N,1,prob=prob)
      if (tm < bt[i+1] & BD < (N+1)){
        tm_ext = rexp(1,mu0)
      }
    }
    i=i+1
    #print(paste('nbranch=',length(bt),'iter',i))
  }
  return(list(t=bt,n=Nu,E=E))
}


