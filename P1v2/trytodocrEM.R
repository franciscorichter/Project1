#generate phylogenetic tree under dd model
phyl <- function(tt=15, lambda0=0.8, mu0=0.1, K=40, draw=TRUE, model="dd",printEv=FALSE,seed=1){
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
    if(model == "dd"){  #diversity-dependence model
      lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
      mu = mu0
      lambda = rep(lambda,N)
      mu = rep(mu,N)
    }
    if(model == 'cr'){ #constant-rate model
      lambda = rep(lambda0,N)
      mu = rep(mu0,N)
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
      i<-i+1
    }
  }
  vals = data.frame(time=cumsum(Tm),n=n)
  newick = compphyl(newi=newick,identf=identf,sumt=sumt)
  newick = read.tree(text=newick)
  treeD = list(t=Tm, E=E, r=reboot, i=i, n=n, vals=vals, newick=newick)
}


rec_tree_cr <- function(bt, pars, seed = F, ct){
  if(seed) set.seed(seed)
  cum_bt = c(0,cumsum(bt))
  n_bt = length(bt)
  E = rep(1,n_bt)
  N = 1:(n_bt+1)
  o_bt = bt # original bt
  lambda = pars[1]
  mu = pars[2]
  for (i in 1:n_bt){
    period <- o_bt[i]
    n = rpois(1,lambda*N[i])
    if (n != 0){
      speciations = sort(runif(n, min=0, max=1))
      intervals = diff(c(0,speciations))*period
      speciations <- cum_bt[i] + cumsum(intervals)
      extinctions = rexp(n,mu)
      extinctions <- extinctions + cum_bt[i] + speciations
      acept = extinctions < ct
      set = 1:n
      for (j in set[acept]){
        up = update_tree_cr(bt=bt, t_spe=speciations[j], t_ext=extinctions[j], E=E, Nu=N, ct=ct)
        E = up$E
        N = up$Nu
        bt = up$bt
      }
    }
  }
  return(list(bt = bt,N = N,E = E))
}

update_tree_cr <- function(bt, t_spe, t_ext, ct, last_bt, E, Nu){
  #adding speciation
  bt = c(bt,ct)
  k = length(bt[cumsum(bt)<t_spe])
  bt = c(bt[0:k], t_spe-sum(bt[0:k]), bt[k+1]-(t_spe-sum(bt[0:k])), bt[(k+2):length(bt)])
  E = c(E[0:k],1,E[(k+1):length(E)])
  Nu = c(Nu[0:k],Nu[k:length(Nu)]+1)
  #adding extinction
  k = length(bt[cumsum(bt)<t_ext])
  E = c(E[0:k],0,E[(k+1):length(E)])
  Nu = c(Nu[0:k], Nu[k:length(Nu)] - 1)
  bt = c(bt[0:k], t_ext-sum(bt[0:k]), bt[k+1]-(t_ext-sum(bt[0:k])), bt[(k+2):length(bt)])
  return(list(bt=bt,E=E,Nu=Nu))
}


########
reconst_tree <- function(bt,pars,model="dd",seed = F){
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


  sl = paste(letters[1],letters,":0",sep="")
  for (i in 2:26){
    ll = paste(letters[i],letters,":0",sep="")
    sl = c(sl,ll)
  }

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


