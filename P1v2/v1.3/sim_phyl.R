### Phylogenetic tree simulation
sim_phyl <- function(ct=15, lambda0=0.8, mu0=0.1, K=40, draw=TRUE, model="dd",printEv=FALSE,seed=7){
  ## Set up
  set.seed(seed)
  i = 1
  N = 1 # Starting number of species
  Tm = NULL # Waiting times
  E = NULL # vector with 0 if extinction and 1 if speciation
  n = NULL # vector with number of species at time t_i
  sumt = 0 # time in simulation
  reboot = 0 # this is in case we want to check how many reboots the simulation had.
  newick = paste(sl[1],";",sep="")  # Newick tree
  identf = data.frame(Spec="aa",Time=0) # Labels of species
  while (sumt<ct){
    if(model == "dd"){  # diversity-dependence model
      lambda = max(0,lambda0 - (lambda0-mu0)*N/K)
      mu = mu0
      lambda = rep(lambda,N)
      mu = rep(mu,N)
    }
    if(model == 'cr'){ # constant-rate model
      lambda = rep(lambda0,N)
      mu = rep(mu0,N)
    }
    s = sum(lambda)+sum(mu)
    if (s == 0){break}
    tm = rexp(1,s)  # waiting time of iteration i
    if(tm+sumt>ct){break}
    sumt = tm + sumt
    prob = c(lambda,mu)/s  # Probability of extinctions and speciations
    BD = sample(2*N,1,prob=prob)  # speciation/extinction & identification of the species.
    n[i] = N
    if(BD > N){   # Extinction
      E[i] = 0
      ## for newick output
      species = identf[BD-N,1]
      ind = regexpr(species,newick)[1] + 2
      atm = sumt-identf[which(identf[,1]==species),2]
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
      E = NULL # vector with 0 if extinction and 1 if speciation
      n = NULL # vector with number of species at time t_i
      newick = paste(sl[1],";",sep="")  # Newick tree
      identf = data.frame(Spec="aa",Time=0)
    }else { # Otherwise, update values and go to next iteration
      Tm[i] = tm
      i<-i+1
    }
  }
  newick = compphyl(newi=newick,identf=identf,sumt=sumt)
  newick = read.tree(text=newick)
  return(list(t=Tm, E=E, i=i, n=n, newick=newick, br = cumsum(Tm)))
}
