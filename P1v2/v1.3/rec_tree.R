rec_tree <- function(wt, pars, ct){
  lambda0 = pars[1]
  mu0 = pars[3]
  K = (lambda0-mu0)/pars[2]
  n = 1:length(wt)
  i = 1
  fake = FALSE
  while(i < length(wt)){
    N = n[i]
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
    s = sum(lambda)
    if(fake){ # in case there was an speciation but not extinction previously
      cwt = cwt - t_spe
      cbt = cbt + t_spe
    }
    else{
      cwt = wt[i]
      cbt = cumsum(wt)[i]
    }
    t_spe = rexp(1,s)
    if (t_spe < cwt){
      t_ext = rexp(1,mu0) # this is not as general as trees with trait-dependance species yet,
      up = update_tree(wt=wt,t_spe = t_spe, t_ext = t_ext, E = ...)
    }
  }

}
