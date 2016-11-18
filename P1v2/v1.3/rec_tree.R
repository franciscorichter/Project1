rec_tree <- function(wt, pars=c(0.8,0.0175,0.1), model='dd'){
  lambda0 = pars[1]
  mu0 = pars[3]
  K = (lambda0-mu0)/pars[2]
  n = 1:length(wt)
  i = 1
  E = rep(1,(length(wt)-1))
  fake = FALSE
  ct = sum(wt)
  prob = 1
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
    if(s==0){
      #print('s=0')
      break
    }
    if(fake){ # in case there was an speciation but not extinction previously
      cwt = cwt - t_spe
      cbt = cbt + t_spe
      #print(paste('fake',cwt))
    }
    else{
      cwt = wt[i]
      cbt = cumsum(wt)[i] - cwt
      #print(paste('nofake',cwt))
    }
    t_spe = rexp(1,s)
    if (t_spe < cwt){
      t_ext = rexp(1,mu0) # this is not as general as trees with trait-dependance species yet,
      t_ext = cbt + t_spe + t_ext
      rt =  ct - #remaining time
      prob = prob + log(1-exp(s*cwt))
      if (t_ext < ct){
        up = update_tree(wt=wt,t_spe = (cbt + t_spe), t_ext = t_ext, E = E, n = n)
        E = up$E
        n = up$n
        wt = up$wt
        fake = FALSE
        prob = prob + log(1-exp(mu0*(ct-t_spe-cbt)))
      }else{
        prob = prob + mu0*(ct-t_spe-cbt)
        fake = TRUE
        i = i-1
      }
    }else{
      fake = FALSE
      prob = prob + s*cwt
    }
    i = i+1
  }
  return(list(wt=wt,E=E,n=n,prob=prob))

}
