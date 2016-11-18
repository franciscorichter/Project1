sl = paste(letters[1],letters,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],letters,":0",sep="")
  sl = c(sl,ll)
}

compphyl <- function(newi,identf,ct){
  #set to extant species to the present time
  identf[,1] = as.character(identf[,1])
  identf[,2] = ct-identf[,2]
  for(i in 1:length(identf[,1])){
    ind = regexpr(identf[i,1],newi)[1] + 2
    newi = paste(substr(newi,1,ind),as.character(identf[i,2]),substring(newi,ind+2),sep="")
  }
  return(newi)
}

# TODO: add newick output to update_tree
# add warning message when spe or ext is begond present time
update_tree <- function(wt, t_spe, t_ext, E, n){
  #adding speciation
  ct = sum(wt)
  k = length(wt[cumsum(wt)<t_spe])
   if(k<length(E)){
     lastbitE = E[(k+1):length(E)]
     lastbitN = c(n[k+1],n[(k+1):length(n)]+1)
     lastbitt = wt[(k+2):length(wt)]
   }else{
     lastbitE = NULL
     lastbitN = c(n[k+1],n[k+1]+1)#n[k]+1
     lastbitt = NULL
   }
  wt = c(wt[0:k], t_spe-sum(wt[0:k]), wt[k+1]-(t_spe-sum(wt[0:k])), lastbitt)
  E = c(E[0:k],1,lastbitE)
  n = c(n[0:k],lastbitN)
  #adding extinction
  k = length(wt[cumsum(wt)<t_ext])
  #print(k)
  if(k<length(E)){
    lastbitE = E[(k+1):length(E)]
    lastbitN = c(n[k+1],n[(k+1):length(n)]-1)
    lastbitt = wt[(k+2):length(wt)]
  }else{
    lastbitE = NULL
    lastbitN = c(n[k+1],n[(k+1):length(n)]-1)#n[k]-1
    lastbitt = NULL
  }
  E = c(E[0:k],0,lastbitE)
  n = c(n[0:k],lastbitN)
  wt = c(wt[0:k], t_ext-sum(wt[0:k]), wt[k+1]-(t_ext-sum(wt[0:k])), lastbitt)
  #print(n)
  return(list(wt=wt,E=E,n=n))
}



###
llik_st = function(pars,setoftrees,impsam = F,weith=0){
  m = length(setoftrees)
  l = NULL
  for(i in 1:m){
    s = setoftrees[[i]]
    if(impsam){
      l[i] = llik(b=pars,n=s$n,E=s$E,t=s$wt)#*llik(b=pars,n=s$n,E=s$E,t=s$wt)/weith
    }else{
      l[i] = llik(b=pars,n=s$n,E=s$E,t=s$wt)
    }
  }
  L = sum(l)
  return(L)
}


par_est_vis <- function(P,par,PR){
  #P is the recost estim values
  # par is the parameter you want to use
  #PR is the real estim
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


