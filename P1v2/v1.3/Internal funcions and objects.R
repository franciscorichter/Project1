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
update_tree <- function(wt, t_spe, t_ext, E, Nu){
  #adding speciation
  ct = sum(wt)
  k = length(wt[cumsum(wt)<t_spe])
  wt = c(wt[0:k], t_spe-sum(wt[0:k]), wt[k+1]-(t_spe-sum(wt[0:k])), wt[(k+2):length(wt)])
  E = c(E[0:k],1,E[(k+1):length(E)])
  Nu = c(Nu[0:(k+1)],Nu[(k+1):length(Nu)]+1)
  print(Nu)
  #adding extinction
  k = length(wt[cumsum(wt)<t_ext])
  print(k)
  if(k<length(E)){
    lastbitE = E[(k+1):length(E)]
    lastbitN = c(Nu[k+1],Nu[(k+1):length(Nu)]-1)
    lastbitt = wt[(k+2):length(wt)]
  }else{
    lastbitE = NULL
    lastbitN = Nu[k]-1
    lastbitt = NULL
  }
  E = c(E[0:k],0,lastbitE)
  Nu = c(Nu[0:k],lastbitN)
  wt = c(wt[0:k], t_ext-sum(wt[0:k]), wt[k+1]-(t_ext-sum(wt[0:k])), lastbitt)
  print(Nu)
  return(list(wt=wt,E=E,Nu=Nu))
}
