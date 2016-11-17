sl = paste(letters[1],letters,":0",sep="")
for (i in 2:26){
  ll = paste(letters[i],letters,":0",sep="")
  sl = c(sl,ll)
}

compphyl <- function(newi,identf,sumt){
  #set to extant species to the present time
  identf[,1] = as.character(identf[,1])
  identf[,2] = sumt-identf[,2]
  for(i in 1:length(identf[,1])){
    ind = regexpr(identf[i,1],newi)[1] + 2
    newi = paste(substr(newi,1,ind),as.character(identf[i,2]),substring(newi,ind+2),sep="")
  }
  return(newi)
}



update_tree <- function(bt, t_spe, t_ext, ct, last_bt, E, Nu){
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
