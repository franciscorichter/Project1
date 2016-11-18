phylo2p <- function(tree){
  # to map newick trees into ther xxxx format
  ltt = ltt.plot.coords(tree)
  t = diff(ltt[,1])
  ltt = ltt[-1,]
  n = ltt[,2]
  E = diff(n)
  n = n[-length(n)]
  E[E==-1]=0
  return(list(t=t,E=E,n=n))
}
