#TODO: create the general one
llik = function(b,n,E,t){
  sigma = n*(b[1]-b[2]*n + b[3]) #n-dimentional
  rho = pmax(b[1]*E-b[2]*n*E+b[3]*(1-E),0)
  l = -sum(-sigma*t+log(rho))
  if(min(b)<0){l = Inf}
  return(l)
}
