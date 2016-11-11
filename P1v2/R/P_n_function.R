P_2 <- function(n=2,t_0=0.7,t_1=0.1,par=c(0.8,0.0175,0.1)){
  lambda = par[1]-par[2]*1:n
  mu = rep(par[3],n)
  # here it shoul be a loop where I can create n functions instead of one, something like
  #f_i = function(x){exp(x*(lambda[i]-lambda[1]))}
  #g_i = function(x){dexp(,mu[i])}
  f_1 = function(x){exp(x*(lambda[2]-lambda[1]))}
  g_1 = function(x){dexp(x,mu[1])}
  #this should be a generic function like f = function(x) {f_1(x)*g_1(x)*f_2(x)*g_2(x)*...*f_n(x)*g_n(x)}
  f = function(x) {f_1(x)*g_1(x)}
  result = integrate(f = f, lower = 0, upper = t_0)$value
  return(result)
}

## list of functions
#x <- sapply(1:10, function(i) function(x) exp(x*i) ); x[[2]](3)

P_3 <- function(n=3,t_0=0.7,t_1=0.1,par=c(0.8,0.0175,0.1)){
  lambda = par[1]-par[2]*1:n
  mu = rep(par[3],n)
  # here it shoul be a loop where I can create n functions instead of one, something like
  #f_i = function(x){exp(x*(lambda[i]-lambda[1]))}
  #g_i = function(x){dexp(,mu[i])}
  f_1 = function(x1){exp(x1*(lambda[2]-lambda[1]))}
  g_1 = function(x1){dexp(x1,mu[1])}
  f_2 = function(x2){exp(x2*(lambda[3]-lambda[1]))}
  g_2 = function(x1,x2){dexp(sum(x1,x2),mu[2])}
  #this should be a generic function like f = function(x) {f_1(x)*g_1(x)*f_2(x)*g_2(x)*...*f_n(x)*g_n(x)}
  f_comp = function(x2) {f_1(x1)*g_1(x1)*f_2(x2)*g_2(x1,x2)}
  #here a nested integral should be created, I guess with a loop?
  f_i1 = function(x1) {integrate(f_comp,upper = t_0, lower = x1)}
  result = integrate(f = f_1, lower = 0, upper = t_0)$value
  return(result)
}

## example for stackoverflow
  t_0 = 15
  mu = 0.1
  lambda = 0.8
  f = function(x1,x2) exp(mu*(x1+x2))*dexp(log(lambda)*(x1+x2))
  f_comp = function(x2) f(x1,x2)
  f_1 = function(x1) {integrate(f_comp,upper = t_0, lower = x1)}
  result = integrate(f = f_1, lower = 0, upper = t_0)$value
  result

integrate(function(x1) {
  sapply(x1, function(x1){
    integrate(function(x2)  exp(mu*(x1+x2))*pexp(log(lambda)*(x1+x2)), lower = x1, upper = t_0)$value
  })
}, 0, t_0)



integrate(function(y) {
       sapply(y, function(y) {
           integrate(function(x) tan(x + y), -.5, y)$value
         })
     }, 0, .5)




#####

if g(x, y) can only take scalar input, we can do

> integrate(function(y) {
  +   sapply(y, function(y) {
    +     integrate(function(x) {
      +       sapply(x, function(x) tan(x + y))
      +     }, -.5, y)$value
    +   })
  + }, 0, .5)
0.07238855 with absolute error < 1.1e-15
