library(mvtnorm)

#### Common Copulas ##
Independence = Copula(density = function(u,v){1})
AliMikhailHaq = Copula(copula = function(u,v){
  theta = -0.99
  u*v/(1-theta*(1-u)*(1-v))})
Clayton = Copula(copula = function(u,v){
  theta  = 3

  A = u^(-theta) + v^(-theta) - 1

  (A*(A>0))^(-1/theta)

})
Normal = Copula(copula = function(u,v){
  apply(X = cbind(qnorm(u),qnorm(v)),MARGIN = 1,FUN =  pmvnorm, lower = -Inf, sigma = cbind(c(1,0.9),c(0.9,1)))
})
