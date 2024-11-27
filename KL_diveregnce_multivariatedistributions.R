library(emdbook) ##for multivariate normal distribution
library(cubature) ##for 4 dimensional integral 
library(compositions) ##for log normal 

start_time<- Sys.time()
print(start_time)



##kullback leibler divergence between 4 dimensional log normal and normal distributed data
f_multi_lognormal<- function(X) {
  x <- c(X[1], X[2], X[3],X[4])
  dlnorm.rplus(x,mu,Sigma) }

#f_multi_lognormal(c(1,2,1,3))

f_multi_normal<- function(X) {
  x <- c(X[1], X[2], X[3],X[4])
  dmvnorm(x, c(5,4,4,5), Sigma, log = FALSE, tol = 1e-06) }

#f_multi_normal(c(1,2,1,3))

##differences between log normal- normal
KL_fun<- function(X) {
  f_multi_lognormal(X)* log(f_multi_lognormal(X)/f_multi_normal(X)) }
KL_fun(c(2,1,1,3))
int1<- adaptIntegrate(KL_fun, lowerLimit = c(5,4,4,5), upperLimit = c(40,30,30,30))
KL<- int1$integral

##differences between normal- log normal
KL_fun2<- function(X) {
  f_multi_normal(X)* log(f_multi_normal(X)/f_multi_lognormal(X)) }
KL_fun2(c(2,1,1,3))
int2<- adaptIntegrate(KL_fun2, lowerLimit = c(5,4,4,5), upperLimit = c(40,30,30,30))
KL2<- int2$integral


adaptIntegrate(KL_fun2, lowerLimit = c(5,4,4,5), upperLimit = c(40,30,30,30))



end_time<- Sys.time()
print(end_time)