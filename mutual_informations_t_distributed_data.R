##mutual informations in t distributed data

library(cubature) ##for integral calculation
library(extraDistr) #for bivariate dist
library(mvtnorm)
library(MASS)#for generate multinormal data
library(dplyr)
library(LaplacesDemon)
library(pracma)
library(compositions) ##multiple log normal
#library(Metrics) ##FOR bias
library(devtools) ## local functiom
library(akima) ##for approachfun2d
library(KernSmooth) ##BKDE2D
library(FNN) ##to use mutinfo function

start_time<- Sys.time()
print(start_time)

##mutual informations
#mutual_informations<- list()

##generate data rho = medium
#for (j in 5:6) {
  
  #set.seed(j)
  nu<-3
  mu<- c(1.2, 2, 1.5, 3)
  Sigma<- rbind(c(1.0780992, -0.3676667,  0.2402575, -0.4867786),
                c(-0.3676667,  1.0623604, -0.4322828,  0.3046945),
                c(0.2402575 ,-0.4322828 , 1.0594279, -0.2495095),
                c(-0.4867786 , 0.3046945, -0.2495095,  1.0704080))
  
  data_t <- rmvt(1e5, mu, Sigma, df = nu)
  T0<-data_t[,1]
  T1<-data_t[,2]
  S0<-data_t[,3]
  S1<-data_t[,4]
  
  deltaT <- T1-T0
  deltaS <- S1-S0
  
  rho<- cor(deltaS,deltaT)
  
  ##new formula
  zeta_fun <- function(nu) 2*log( (gamma(nu/2)*gamma(1/2))/(sqrt(pi)*gamma( (1+nu)/2 ))*sqrt(nu/2))- (nu/2)*psi(nu/2)+(nu+1)*psi((1+nu)/2)- ((nu+2)/2)*psi((nu+2)/2)
  zeta_fun(3)
  
  ##marginal functions of delta T and delta S
  fT_1<- function(x) dst(x, mu=mean(deltaT), sigma=sd(deltaT), nu=3, log=FALSE)
  fS_1<- function(y) dst(y, mu=mean(deltaS), sigma=sd(deltaS), nu=3, log=FALSE)
  
  cov_matrix<- matrix(c(var(deltaT), cov(deltaT, deltaS),cov(deltaT, deltaS),var(deltaS)), nrow = 2)
  inv_cov<- inv(cov_matrix)
  
  ##joint functions of delta T and delta S
  p<-2
  fxy_1 <- function(x) {
    gamma( (nu+p)/2)/ (gamma(nu/2)*pi*nu*sqrt(det(cov_matrix) )) * (1+ (1/nu)*( (c(x[1],x[2])-(c(mean(deltaT),mean(deltaS) )))%*%inv_cov %*% ((c(x[1],x[2])-(c(mean(deltaT),mean(deltaS) ))))  ) )^( -2.5) 
  }
  
  ##Kullback-Leibler function 
  KLfun<- function(x)   fxy_1(c(x[1],x[2]))*  ( log(  (fxy_1(c(x[1],x[2]))) /  (fT_1(x[1])*fS_1(x[2])) ))
  
  ##mutual information calculations
  mutual_info_KL_full<-cubature::cubintegrate(f = KLfun, lower = c(-Inf,-Inf), upper = c(Inf,Inf),  method = "pcubature")
  integral_error<- mutual_info_KL_full$error
  mutual_info_KL_full<- mutual_info_KL_full$integral
  
  ##copula
  a<- rank(deltaT) / (length(deltaT) + 1)
  b<- rank(deltaS) / (length(deltaS) + 1)
  temp_matrix<- matrix(c(a,b), ncol = 2)
  fit_kdecop <- kdecop(temp_matrix, knots = 151)
  mutual_info_copula <- dep_measures(fit_kdecop, n_qmc = 1e5, measures = "minfo")
  
  ##from the FNN package 
  mutual_info_package <- mutinfo(deltaT, deltaS, k = 10, direct = TRUE)
  mutual_info_package_based_entropy <- mutinfo(deltaT, deltaS, k = 10, direct = FALSE)
  mutual_info_package_based_entropy<- mean(mutual_info_package_based_entropy)
  
  
  ##mutual infor from the formula for t distribution
  mutual_info_normal<- -0.5*log(1-rho^2)
  mutual_info_formula<- zeta_fun(nu)+mutual_info_normal
  
  
  all_infos<- cbind(mutual_info_normal, zeta_fun(nu),  mutual_info_formula, mutual_info_KL_full,integral_error, mutual_info_package_based_entropy, mutual_info_package, mutual_info_copula)
  
  mutual_informations_all<- cbind(mutual_info_normal, mutual_info_formula, mutual_info_KL_full,mutual_info_package_based_entropy, mutual_info_package, mutual_info_copula)
  rh_2_all<- 1- exp(-2*mutual_informations_all)
  
  write(rh_2_all, file="t_dist_rh_2.txt", append=TRUE) 
  write(all_infos, file="t_dist_mutualinfos.txt", append=TRUE) 

end_time<- Sys.time()
print(end_time)


##read file
mutual_infos<-read.table(file= '/Users/gokcedeliorman/t_dist_mutualinfos.txt', header = FALSE, sep = "\t", dec = ".")
rh_2_t_dist<-read.table(file= '/Users/gokcedeliorman/t_dist_rh_2.txt', header = FALSE, sep = "\t", dec = ".")


