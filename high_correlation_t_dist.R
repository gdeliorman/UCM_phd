###

library(extraDistr) #for bivariate dist
library(mvtnorm)
library(MASS) #for generate multinormal data
library(dplyr)
library(LaplacesDemon)
library(pracma)
library(compositions) ##multivariate log normal
library(Metrics) ##FOR Bias
library(devtools) ## local function
library(ks) ###kernel package
library(Surrogate)
library(compositions) ##generate multi log normal
library(kdensity)
library(akima) ##for approachfun2d
library(ADGofTest)
library(KernSmooth) ##BKDE2D
library(kdecopula)
library(FNN)
library(tmvtnorm)
library(fMultivar)
library(fGarch)
library(VGAM)
library(matlib)


start_time<- Sys.time()
print(start_time)

# Set the degrees of freedom parameter for the multivariate t-distribution.
nu <- 3
# Define the zeta function as defined in the paper draft.
zeta <- function(x) {
  2 * log(gamma(nu / 2) / gamma((1 + nu) / 2) * sqrt(nu / 2)) - (1 + nu) *
    psi((1 + nu) / 2) + (1 + nu) * psi(nu / 2) + ((2 + nu) / nu)
}  

zeta(3)

# Set the seed for reproducibility.
set.seed(1)
# Mean vector and Sigma matrix for the multivariate t-distribution.
mu<- c(0.3, 1, 1.5, 1)
Sigma<-rbind(c(544.3285, 0, 300.9816, 0), c(0,550.6597,0, 304.4222), c(300.9816,0, 180.6831, 0), c(0,304.4222, 0,180.9433))


# Generate 1e5 samples from the multivariate t-distribution. We are generating a
# large number of samples because we want to exclude the influence of sampling
# variability. 
data_t <- rmvt(1e5, mu, Sigma, df = nu)

T0 <- data_t[, 1]
T1 <- data_t[, 2]
S0 <- data_t[, 3]
S1 <- data_t[, 4]

# Compute the sampled individual causal effects.
deltaT <- T1-T0
deltaS <- S1-S0

Rho<-cor(deltaT, deltaS)
A = matrix(c(-1, 1, 0, 0,
             0, 0, -1, 1), ncol = 4, byrow = TRUE)
Sigma_Delta = A %*% Sigma %*% t(A)
rho_Delta = Sigma_Delta[1, 2] / sqrt(Sigma_Delta[1, 1] * Sigma_Delta[2, 2])


##RH^2 from copulas
# version of this method is implemented in the Surrogate package.
a = rank(deltaS) / (length(deltaS) + 1)
b = rank(deltaT) / (length(deltaT) + 1)
temp_matrix = matrix(data = c(a, b), ncol = 2)
fit <- kdecop(temp_matrix)
minfo_numerical_1 <- dep_measures(fit, n_qmc = 1e5, measures = "minfo")
minfo_numerical_2 <- mutinfo(deltaT, deltaS, k = 10, direct = TRUE)

RH_square_copula_1 <- 1-exp(-2*minfo_numerical_1)
RH_square_copula_2 <- 1-exp(-2*minfo_numerical_2)


# Compute the mutual information based on the formula in the paper draft.
minfo_formula <- (-0.5 * log(1 - Rho ** 2)) + zeta(nu)
RH_2_t_dist<- 1-exp(-2*minfo_formula)

1-(1-Rho^2)*exp(-2*zeta(3))


# Print the three mutual information computed by the three different approaches.
# The numerical methods agree with each other. The value calculated from the
# formula does not agree with the numerical methods.
print(c(minfo_formula, minfo_numerical_1, minfo_numerical_2))

# To examine whether the formula is correct, we will compute the mutual
# information through numerical integration and using the "true" copula density.

# Define copula that corresponds to distribution of the individual causal
# effects.
t_copula = copula::ellipCopula(
  family = "t",
  param = rho_Delta,
  dim = 2,
  dispstr = "ex",
  df = nu
)
t_copula

# Define the integrand in the definition of the mutual information. We use two
# "versions": (i) integrand on the copula scale, and (ii) integrand after
# transforming the margins to standard normal distribution. The latter basically
# amounts to using the substition method for solving the integral.
minfo_integrand = function(u) {
  d = copula::dCopula(copula = t_copula, u = t(u), log = FALSE)
  integrand = integrand = ifelse(d == 0, 0, d * log(d))
  return(matrix(integrand, ncol = 1))
}

minfo_integrand_substitution = function(z) {
  u = pnorm(z)
  d = copula::dCopula(copula = t_copula, u = t(u), log = FALSE)
  integrand = ifelse(d == 0, 0, d * log(d) * dnorm(z[1, ]) * dnorm(z[2, ]))
  return(matrix(integrand, ncol = 1))
}

# We compute the integrals numerically using different methods. To minimize the
# probability that we get misleading results because one intergration method
# fails, we use all methods that are available in the cubature R-package.

# The following numerical integration methods are deterministic. They all give
# the same results.
minfo_hcubature_1 = cubature::cubintegrate(
  f = minfo_integrand, 
  lower = c(0, 0),
  upper = c(1, 1), 
  method = "hcubature",
  nVec = 1024
)
minfo_hcubature_1
RH_square_hcubature_1 <- 1-exp(-2*minfo_hcubature_1$integral)


minfo_cuhre_1 = cubature::cubintegrate(
  f = minfo_integrand, 
  lower = c(0, 0),
  upper = c(1, 1), 
  method = "cuhre",
  nVec = 1024
)
minfo_cuhre_1
RH_square_minfo_cuhre_1 <- 1-exp(-2*minfo_cuhre_1$integral)



# The following integration methods have a stochastic component. They give the
# same results, but some methods give indicate that there may be an issue. Some
# methods return a high probability that the estimated integration error is not
# reliable.
minfo_vegas_1 = cubature::cubintegrate(
  f = minfo_integrand, 
  lower = c(0, 0),
  upper = c(1, 1), 
  method = "vegas",
  nVec = 1024, maxEval = 5e6
)
minfo_vegas_1
RH_square_minfo_vegas_1<- 1-exp(-2*minfo_vegas_1$integral)

minfo_suave_1 = cubature::cubintegrate(
  f = minfo_integrand, 
  lower = c(0, 0),
  upper = c(1, 1), 
  method = "suave",
  nVec = 1024, maxEval = 5e6
)
minfo_suave_1
RH_square_minfo_suave_1<- 1-exp(-2*minfo_suave_1$integral)


minfo_divonne_1 = cubature::cubintegrate(
  f = minfo_integrand, 
  lower = c(0, 0),
  upper = c(1, 1), 
  method = "divonne",
  nVec = 1024, maxEval = 5e6
)
minfo_divonne_1
RH_square_minfo_divonne_1<- 1-exp(-2*minfo_divonne_1$integral)



##parameter estimation of deltaT and deltaS
sstdFit(deltaT)
sstdFit(deltaS)
fitdistr(deltaT, "t")


##marginals 
fT_1<- function(x) dst(x, mu=mean(deltaT), sigma=sd(deltaT), nu=3, log=FALSE)
plot(fT_1, xlim=c(-500,500))
fS_1<- function(y) dst(y, mu=mean(deltaS), sigma=sd(deltaS), nu=3, log=FALSE)
plot(fS_1, xlim=c(-500,500))


##find min and max
fT_bound<-matrix(NA, ncol=2)
for (i in -300:300) {
  
  if (fT_1(i)>=0.0001)
  {
    
    c<-cbind(i, fT_1(i))
    fT_bound<-rbind(c,fT_bound)
  }
}

###finding bound for fT and fS
fS_bound<-matrix(NA, ncol=2)
for (i in -300:300) {
  if (fS_1(i)>=0.0001)
  {
    
    c<-cbind(i,fS_1(i))
    fS_bound<-rbind(c,fS_bound)
    
  }
}


xmax<- max(na.omit(fT_bound[,1]))
xmin<- min(na.omit(fT_bound[,1]))
ymax<- max(na.omit(fS_bound[,1]))
ymin<- min(na.omit(fS_bound[,1]))



##Bivariate t-dist functions 
##kernel density estimation
est1 <- bkde2D(cbind(deltaT, deltaS), bandwidth=c(bw.nrd0(deltaT),bw.nrd0(deltaS)) , gridsize = c(512L, 512L), range.x = list(c(xmin,xmax), c(ymin,ymax)))

approxfun21 <- function(x=est1$x1, y=est1$x2, z=est1$fhat, method = "linear") {
  function(xp, yp) interp2(est1$x1, est1$x2, est1$fhat, xp, yp, method) }
fTS_1 <- approxfun21(x=est1$x1, y=est1$x2, z=est1$fhat)

##kernel density estimation 2
est2<- kde2d(deltaT, deltaS, n = 512, lims = c(xmin, xmax, ymin, ymax))   # from MASS package

approxfun22 <- function(x=est2$x, y=est2$y, z=est2$z, method = "linear") {
  function(xp, yp) interp2(est2$x, est2$y, est2$z, xp, yp, method) }
fTS_2 <- approxfun22(x=est2$x, y=est2$y, z=est2$z)

##they have so close results
fTS_1(1,2)
fTS_2(1,2)

## ##Kullback-Leibler divergence 
KL_fun_1<- function(x,y)   fTS_1(x,y)*  ( log(  (fTS_1(x,y)) /  (fT_1(x)*fS_1(y)) ))


###bivariate t distribution function
delta_TS<- matrix( c(as.numeric(deltaT), as.numeric(deltaS)), ncol=2 )
mu_2 <- colMeans(delta_TS)
Sigma_2<- var(delta_TS)

fTS_3<- function(x,y) (( gamma((nu+2)/2)  ) / ( gamma(nu/2)*pi*nu* sqrt(det(Sigma_2))  ) )*(1+ ( (1/nu)*(( matrix( c(x,y), ncol=2) )-mu_2) %*% inv(Sigma_2) %*% ( t(matrix(c(x,y), ncol=2))-mu_2) ))^(-2.5)
KL_fun_3<- function(x,y)   fTS_3(x,y)*  ( log(  (fTS_3(x,y)) /  (fT_1(x)*fS_1(y)) ))


##integration 
ICA_t<- matrix(NA, ncol = 1)
mut_kernel<- matrix(NA, ncol = 1)
for (k in 1:100) {
  point_xx<- runif(n=10000, min=xmin, max=xmax)
  point_yy<- runif(n=10000, min=ymin, max=ymax)
  full_value<- (KL_fun_1(point_xx, point_yy))
  elimination<- which(fTS_1(point_xx, point_yy)< 10^(-6) & fT_1(point_xx)*fS_1(point_yy)>= 10^(-6) )
  full_value[c(elimination)]<-0
  ele<-which(is.na(full_value))
  if(isempty(ele) ==TRUE){full_value<- full_value} else {full_value<- full_value[-c(ele)] }
  int1<-mean(na.omit(full_value))*(xmax-xmin)*(ymax-ymin)
  mut_kernel[k]<- int1
  ICA_tt<-1-exp(-2*int1)
  ICA_t[k]<- ICA_tt
}

RH_square_kernel<- mean(ICA_t)
minfo_kernel<- mean(mut_kernel)


##integral for KL_fun3 (it should work but it doesnt work)
#full_value_2<- (KL_fun_3(point_xx, point_yy))
#elimination_2<- which(fTS_3(point_xx, point_yy)< 10^(-6) & fT_1(point_xx)*fS_1(point_yy)>= 10^(-6) )
#full_value_2[c(elimination_2)]<-0
#ele_2<-which(is.na(full_value_2))
#if(isempty(ele) ==TRUE){full_value_2<- full_value_2} else {full_value<- full_value_2[-c(ele_2)] }
#int2<-mean(na.omit(full_value_2))*(xmax-xmin)*(ymax-ymin)
#mut_t<- int2
#ICA_t<-1-exp(-2*int2)


print(c(RH_2_t_dist, RH_square_kernel, RH_square_copula_1, RH_square_copula_2, RH_square_minfo_cuhre_1,
        RH_square_hcubature_1, RH_square_minfo_divonne_1, RH_square_minfo_vegas_1))


print(c(minfo_formula, minfo_kernel, minfo_numerical_1, minfo_numerical_2, 
        minfo_hcubature_1$integral, minfo_divonne_1$integral, minfo_vegas_1$integral))


end_time<- Sys.time()
print(end_time)




