##Example of Kullback-Leibler Divergence calculation between fx and fy with formulas and integral  

##PAPER 1: Renyi Divergence Measures for Commonly Used Univariate Continuous Distributions (M. Gil) in page 8
##PAPER 2: Renyi Divergence and Kullback-Leibler Divergence (Tim van Erven) page 6

##example parameters 
sigma_x<- 50
sigma_y<- 40
sd_x<- sqrt(sigma_x)
sd_y<-sqrt(sigma_y)
beta<- 5
alpha<- 3


##functions from normal distribution
fx<-function(x) dnorm(x,mean = beta,sd=sqrt(sigma_x))
fy<-function(y) dnorm(y,mean = alpha,sd=sqrt(sigma_y))

##find integral limits
fT_bound<-matrix(NA, ncol=2)
for (k in -300:300) {
  if (fx(k)>=0.00001)
  {
    
    c<-cbind(k,fx(k))
    fT_bound<-rbind(c,fT_bound)
    
  }
}


xmax<- max(na.omit(fT_bound[,1]))
xmin<- min(na.omit(fT_bound[,1]))

##KL calculation using integral
KL_function<- function(x) fx(x)* log( (fx(x)/fy(x)) )
KL_integral<-integral(KL_function, xmin,xmax, reltol = 1e-10)
#KL2<-integral(KL_function, -100,100, reltol = 1e-10)
#KL3<-integral(KL_function, -200,200, reltol = 1e-10)

##KL calculation using R packages
f1<- rnorm(100000, mean = beta, sd = sd_x)
f2<- rnorm(100000, mean = alpha, sd = sd_y)

KL_package<- mean(KL.divergence(f1, f2, k=10))

##KL calculation using Gil and Tim Ervan formulas
KL_formula_Gil<- (1/(2*sigma_y))*( (beta-alpha)^2+sigma_x-sigma_y )+ log(sd_y/sd_x)
KL_formula_Tim_and_Erven<- 0.5*( ((beta-alpha)^2/sigma_y) + log(sigma_y/sigma_x)+ (sigma_x/sigma_y)-1 ) 
  
  
  

