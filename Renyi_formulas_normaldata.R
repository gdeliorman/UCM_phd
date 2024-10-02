##Example of Renyi Divergence calculation between fx and fy with formulas and integral  

##PAPER 1: Renyi Divergence Measures for Commonly Used Univariate Continuous Distributions (M. Gil) in page 8
##PAPER 2: Renyi Divergence and Kullback-Leibler Divergence (Tim van Erven) page 6

##example parameters 
sigma_x<- 50
sigma_y<- 40
sd_x<- sqrt(sigma_x)
sd_y<-sqrt(sigma_y)
beta<- 5
alpha<- 3

##alpha>0
alpha_k<- 1.5
#alpha_k<- 0.5
#alpha_k<- 2


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

##Renyi calculation using integral
Renyi_function<- function(x) fx(x)^(alpha_k) * ( (fy(x)) ^(1-alpha_k))
Renyi_integral<- (1/(alpha_k-1))*log(integral(Renyi_function, xmin,xmax, reltol = 1e-10))



##KL calculation using Gil and Tim Ervan formulas
sigma_alpha<- alpha_k*sigma_y+ sigma_x*(1-alpha_k) ##it should be positive


first_part<-    (alpha_k*(beta-alpha)^2)/(2*sigma_alpha)
second_part<-  (1/(1-alpha_k))*log(sqrt(sigma_alpha)/ ( (sd_y^alpha_k)*(sd_x)^(1-alpha_k) ) )
second_part_2<- (1/ (2*(alpha_k-1)))*log( (sigma_alpha)/ ( (sigma_y^alpha_k)*(sigma_x)^(1-alpha_k) ) )

Renyi_formula_Gil<- first_part-second_part_2
Renyi_formula_Tim_and_Erven<- first_part+second_part




