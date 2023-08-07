##true value of ICA in log normal case

library(extraDistr) #for bivariate dist
library(mvtnorm)
library(MASS)#for generate multinormal data
library(dplyr)
library(LaplacesDemon)
library(pracma)
library(compositions) ##multiple log normal
library(Metrics) ##FOR BÄ°AS
library(devtools) ## local functiom
library(ks) ###kernel package
library(Surrogate)
library(compositions) ##generate multi log normal
library(kdensity)
library(akima) ##for approachfun2d
library(ADGofTest)
library(KernSmooth) ##BKDE2D


print("medium log-normal works")
start_time<- Sys.time()
print(start_time)

for(i in 1:10000) {

###SİMULATİON
  print(i)
mu<- c(1.2, 2, 1.5, 3)
sigma<- rbind(c(1.0780992, -0.3676667,  0.2402575, -0.4867786),
              c(-0.3676667,  1.0623604, -0.4322828,  0.3046945),
              c(0.2402575 ,-0.4322828 , 1.0594279, -0.2495095),
              c(-0.4867786 , 0.3046945, -0.2495095,  1.0704080))

data_log<-rlnorm.rplus(10000,mu,sigma)
T0<-data_log[,1]
T1<-data_log[,2]
S0<-data_log[,3]
S1<-data_log[,4]
#treatment<-rbinom(10000,1,0.5)

#data<-data.frame(T0,T1,S0,S1, treatment)
cor(T0,S0)
cor(T1,S1)

deltaT<- T1-T0
deltaS<- S1-S0

cor(deltaT, deltaS)


try(kdeT<- kdensity(deltaT, start = "gaussian", kernel= "gaussian", bw=bw.nrd0(deltaT)), silent = TRUE)
try(kdeS<- kdensity(deltaS, start = "gaussian", kernel= "gaussian", bw=bw.nrd0(deltaS)), silent = TRUE)

if (exists("kdeS")== TRUE | exists("kdeT")== TRUE) {
  
  #kdeT<- kdensity(deltaT, start = "gaussian", kernel= "gaussian", bw=bw.nrd0(deltaT))
  #kdeS<- kdensity(deltaS, start = "gaussian", kernel= "gaussian", bw=bw.nrd0(deltaS))
  
  fT<- kdeT
  fS<- kdeS
  
  fT_bound<-matrix(NA, ncol=2)
  for (i in -100:100) {
    
    if (fT(i)>=0.001)
    {
      
      c<-cbind(i, fT(i))
      fT_bound<-rbind(c,fT_bound)
    }
  }
  
  ###finding boundaries fS
  fS_bound<-matrix(NA, ncol=2)
  for (i in -100:100) {
    if (fS(i)>=0.001)
    {
      
      c<-cbind(i,fS(i))
      fS_bound<-rbind(c,fS_bound)
      
    }
  }
  
  xmax<- max(na.omit(fT_bound[,1]))
  xmin<- min(na.omit(fT_bound[,1]))
  ymax<- max(na.omit(fS_bound[,1]))
  ymin<- min(na.omit(fS_bound[,1]))
  
  
  ##bivariate kernel
  est1 <- bkde2D(cbind(deltaT, deltaS), bandwidth=c(bw.nrd0(deltaT),bw.nrd0(deltaS)) , gridsize = c(1000L, 1000L), range.x = list(c(xmin,xmax), c(ymin,ymax)))
  
  approxfun21 <- function(x=est1$x1, y=est1$x2, z=est1$fhat, method = "linear") {
    function(xp, yp) interp2(est1$x1, est1$x2, est1$fhat, xp, yp, method) }
  fTS <- approxfun21(x=est1$x1, y=est1$x2, z=est1$fhat)
  KLfun<- function(x,y)   fTS(x,y)*  ( log(  (fTS(x,y)) /  (fT(x)*fS(y)) ))
  
  
  point_xx<- runif(n=10000, min=xmin, max=xmax)
  point_yy<- runif(n=10000, min=ymin, max=ymax)
  elimination2<- which(fTS(point_xx, point_yy)< 10^(-6) & fT(point_xx)*fS(point_yy)>= 10^(-6) )
  
  fullvalue1<- (KLfun(point_xx, point_yy))
  ele<-which(is.na(fullvalue1))
  fullvalue1[c(elimination2)]<-0
  
  if(isempty(ele) ==TRUE){
    fullvalue1<- fullvalue1
    
  } else {fullvalue1<- fullvalue1[-c(ele)] }
  
  int1<-mean(na.omit(fullvalue1))*(xmax-xmin)*(ymax-ymin)
  ICA_L<-1-exp(-2*int1)
  
  write(ICA_L, file="log_normal_ICA_L_true_value_2august.txt", append=TRUE) }}

##read text file

#true_ICA<- read.delim(file= '/Users/gdy/Desktop/2022-2023PHD/log_normal_ICA_L_true_value_2august.txt', header = FALSE, sep = "\t", dec = ".")
true_ICA<- read.delim(file= '/Users/gdy/Desktop/phd_second_paper/true_ICA_log_normal_ICA_L_true_value.txt', header = FALSE, sep = "\t", dec = ".")

true_ICA<- as.numeric(as.matrix(true_ICA))


mean(true_ICA)
hist(true_ICA)
plot(density(true_ICA))
plot(true_ICA)

