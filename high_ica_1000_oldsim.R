##libraries
#library(Surrogate) #for ICA
library(extraDistr) #for bivariate dist
library(MASS)#for generate multinormal data
library(dplyr)
library(LaplacesDemon)
library(pracma)  ##for integral


##n=1000
#ica: 0.941785
start_time<- Sys.time()
print(start_time)
print("n1000 high normal")
highmetricn1000<-matrix(c(-0.40, 0.963606, -0.50, -0.40, 0.963032, -0.50), ncol = 6)


ralpha01<-0.1
ralpha05<-0.5
ralpha125<-1.25
ralpha15<-1.5
ralpha2<-2
ralpha25<-2.5
ralpha3<-3

the_last_ICAS<- list()
last_ICA_alpha01<-list()
last_ICA_alpha05<-list()
last_ICA_alpha125<-list()
last_ICA_alpha15<-list()
last_ICA_alpha2<-list()
last_ICA_alpha25<-list()
last_ICA_alpha3<-list()
last_ICAformul1high<- list()

min_last_ICAS<- list()
min_last_ICA_alpha01<-list()
min_last_ICA_alpha05<-list()
min_last_ICA_alpha125<-list()
min_last_ICA_alpha15<-list()
min_last_ICA_alpha2<-list()
min_last_ICA_alpha25<-list()
min_last_ICA_alpha3<-list()
min_last_ICAformul1high<- list()

max_last_ICAS<- list()
max_last_ICA_alpha01<-list()
max_last_ICA_alpha05<-list()
max_last_ICA_alpha125<-list()
max_last_ICA_alpha15<-list()
max_last_ICA_alpha2<-list()
max_last_ICA_alpha25<-list()
max_last_ICA_alpha3<-list()
max_last_ICAformul1high<- list()

range_last_ICAS<- list()
range_last_ICA_alpha01<-list()
range_last_ICA_alpha05<-list()
range_last_ICA_alpha125<-list()
range_last_ICA_alpha15<-list()
range_last_ICA_alpha2<-list()
range_last_ICA_alpha25<-list()
range_last_ICA_alpha3<-list()
range_last_ICAformul1high<- list()

median_last_ICAS<- list()
median_last_ICA_alpha01<-list()
median_last_ICA_alpha05<-list()
median_last_ICA_alpha125<-list()
median_last_ICA_alpha15<-list()
median_last_ICA_alpha2<-list()
median_last_ICA_alpha25<-list()
median_last_ICA_alpha3<-list()
median_last_ICAformul1high<- list()

sd_last_ICAS<- list()
sd_last_ICA_alpha01<-list()
sd_last_ICA_alpha05<-list()
sd_last_ICA_alpha125<-list()
sd_last_ICA_alpha15<-list()
sd_last_ICA_alpha2<-list()
sd_last_ICA_alpha25<-list()
sd_last_ICA_alpha3<-list()
sd_last_ICAformul1high<- list()

for (N in 1:1000) {
  set.seed(N)
  ###10000 Data generation
  ICAformul1high<-list()
  sigmadeltashigh<-list()
  sigmadeltathigh<-list()
  truefalseresulthigh<-list()
  
  alphas<-list()
  betas<-list()
  
  for(m in 1:1000) {
    
    ##Data generation x times (unobserved data)
    sigma1<-matrix(c(544.3285, 0, 300.9816, 0, 
                     0,550.6597,0, 304.4222, 
                     300.9816,0, 180.6831, 0, 
                     0,304.4222, 0,180.9433), nrow=4, ncol=4)
    
    
    colnames(sigma1) = c("T0", "T1", "S0", "S1")
    rownames(sigma1) = c("T0", "T1", "S0", "S1")
    is.positive.definite(sigma1)
    # create the mean vector (choose from real cases)
    mu<-c(-11.47191,-16.0151,-6.646067,-8.938326)
    # generate the multivariate normal distribution
    df<-as.data.frame(mvrnorm(n=1000, mu=mu, Sigma=sigma1), empirical = TRUE)
    
    T0<-df[,1]
    T1<-df[,2]
    S0<-df[,3]
    S1<-df[,4]
    treatment<-rbinom(1000,1,0.5)
    data<-data.frame(T0,T1,S0,S1, treatment)
    cor(T0,S0)
    cor(T1,S1)
    var(T0,T0)
    
    ##select 2 coloumn for observed data
    newdata<-matrix(NA,nrow= 1000, ncol=3)
    for (i in 1:1000) {
      if(data$treatment[i]==0)
      {
        newdata[i,]<-c(data$T0[i], data$S0[i], data$treatment[i])
      }
      else
        newdata[i,]<-cbind(data$T1[i], data$S1[i], data$treatment[i]) }
    
    newdata<-as.data.frame(newdata)
    placebo<-newdata[newdata$V3 == "0",]
    experimental<-newdata[newdata$V3  == "1",]
    
    ##basic measures
    Tr<-newdata$V1
    S<-newdata$V2
    T0<- placebo$V1
    T1<- experimental$V1
    S0<- placebo$V2
    S1<- experimental$V2
    var(Tr,Tr)
    ##means
    mT0<-mean(T0)
    mT1<-mean(T1)
    mS0<-mean(S0)
    mS1<-mean(S1)
    
    
    
    beta<-mT1-mT0
    alpha<-mS1-mS0
    
    betas<-append(beta, betas)
    alphas<-append(alpha, alphas)
    
    ##variances
    vT0T0<-var(T0,T0)
    vT1T1<-var(T1,T1)
    vS1S1<-var(S1,S1)
    vS0S0<-var(S0,S0)
    vT1S1<-var(T1,S1)
    vT0S0<-var(T0,S0)
    
    ##correlations
    cT0S0<-cor(T0,S0)
    cT1S1<-cor(T1,S1)
    cTS<-cor(Tr,S)
    
    Sigma_matrix<- matrix(NA, nrow=4, ncol=4)
    Sigma_matrix[1,1]<- vT0T0
    Sigma_matrix[2,2]<- vT1T1
    Sigma_matrix[3,3]<- vS0S0
    Sigma_matrix[4,4]<- vS1S1
    Sigma_matrix[1,3]<-Sigma_matrix[3,1]<- vT0S0
    Sigma_matrix[2,4]<-Sigma_matrix[4,2]<- vT1S1
    
    Sigma_matrix[1,2]<-Sigma_matrix[2,1]<- highmetricn1000[,1]*sd(T0)*sd(T1)
    Sigma_matrix[1,4]<-Sigma_matrix[4,1]<- highmetricn1000[,3]*sd(T0)*sd(S1)
    Sigma_matrix[2,3]<-Sigma_matrix[3,2]<- highmetricn1000[,4]*sd(T1)*sd(S0)
    Sigma_matrix[3,4]<-Sigma_matrix[4,3]<- highmetricn1000[,6]*sd(S0)*sd(S1)
    
    if(is.positive.definite(Sigma_matrix)==TRUE){
      
      payhigh<-(sqrt(vT0T0*vS0S0)*cT0S0)+(sqrt(vT1T1*vS1S1)*cT1S1)-(sqrt(vT1T1*vS0S0)*highmetricn1000[,4])-(sqrt(vT0T0*vS1S1)*highmetricn1000[,3])
      paydahigh<-sqrt((vT0T0+vT1T1-2*sqrt(vT0T0*vT1T1)*highmetricn1000[,1])*(vS0S0+vS1S1-2*sqrt(vS0S0*vS1S1)*highmetricn1000[,6]))
      ICAformulhigh<-payhigh/paydahigh
      
      ICAformul1high<- append(ICAformulhigh, ICAformul1high)
      ICAformul1high<-as.numeric(as.matrix(ICAformul1high))
      sigmadelta_shigh<-var(S0,S0)+var(S1,S1)-2*highmetricn1000[,6]*sqrt(var(S0,S0)*var(S1,S1))
      sigmadeltashigh<- append(sigmadelta_shigh,sigmadeltashigh)
      sigmadelta_thigh<-var(T0,T0)+var(T1,T1)-2*highmetricn1000[,1]*sqrt(var(T0,T0)*var(T1,T1))
      sigmadeltathigh<- append(sigmadelta_thigh,sigmadeltathigh)
      cormatrixhigh<- matrix(c(1, highmetricn1000[,1], cT0S0,highmetricn1000[,3],
                               highmetricn1000[,1],1,highmetricn1000[,4],cT1S1,
                               cT0S0,highmetricn1000[,4],1,highmetricn1000[,6],
                               highmetricn1000[,3],cT1S1,highmetricn1000[,6], 1 ), nrow = 4)
      
      
      alphas<-as.numeric(as.matrix(alphas))
      betas<-as.numeric(as.matrix(betas))
      
      sigmadeltashigh<-as.numeric(as.matrix(sigmadeltashigh))
      sigmadeltathigh<-as.numeric(as.matrix(sigmadeltathigh))
      
      highmetrics1000<- cbind( ICAformul1high, sigmadeltathigh, sigmadeltashigh, betas, alphas)
      highmetrics1000<- as.matrix(highmetrics1000) }}
  
  
  
  ##high ICA
  ###KL and Renyi high ICA n=1000, 1000 times
  
  ICA_KL<-list()
  ICA_alpha01<-list()
  ICA_alpha05<-list()
  ICA_alpha125<-list()
  ICA_alpha15<-list()
  ICA_alpha2<-list()
  ICA_alpha25<-list()
  ICA_alpha3<-list()
  
  
  for (l in 1:nrow(highmetrics1000) ) {
    
    fx<-function(x) dnorm(x,mean = highmetrics1000[l,4],sd=sqrt( as.numeric(highmetrics1000[l,2])))
    fy<-function(y) dnorm(y,mean =highmetrics1000[l,5],sd=sqrt( as.numeric(highmetrics1000[l,3])))
    fxy<- function(x,y) dbvnorm(x,y, mean1 = as.numeric(highmetrics1000[l,4]), mean2 = as.numeric(highmetrics1000[l,5]), sd1 =sqrt( as.numeric(highmetrics1000[l,2])) , sd2 = sqrt(as.numeric(highmetrics1000[l,3])),   cor=as.numeric(highmetrics1000[l,1]))
    
    RY25<- function(x,y) (fxy(x,y)^ralpha25) * ( (fx(x)*fy(y)) ^(1-ralpha25) )
    RY3<- function(x,y) (fxy(x,y)^ralpha3) * ( (fx(x)*fy(y)) ^(1-ralpha3) )
    
    fT_bound<-matrix(NA, ncol=2)
    for (k in -400:400) {
      
      if (fx(k)>=0.00001)
      {
        
        c<-cbind(k, fx(k))
        fT_bound<-rbind(c,fT_bound)
      }
    }
    
    ###finding boundaries fS
    fS_bound<-matrix(NA, ncol=2)
    for (j in -400:400) {
      if (fy(j)>=0.00001)
      {
        
        c<-cbind(j,fy(j))
        fS_bound<-rbind(c,fS_bound)
        
      }
    }
    
    xmax<- max(na.omit(fT_bound[,1]))
    xmin<- min(na.omit(fT_bound[,1]))
    ymax<- max(na.omit(fS_bound[,1]))
    ymin<- min(na.omit(fS_bound[,1]))
    
    rh<- as.numeric(highmetrics1000[l,1])
    ICA_1<-  rh^2
    ICA_KL<- append(ICA_1, ICA_KL)
    
    ICA_01<- 1-(1-rh^2)*(1-(1-ralpha01)^2*rh^2)^(-1/(1-ralpha01))
    ICA_alpha01<-append(ICA_01,ICA_alpha01)
    
    ICA_05<- 1-(1-rh^2)*(1-(1-ralpha05)^2*rh^2)^(-1/(1-ralpha05))
    ICA_alpha05<-append(ICA_05,ICA_alpha05)
    
    ICA_125<- 1-(1-rh^2)*(1-(1-ralpha125)^2*rh^2)^(-1/(1-ralpha125))
    ICA_alpha125<-append(ICA_125,ICA_alpha125)
    
    ICA_15<- 1-(1-rh^2)*(1-(1-ralpha15)^2*rh^2)^(-1/(1-ralpha15))
    ICA_alpha15<-append(ICA_15,ICA_alpha15)
    
    ICA_2<- 1-(1-rh^2)*(1-(1-ralpha2)^2*rh^2)^(-1/(1-ralpha2))
    ICA_alpha2<-append(ICA_2,ICA_alpha2)
    
    R_25<-integral2(RY25,xmin,xmax, ymin, ymax,reltol = 1e-10)
    R25<- log(R_25$Q)* (1/ (ralpha25-1))
    ICA_25<-1-exp(-2*R25)
    ICA_alpha25<-append(ICA_25,ICA_alpha25)
    
    R_3<-integral2(RY3,xmin,xmax, ymin, ymax,reltol = 1e-10)
    R3<- log(R_3$Q)* (1/ (ralpha3-1))
    ICA_3<-1-exp(-2*R3)
    ICA_alpha3<-append(ICA_3,ICA_alpha3) }
  
  
  
  the_last_ICAS<- append(the_last_ICAS, mean(na.omit(as.numeric(as.matrix(ICA_KL)))))
  last_ICA_alpha01<- append(last_ICA_alpha01, mean(na.omit(as.numeric(as.matrix(ICA_alpha01)))))
  last_ICA_alpha05<- append(last_ICA_alpha05, mean(na.omit(as.numeric(as.matrix(ICA_alpha05)))))
  last_ICA_alpha125<- append(last_ICA_alpha125, mean(na.omit(as.numeric(as.matrix(ICA_alpha125)))))
  last_ICA_alpha15<- append(last_ICA_alpha15, mean(na.omit(as.numeric(as.matrix(ICA_alpha15)))) )
  last_ICA_alpha2<- append(last_ICA_alpha2, mean(na.omit(as.numeric(as.matrix(ICA_alpha2)))) )
  last_ICA_alpha25<- append(last_ICA_alpha25, mean(na.omit(as.numeric(as.matrix(ICA_alpha25)))) )
  last_ICA_alpha3<- append(last_ICA_alpha3, mean(na.omit(as.numeric(as.matrix(ICA_alpha3)))) )
  last_ICAformul1high<-append( last_ICAformul1high ,mean(as.numeric(as.matrix(ICAformul1high))))
  
  min_last_ICAS<- append(min_last_ICAS, min(na.omit(as.numeric(as.matrix(ICA_KL)))))
  min_last_ICA_alpha01<- append(min_last_ICA_alpha01, min(na.omit(as.numeric(as.matrix(ICA_alpha01)))))
  min_last_ICA_alpha05<- append(min_last_ICA_alpha05, min(na.omit(as.numeric(as.matrix(ICA_alpha05)))))
  min_last_ICA_alpha125<- append(min_last_ICA_alpha125, min(na.omit(as.numeric(as.matrix(ICA_alpha125)))))
  min_last_ICA_alpha15<- append(min_last_ICA_alpha15, min(na.omit(as.numeric(as.matrix(ICA_alpha15)))) )
  min_last_ICA_alpha2<- append(min_last_ICA_alpha2, min(na.omit(as.numeric(as.matrix(ICA_alpha2)))) )
  min_last_ICA_alpha25<- append(min_last_ICA_alpha25, min(na.omit(as.numeric(as.matrix(ICA_alpha25)))) )
  min_last_ICA_alpha3<- append(min_last_ICA_alpha3, min(na.omit(as.numeric(as.matrix(ICA_alpha3)))) )
  
  max_last_ICAS<-        append(max_last_ICAS, max(na.omit(as.numeric(as.matrix(ICA_KL)))))
  max_last_ICA_alpha01<- append(max_last_ICA_alpha01, max(na.omit(as.numeric(as.matrix(ICA_alpha01)))))
  max_last_ICA_alpha05<- append(max_last_ICA_alpha05, max(na.omit(as.numeric(as.matrix(ICA_alpha05)))))
  max_last_ICA_alpha125<-append(max_last_ICA_alpha125, max(na.omit(as.numeric(as.matrix(ICA_alpha125)))))
  max_last_ICA_alpha15<- append(max_last_ICA_alpha15, max(na.omit(as.numeric(as.matrix(ICA_alpha15)))) )
  max_last_ICA_alpha2<-  append(max_last_ICA_alpha2, max(na.omit(as.numeric(as.matrix(ICA_alpha2)))) )
  max_last_ICA_alpha25<- append(max_last_ICA_alpha25, max(na.omit(as.numeric(as.matrix(ICA_alpha25)))) )
  max_last_ICA_alpha3<-  append(max_last_ICA_alpha3, max(na.omit(as.numeric(as.matrix(ICA_alpha3)))) )
  
  median_last_ICAS<-append(median_last_ICAS, median(na.omit(as.numeric(as.matrix(ICA_KL)))))
  median_last_ICA_alpha01<- append(median_last_ICA_alpha01, median(na.omit(as.numeric(as.matrix(ICA_alpha01)))))
  median_last_ICA_alpha05<- append(median_last_ICA_alpha05, median(na.omit(as.numeric(as.matrix(ICA_alpha05)))))
  median_last_ICA_alpha125<-append(median_last_ICA_alpha125, median(na.omit(as.numeric(as.matrix(ICA_alpha125)))))
  median_last_ICA_alpha15<- append(median_last_ICA_alpha15, median(na.omit(as.numeric(as.matrix(ICA_alpha15)))) )
  median_last_ICA_alpha2<-  append(median_last_ICA_alpha2, median(na.omit(as.numeric(as.matrix(ICA_alpha2)))) )
  median_last_ICA_alpha25<- append(median_last_ICA_alpha25, median(na.omit(as.numeric(as.matrix(ICA_alpha25)))) )
  median_last_ICA_alpha3<-  append(median_last_ICA_alpha3, median(na.omit(as.numeric(as.matrix(ICA_alpha3)))) )
  
  sd_last_ICAS <-       append(sd_last_ICAS, sd(na.omit(as.numeric(as.matrix(ICA_KL)))))
  sd_last_ICA_alpha01<- append(sd_last_ICA_alpha01, sd(na.omit(as.numeric(as.matrix(ICA_alpha01)))))
  sd_last_ICA_alpha05<- append(sd_last_ICA_alpha05, sd(na.omit(as.numeric(as.matrix(ICA_alpha05)))))
  sd_last_ICA_alpha125<-append(sd_last_ICA_alpha125,sd(na.omit(as.numeric(as.matrix(ICA_alpha125)))))
  sd_last_ICA_alpha15<- append(sd_last_ICA_alpha15, sd(na.omit(as.numeric(as.matrix(ICA_alpha15)))) )
  sd_last_ICA_alpha2<-  append(sd_last_ICA_alpha2,  sd(na.omit(as.numeric(as.matrix(ICA_alpha2)))) )
  sd_last_ICA_alpha25<- append(sd_last_ICA_alpha25, sd(na.omit(as.numeric(as.matrix(ICA_alpha25)))) )
  sd_last_ICA_alpha3<-  append(sd_last_ICA_alpha3,  sd(na.omit(as.numeric(as.matrix(ICA_alpha3)))) )}

end_time<- Sys.time()
print(end_time)

###find ICA and ICA_alphas histograms, mean, sd, max, min, range
ICA_KL<-as.numeric(as.matrix(the_last_ICAS))
ICAformul1high<- as.numeric(as.matrix( last_ICAformul1high))
last_ICAformul1high<- ICAformul1high*ICAformul1high
ICA_alpha01<- as.numeric(as.matrix(last_ICA_alpha01))
ICA_alpha05<- as.numeric(as.matrix(last_ICA_alpha05))
ICA_alpha125<- as.numeric(as.matrix(last_ICA_alpha125))
ICA_alpha15<- as.numeric(as.matrix(last_ICA_alpha15))
ICA_alpha2<- as.numeric(as.matrix(last_ICA_alpha2))
ICA_alpha25<- as.numeric(as.matrix(last_ICA_alpha25))
ICA_alpha3<- as.numeric(as.matrix(last_ICA_alpha3))

min_ICA_KL<-as.numeric(as.matrix(min_last_ICAS))
min_ICA_alpha01<- as.numeric(as.matrix(min_last_ICA_alpha01))
min_ICA_alpha05<- as.numeric(as.matrix(min_last_ICA_alpha05))
min_ICA_alpha125<- as.numeric(as.matrix(min_last_ICA_alpha125))
min_ICA_alpha15<- as.numeric(as.matrix(min_last_ICA_alpha15))
min_ICA_alpha2<- as.numeric(as.matrix(min_last_ICA_alpha2))
min_ICA_alpha25<- as.numeric(as.matrix(min_last_ICA_alpha25))
min_ICA_alpha3<- as.numeric(as.matrix(min_last_ICA_alpha3))

max_ICA_KL<-as.numeric(as.matrix(max_last_ICAS))
max_ICA_alpha01<- as.numeric(as.matrix(max_last_ICA_alpha01))
max_ICA_alpha05<- as.numeric(as.matrix(max_last_ICA_alpha05))
max_ICA_alpha125<- as.numeric(as.matrix(max_last_ICA_alpha125))
max_ICA_alpha15<- as.numeric(as.matrix(max_last_ICA_alpha15))
max_ICA_alpha2<- as.numeric(as.matrix(max_last_ICA_alpha2))
max_ICA_alpha25<- as.numeric(as.matrix(max_last_ICA_alpha25))
max_ICA_alpha3<- as.numeric(as.matrix(max_last_ICA_alpha3))

sd_ICA_KL<-as.numeric(as.matrix(sd_last_ICAS))
sd_ICA_alpha01<- as.numeric(as.matrix(sd_last_ICA_alpha01))
sd_ICA_alpha05<- as.numeric(as.matrix(sd_last_ICA_alpha05))
sd_ICA_alpha125<- as.numeric(as.matrix(sd_last_ICA_alpha125))
sd_ICA_alpha15<- as.numeric(as.matrix(sd_last_ICA_alpha15))
sd_ICA_alpha2<- as.numeric(as.matrix(sd_last_ICA_alpha2))
sd_ICA_alpha25<- as.numeric(as.matrix(sd_last_ICA_alpha25))
sd_ICA_alpha3<- as.numeric(as.matrix(sd_last_ICA_alpha3))

median_ICA_KL<-as.numeric(as.matrix(median_last_ICAS))
median_ICA_alpha01<- as.numeric(as.matrix(median_last_ICA_alpha01))
median_ICA_alpha05<- as.numeric(as.matrix(median_last_ICA_alpha05))
median_ICA_alpha125<- as.numeric(as.matrix(median_last_ICA_alpha125))
median_ICA_alpha15<- as.numeric(as.matrix(median_last_ICA_alpha15))
median_ICA_alpha2<- as.numeric(as.matrix(median_last_ICA_alpha2))
median_ICA_alpha25<- as.numeric(as.matrix(median_last_ICA_alpha25))
median_ICA_alpha3<- as.numeric(as.matrix(median_last_ICA_alpha3))


print(paste("high ica 1000"))

##mean of mean 
print(paste("ICA_KL mean", "=",  mean(ICA_KL)))
print(paste("ICA_alpha01 mean", "=",  mean(ICA_alpha01)))
print(paste("ICA_alpha05 mean", "=",  mean(ICA_alpha05)))
print(paste("ICA_alpha125 mean", "=",  mean(ICA_alpha125)))
print(paste("ICA_alpha15 mean", "=",  mean(ICA_alpha15)))
print(paste("ICA_alpha2 mean", "=",  mean(ICA_alpha2)))
print(paste("ICA_alpha25 mean", "=",  mean(ICA_alpha25)))
print(paste("ICA_alpha3 mean", "=",  mean(ICA_alpha3)))


##mean of the minimums
print(paste("ICA_KL mean of min", "=",  mean(min_ICA_KL)))
print(paste("ICA_alpha01 mean of min", "=",  mean(min_ICA_alpha01)))
print(paste("ICA_alpha05 mean of min", "=",  mean(min_ICA_alpha05)))
print(paste("ICA_alpha125 mean of min", "=",  mean(min_ICA_alpha125)))
print(paste("ICA_alpha15 mean of min", "=",  mean(min_ICA_alpha15)))
print(paste("ICA_alpha2  mean of min", "=",  mean(min_ICA_alpha2)))
print(paste("ICA_alpha25 mean of min", "=",  mean(min_ICA_alpha25)))
print(paste("ICA_alpha3 mean of min", "=",  mean(min_ICA_alpha3)))

##mean of the maksimum
print(paste("ICA_KL mean of max", "=", mean(max_ICA_KL)))
print(paste("ICA_alpha01 mean of max", "=",  mean(max_ICA_alpha01)))
print(paste("ICA_alpha05 mean of max", "=",  mean(max_ICA_alpha05)))
print(paste("ICA_alpha125 mean of max", "=",  mean(max_ICA_alpha125)))
print(paste("ICA_alpha15 mean of max", "=",  mean(max_ICA_alpha15)))
print(paste("ICA_alpha2  mean of max", "=",  mean(max_ICA_alpha2)))
print(paste("ICA_alpha25 mean of max", "=",  mean(max_ICA_alpha25)))
print(paste("ICA_alpha3 mean of max", "=",  mean(max_ICA_alpha3)))

##mean of the sd
print(paste("ICA_KL mean of sd", "=", mean(sd_ICA_KL)))
print(paste("ICA_alpha01 mean of sd", "=",  mean(sd_ICA_alpha01)))
print(paste("ICA_alpha05 mean of sd", "=",  mean(sd_ICA_alpha05)))
print(paste("ICA_alpha125 mean of sd", "=",  mean(sd_ICA_alpha125)))
print(paste("ICA_alpha15 mean of sd", "=",  mean(sd_ICA_alpha15)))
print(paste("ICA_alpha2  mean of sd", "=",  mean(sd_ICA_alpha2)))
print(paste("ICA_alpha25 mean of sd", "=",  mean(sd_ICA_alpha25)))
print(paste("ICA_alpha3 mean of sd", "=",  mean(sd_ICA_alpha3)))


#range
print(paste("ICA_KL range", "=", mean(max_ICA_KL)-mean(min_ICA_KL)  ))
print(paste("ICA_alpha01 range", "=",  mean(max_ICA_alpha01)-mean(min_ICA_alpha01)  ))
print(paste("ICA_alpha05 range", "=",  mean(max_ICA_alpha05)-mean(min_ICA_alpha05)  ))
print(paste("ICA_alpha125 range", "=",  mean(max_ICA_alpha125)-mean(min_ICA_alpha125)  ))
print(paste("ICA_alpha15 range", "=",  mean(max_ICA_alpha15)-mean(min_ICA_alpha15)  ))
print(paste("ICA_alpha2  range", "=", mean(max_ICA_alpha2)-mean(min_ICA_alpha2)  ))
print(paste("ICA_alpha25 range", "=",  mean(max_ICA_alpha25)-mean(min_ICA_alpha25)  ))
print(paste("ICA_alpha3 range", "=",  mean(max_ICA_alpha3)-mean(min_ICA_alpha3)  ))

#median
print(paste("ICA_KL mean of median", "=", mean(median_ICA_KL)))
print(paste("ICA_alpha01 mean of median", "=",  mean(median_ICA_alpha01)))
print(paste("ICA_alpha05 mean of median", "=",  mean(median_ICA_alpha05)))
print(paste("ICA_alpha125 mean of median", "=",  mean(median_ICA_alpha125)))
print(paste("ICA_alpha15 mean of median", "=",  mean(median_ICA_alpha15)))
print(paste("ICA_alpha2  mean of median", "=",  mean(median_ICA_alpha2)))
print(paste("ICA_alpha25 mean of median", "=",  mean(median_ICA_alpha25)))
print(paste("ICA_alpha3 mean of median", "=",  mean(median_ICA_alpha3)))