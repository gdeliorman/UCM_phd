##Numerical examples to calculate ICA for normal case

###packages
library(extraDistr) #for bivariate dist
library(mvtnorm) ## to generate multivariate data
library(MASS)#to generate multinormal data
library(dplyr) ##
library(LaplacesDemon)
library(pracma) ##
library(compositions) ##multiple log normal
library(Metrics) ##FOR BiaS
library(devtools) ## local functiom
library(ks) ###kernel package
library(Surrogate)
library(compositions) ##generate multi log normal
library(kdensity) ##
library(akima) #
library(KernSmooth)


##high parameters to generate 4-dimensional data 
sigma<-rbind(c(544.3285, 0, 300.9816, 0), c(0,550.6597,0, 304.4222), c(300.9816,0, 180.6831, 0), c(0,304.4222, 0,180.9433))
mu<-c(-11.47191,-16.0151,-6.646067,-8.938326) 

##generate data from normal distribution
data_log2<-as.data.frame(mvrnorm(n=1e5, mu=mu, Sigma=sigma), empirical = TRUE)


T0<-data_log2[,1]
T1<-data_log2[,2]
S0<-data_log2[,3]
S1<-data_log2[,4]
treatment<-rbinom(1e5,1,0.5)

##unobserved data
deltaT<- T1-T0
deltaS<- S1-S0

unobserved_data<-data.frame(T0,T1,S0,S1, treatment)
rho<- cor(deltaT, deltaS)
##first step: use unobserved data
##ICA: correlation square suggested by Alonso 2015 for normal case
ICA_rho<- rho*rho


##use joint and marginal distribution of DeltaT and DeltaS
fx<-function(x) dnorm(x,mean = mean(deltaT),sd=sd(deltaT) )
fy<-function(y) dnorm(y,mean =mean(deltaS),sd=sd(deltaS))
fxy<- function(x,y) dbvnorm(x,y, mean1 = mean(deltaT), mean2 = mean(deltaS), sd1 =sd(deltaT) , sd2 = sd(deltaS), cor=cor(deltaT, deltaS))
KLfun<- function(x,y)  fxy(x,y)* ( log(  (fxy(x,y)) / (fx(x)*fy(y)) ))


##use kernel estimation of DeltaT and DeltaS
kdeT<- kdensity(deltaT, start = "gaussian", kernel= "gaussian", bw=dpih(deltaT))
kdeS<- kdensity(deltaS, start = "gaussian", kernel= "gaussian", bw=dpih(deltaS))
fT<- kdeT
fS<- kdeS

###bivariate kernel
est <- bkde2D(cbind(deltaT, deltaS), bandwidth=c(dpih(deltaT),dpih(deltaS)), range.x = list(c(xmin,xmax), c(ymin,ymax)) , gridsize = c(2000L, 2000L)) 

approxfun2 <- function(x=est$x1, y=est$x2, z=est$fhat, method = "linear") {
  function(xp, yp) interp2(est$x1, est$x2, est$fhat, xp, yp, method) }
fTS <- approxfun2(x=est$x1, y=est$x2, z=est$fhat)

##we can test some points
fxy(-1,-3)
fTS(-1,-3)
fxy(0,0)
fTS(0,0)

##plot of DeltaT and DeltaS: real distribution vs kernel estimate 
plot(fx, xlim=c(-100,100), main="T1-T0")
lines(kdeT, col="2")
legend("topleft", legend = paste(c("normal","kernel") ), col = 1:2, pch = 19, bty = "n")

plot(fy, xlim=c(-100,100), main="S1-S0")
lines(kdeS, col="2")
legend("topleft", legend = paste(c("normal","kernel") ), col = 1:2, pch = 19, bty = "n")


##to find bound (max and min values of DeltaT and DeltaS)
fT_bound<-matrix(NA, ncol=2)
for (i in -500:500) {
  
  if(fT(i)== 0.000000e+00) {
    
  }
  
  else if (fT(i)<0.001)
  {
    
    c<-cbind(i, fT(i))
    fT_bound<-rbind(c,fT_bound)
  }
}

###finding boundaries fS
fS_bound<-matrix(NA, ncol=2)
for (i in -2000:1000) {
  if(fS(i)== 0.000000e+00) {
    
  }
  
  else if (fS(i)<0.001)
  {
    
    c<-cbind(i,fS(i))
    fS_bound<-rbind(c,fS_bound)
    
  }
}


###set xmin ymin xmax ymax
if (all(fT_bound[,1] >= 0) ){
  ##write function
} else if (all(fT_bound[,1] <= 0)){
}  else {
  xmin= (max(fT_bound[(which(fT_bound[,1]<0))]))
}

if (all(fS_bound[,1] >= 0) ){
  ##write function
} else if (all(na.omit(fS_bound[,1]) <= 0)){
}  else {
  ymin= (max(fS_bound[(which(na.omit(fS_bound[,1])<0))]))
}


if (all(fT_bound[,1] >= 0) ){
  xmax= max(as.data.frame(na.omit(fT_bound[,1]))) ##write function
} else if (all(na.omit(fT_bound[,1]) <= 0)){
}  else {
  xmax<-(min(fT_bound[(which(fT_bound[,1]>0))]))
}

if (all(fS_bound[,1] >= 0) ){
  ymax= max(as.data.frame(na.omit(fS_bound[,1]))) ##write function
} else if (all(na.omit(fS_bound[,1]) <= 0)){
}  else {
  ymax<-(min(fS_bound[(which(na.omit(fS_bound[,1])>0))]))
}


###ICA calculation with Information Correlation Coefficient (ICC) using normal functions
KLfor<-integral2(KLfun, xmin, xmax, ymin, ymax,  reltol = 1e-10) 
ICA_normal<-1-exp(-2*KLfor$Q)

##it can be also calculated the following numerical method:##grids for integral
point_xx<- runif(n=100000, min=xmin, max=xmax)
point_yy<- runif(n=100000, min=ymin, max=ymax)

###elimination and asign 0
elimination1<- which(fxy(point_xx, point_yy)< 10^(-6) & fx(point_xx)*fy(point_yy)< 10^(-6))
elimination2<- which(fxy(point_xx, point_yy)< 10^(-6) & fx(point_xx)*fy(point_yy)>= 10^(-6) ) 

##ICA_1=ICA_KL
KLfunxy<- function(x,y)  fxy(x,y)* ( log(  (fxy(x,y)) / (fx(x)*fy(y)) ))
fullvalue<- (KLfunxy(point_xx, point_yy))

#fullvalue[c(elimination1,elimination2)]<-0
fullvalue[c(elimination2)]<-0

xmin<- min(point_xx)
xmax<- max(point_xx)
ymin<- min(point_yy)
ymax<- max(point_yy)

int<-mean(fullvalue)*(xmax-xmin)*(ymax-ymin)
ICA_KL<-1-exp(-2*mean(fullvalue)*(xmax-xmin)*(ymax-ymin)) 

##ICA calculation using mutual information directly from mutinfo() function
mutual_info<- mutinfo(deltaT, deltaS, k=10, direct=FALSE)
mutual_info<- mean(mutual_info)
ICA_mutual_info<- 1-exp(-2*mutual_info)


###ICA calculation with Information Correlation Coefficient (ICC) using kernel functions
###elimination and asign 0
elimination1<- which(fTS(point_xx, point_yy)< 10^(-6) & fT(point_xx)*fS(point_yy)< 10^(-6))
elimination2<- which(fTS(point_xx, point_yy)< 10^(-6) & fT(point_xx)*fS(point_yy)>= 10^(-6) ) 

##ICA_1=ICA_KL
KLfunTS<- function(x,y)  fTS(x,y)* ( log(  (fTS(x,y)) / (fT(x)*fS(y)) ))
fullvalue<- (KLfunTS(point_xx, point_yy))

fullvalue[c(elimination1,elimination2)]<-0
int<-mean(na.omit(fullvalue))*(xmax-xmin)*(ymax-ymin)
ICA_kernel<-1-exp(-2*mean(fullvalue)*(xmax-xmin)*(ymax-ymin)) 

##second step: use observed data
##ICA: sensitivity analysis suggested by Alonso (2015) for normal case ICA.ContCont()

##select 2 column for observed data
newdata<-matrix(NA,nrow= 1e5, ncol=3)
for (i in 1:100000) {
  if(unobserved_data$treatment[i]==0)
  {
    newdata[i,]<-c(unobserved_data$T0[i], unobserved_data$S0[i], unobserved_data$treatment[i])
  }
  else
    newdata[i,]<-cbind(unobserved_data$T1[i], unobserved_data$S1[i], unobserved_data$treatment[i])
}

observed_data<-as.data.frame(newdata)
placebo<-observed_data[observed_data$V3 == "0",]
experimental<-observed_data[observed_data$V3== "1",]

#estiamble parameters from data
Tr<-observed_data$V1
S<-observed_data$V2
T0<- placebo$V1
T1<- experimental$V1
S0<- placebo$V2
S1<- experimental$V2

beta<-mean(T1)-mean(T0)
alpha<-mean(S1)-mean(S0)

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

##find unestiamble cor and variances, acoording to grid matrix G={-1,...,1}
realICA<- ICA.ContCont(cor(T0,S0), cor(T1,S1), var(T0,T0), var(T1,T1), var(S0,S0),var(S1,S1), T0T1=seq(-1, 1, by=.1), 
                       T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))

##desc. statistics
mean(realICA$ICA*realICA$ICA)
median(realICA$ICA*realICA$ICA)
max(realICA$ICA*realICA$ICA)
min(realICA$ICA*realICA$ICA)
sd(realICA$ICA*realICA$ICA)
hist(realICA$ICA*realICA$ICA, main="Histogram of ICA", xlab="")

##ICA_1: sensitivity analysis and use normal joint and marginal functions of DeltaT and DeltaS
rhodelta<-realICA$ICA
Ica<-rhodelta*rhodelta
pdf<-cbind(realICA$Pos.Def,rhodelta,Ica)


sigmadelta_ss<-list()
sigmadelta_tt<-list()
truefalseresult<-list()
ICAformul1<-list()

for (i in 1:nrow(ICA2$Pos.Def) )
  
{
  
  ### ICA values calculated agian for double check
  pay<-(sqrt(vT0T0*vS0S0)*cT0S0)+(sqrt(vT1T1*vS1S1)*cT1S1)-(sqrt(vT1T1*vS0S0)*ICA2$Pos.Def[i,4])-(sqrt(vT0T0*vS1S1)*ICA2$Pos.Def[i,3])
  payda<-sqrt((vT0T0+vT1T1-2*sqrt(vT0T0*vT1T1)*ICA2$Pos.Def[i,1])*(vS0S0+vS1S1-2*sqrt(vS0S0*vS1S1)*ICA2$Pos.Def[i,6]))
  ICAformul<-pay/payda
  
  ICAformul1<- append(ICAformul, ICAformul1)
  ICAformul1<-as.matrix(ICAformul1)
  ICAformul1<-as.numeric(ICAformul1)
  
  sigmadelta_s<-var(S0,S0)+var(S1,S1)- (2*ICA2$Pos.Def[i,6]*sqrt(var(S0,S0)*var(S1,S1)))
  sigmadelta_ss<- append(sigmadelta_s,sigmadelta_ss)
  sigmadelta_t<-var(T0,T0)+var(T1,T1)-2*ICA2$Pos.Def[i,1]*sqrt(var(T0,T0)*var(T1,T1))
  sigmadelta_tt<- append(sigmadelta_t,sigmadelta_tt)
  cormatrix<- matrix(c(1, ICA2$Pos.Def[i,1], ICA2$Pos.Def[i,2],ICA2$Pos.Def[i,3],
                       ICA2$Pos.Def[i,1],1,ICA2$Pos.Def[i,4],ICA2$Pos.Def[i,5],
                       ICA2$Pos.Def[i,2],ICA2$Pos.Def[i,4],1,ICA2$Pos.Def[i,6],
                       ICA2$Pos.Def[i,3],ICA2$Pos.Def[i,5],ICA2$Pos.Def[i,6], 1 ), nrow = 4)
  result<- is.positive.definite(cormatrix)
  truefalseresult<-append(result,truefalseresult)
}

###combine all results and check
ICAformul1<-matrix(c(ICAformul1), ncol = 1)
df2=ICAformul1[order(nrow(ICAformul1):1),] ##for order
pdf<- cbind(pdf,df2)

sigmadelta_tt<-as.matrix(sigmadelta_tt)
sigmadelta_ss<-as.matrix(sigmadelta_ss)

sigmadeltass<-cbind(sigmadelta_tt, sigmadelta_ss)
sigmadelta<-sigmadeltass[order(nrow(sigmadeltass):1),]

ICAS<-cbind(ICA2$Pos.Def, df2)
matrixforr<- cbind(sigmadelta_tt,   sigmadelta_ss,  ICAformul1, truefalseresult)
matrixforr<- matrixforr[order(nrow(matrixforr):1),]
fullmatrix<-cbind(pdf, sigmadelta)
print(fullmatrix)

##for using above results KL calculate
ICA_KL<-list()

for (i in 1:nrow(fullmatrix)  )
  
{
  
  fx<-function(x) dnorm(x,mean = beta,sd=sqrt( as.numeric(fullmatrix[i,10] )))
  fy<-function(y) dnorm(y,mean=alpha,sd=sqrt( as.numeric(fullmatrix[i,11])))
  fxy<- function(x,y) dbvnorm(x,y, mean1 = beta, mean2 = alpha, sd1 =sqrt( as.numeric(fullmatrix[i,10])) , sd2 = sqrt(as.numeric(fullmatrix[i,11])), cor =as.numeric(fullmatrix[i,7]))
  KLfun<- function(x,y)  fxy(x,y)* ( log(  (fxy(x,y)) / (fx(x)*fy(y)) ))
  
  
  ## bounds
  fT_bound<-matrix(NA, ncol=2)
  for (k in -200:200) {
    
    if (fx(k)>=0.0001)
    {
      
      c<-cbind(k, fx(k))
      fT_bound<-rbind(c,fT_bound)
    }
  }
  
  ###finding boundaries fS
  fS_bound<-matrix(NA, ncol=2)
  for (j in -200:200) {
    if (fy(j)>=0.0001)
    {
      
      c<-cbind(j,fy(j))
      fS_bound<-rbind(c,fS_bound)
      
    }
  }
  
  xmax<- max(na.omit(fT_bound[,1]))
  xmin<- min(na.omit(fT_bound[,1]))
  ymax<- max(na.omit(fS_bound[,1]))
  ymin<- min(na.omit(fS_bound[,1]))
  
  ##numeric integral
  
  point_xx<- runif(n=1000, min=xmin, max=xmax)
  point_yy<- runif(n=1000, min=ymin, max=ymax)
  
  ###elimination and asign 0
  elimination1<- which(fxy(point_xx, point_yy)< 10^(-6) & fx(point_xx)*fy(point_yy)< 10^(-6))
  elimination2<- which(fxy(point_xx, point_yy)< 10^(-6) & fx(point_xx)*fy(point_yy)>= 10^(-6) ) 
  elimination3<- which(fxy(point_xx, point_yy)>=10^(-6) & (fx(point_xx)*fy(point_yy)< 10^(-6) ))
  
  ##fullvalues
  fullvalue_KL<- (KLfun(point_xx, point_yy))
  fullvalue_KL[c(elimination2, elimination1)]<-0
  
  
  if(isempty(elimination3) ==TRUE ){
    
  } else {fullvalue_KL<- fullvalue_KL[-c(elimination3)] }
  
  try(R1<-integral2(KLfun, xmin, xmax, ymin, ymax, reltol = 1e-10))
  ICA_1<-1-exp(-2*R1$Q) 
  write(ICA_1, file="ICA_KL1.txt", append=TRUE)  
  
  
  int_KL<-mean(na.omit(fullvalue_KL))*(xmax-xmin)*(ymax-ymin)
  ICA_KL1<-1-exp(-2*int_KL) 
  write(ICA_KL1, file=" ICA_KL2.txt", append=TRUE)    }

##read text files
ICA_11<-read.delim(file= '/Users/gokcedeliorman/ICA_KL1.txt', header = FALSE, sep = "\t", dec = ".")
ICA_12<-read.delim(file= '/Users/gokcedeliorman/ ICA_KL2.txt', header = FALSE, sep = "\t", dec = ".")

ICA_11<- as.numeric(as.matrix(ICA_11))
ICA_12<- as.numeric(as.matrix(ICA_12))

mean(ICA_11)
mean(ICA_12)

hist(ICA_11, main = "Histogram of ICA_1", xlab="")
sd(ICA_11)
max(ICA_11)
min(ICA_11)
median(ICA_11)

hist(ICA_12)
sd(ICA_12)
max(ICA_12)
min(ICA_12)


##way 3: generate new mormal data with the parameters from pdf matrices then use mutinfo() function

for (i in 1:nrow(fullmatrix)  )
  
{
  
  Sigma_22<- matrix(NA, nrow = 2, ncol = 2)
  Sigma_22[1,1]<- as.numeric(fullmatrix[i, 10])
  Sigma_22[2,2]<- as.numeric(fullmatrix[i, 11])
  Sigma_22[1,2]<-   Sigma_22[2,1]<- fullmatrix[i, 9]* sqrt(as.numeric(fullmatrix[i, 10]))* sqrt(as.numeric(fullmatrix[i, 11]))
  
  data_log22<-as.data.frame(mvrnorm(n=1e5, mu=c(beta, alpha), Sigma=Sigma_22), empirical = TRUE)
  delta_T2<-  data_log22[,1]
  delta_S2<-  data_log22[,2]
  
  mut_info2<- mutinfo(delta_S2, delta_T2, k=10, direct=FALSE)
  mut_info2<- mean(mut_info2)
  ICA_mutinfo2<- 1-exp(-2*mut_info2)
  
  write(ICA_mutinfo2, file="ICA_mutinfo.txt", append=TRUE)    }

##read this file 
ICA_mutinfo2<-read.delim(file= '/Users/gokcedeliorman/ICA_mutinfo.txt', header = FALSE, sep = "\t", dec = ".")
ICA_mutinfo2<- as.numeric(as.matrix(ICA_mutinfo2))
hist(ICA_mutinfo2, xlab="", main=expression(ICA[mutinfo]))
mean(ICA_mutinfo2)
sd(ICA_mutinfo2)
max(ICA_mutinfo2)
min(ICA_mutinfo2)
median(ICA_mutinfo2)

##way 4: use kernel density estimation fucntions instead of normal functions:
point_xx<- runif(n=10000, min=xmin, max=xmax)
point_yy<- runif(n=10000, min=ymin, max=ymax)

i<-343
for (i in 1:nrow(fullmatrix)  )
  
{
  
  Sigma_22<- matrix(NA, nrow = 2, ncol = 2)
  Sigma_22[1,1]<- as.numeric(fullmatrix[i, 10])
  Sigma_22[2,2]<- as.numeric(fullmatrix[i, 11])
  Sigma_22[1,2]<- Sigma_22[2,1]<- fullmatrix[i, 9]* sqrt(as.numeric(fullmatrix[i, 10]))* sqrt(as.numeric(fullmatrix[i, 11]))
  
  data_log22<-as.data.frame(mvrnorm(n=1e5, mu=c(beta, alpha), Sigma=Sigma_22), empirical = TRUE)
  delta_T2<-  data_log22[,1]
  delta_S2<-  data_log22[,2]
  
  kdeT2<- kdensity(deltaT, start = "gaussian", kernel= "gaussian", bw=dpih(delta_T2))
  kdeS2<- kdensity(deltaS, start = "gaussian", kernel= "gaussian", bw=dpih(delta_S2))
  fT2<- kdeT2
  fS2<- kdeS2
  
  ## bounds
  fT_bound2<-matrix(NA, ncol=2)
  for (k in -200:200) {
    
    if (fT2(k)>=0.0001)
    {
      
      c<-cbind(k, fT2(k))
      fT_bound2<-rbind(c,fT_bound2)
    }
  }
  
  ###finding boundaries fS
  fS_bound2<-matrix(NA, ncol=2)
  for (j in -200:200) {
    if (fS2(j)>=0.0001)
    {
      
      c<-cbind(j,fS2(j))
      fS_bound2<-rbind(c,fS_bound2)
      
    }
  }
  
  xmax<- max(na.omit(fT_bound2[,1]))
  xmin<- min(na.omit(fT_bound2[,1]))
  ymax<- max(na.omit(fS_bound2[,1]))
  ymin<- min(na.omit(fS_bound2[,1]))
  
  
  fT_bound2<-matrix(NA, ncol=2)
  for (j in -500:500) {
    
    if(fT2(j)== 0.000000e+00) {
      
    }
    
    else if (fT2(j)<0.001)
    {
      
      c<-cbind(j, fT2(j))
      fT_bound2<-rbind(c,fT_bound2)
    }
  }
  
  ###finding boundaries fS
  fS_bound2<-matrix(NA, ncol=2)
  for (k in -2000:1000) {
    if(fS2(k)== 0.000000e+00) {
      
    }
    
    else if (fS2(k)<0.001)
    {
      
      c<-cbind(k,fS2(k))
      fS_bound2<-rbind(c,fS_bound2)
      
    }
  }
  
  
  ###set xmin ymin xmax ymax
  if (all(fT_bound2[,1] >= 0) ){
    ##write function
  } else if (all(fT_bound2[,1] <= 0)){
  }  else {
    xmin= (max(fT_bound2[(which(fT_bound2[,1]<0))]))
  }
  
  if (all(fS_bound2[,1] >= 0) ){
    ##write function
  } else if (all(na.omit(fS_bound2[,1]) <= 0)){
  }  else {
    ymin= (max(fS_bound2[(which(na.omit(fS_bound2[,1])<0))]))
  }
  
  
  if (all(fT_bound2[,1] >= 0) ){
    xmax= max(as.data.frame(na.omit(fT_bound2[,1]))) ##write function
  } else if (all(na.omit(fT_bound2[,1]) <= 0)){
  }  else {
    xmax<-(min(fT_bound2[(which(fT_bound2[,1]>0))]))
  }
  
  if (all(fS_bound2[,1] >= 0) ){
    ymax= max(as.data.frame(na.omit(fS_bound2[,1]))) ##write function
  } else if (all(na.omit(fS_bound2[,1]) <= 0)){
  }  else {
    ymax<-(min(fS_bound2[(which(na.omit(fS_bound2[,1])>0))]))
  }
  
  ###bivariate kernel
  est2 <- bkde2D(cbind(delta_T2, delta_S2), bandwidth=c(dpih(delta_T2),dpih(delta_S2)), range.x = list(c(xmin,xmax), c(ymin,ymax)) , gridsize = c(2000L, 2000L)) 
  
  
  approxfun2 <- function(x=est2$x1, y=est2$x2, z=est2$fhat, method = "linear") {
    function(xp, yp) interp2(est2$x1, est2$x2, est2$fhat, xp, yp, method) }
  fTS2 <- approxfun2(x=est2$x1, y=est2$x2, z=est2$fhat)
  
  ###elimination and asign 0
  elimination1<- which(fTS2(point_xx, point_yy)< 10^(-6) & fT2(point_xx)*fS2(point_yy)< 10^(-6))
  elimination2<- which(fTS2(point_xx, point_yy)< 10^(-6) & fT2(point_xx)*fS2(point_yy)>= 10^(-6) ) 

  KLfunTS2<- function(x,y)  fTS2(x,y)* ( log(  (fTS2(x,y)) / (fT2(x)*fS2(y)) ))
  fullvalue2<- (KLfunTS2(point_xx, point_yy))
  
  fullvalue2[c(elimination1,elimination2)]<-0
  int<-mean(na.omit(fullvalue_2))*(xmax-xmin)*(ymax-ymin)
  ICA_kernel2<-1-exp(-2*mean(na.omit(fullvalue_2))*(xmax-xmin)*(ymax-ymin)) 
  
  write(ICA_kernel2, file="ICA_kernel2.txt", append=TRUE)    }


##read this file 
ICA_kernel2<-read.delim(file= '/Users/gokcedeliorman/ICA_kernel2.txt', header = FALSE, sep = "\t", dec = ".")
ICA_kernel2<- as.numeric(as.matrix(ICA_kernel2))
hist(ICA_kernel2, xlab="", main=expression(ICA[kernel]))
mean(ICA_kernel2)
sd(ICA_kernel2)
max(ICA_kernel2)
min(ICA_kernel2)
median(ICA_kernel2)



