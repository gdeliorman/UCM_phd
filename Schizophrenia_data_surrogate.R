install.packages("Surrogate")
install.packages("MASS")
install.packages("pracma")
install.packages("matrixcalc")
install.packages("extraDistr")


library(Surrogate)
library(MASS) ##for generate multinormal data
library(pracma) ##for integral
library(matrixcalc) ##for positive definite
library(extraDistr) ##for dvnorm

data("Schizo")

##remove NA
which(is.na(Schizo$PANSS)=='TRUE')

Schizo<- Schizo[-c(405,  705, 1358, 1719, 2111),]
placebo<-Schizo[Schizo$Treat == "-1",]
experimental<-Schizo[Schizo$Treat == "1",]

Tr<-Schizo$PANSS
S<-Schizo$BPRS

T0<- placebo$PANSS
T1<- experimental$PANSS
S0<- placebo$BPRS
S1<- experimental$BPRS

##means
mT0<-mean(T0)
mT1<-mean(T1)
mS0<-mean(S0)
mS1<-mean(S1)

beta<-mT1-mT0
alpha<- mS1-mS0
muu<- c(mT0, mT1, mS0, mS1)
mudelta<-c(beta,alpha)

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


sigma<- matrix( c(vT0T0, 0, vT0S0,0,
                  0,vT1T1,0,vT1S1,
                  vT0S0,0,vS0S0,0,
                  0,vT1S1,0,vS1S1), nrow = 4, ncol = 4)
A<-matrix(c(-1,0,1,0,0,-1,0,1),  ncol=4,nrow = 2) ##corresponding contrast matrix
sigmadelta<- A%*%sigma%*% t(A)

sigmadelta2<- matrix(c(1094.9882,602.108,
                       602.108,361.6264), nrow = 2)


ICA<-ICA.ContCont(cT0S0, cT1S1, vT0T0, vT1T1, vS0S0, vS1S1, T0T1=seq(-1, 1, by=.05), 
                  T0S1=seq(-1, 1, by=.05), T1S0=seq(-1, 1, by=.05), S0S1=seq(-1, 1, by=.05))

ICA2<-ICA.ContCont(cT0S0, cT1S1, vT0T0, vT1T1, vS0S0, vS1S1, T0T1=seq(-1, 1, by=.1), 
                   T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))

mean(ICA2$ICA)
##ICA features
rhodelta<-ICA2$ICA
Ica<-rhodelta*rhodelta

##prints
print( paste("mean ICA=",mean(ICA2$ICA)))
print( paste("median ICA=",median(ICA2$ICA)))
print( paste("maximum ICA=",max(ICA2$ICA)))
print( paste("minimum ICA=",min(ICA2$ICA)))
print( paste("sd ICA=",sd(ICA2$ICA)))

##positive definite matrices 
pdf<-cbind(ICA2$Pos.Def,rhodelta,Ica)



fx<-function(x) dnorm(x,mean = beta,sd=sqrt(1094.9882))
fy<-function(y) dnorm(y,mean =alpha,sd=sqrt(361.6264))
fxy<- function(x,y) dbvnorm(x,y, mean1 = beta, mean2 = alpha, sd1 =sqrt(1094.9882) , sd2 = sqrt(361.6264), cor = mean(rhodelta))



##PLOT Delta T and Delta S
plot(fx, xlim = c(-150,150), ylab = bquote(f(Delta[T])))
abline(v=c(-106,97), col=c("red", "red"), lty=c(3,3), lwd=c(2, 2))
text(-132,0.002, expression(f(Delta[T])<0.0001 ) )
text(130,0.002, expression(f(Delta[T])<0.0001 ) )

plot(fy, xlim = c(-100,100), ylab = bquote(f(Delta[S])))
abline(v=c(-64,59), col=c("red", "red"), lty=c(3,3), lwd=c(2, 2))
text(-85,0.005, expression(f(Delta[S])<0.0001 ) )
text(85,0.005, expression(f(Delta[S])<0.0001 ) )


##DIVERGENCES
ralpha_05<-0.5
ralpha_125<-1.25
ralpha_15<-1.5
ralpha_2<-2

##FUNCTIONS
KLfun<- function(x,y) fxy(x,y)* ( log( (fxy(x,y)) / (fx(x)*fy(y)) ))
RY_05<- function(x,y) fxy(x,y)^ralpha_05 * ( (fx(x)*fy(y)) ^(1-ralpha_05))
RY_125<- function(x,y) fxy(x,y)^ralpha_125 * ( (fx(x)*fy(y))^(1-ralpha_125))
RY_15<- function(x,y) fxy(x,y)^ralpha_15 * ( (fx(x)*fy(y)) ^(1-ralpha_15))
RY_2<- function(x,y) fxy(x,y)^ralpha_2 * ( (fx(x)*fy(y)) ^(1-ralpha_2))



##INTEGRALS with random limits
KL<-integral2(KLfun, -135, 135, -135, 135, reltol = 1e-10)
1-exp(-2*KL$Q)

R_05<-integral2(RY_05, -135, 135, -135, 135, reltol = 1e-10)
R05<-(1/ (ralpha_05-1))* log(R_05$Q)
1-exp(-2*R05)


R_125<-integral2(RY_125, -135, 135, -135, 135, reltol = 1e-10)
R125<-( 1/(ralpha_125-1))* log(R_125$Q)
1-exp(-2*R125)


R_15<-integral2(RY_15, -135, 135, -135, 135, reltol = 1e-10)
R15<- log(R_15$Q)* (1/ (ralpha_15-1))
1-exp(-2*R15)


R_2<-integral2(RY_2,-135, 135, -135, 135, reltol = 1e-10)
R2<- log(R_2$Q)* (1/ (ralpha_2-1))
1-exp(-2*R2)



##according to pdf matrix: sigma_delta_T and sigma_delta_S and check ICA again
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
  
  sigmadelta_s<-var(S0,S0)+var(S1,S1)-2*ICA2$Pos.Def[i,6]*sqrt(var(S0,S0)*var(S1,S1))
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

ICAformul1<-matrix(c(ICAformul1), ncol = 1)
df2=ICAformul1[order(nrow(ICAformul1):1),]
fullmatrix<-cbind(as.numeric(df2),as.numeric(sigmadelta_tt), as.numeric(sigmadelta_ss))

ralpha05<-0.5
ralpha125<-1.25
ralpha15<-1.5
ralpha2<-2


ICA_KL<-list()
ICA_alpha05<-list()
ICA_alpha125<-list()
ICA_alpha15<-list()
ICA_alpha2<-list()

for (i in 1:nrow(fullmatrix) ) {
  
  fx<-function(x) dnorm(x,mean = betas,sd=sqrt( as.numeric(fullmatrix[i,2])))
  fy<-function(y) dnorm(y,mean =alphas,sd=sqrt( as.numeric(fullmatrix[i,3])))
  fxy<- function(x,y) dbvnorm(x,y, mean1 = betas, mean2 = alphas, sd1 =sqrt( as.numeric(fullmatrix[i,2])) , sd2 = sqrt(as.numeric(fullmatrix[i,3])),   cor=fullmatrix[i,1])
  KLfun<- function(x,y)  fxy(x,y)* ( log(  (fxy(x,y)) / (fx(x)*fy(y)) ))
  
  RY05<- function(x,y) fxy(x,y)^ralpha05 * ( (fx(x)*fy(y)) ^(1-ralpha05))
  RY125<- function(x,y) (fxy(x,y)^ralpha125) * ( (fx(x)*fy(y)) ^(1-ralpha125))
  RY15<- function(x,y) (fxy(x,y)^ralpha15) * ( (fx(x)*fy(y)) ^(1-ralpha15) )
  RY2<- function(x,y) (fxy(x,y)^ralpha2) * ( (fx(x)*fy(y)) ^(1-ralpha2) )
  
  
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
  for (j in -100:100) {
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
  
  
  point_xx<- runif(n=10000, min=xmin, max=xmax)
  point_yy<- runif(n=10000, min=ymin, max=ymax)
  
  ###elimination and asign 0
  elimination1<- which(fxy(point_xx, point_yy)< 10^(-6) & fx(point_xx)*fy(point_yy)< 10^(-6))
  elimination2<- which(fxy(point_xx, point_yy)< 10^(-6) & fx(point_xx)*fy(point_yy)>= 10^(-6) ) 
  elimination3<- which(fxy(point_xx, point_yy)>=10^(-6) & (fx(point_xx)*fy(point_yy)< 10^(-6) ))
  
  ##fullvalues
  fullvalue_KL<- (KLfun(point_xx, point_yy))
  fullvalue_05<- ( RY05(point_xx, point_yy))
  fullvalue_125<- ( RY125(point_xx, point_yy))
  fullvalue_15<- ( RY15(point_xx, point_yy))
  fullvalue_2<- ( RY2(point_xx, point_yy))
  
  ele_KL<-which(is.na(fullvalue_KL))
  ele_05<-which(is.na(fullvalue_05))
  ele_125<-which(is.na(fullvalue_125))
  ele_15<-which(is.na(fullvalue_15))
  ele_2<-which(is.na(fullvalue_2))
  
  fullvalue_KL[c(elimination2, elimination1)]<-0
  fullvalue_05[c(elimination2, elimination1)]<-0
  fullvalue_125[c(elimination2, elimination1)]<-0
  fullvalue_15[c(elimination2, elimination1)]<-0
  fullvalue_2[c(elimination2,elimination1)]<-0
  
  
  if(isempty(elimination3) ==TRUE ){
    
  } else {fullvalue_KL<- fullvalue_KL[-c(elimination3)] 
  fullvalue_05<- fullvalue_05[-c(elimination3)]
  fullvalue_125<- fullvalue_125[-c(elimination3)]
  fullvalue_15<- fullvalue_15[-c(elimination3)]
  fullvalue_2<- fullvalue_2[-c(elimination3)] }
  
  int_KL<-mean(na.omit(fullvalue_KL))*(xmax-xmin)*(ymax-ymin)
  ICA_KL1<-1-exp(-2*int_KL) 
  write(ICA_KL1, file="all_ICAKL1_normal.txt", append=TRUE)  

  int_R05<- mean(na.omit(fullvalue_05))*(xmax-xmin)*(ymax-ymin)
  R05for<- (1/ (ralpha05-1))* log(int_R05) 
  ICA_alpha05for<-1-exp(-2*R05for)  
  write(ICA_alpha05for, file="all_ICA05_normal.txt", append=TRUE)  

  int_R125<- mean(na.omit(fullvalue_125))*(xmax-xmin)*(ymax-ymin)
  R125for<- (1/(ralpha125-1))* log(int_R125)
  ICA_alpha125for<-1-exp(-2*R125for) 
  write(ICA_alpha125for, file="all_ICA125_normal.txt", append=TRUE)  

  int_R15<- mean(na.omit(fullvalue_15))*(xmax-xmin)*(ymax-ymin)
  R15for<- (1/(ralpha15-1))* log(int_R15)
  ICA_alpha15for<-1-exp(-2*R15for)
  write(ICA_alpha15for, file="all_ICA15_normal.txt", append=TRUE)  

  
  int_R2<- mean(na.omit(fullvalue_2))*(xmax-xmin)*(ymax-ymin)
  R2for<- (1/(ralpha2-1))* log(int_R2)
  ICA_alpha2for<-1-exp(-2*R2for) 
  write(ICA_alpha2for, file="all_ICA2_normal.txt", append=TRUE)
  
  
  ##second way integral example
  R_2<-integral2(RY_2, xmin, xmax, ymin, ymax, reltol = 1e-10)
  R2<- log(R_2$Q)* (1/ (ralpha_2-1))
  1-exp(-2*R2)  }

  
  
