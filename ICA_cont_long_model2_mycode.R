#Setting 1: BPRS as a Surrogate for PANSS
##Linear model with random intercept

library(Surrogate)
library(matrixcalc) ##to check positive definite matrices

##from the sas joint model V
Sigma<- V<- matrix(c(0.03804 , NA,0.03401 ,NA,
                     NA, 0.04099, NA, 0.03640,
                     0.03401, NA,0.03485 , NA,
                     NA, 0.03640, NA, 0.03653), nrow=4, ncol=4)


##estimable corr
corr<- cov2cor(Sigma)


D_matrix<- matrix(c( 0.04034, NA, 0.03497,NA,
                     NA, 0.03470, NA,0.03412 ,
                     0.03497, NA,0.03492 , NA,
                     NA, 0.03412, NA,0.03800 ), nrow=4, ncol=4)

d_corr<- cov2cor(D_matrix)

ICA_1<- ICA.ContCont(T0S0=corr[1,3], T1S1=corr[2,4], T0T0=Sigma[1,1], T1T1=Sigma[2,2], S0S0=Sigma[3,3], S1S1=Sigma[4,4], T0T1=seq(-1, 1, by=.1), 
                     T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))


ICA_1<- ICA.ContCont(T0S0=corr[1,3], T1S1=corr[2,4], T0T0=Sigma[1,1], T1T1=Sigma[2,2], S0S0=Sigma[3,3], S1S1=Sigma[4,4], T0T1=seq(-1, 1, by=.2), 
                     T0S1=seq(-1, 1, by=.2), T1S0=seq(-1, 1, by=.2), S0S1=seq(-1, 1, by=.2))
rho_delta_1<- ICA_1$ICA
pos_def_1<-ICA_1$Pos.Def


ICA_2<- ICA.ContCont(T0S0=d_corr[1,3], T1S1=d_corr[2,4], T0T0=D_matrix[1,1], T1T1=D_matrix[2,2], S0S0=D_matrix[3,3], S1S1=D_matrix[4,4], T0T1=seq(-1, 1, by=.1), 
                     T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))

ICA_2<- ICA.ContCont(T0S0=d_corr[1,3], T1S1=d_corr[2,4], T0T0=D_matrix[1,1], T1T1=D_matrix[2,2], S0S0=D_matrix[3,3], S1S1=D_matrix[4,4], T0T1=seq(-1, 1, by=.2), 
                     T0S1=seq(-1, 1, by=.2), T1S0=seq(-1, 1, by=.2), S0S1=seq(-1, 1, by=.2))

pos_def_2<-ICA_2$Pos.Def
rho_delta_2<- ICA_2$ICA

p<-6
AR<- 0.6223
alpha<- (p-(p*AR)+2*AR)/(1+AR)
Q<- matrix(c(-1, 0,1,0,
             0, -1,0,1), nrow=2, ncol=4)
I6<- diag(6)

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}
r<-as.matrix(ar1_cor(6, AR), nrow=6, ncol=6)
p1 <- matrix(c(1,1,1,1,1,1), nrow=1, ncol=6)



#RH_square<- data.frame()
for (i in 1:length(rho_delta_1)) {
  for (j in 1:length(rho_delta_2)) {
    
    ##covariances
    cov_sigma_T0T1<-  pos_def_1[i,1]* sqrt(Sigma[1,1]*Sigma[2,2])
    cov_sigma_T0S1<-  pos_def_1[i,3]* sqrt(Sigma[1,1]*Sigma[4,4])
    cov_sigma_T1S0<-  pos_def_1[i,4]* sqrt(Sigma[2,2]*Sigma[3,3])
    cov_sigma_S0S1<-  pos_def_1[i,6]* sqrt(Sigma[3,3]*Sigma[4,4])
    
    
    cov_d_T0T1<-  pos_def_2[j,1]* sqrt(D_matrix[1,1]*D_matrix[2,2])
    cov_d_T0S1<-  pos_def_2[j,3]* sqrt(D_matrix[1,1]*D_matrix[4,4])
    cov_d_T1S0<-  pos_def_2[j,4]* sqrt(D_matrix[2,2]*D_matrix[3,3])
    cov_d_S0S1<-  pos_def_2[j,6]* sqrt(D_matrix[3,3]*D_matrix[4,4])
    
    
    Sigma[1,2]<- Sigma[2,1]<- cov_sigma_T0T1
    Sigma[1,4]<- Sigma[4,1]<- cov_sigma_T0S1
    Sigma[3,2]<- Sigma[2,3]<- cov_sigma_T1S0
    Sigma[3,4]<- Sigma[4,3]<- cov_sigma_S0S1
    
    
    
    V_1<- sigma_u <- matrix (c((Sigma[1,1])+(Sigma[2,2])-2*cov_sigma_T0T1,
                               (Sigma[1,3])+(Sigma[2,4])- cov_sigma_T1S0-cov_sigma_T0S1,
                               (Sigma[1,3])+(Sigma[2,4])- cov_sigma_T1S0-cov_sigma_T0S1,
                               (Sigma[3,3])+(Sigma[4,4])-2*cov_sigma_S0S1), nrow=2, ncol=2)
    
    
    
    
    D_matrix[1,2]<- D_matrix[2,1]<- cov_d_T0T1
    D_matrix[1,4]<- D_matrix[4,1]<- cov_d_T0S1
    D_matrix[3,2]<- D_matrix[2,3]<- cov_d_T1S0
    D_matrix[3,4]<- D_matrix[4,3]<- cov_d_S0S1
    
    
    
    D_1<- Q%*% D_matrix %*% t(Q)
    D_1<- round(D_1, digits = 6)
    
    u11 <- sigma_u[1, 1]
    u12 <- sigma_u[1, 2]
    u22 <- sigma_u[2, 2]
    d11 <- D_1[1, 1]
    d12 <- D_1[1, 2]
    d22 <- D_1[2, 2]
    
    
    rho_u <- u12 / sqrt(u11 * u22)
    rho_u2 <- rho_u^2
    
    denom_ud <- (u11 + alpha * d11) * (u22 + alpha * d22)
    
    rho_ud <- (u12 + alpha * d12) / sqrt(denom_ud)
    rho_ud2 <- rho_ud^2
    
    
    if (is.positive.definite(V_1)==TRUE & is.positive.definite(D_1)==TRUE ) {
      
      sigma_delta<- kronecker(D_1, t(p1)%*%(p1))+  kronecker(V_1, r)
      rh_square_check<- 1- (det(sigma_delta)/(det(sigma_delta[1:6, 1:6])*det(sigma_delta[7:12,7:12])))
      rh_square_check2 <- 1 - (1 - rho_u2)^(p - 1) * (1 - rho_ud2)
      write(as.numeric(rh_square_check),"~/Downloads/RH_square_log_set1_model_2_25march_det.txt",  append=TRUE, sep = "\t")
      write(as.numeric(rh_square_check2),"~/Downloads/RH_square_log_set1_model_2_25march_close.txt",  append=TRUE, sep = "\t")
      
    }
  }
}


RH_square_check_det<-read.delim(file= '~/Downloads/RH_square_log_set1_model_2_25march_det.txt', header = FALSE, sep = "\t", dec = ".")
RH_square_check2_close<-read.delim(file= '~/Downloads/RH_square_log_set1_model_2_25march_close.txt', header = FALSE, sep = "\t", dec = ".")




summary(RH_square_check_det$V1)
summary(RH_square_check2_close$V1)

hist(RH_square_check_det$V1)
hist(RH_square_check2_close$V1)

min(RH_square_check_det$V1)
min(RH_square_check2_close$V1)


median(RH_square_check_det$V1)
median(RH_square_check2_close$V1)


sd(RH_square_check_det$V1)
sd(RH_square_check2_close$V1)

mean(RH_square_check_det$V1)
mean(RH_square_check2_close$V1)


min(RH_square_check$V1)
quantile(RH_square_check$V1)
mean(RH_square_check$V1)
median(RH_square_check$V1)
length(RH_square_check$V1)

#hist( as.numeric(RH_square_check$V1), main="", xlab=bquote("R"[Lambda]^2),  breaks = 15, labels = TRUE,ylim = c(0,1200000), xlim=c(0,1))
hist( as.numeric(RH_square_check$V1), main="", xlab=bquote("R"[Lambda]^2),  breaks = 15, ylim = c(0,1200000), xlim=c(0,1))

hist( as.numeric(RH_square_check$V1), main="", xlab=bquote("R"[Lambda]^2),  breaks = 15, ylim = c(0,1180000), xlim=c(0,1))

grid()
box()
