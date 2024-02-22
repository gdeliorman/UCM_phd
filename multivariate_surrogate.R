##multivariate surrogate 

install.packages("Surrogate")
library(Surrogate)

data("ARMD.MultS")
data_multi<- ARMD.MultS

placebo<-data_multi[data_multi$Treat == "-1",]
experimental<-data_multi[data_multi$Treat == "1",]

Tr<- data_multi$Diff52
S1<-data_multi$Diff4
S2<-data_multi$Diff12
S3<-data_multi$Diff24


T0<- placebo$Diff52
T1<- experimental$Diff52
S10<- placebo$Diff4
S11<- experimental$Diff4
S20<- placebo$Diff12
S21<- experimental$Diff12
S30<- placebo$Diff24
S31<- experimental$Diff24

##correlastions
cor(cbind(T0, S10, S20, S30))
cor(cbind(T1, S11, S21, S31))

##multiple adjusted association:

##MODELS
Model_T <- lm(Diff52~Treat, data=ARMD.MultS)
summary(Model_T)
Model_S1 <- lm(Diff4~Treat, data=ARMD.MultS)
summary(Model_S1)
Model_S2 <- lm(Diff12~Treat, data=ARMD.MultS)
summary(Model_S2)
Model_S3 <- lm(Diff24~Treat, data=ARMD.MultS)
summary(Model_S3)


##residuals
Res_T <- residuals(lm(Diff52~Treat, data=ARMD.MultS))
Res_S1 <- residuals(lm(Diff4~Treat, data=ARMD.MultS))
Res_S2 <- residuals(lm(Diff12~Treat, data=ARMD.MultS))
Res_S3 <- residuals(lm(Diff24~Treat, data=ARMD.MultS))
Residuals <- cbind(Res_T, Res_S1, Res_S2, Res_S3)

#Make covariance matrix of residuals
Sigma_gamma <- cov(Residuals)

##all covariances
dims<- dim(Sigma_gamma)[1] 
sigma_TT <- matrix(data = Sigma_gamma[1,1], nrow = 1)
sigma_ST <- matrix(data = Sigma_gamma[2:dims,1], nrow = (dims-1))
sigma_SS <- matrix(data = Sigma_gamma[2:dims,2:dims], nrow = (dims-1))

Gamma.Delta <- as.numeric((t(sigma_ST) %*% solve(sigma_SS) %*% sigma_ST) / sigma_TT)


# Compute the multiple surrogate adjusted association
Result <- AA.MultS(Sigma_gamma = Sigma_gamma, N = 188, Alpha = .05)

# Explore results
summary(Result)

N<-188
Alpha<-0.05
sd_val <- sqrt((4*Gamma.Delta*(1-Gamma.Delta)^2)/(N-3))
lb <- max(0, as.numeric(as.numeric(Gamma.Delta) + qnorm(Alpha/2) * (sd_val)))
ub <- min(1, as.numeric(Gamma.Delta + qnorm(1-Alpha/2)*(sd_val)))
Gamma.Delta_Results <- data.frame(cbind(Gamma.Delta, sd_val, lb, ub), stringsAsFactors = TRUE)
colnames(Gamma.Delta_Results) <- c("Multivariate AA", "Standard Error", "CI lower limit", "CI upper limit")
rownames(Gamma.Delta_Results) <- c(" ") 


###ICA multiple surrogate
all_possible_cors<- seq(from=-1, to =1, by=0.1)

delta_T<- matrix(NA, ncol=1, nrow = length(all_possible_cors) )
for (i in 1:length(all_possible_cors)) {
  delta_T[i,1]<-var(T0)+var(T1)+2*all_possible_cors[i]*sqrt(var(T0)+var(T1))
}

mean(delta_T)
max(delta_T)
min(delta_T)

Sigma_matrix_0<- var(cbind(T0, S10, S20, S30))
Sigma_matrix_1<- var(cbind(T1, S11, S21, S31))

Sigma_matrix<- matrix(NA, nrow = 8, ncol = 8)
Sigma_matrix[1,1]<- Sigma_matrix_0[1,1]
Sigma_matrix[2,2]<- Sigma_matrix_1[1,1]
Sigma_matrix[3,3]<- Sigma_matrix_0[2,2]
Sigma_matrix[4,4]<- Sigma_matrix_1[2,2]
Sigma_matrix[5,5]<- Sigma_matrix_0[3,3]
Sigma_matrix[6,6]<- Sigma_matrix_1[3,3]
Sigma_matrix[7,7]<- Sigma_matrix_0[4,4]
Sigma_matrix[8,8]<- Sigma_matrix_1[4,4]

Sigma_matrix[1,3]<-Sigma_matrix[3,1]<- Sigma_matrix_0[2,1]
Sigma_matrix[1,5]<-Sigma_matrix[5,1]<- Sigma_matrix_0[3,1]
Sigma_matrix[1,7]<-Sigma_matrix[7,1]<- Sigma_matrix_0[4,1]

Sigma_matrix[2,4]<-Sigma_matrix[4,2]<- Sigma_matrix_1[2,1]
Sigma_matrix[2,6]<-Sigma_matrix[6,2]<- Sigma_matrix_1[3,1]
Sigma_matrix[2,8]<-Sigma_matrix[8,2]<- Sigma_matrix_1[4,1]

Sigma_matrix[3,5]<-Sigma_matrix[5,3]<- Sigma_matrix_0[2,3]
Sigma_matrix[3,7]<-Sigma_matrix[7,3]<- Sigma_matrix_0[2,4]
Sigma_matrix[5,7]<-Sigma_matrix[7,5]<- Sigma_matrix_0[3,4]

Sigma_matrix[4,6]<-Sigma_matrix[6,4]<- Sigma_matrix_1[2,3]
Sigma_matrix[4,8]<-Sigma_matrix[8,4]<- Sigma_matrix_1[2,4]
Sigma_matrix[6,8]<-Sigma_matrix[8,6]<- Sigma_matrix_1[3,4]

Sigma<-Sigma_matrix

ICA_multi<- ICA.ContCont.MultS(M = 500, 188, Sigma, 
                   G = seq(from=-1, to=1, by = .1),
                   Seed=c(15), Show.Progress=TRUE)


hist(ICA_multi$R2_H)
ICA_multi$Corr.R2_H

ICA_multi$Call

ICA_multi_PC<- ICA.ContCont.MultS.PC(M=500,N,Sigma,Seed=5,Show.Progress=TRUE)
ICA_multi_PC$R2_H
hist(ICA_multi_PC$R2_H)
ICA_multi_PC$Corr.R2_H
hist(ICA_multi_PC$Corr.R2_H)
ICA_multi_PC$Lower.Dig.Corrs.All
