library(cateSurvival)

# Functions to generate data for transport:

# make W different for the two sites:
n=1e6
W = rnorm(n)
f_S = function(W) plogis(W +.7)
P_SW = f_S(W)
S = rbinom(n,1,P_SW)
mean(S)

# make a pscore model
f_A = function(S,W) plogis(-.2*S-.2*W+.7)
pscores = f_A(S,W)
max(pscores)
min(pscores)
hist(pscores, 200)
A = rbinom(n, 1, pscores)

# make a intermediate confounder model
f_Z = function(A,S,W) plogis(.1*S-.4*W+2*A-.8)
pzscores = f_Z(A,S,W)
hist(pzscores,200)
Z = rbinom(n, 1, pzscores)

# make an M model according to the restrictions
f_M = function(Z,W,S) plogis(-.14*S+.25*W+2*Z+.6)
Mscores = f_M(Z,W,S)
hist(Mscores, 200)
M = rbinom(n, 1, Mscores)

# make a Y model according to the restrictions
f_Y = function(M,Z,W) plogis(1*M+.25*W+1*Z-1)
Yscores = f_Y(M,Z,W)
Y = rbinom(n, 1, Yscores)
hist(Yscores, 200)

# pack these functions into a DGP
gendata = function(n, f_S, f_A, f_Z, f_M, f_Y) {
  W = rnorm(n)
  P_SW = f_S(W)
  S = rbinom(n,1,P_SW)

  # make a pscore model
  pscores = f_A(S,W)
  A = rbinom(n, 1, pscores)
  
  # make a intermediate confounder model
  pzscores = f_Z(A,S,W)
  Z = rbinom(n, 1, pzscores)
  
  # make an M model according to the restrictions
  Mscores = f_M(Z,W,S)
  M = rbinom(n, 1, Mscores)
  
  # make a Y model according to the restrictions
  Yscores = f_Y(M,Z,W)
  Y = rbinom(n, 1, Yscores)
  
  return(data.frame(W = W, S = S, A = A, Z = Z, M = M, Y = Y))
}

data_pop = gendata(5*1e6, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
mean(data_pop$Y)

data = gendata(1000, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)

# Let's generate the truth!!!
# figure out the chosen formula for M

Mfit = glm(M ~ Z + W, data = data[data$S==1, ], family = binomial())
Zfit = glm(Z ~ A + W + S, data = data, family = binomial())

f_M1 = function(Z, W, S) {
  # W = W[1:10]
  nn = length(W)
  dataM1 = data.frame(Z = rep(1,nn), W = W)
  predM1 = predict(Mfit, newdata = dataM1, type = 'response')
  dataM0 = data.frame(Z = rep(0,nn), W = W)
  predM0 = predict(Mfit, newdata = dataM0, type = 'response')
  dataZ = data.frame(A = rep(1,nn), W = W, S = rep(1, nn))
  predZ = predict(Zfit, newdata = dataZ, type = 'response')
  gM = predM1*predZ + predM0*(1 - predZ)
  return(gM)
}  

f_M0 = function(Z, W, S) {
  # W = W[1:10]
  nn = length(W)
  dataM1 = data.frame(Z = rep(1,nn), W = W)
  predM1 = predict(Mfit, newdata = dataM1, type = 'response')
  dataM0 = data.frame(Z = rep(0,nn), W = W)
  predM0 = predict(Mfit, newdata = dataM0, type = 'response')
  dataZ = data.frame(A = rep(0,nn), W = W, S = rep(1, nn))
  predZ = predict(Zfit, newdata = dataZ, type = 'response')
  gM = predM1*predZ + predM0*(1 - predZ)
  return(gM)
}  
# We can now set A to 1, and stochastically intervene for A = 0 to generate M and our outcome for pop

data_popS11 = gendata(2*1e6, f_S = f_S, f_A = ff, f_Z = f_Z, f_M = f_M1, f_Y = f_Y)
true11 = mean(data_popS11$Y[data_popS11$S==0])

data_popS10 = gendata(2*1e6, f_S = f_S, f_A = function(S,W) 1, f_Z = f_Z, f_M = f_M0, f_Y = f_Y)
true10 = mean(data_popS10$Y[data_popS10$S==0])

data_popS01 = gendata(2*1e6, f_S = f_S, f_A = function(S,W) 0, f_Z = f_Z, f_M = f_M1, f_Y = f_Y)
true01 = mean(data_popS01$Y[data_popS01$S==0])

data_popS00 = gendata(2*1e6, f_S = f_S, f_A = function(S,W) 0, f_Z = f_Z, f_M = f_M0, f_Y = f_Y)
true00 = mean(data_popS00$Y[data_popS00$S==0])

c(unlist(lapply(0:1, FUN = function(a_star) fcn.SDE(data, f_M0, f_M1, Mfit, Zfit, 0, a_star))),
unlist(lapply(0:1, FUN = function(a_star) fcn.SDE(data, f_M0, f_M1, Mfit, Zfit, 1, a_star))))

c(true00, true01, true10, true11)

fcn.SDE = function(data, f_M0, f_M1, Mfit, Zfit, a, a_star) {
  # computing a TMLE for each of these parameters
  # start with a regression of Y on the past given S = 1
  Yfit = glm(Y ~ M + Z + W, data = data[data$S==1,], family = binomial)
  
  # We have our Mfit, Zfit from before for S = 1, which determines the stoc int.\
  
  # treatment mechanism
  Afit = glm(A ~ W + S, data = data, family = binomial)
  
  # S mechanism
  Sfit = glm(S ~ W, data = data, family = binomial)
  
  # compute the first clever covariate
  
  if (a_star == 1) gstarM_ps = f_M1(Z=data$Z, W = data$W, S = data$S) else
    gstarM_ps = f_M0(Z=data$Z, W = data$W, S = data$S)
  
  # compute Z preds for S = 0
  dataS0 = data
  dataS0$S = 0
  Z0_ps = predict(Zfit,newdata = dataS0, type = 'response')
  
  
  # compute probs S = 0 given W
  S0preds = 1 - predict(Sfit, type = 'response')
  
  # compute M preds 
  M_ps = predict(Mfit, newdata = data, type = 'response')
  
  # compute Z preds 
  Z_ps = predict(Zfit, type = 'response')
  
  # compute A=1 preds for S = 1
  A_ps = predict(Afit, type = 'response')
  
  #compute prob S = 1 given W
  S1preds = 1 - S0preds
  
  # compute prob S = 0
  PS0 = mean(data$S == 0)
  
  # 1st clever cov
  dataM1 = dataM0 = data
  dataM1$M = 1
  dataM0$M = 0
  Hm = with(data, ((S == 1)*(A == a)*
                     ((M == 1)*gstarM_ps + (M == 0)*(1 - gstarM_ps))*
                     ((Z == 1)*Z0_ps + (Z == 0)*(1 - Z0_ps))*S0preds)/
              (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                 ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                 ((A == 1)*A_ps + (A == 0)*(1 - A_ps))*S1preds*PS0))
  
  Hm1 = with(dataM1, ((S == 1)*(A == a)*
                        ((M == 1)*gstarM_ps + (M == 0)*(1 - gstarM_ps))*
                        ((Z == 1)*Z0_ps + (Z == 0)*(1 - Z0_ps))*S0preds)/
               (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                  ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                  ((A == 1)*A_ps + (A == 0)*(1 - A_ps))*S1preds*PS0))
  
  Hm0 = with(dataM0, ((S == 1)*(A == a)*
                        ((M == 1)*gstarM_ps + (M == 0)*(1 - gstarM_ps))*
                        ((Z == 1)*Z0_ps + (Z == 0)*(1 - Z0_ps))*S0preds)/
               (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                  ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                  ((A == 1)*A_ps + (A == 0)*(1 - A_ps))*S1preds*PS0))
  
  Ypreds = predict(Yfit, newdata = data, type = 'response')
  Qfit = glm(data$Y ~ Hm - 1 + offset(qlogis(Ypreds)), family = binomial)
  eps = Qfit$coefficients
  
  # perform the stochastic intervention on Qstar, fcn of Z, W, S.  These are outcomes
  YpredsM1 = predict(Yfit, newdata = dataM1, type = 'response')
  YpredsM0 = predict(Yfit, newdata = dataM0, type = 'response')
  
  QstarM1 = plogis(qlogis(YpredsM1) + eps*Hm1)
  QstarM0 = plogis(qlogis(YpredsM0) + eps*Hm0)
  QstarM = QstarM1*gstarM_ps + QstarM0*(1 - gstarM_ps)
  
  # regress on Z,W,S
  QZfit = glm(QstarM ~ Z + W + S, data = data, family = 'binomial')
  
  # compute the clever covariate 2
  Apreds0_ps = predict(Afit, type = 'response')
  Apreds0 = data$A*Apreds0_ps + (1 - data$A)*(1 - Apreds0_ps)
  if (a == 1) Hz = (data$A == a)*(data$S == 0)/Apreds0/PS0 else
    Hz = (data$A == a)*(data$S == 0)/(1 - Apreds0)/PS0
  
  QZpreds = predict(QZfit, type = 'response')
  # take the mean over Z given S = 0, A = 1 
  data_a = data
  data_a$A = a
  Zpreds_a = predict(Zfit, newdata = data_a, type = 'response')
  dataZpreds_a = data
  dataZpreds_a$Z = Zpreds_a
  QZ = predict(QZfit, newdata = dataZpreds_a, type = 'response')

  QZfit_tmle = glm(QZpreds ~ Hz - 1 + offset(qlogis(QZ)), family = binomial)
  
  # update QZ 
  if (a == 1) Hza = (data$S==0)/Apreds0_ps else Hza = (data$S==0)/(1 - Apreds0_ps)
  QZstar = plogis(qlogis(QZ) + Hza*QZfit_tmle$coefficients)
  
  # compute the parameter estimate
  return(mean(QZstar[data$S==0]))
}

