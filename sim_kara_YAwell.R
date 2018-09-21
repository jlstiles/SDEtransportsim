# library(cateSurvival)
# devtools::install_github("jlstiles/SDE_transport")
# library(Simulations)
library(SDEtransport)
# Functions to generate data for transport:

# make W different for the two sites:
n=1e6

f_W = function(n) {
  W1 = rbinom(n, 1, 0.5)
  W2 = rbinom(n, 1, 0.4 + 0.2*W1)
  data.frame(W1=W1, W2=W2)
}
W = f_W(n)
f_S = function(W) {
  with(W, plogis(W1*W2 + .2*W2 + .7))
}
P_SW = f_S(W)
S = rbinom(n,1,P_SW)
mean(S)
hist(P_SW, breaks = 200)
# make a pscore model

f_A = function(S,W) {
  df = cbind(S=S, W)
  with(df, plogis(-.6*S-.7*W1-W2 +1))
}

pscores = f_A(S,W)
max(pscores)
min(pscores)
hist(pscores, 200)
A = rbinom(n, 1, pscores)
mean(A)
# make a intermediate confounder model

f_Z = function(A,S,W) {
  df = cbind(S=S, W, A = A)
  # with(df, plogis(.1*S-.4*W1+.3*W2+1*A-.3))
  with(df, plogis(.1*S-.4*W1*S+ 1*W2*A-.3 + A))
}

pzscores = f_Z(A,S,W)
hist(pzscores,200)
Z = rbinom(n, 1, pzscores)

# make an M model according to the restrictions

# make an M model according to the restrictions
f_M = function(Z,W,S) {
  df = cbind(S=S, W, Z = Z)
  with(df, plogis(1*W1*W2 + 1.2*Z*W1 + .2*Z + .1))
  # with(df, plogis(-.14*S + 1*W1-.5*W2 + 1.2*Z +.1))
}
Mscores = f_M(Z,W,S)
hist(Mscores, 200)
M = rbinom(n, 1, Mscores)

# make a Y model according to the restrictions
f_Y = function(M,Z,W) {
  df = cbind(M=M, Z = Z, W)
  # with(df, plogis(1*M*Z + 1.5*W2*Z + 1*Z*W2*W1 - .7*M*W1 + .6))
  # with(df, plogis(1*M  + Z - .68*W - 1))
  # with(df, plogis(1*M  - Z*M - .68*W1^2*Z - .3*Z*W2 - (W2 > .7)*M*Z + (W1 < -.4)*M - 1))
  # with(df, plogis(1*M + 1.5*W + 1*Z - 1))
  # with(df, plogis(1*M -.5*cos(W) +.68*cos(W)*Z - .38*Z*M + (W < -.4)*M - 1))
  # with(df, plogis(1*M + 1.5*W + 1*Z - 1))
  with(df, plogis(1*M + 1.5*W1 - .37*W2 + .3*Z - 1))
  # with(df, plogis(W * M + 3 * cos(W) * Z - 0.4))
  # with(df, plogis(.3*M*W2^2-.3*W1*W2+.2*cos(W1)*Z-.3*Z*W1^2+.3))
}

Yscores = f_Y(M,Z,W)
Y = rbinom(n, 1, Yscores)
hist(Yscores, 200)
min(Y*log(Yscores)+(1-Y)*log(1-Yscores))
mean(Y)
# pack these functions into a DGP
func_list = list(f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
covariates = list(covariates_S = c("W1","W2"),
                  covariates_A = c("S","W1","W2"),
                  covariates_Z = c("S","A","W1","W2"),
                  covariates_M = c("Z","W1","W2"),
                  covariates_Y = c("M","Z","W1","W2"),
                  covariates_QZ = c("S","W1","W2"))


sim_kara = function(n, covariates, truth) {
  data = gendata.SDEtransport(n, 
                              f_W = truth$f_W, 
                              f_S = truth$f_S, 
                              f_A = truth$f_A, 
                              f_Z = truth$f_Z, 
                              f_M = truth$f_M, 
                              f_Y = truth$f_Y)
  SDE_tmle4(data, sl = NULL, covariates= covariates, truth = truth,
            truncate = list(lower =.0001, upper = .9999), glm_only = TRUE,
            B=500)
}

library(parallel)

B = 1000
n=100

res100_YAwell = mclapply(1:B, FUN = function(x) sim_kara(n, covariates, func_list), 
                       mc.cores = getOption("mc.cores", 24L))

save(res100_YAwell, func_list, covariates, file = "results/res100_YAwell.RData")

B = 1000
n=500

res500_YAwell = mclapply(1:B, FUN = function(x) sim_kara(n, covariates, func_list), 
                       mc.cores = getOption("mc.cores", 24L))

save(res500_YAwell, func_list, covariates, file = "results/res500_YAwell.RData")

B = 500
n=5000

res5000_YAwell = mclapply(1:B, FUN = function(x) sim_kara(n, covariates, func_list), 
                        mc.cores = getOption("mc.cores", 24L))

save(res5000_YAwell, func_list, covariates, file = "results/res5000_YAwell.RData")

B = 500
n=5000

res5000_YAwell = mclapply(1:B, FUN = function(x) sim_kara(n, covariates, func_list), 
                          mc.cores = getOption("mc.cores", 24L))

save(res5000_YAwell, func_list, covariates, file = "results/res5000_YAwell1.RData")