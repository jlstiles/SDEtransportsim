# library(cateSurvival)
# devtools::install_github("jlstiles/SDE_transport")
# library(Simulations)
library(SDEtransport)
# Functions to generate data for transport:

# make W different for the two sites:
n=1e6

f_W = function(n) {
  W1 = rbinom(n, 1, 0.5)
  W2 = rbinom(n, 1, 0.4 + 0.2 * W1)
  data.frame(W1 = W1, W2 = W2)
}

W = f_W(n)
f_S = function(W) {
  with(W, plogis(3*W2 - 1))
}

P_SW = f_S(W)
S = rbinom(n,1,P_SW)
mean(S)
hist(P_SW, breaks = 200)
max(P_SW)
min(P_SW)

f_A = function(S,W) {
  rep(.5, length(S))
}

pscores = f_A(S,W)
A = rbinom(n, 1, pscores)
mean(A)
# make an intermediate confounder model


f_Z = function(A,S,W) {
  df = cbind(S=S, W, A = A)
  with(df, plogis(A*log(40) - log(4)*W2 - log(10)*S))
}

pzscores = f_Z(A,S,W)
hist(pzscores,200)
Z = rbinom(n, 1, pzscores)
max(pzscores)
min(pzscores)
# make an M model according to the restrictions

f_M = function(Z,W,S) {
  df = cbind(S=S, W, Z = Z)
  with(df, plogis(-log(3) + log(20)*Z - log(5)*W2 + .1*S))
}
Mscores = f_M(Z,W,S)
hist(Mscores, 200)
M = rbinom(n, 1, Mscores)
mean(M)
max(Mscores)
min(Mscores)

f_Y = function(M,Z,W) {
  df = cbind(M=M, Z = Z, W)
  with(df, plogis(log(1.2)  + log(40)*Z  - log(100)*M - log(1.2)*W2 - log(80)*W2*Z))
}

Yscores = f_Y(M,Z,W)
Y = rbinom(n, 1, Yscores)
hist(Yscores, 200)
min(Y*log(Yscores)+(1-Y)*log(1-Yscores))
mean(Y)
min(Yscores)
max(Yscores)

# pack these functions into a DGP
func_list = list(f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)

forms = list(Sform = "S~W2", Aform = NULL, Zstarform = "Z ~ A+W2+S", Mstarform = "M ~ Z+W2", 
             Yform = "Y ~ Z", QZform = "Qstar_Mg ~ W2 + S")
# 
# covariates = list(covariates_S = c("W1","W2"),
#                   covariates_A = c("S","W1","W2"),
#                   covariates_Z = c("S","A","W1","W2"),
#                   covariates_M = c("Z","W1","W2"),
#                   covariates_Y = c("M","Z","W1","W2"),
#                   covariates_QZ = c("S","W1","W2"))


# this gives CI's for tmle, EE and iptw for SDE and SIE as well as the truths'
# debug(SDE_glm4)
# p = sim_kara(5000, forms, truth = func_list, B=NULL)
# c(p$CI_SDE, p$CI_SDE_1s,p$CI_SDE_iptw, SDE_0 = p$SDE_0, SE_SDE_0 = p$SE_SDE_0)
# 
# c(p$CI_SIE, p$CI_SIE_1s,p$CI_SIE_iptw, SIE_0 = p$SIE_0, SE_SIE_0 = p$SE_SIE_0)

# undebug(get_gstarM_glm)
# IC_info = get_trueIC(5000, truth = func_list, forms = forms)
# max(IC_info$Hm_astar0a1_0)
# max(IC_info$Hm_astar0a0_0)
# max(IC_info$Hm_astar0a0_0)
# 
# max(IC_info$Hz_astar0a1_0)
# max(IC_info$Hz_astar0a0_0)
# max(IC_info$Hz_astar0a0_0)
# 
# hist(IC_info$Hm_astar0a0_0)

sim_kara = function(n, forms, truth, B = NULL) {
  
  data = gendata.SDEtransport(n, 
                              f_W = truth$f_W, 
                              f_S = truth$f_S, 
                              f_A = truth$f_A, 
                              f_Z = truth$f_Z, 
                              f_M = truth$f_M, 
                              f_Y = truth$f_Y)
  SDE_glm4(data, truth = truth,
           truncate = list(lower =.0001, upper = .9999),
           B=B, forms = forms, RCT = 0.5)
}

library(parallel)

B = 1000
n=100

res100_Ymis = mclapply(1:B, FUN = function(x) sim_kara(n=100, forms=forms, truth=func_list, B = NULL), 
                        mc.cores = getOption("mc.cores", 20L))

save(res100_Ymis, func_list, forms, file = "results5/res100_Ymis.RData")

B = 1000
n=500

res500_Ymis = mclapply(1:B, FUN = function(x) sim_kara(n=500, forms=forms, truth=func_list, B = NULL), 
                        mc.cores = getOption("mc.cores", 20L))

save(res500_Ymis, func_list, forms, file = "results5/res500_Ymis.RData")

B = 1000
n=5000

res5000_Ymis = mclapply(1:B, FUN = function(x) sim_kara(n=5000, forms=forms, truth=func_list, B = NULL), 
                         mc.cores = getOption("mc.cores", 20L))

save(res5000_Ymis, func_list, forms, file = "results5/res5000_Ymis.RData")

