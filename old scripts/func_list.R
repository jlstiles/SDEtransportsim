# library(cateSurvival)
# devtools::install_github("jlstiles/SDE_transport")
# library(Simulations)
library(SDEtransport)
# Functions to generate data for transport:

# func_list = sim_kara_results_noboot5$res5000_YMmis$func_list
forms = sim_kara_results_noboot5$res5000_YMmis$forms


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

f_Z = function (A, S, W) {
  df = cbind(S = S, W, A = A)
  with(df, plogis(-1*A - .2*S +1*A*W2 + 1))
}

func_list$f_Z = f_Z
pzscores = f_Z(A,S,W)
hist(pzscores,200)
Z = rbinom(n, 1, pzscores)
max(pzscores)
min(pzscores)
# make an M model according to the restrictions

f_M = function(Z,W,S) {
  df = cbind(S=S, W, Z = Z)
  with(df, plogis(1*Z - 3 + 2*Z*W2))
}
# func_list$f_M = f_M

Mscores = f_M(Z,W,S)
hist(Mscores, 200)
M = rbinom(n, 1, Mscores)
mean(M)
max(Mscores)
min(Mscores)

# f_Y = func_list$f_Y

f_Y = function (M, Z, W) {
  df = cbind(M = M, Z = Z, W)
  with(df, plogis(log(1.2) + 1*Z*M - .2*W2 - 1*W2*Z + 1.5*M))
}

Yscores = f_Y(M,Z,W)
Y = rbinom(n, 1, Yscores)
hist(Yscores, 200)
min(Y*log(Yscores)+(1-Y)*log(1-Yscores))
mean(Y)
min(Yscores)
max(Yscores)

func_list = list(f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)

func_list = func_formsYM$func_listYMmis
forms = func_formsYM$formsYMmis


func_formsYM$func_listYMmis$f_M

ic = get_trueIC(100000,truth = func_list, forms = forms)

max(ic$Hm_astar0a1_0)
max(ic$Hm_astar0a0_0)
max(ic$Hm_astar1a1_0)
max(ic$Hz_astar0a1_0)
max(ic$Hz_astar0a0_0)
max(ic$Hz_astar1a1_0)

f_M = function (Z, W, S) 
{
  df = cbind(S=S, W, Z = Z)
  with(df, plogis(1*Z + 6*W2*Z - 2*W2 - 2))
}

f_Y = function (M, Z, W) 
{
  df = cbind(M = M, Z = Z, W)
  with(df, plogis(log(1.2) + log(40) * Z - log(30) * M - log(1.2) * 
                    W2 - log(40) * W2 * Z))
}

func_list = sim_kara_resultsYZmis$res100_well$func_list
forms = sim_kara_resultsYZmis$res100_YZmis$covariates

p = sim_kara(100000, forms = forms, truth = func_list, B = NULL)
# p = sim_karaP(5000, forms = forms, truth = func_list, B = NULL, PS0 = PS0)
c(p$CI_SDE, p$CI_SDE_1s, p$CI_SDE_iptw, p$SDE_0)
(p$CI_SDE[3] - p$CI_SDE[2])/(2*1.96*p$SE_SDE_0)
c(p$CI_SIE, p$CI_SIE_1s, p$CI_SIE_iptw, p$SIE_0)
(p$CI_SIE[3] - p$CI_SIE[2])/(2*1.96*p$SE_SIE_0)

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
# 
# get_trueIC = function(n,truth, forms) {
#   data = gendata.SDEtransport(n, 
#                               f_W = truth$f_W, 
#                               f_S = truth$f_S, 
#                               f_A = truth$f_A, 
#                               f_Z = truth$f_Z, 
#                               f_M = truth$f_M, 
#                               f_Y = truth$f_Y)
#   get_gstarM_glm(data = data, truth = truth, forms = forms)
# }
