# library(cateSurvival)
devtools::install_github("jlstiles/SDE_transport")
library(SDEtransport)
# Functions to generate data for transport:

f_W = function(n) {
  W1 = rbinom(n, 1, 0.5)
  W2 = rbinom(n, 1, 0.4 + 0.2 * W1)
  data.frame(W1 = W1, W2 = W2)
}

W = f_W(n)
f_S = function(W) {
  with(W, plogis(6*W1+W2 - 3))
}
P_SW = f_S(W)

# make a pscore model

f_A = function(S,W) {
  df = cbind(S = S, W)
  with(df, plogis(-0.6 * S - 0.7 * W1 - W2 + 1))
}

pscores = f_A(S,W)
# make a intermediate confounder model

f_Z = function(A,S,W) {
  df = cbind(S=S, W, A = A)
  with(df, plogis(2 * S - 2 * W1 + 0.3 * W2 + 6* A - 2))
}

pzscores = f_Z(A,S,W)

# make an M model according to the restrictions

f_M = function(Z,W,S) {
  df = cbind(S=S, W, Z = Z)
  # with(df, plogis(5*Z*S*W1 - 4))
  with(df, plogis(-.14*S - 1*W1 + .5*W2 + 1.2*Z +.1))
}
Mscores = f_M(Z,W,S)

# make a Y model according to the restrictions
f_Y = function(M,Z,W) {
  df = cbind(M=M, Z = Z, W)
  with(df, plogis(6 * M * Z - 3))
}

Yscores = f_Y(M,Z,W)

# pack these functions into a DGP
func_list = list(f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)

covariates = list(covariates_S = c("W1","W2"),
                  covariates_A = c("S","W1", "W2"),
                  covariates_Z = c("S","A","W1","W2"),
                  covariates_M = c("Z","W1","W2"),
                  covariates_Y = c("M","Z","W1","W2"),
                  covariates_QZ = c("S","W1","W2"))

# make the superlearner using sl3
lglm = make_learner(Lrnr_glm)
lmean = make_learner(Lrnr_mean)
# lxgboost = make_learner(Lrnr_xgboost)
lrnr_stack = make_learner(Stack, list(lglm, lmean))
metalearner = make_learner(Lrnr_nnls)

# define the superlearner
sl <- Lrnr_sl$new(learners = lrnr_stack,
                  metalearner = metalearner)


n = 5e2
data = gendata.SDEtransport(n, 
                            f_W = truth$f_W, 
                            f_S = truth$f_S, 
                            f_A = truth$f_A, 
                            f_Z = truth$f_Z, 
                            f_M = truth$f_M, 
                            f_Y = truth$f_Y)
result = SDE_tmle(data, sl = sl, covariates= covariates, truth = NULL,
                   truncate = list(lower =.0001, upper = .9999), glm_only = FALSE,
                   B=NULL)
result 

