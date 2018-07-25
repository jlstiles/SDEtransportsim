devtools::install_github("jeremyrcoyle/sl3")
devtools::install_github("jlstiles/SDE_transport")

n=1e6
W = rnorm(n)
f_S = function(W) plogis(W +.7)
P_SW = f_S(W)
S = rbinom(n,1,P_SW)
mean(S)

# make a pscore model
f_A = function(S,W) plogis(-.7*S-1*W +.7)
pscores = f_A(S,W)
max(pscores)
min(pscores)
hist(pscores, 200)
A = rbinom(n, 1, pscores)

# make a intermediate confounder model
f_Z = function(A,S,W) plogis(.1*S-1*W+2*A-.8)
pzscores = f_Z(A,S,W)
hist(pzscores,200)
Z = rbinom(n, 1, pzscores)

# make an M model according to the restrictions
f_M = function(Z,W,S) plogis(-.14*S + 1*W + 2*Z +.6)
Mscores = f_M(Z,W,S)
hist(Mscores, 200)
M = rbinom(n, 1, Mscores)

# make a Y model according to the restrictions
f_Y = function(M,Z,W) plogis(1*M + 1.5*W + 1*Z - 1)
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

# define learners in sl3
lglm = make_learner(Lrnr_glm)
lmean = make_learner(Lrnr_mean)
lbayesglm = make_learner(Lrnr_bayesglm)
lxgboost = make_learner(Lrnr_xgboost)
lrnr_stack = make_learner(Stack, list(lglm, lmean, lxgboost))
# lrnr_stack = make_learner(Stack, list(lglm))
# lrnr_stack = make_learner(Stack, list(lmean))
metalearner = make_learner(Lrnr_nnls)

# define the superlearner
sl <- Lrnr_sl$new(learners = lrnr_stack,
                  metalearner = metalearner)

# declare covariates for the learners
covariates=list(covariates_Mstar = c("W", "Z"),
                covariates_Y = c("M", "Z", "W"),
                covariates_M = c("S", "W", "Z"),
                covariates_Z = c("W", "A", "S"),
                covariates_A = c("S","W"),
                covariates_S = c("W"))

# generate the data
n = 1e3
data = gendata(n, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)

# run the tmle
# a is the intervention, a_star is for the stochastic intervention.  For now, the stochastic
# intervention is defined on S = 1 only for both M and Z but can add options easily for that

res = SDE_tmle(data = data, a = 1, a_star = 0, sl = sl, covariates = covariates)
# tmle est
res$est
# mle gcomp
res$est_mle
# IC mean for tmle
mean(res$IC)
# sl coefficients for the learners
res$SL_coef
