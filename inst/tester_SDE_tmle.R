devtools::install_github("jeremyrcoyle/sl3")
devtools::install_github("jlstiles/SDE_transport")
library("SDEtransport")
library(sl3)

# data generating process for 1-d W
f_W = function(n) rnorm(n)
f_S = function(W) plogis(W +.7)
# make a pscore model
f_A = function(S,W) plogis(-.7*S-1*W +.7)
# make a intermediate confounder model
f_Z = function(A,S,W) plogis(.1*S-1*W+2*A-.8)
# make an M model according to the restrictions
f_M = function(Z,W,S) plogis(-.14*S + 1*W + 2*Z +.6)
# make a Y model according to the restrictions
f_Y = function(M,Z,W) plogis(1*M + 1.5*W + 1*Z - 1)

# generate n random samples
n = 1e3
set.seed(1)
data = gendata.SDEtransport(n, f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)

# define learners in sl3
lglm = make_learner(Lrnr_glm)
lmean = make_learner(Lrnr_mean)
lxgboost = make_learner(Lrnr_xgboost)
lrnr_stack = make_learner(Stack, list(lglm, lmean, lxgboost))
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



# run the tmle
# a is the intervention, a_star is for the stochastic intervention.  For now, the stochastic
# intervention is defined on S = 1 only for both M and Z but can add options easily for that.
# There is a truth option for testing purposes to compute the true data adaptive parameter

a = 1
a_star = 1
res = SDE_tmle(data = data, a = a, a_star = a_star, sl = sl, covariates = covariates,
               truth = list(f_S = f_S, f_Z = f_Z, f_Y = f_Y))

# tmle est
res$CI
# mle gcomp
res$est_mle
# IC mean for tmle
mean(res$IC)
# sl coefficients for the learners
res$SL_coef
# true parameter value
res$Psi_0

