# devtools::install_github("jeremyrcoyle/sl3")
# devtools::install_github("jlstiles/SDE_transport")
library("SDEtransport")
library(sl3)

#Set up data generating process:
# data generating process for 1-d W
f_W = function(n) {
  W1 = rnorm(n)
  gender = rbinom(n, 1, 0.72)
  data.frame(W1=W1, gender=gender)
}
f_S = function(W) {
with(W, plogis(W1 - gender + 0.3))
}
# make a pscore model
f_A = function(S,W) {
  df = as.data.frame(cbind(S=S, W))
  with(df, plogis(-.6*S-.7*W1 + 0.1*gender +.17))
}
# make a intermediate confounder model
f_Z = function(A,S,W) {
  df = as.data.frame(cbind(S=S, W, A = A))
  with(df, plogis(.1*S-.4*W1 - gender + 1*A-.3))
}
# make an M model according to the restrictions
f_M = function(Z,W,S) {
  df = as.data.frame(cbind(Z=Z, W, S=S))
  with(df, plogis(-.14*S + 1*W1 - gender + 1.2*Z +.1))
}
# make a Y model according to the restrictions, main terms linear logistic reg.
# plug-in is biased and not robust like tmle
f_Y = function(M,Z,W) {
  df = as.data.frame(cbind(M = M, Z = Z, W))
  with(df, plogis(gender*M + 3*cos(W1)*Z-.4))
}

# generate n random samples
n = 1e2
set.seed(1)
data = gendata.SDEtransport(n, f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
mean(data$S==0)

# define the forms for all the regressions
forms=list(Sform = formula("S ~ W1"),
                Aform = formula("A ~ S + W1"),
                Zstarform = formula("Z ~ S + A + W1"),
                QZform = formula("Qstar_Mg ~ Z + S + W1"),
                Mstarform = formula("M ~ Z + W1"),
                Yform = formula("Y ~ M + Z + W1"))


W = data[,grep("W", colnames(data))]

data = list(W=W, S=data$S, A=data$A, Z=data$Z, M=data$M, Y=data$Y)
# debug(SDE_tmle_glm)

# run the tmle, specifying 2 bootstraps for bootstrap inference but obviously 2 is way too few.
testres = SDE_tmle_glm(data=data, truth = NULL, truncate = list(lower =.0001, upper = .9999), 
                        B = 2, forms, RCT = NULL) 

testres

# run the tmle, with a known data generating process, in case of simulation
truth = list(f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)

testres1 = SDE_tmle_glm(data, truth = truth, truncate = list(lower =.0001, upper = .9999), 
                       B = 2, forms, RCT = NULL) 

testres1

