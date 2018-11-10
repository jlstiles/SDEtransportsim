# devtools::install_github("jeremyrcoyle/sl3")
devtools::install_github("jlstiles/SDE_transport")
library("SDEtransport")
library(sl3)

#Set up data generating process:
# data generating process for 1-d W
f_W = function(n) {
  W1 = rnorm(n)
  W2 = rbinom(n, 1, 0.72)
  data.frame(W1=W1, W2=W2)
}

# make a pscore model
f_A = function(W) {
  with(W, plogis(-.7*W1 + 0.1*W2 +.17))
}
# make a intermediate confounder model
f_Z = function(A,W) {
  df = as.data.frame(cbind(W, A = A))
  with(df, plogis(.4*W1 - W2 + 1*A-.3))
}
# make an M model according to the restrictions
f_M = function(Z,W) {
  df = as.data.frame(cbind(Z=Z, W))
  with(df, plogis(1*W1 - W2 + 1.2*Z +.1))
}
# make a Y model according to the restrictions, main terms linear logistic reg.
# plug-in is biased and not robust like tmle
f_Y = function(M,Z,W) {
  df = as.data.frame(cbind(M = M, Z = Z, W))
  with(df, plogis(W2*M + 3*cos(W1)*Z-.4))
}

# generate n random samples
n = 1e3
set.seed(1)
data = gendata.SDEtransport_alt(n, f_W = f_W, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
W = data[,grep("W", colnames(data))]

data = list(W=W, A=data$A, Z=data$Z, M=data$M, Y=data$Y)

# define the forms for all the regressions
forms=list(Aform = formula("A ~ W2 + W1"),
           Zstarform = formula("Z ~  A + W1"),
           QZform = formula("Qstar_Mg ~ Z + W2 + W1"),
           Mstarform = formula("M ~ Z + W1 +W2"),
           Yform = formula("Y ~ M + Z + W1 + W2"))

# run the tmle, no bootstrapping inference (bad for lasso anyway)
testres = SDE_tmle_lasso(data=data, truth = NULL, truncate = list(lower =.0001, upper = .9999), 
                        B = NULL, forms, RCT = NULL) 


testres
