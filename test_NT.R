# devtools::install_github("jeremyrcoyle/sl3")
# devtools::install_github("jlstiles/SDE_transport")
library("SDEtransportsim")
library("parallel")
library("doParallel")

#Set up data generating process:
# data(example_dgp)
# data = data_example

truth_NT = list(
  f_W = function(n) {
    W1 = rnorm(n)
    W2 = rnorm(n)
    data.frame(W1=W1, W2=W2)
  }
  
  ,f_A = function(W) {
    if (is.vector(W)) nn = length(W) else nn = nrow(W)
    rep(.5, nn)
  }
  
  ,f_Z = function(A,W) {
    df = as.data.frame(cbind(W, A = A))
    with(df, plogis(.4*W1 -.6*W2 + 1*A + .3))
  }
  
  ,f_M = function(Z,W) {
    df = as.data.frame(cbind(Z=Z, W))
    with(df, plogis(.41*W1 +.16*W2 - 3*Z+ 1.6))
  }
  
  ,f_Y = function(M,Z,W) {
    df = as.data.frame(cbind(M = M, Z = Z, W))
    with(df, plogis(W2 + W1 - Z + 2*M - .4))
  }
)

# be sure of no bad positivity
n=100000
W = truth_NT$f_W(n)
A = rbinom(n, 1, .5)
Z = rbinom(n,1,truth_NT$f_Z(A,W))
hist(truth_NT$f_M(Z,W), breaks = 200)
M = rbinom(n,1,truth_NT$f_M(Z,W))
# forms for non-transporting
formsNT=list(Aform = formula("A ~ W1 + W2"),
             Zstarform = formula("Z ~  A + W1 + W2"),
             QZform = formula("Qstar_Mg ~ W1 + W2"),
             Mstarform = formula("M ~ Z + W1 + W2"),
             Yform = formula("Y ~ M + Z + W1 + W2"))

Wnames = c("W1", "W2")
Wnamesalways = c("W1")

res = lapply(1:334, FUN = function(x){
  p = sample(1:1e6, 1)
  set.seed(p)
  data = gendata.SDE(1000, truth_NT$f_W, truth_NT$f_A, 
                     truth_NT$f_Z, truth_NT$f_M, truth_NT$f_Y)
  data$weights = rep(1,nrow(data))
  res1 = try(suppressWarnings(SDE_tmle_lasso(data, formsNT, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                                             transport = FALSE, pooled = FALSE, gstar_S = 0, truth = truth_NT))) 
  if (class(res1)[1]=="try-error") {
    res = rep(0,8)
  } else {
    res = c(res1$CI_SDE, res1$SDE0, res1$CI_SIE, res1$SIE0)
  }
  return(list(res, p, data))
})

save(res, file = "res2.RData")



