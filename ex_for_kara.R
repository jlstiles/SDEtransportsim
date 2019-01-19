# devtools::install_github("jeremyrcoyle/sl3")
# devtools::install_github("jlstiles/SDE_transport")
library("SDEtransportsim")
library("parallel")
library("doParallel")

#Set up data generating process:
# data(example_dgp)
# data = data_example

truth = list(f_W = function(n) {
  W1 = rnorm(n)
  W2 = rnorm(n)
  data.frame(W1=W1, W2=W2)
}
,f_S = function(W) {
  with(W, plogis(-.37*W1 + W2-.17))
}

,f_A = function(S,W) {
  rep(.5, length(S))
}

,f_Z = function(A,S,W) {
  df = as.data.frame(cbind(W, A = A))
  with(df, plogis(.3*S*(.4*W1 - W2 + A) +.4*W1 -.6*W2 + 1*A -.5*S + .3))
}

,f_M = function(Z,W,S) {
  df = as.data.frame(cbind(Z=Z, W, S))
  with(df, plogis(.4*S*(1*W1 +.6*W2 - 1*Z) + .41*W1 +.16*W2 - 3*Z + .51*S + 1.6))
}

,f_Y = function(M,Z,W) {
  df = as.data.frame(cbind(M = M, Z = Z, W))
  with(df, plogis(W2 + W1 - Z + 2*M - .4))
}
)

truth_NT = list(f_W = function(n) {
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
W = truth$f_W(n)
hist(truth$f_S(W))
S = rbinom(n,1,truth$f_S(W))
A = rbinom(n, 1, .5)
hist(truth$f_Z(A,S,W), breaks = 200)
Z = rbinom(n,1,truth$f_Z(A,S,W))
hist(truth$f_M(Z,W,S), breaks = 200)
M = rbinom(n,1,truth$f_M(Z,W,S))
# forms for non-transporting
formsNT=list(Aform = formula("A ~ W1 + W2"),
             Zstarform = formula("Z ~  A + W1 + W2"),
             QZform = formula("Qstar_Mg ~ W1 + W2"),
             Mstarform = formula("M ~ Z + W1 + W2"),
             Yform = formula("Y ~ M + Z + W1 + W2"))

# forms for transporting without a pooled M
formsNP=list(Sform = formula("S ~  W1 + W2"),
           Aform = formula("A ~ W1 + W2 + S"),
           Zform = formula("Z ~  S:(A + W1 + W2) + A + W1 + W2 + S"),
           Zstarform = formula("Z ~  A + W1 + W2"),
           QZform = formula("Qstar_Mg ~ W1 + W2 + S"),
           Mform = formula("M ~ S:(Z + W1 + W2) + Z + W1 + W2 + S"),
           Mstarform = formula("M ~ Z + W1 + W2"),
           Yform = formula("Y ~ M + Z + W1 + W2"))

# forms for transporting with a pooled M
forms=list(Sform = formula("S ~  W1 + W2"),
           Aform = formula("A ~ W2 + W1 + S"),
           Zstarform = formula("Z ~  S:(A + W1 + W2) + A + W1 + W2 + S"),
           QZform = formula("Qstar_Mg ~ S:(W2 + W1) + W2 + W1 + S"),
           Mstarform = formula("M ~ S:(Z + W1 + W2) + Z + W1 + W2 + S"),
           Yform = formula("Y ~ M + Z + W1 + W2"))

Wnames = c("W1", "W2")
Wnamesalways = c("W1")


check_sim = function(n) {
  data = gendata.SDEtransport(n, truth$f_W, truth$f_S, truth$f_A, truth$f_Z, truth$f_M, truth$f_Y)
  
  data$weights = rep(1,nrow(data))
  
  test = SDE_tmle_lasso(data, formsNP, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                        transport = TRUE, pooled = FALSE, gstar_S = 0, truth = truth) 
  
  
  dataNT = data[data$S==0,]
  dataNT$S = NULL
  
  testNT = SDE_tmle_lasso(dataNT, formsNT, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                          transport = FALSE, pooled = FALSE, gstar_S = 0, truth = truth_NT) 
  
  return(c(test$CI_SDE, test$SDE0, test$CI_SIE, test$SIE0, testNT$CI_SDE,
           testNT$SDE0, testNT$CI_SIE, testNT$SIE0))
}

check_truth = check_sim(100000)
# SDE for transported from S =1 to S = 0 using S = 0 mechanism for gstar
check_truth[c(1,4)]
# SIE for transported from S =1 to S = 0 using S = 0 mechanism for gstar
check_truth[c(5,8)]
# SDE for transported for S = 0 mechanism for gstar
check_truth[c(9,12)]
# SIR for transported S = 0 mechanism for gstar
check_truth[c(13,16)]
