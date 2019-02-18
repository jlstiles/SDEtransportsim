# devtools::install_github("jeremyrcoyle/sl3")
# devtools::install_github("jlstiles/SDE_transport")
library("SDEtransportsim")
library("parallel")
library("doParallel")

#Set up data generating process:
# data(example_dgp)
# data = data_example

truth = list(
  f_W = function(n) {
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


check_sim = function(n, truth, truth_NT, forms, pooled, gstar_S) {
  n=1e3
  truth_NT = truth_NT
  pooled = TRUE
  gstar_S = 1
  data = gendata.SDEtransport(n, truth$f_W, truth$f_S, truth$f_A, truth$f_Z, truth$f_M, truth$f_Y)
  
  data$weights = rep(1,nrow(data))
  
  if (is.null(truth_NT)) truth = NULL
  test = SDE_tmle_lasso(data, forms, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                        transport = TRUE, pooled = pooled, gstar_S = gstar_S, truth = truth) 
  
  test1 = SDE_tmle_lasso_eff(data, forms, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                        transport = TRUE, pooled = pooled, gstar_S = gstar_S, truth = truth) 
  
  dataG = data
  dataG$weights = NULL
  res3 = SDE_glm4(dataG, truth = truth, truncate = list(lower =.0001, upper = .9999), 
                  B = NULL, formsNP, RCT = 0.5) 
  
  res4 = SDE_glm_eff(dataG, truth = truth, truncate = list(lower =.0001, upper = .9999), 
                     B = NULL, formsNP, RCT = 0.5) 
  res5 = SDE_glm_eff1(dataG, truth = truth, truncate = list(lower =.0001, upper = .9999), 
                     B = NULL, formsNP, RCT = 0.5) 
  
  dataNT = data[data$S==1,]
  dataNT$S = NULL
  
  testNT = SDE_tmle_lasso(dataNT, formsNT, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                          transport = FALSE, pooled = FALSE, gstar_S = 1, truth = truth_NT) 
  
  testNT1 = SDE_tmle_lasso_eff(dataNT, formsNT, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                          transport = FALSE, pooled = FALSE, gstar_S = 1, truth = truth_NT) 
  
  res = rbind(c(test$CI_SDE, test$SDE0), c(test1$CI_SDE, test1$SDE0),
              c(test$CI_SIE, test$SIE0), c(test1$CI_SIE, test1$SIE0),
              c(testNT$CI_SDE, testNT$SDE0), c(testNT1$CI_SDE,testNT1$SDE0),
              c(testNT$CI_SIE, testNT$SIE0), c(testNT1$CI_SIE, testNT1$SIE0),
              c(res3$CI_SDE, res3$SDE0), c(res4$CI_SDE, res4$SDE0),
              c(res3$CI_SIE, res3$SIE0), c(res4$CI_SIE, res4$SIE0),
              c(res5$CI_SIE, res5$SIE0))
  colnames(res) = c("est", "left", "right","truth")
  rownames(res) = c("SDE_lasso_trans", "SDE_lasso_trans_eff", "SIE_lasso_trans", "SIE_lasso_trans_eff",
                    "SDE_lasso", "SDE_lasso_eff", "SIE_lasso", "SIE_lasso_eff",
                    "SDE_glm_trans", "SDE_glm_trans_eff", "SIE_glm_trans", "SIE_glm_trans_eff")
  
  return(res)
}


res3$SDE_0
res3$SIE_0

res3$SE_SDE_0
res3$SE_SIE_0
res3$SE_SDE_eff0
res3$SE_SIE_eff0

check_truth = check_sim(1e5, truth = truth, truth_NT=truth_NT, forms = formsNP, pooled = FALSE, gstar_S = 0)
check_truth
check_truth[,3] - check_truth[,2]

check_truth = check_sim(1e5, truth, truth_NT=truth_NT, forms = forms, pooled = TRUE, gstar_S = 0)
check_truth
check_truth[,3] - check_truth[,2]

check_truth = check_sim(1e5, truth, truth_NT=truth_NT, forms = formsNP, pooled = FALSE, gstar_S = 1)
check_truth
check_truth[,3] - check_truth[,2]

check_truth = check_sim(1e5, truth, truth_NT=truth_NT, forms = formsNP, pooled = TRUE, gstar_S = 1)
check_truth
check_truth[,3] - check_truth[,2]

