# devtools::install_github("jeremyrcoyle/sl3")
# devtools::install_github("jlstiles/SDE_transport")
library("SDEtransportsim")
library("parallel")
library("doParallel")

#Set up data generating process:
# data(example_dgp)
# data = data_example
f_W = function(n) {
  W1 = rnorm(n)
  W2 = rnorm(n)
  data.frame(W1=W1, W2=W2)
}

f_S = function(W) {
  with(W, plogis(-.37*W1 + W2-.17))
}

# make a pscore model
f_A = function(S,W) {
  rep(.5, length(S))
}
# make a intermediate confounder model
f_Z = function(A,S,W) {
  df = as.data.frame(cbind(W, A = A))
  with(df, plogis(.3*S*(.4*W1 - W2 + A) +.4*W1 -.6*W2 + 1*A -.5*S + .3))
}
# make an M model according to the restrictions
f_M = function(Z,W,S) {
  df = as.data.frame(cbind(Z=Z, W, S))
  with(df, plogis(S*(1*W1 +.6*W2 - 1*Z) + 1*W1 +.6*W2 - 1*Z+ 1*S + .1))
}
# make a Y model according to the restrictions, main terms linear logistic reg.
# plug-in is biased and not robust like tmle
f_Y = function(M,Z,W) {
  df = as.data.frame(cbind(M = M, Z = Z, W))
  with(df, plogis(W2 + W1 - Z + M - .4))
}

n=100000
data = gendata.SDEtransport(n, f_W, f_S, f_A, f_Z, f_M, f_Y)
truth = list(f_W = f_W, f_S=f_S, f_A=f_A, f_Z=f_Z, f_M=f_M, f_Y=f_Y)

# forms for non-transporting
formsNT=list(Aform = formula("A ~ W1 + W2"),
             Zstarform = formula("Z ~  A + W1 + W2"),
             QZform = formula("Qstar_Mg ~ Z + W1 + W2"),
             Mstarform = formula("M ~ Z + W1 + W2"),
             Yform = formula("Y ~ M + Z + W1 + W2"))

# forms for transporting without a pooled M
formsNP=list(Sform = formula("S ~  W1 + W2"),
           Aform = formula("A ~ W1 + W2 + S"),
           Zform = formula("Z ~  S:(A + W1 + W2) + A + W1 + W2 + S"),
           Zstarform = formula("Z ~  A + W1 + W2"),
           QZform = formula("Qstar_Mg ~ Z + W1 + W2 + S"),
           Mform = formula("M ~ S:(Z + W1 + W2) + Z + W1 + W2 + S"),
           Mstarform = formula("M ~ Z + W1 + W2"),
           Yform = formula("Y ~ M + Z + W1 + W2"))

# forms for transporting with a pooled M
forms=list(Sform = formula("S ~  W1 + W2"),
           Aform = formula("A ~ W2 + W1 + S"),
           Zstarform = formula("Z ~  S:(A + W1 + W2) + A + W1 + W2 + S"),
           QZform = formula("Qstar_Mg ~ S:(Z + W2 + W1) + Z + W2 + W1 + S"),
           Mstarform = formula("M ~ S:(Z + W1 + W2) + Z + W1 + W2 + S"),
           Yform = formula("Y ~ M + Z + W1 + W2"))

Wnames = c("W1", "W2")
Wnamesalways = c("W1")
# run the tmle, no bootstrapping inference (bad for lasso anyway)
data$weights = rep(1,nrow(data))
head(data)

test = SDE_tmle_lasso(data, forms, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
               transport = TRUE, pooled = TRUE, gstar_S = 1, truth = truth) 

test$CI_SDE 
test$SDE0

test$CI_SIE 
test$SIE0

undebug(SDE_tmle_lasso)
undebug(get_gstarM_lasso)
debug(get.mediation.initdata_lasso)
debug(mediation.step1_lasso)
debug(mediation.step2_lasso)





i=1
res = data.frame(matrix(rep(NA,36),nrow = 6))
colnames(res) = c("SDE", "left", "right","SIE", "left", "right")
rownames(res) = c("Trans, pooled, gstar_S=0", 
                  "Trans, pooled, gstar_S=1", 
                  "Trans, not pooled, gstar_S=0", 
                  "Trans, not pooled, gstar_S=1",
                  "Not trans, gstar_S=0", 
                  "Not trans, gstar_S=1")

for (pooled in c(FALSE, TRUE)) {
  for (gstar_S in c(0,1)) {
    if (pooled) forms = forms else forms = formsNP
    ans = SDE_tmle_lasso(data, forms, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                         B = NULL, transport = TRUE, pooled = pooled, gstar_S = gstar_S) 
    res[i,]=c(ans$CI_SDE,ans$CI_SIE)
    i = i + 1
  }
}

for (gstar_S in c(0,1)) {
  ans = SDE_tmle_lasso(data, formsNT, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                       B = NULL, transport = FALSE, pooled = pooled, gstar_S = gstar_S) 
  res[i,]=c(ans$CI_SDE,ans$CI_SIE)
  i = i + 1
}

res



