# uncomment these if testing the lasso functions without loading the package
# library("parallel")
# library("doParallel")
# library("glmnet")

#Set up data generating process:
data(example_dgp)
data = data_example

# forms for non-transporting contains no S variable
formsNT=list(Aform = formula("A ~ W2 + W1"),
             Zstarform = formula("Z ~  A + W1 + W2"),
             QZform = formula("Qstar_Mg ~ Z + W2 + W1"),
             Mstarform = formula("M ~ Z + W1 + W2"),
             Yform = formula("Y ~ M + Z + W1 + W2"))

# forms for transporting without a pooled (over S) M or Z fit
formsNP=list(Sform = formula("S ~  W1 + W2"),
           Aform = formula("A ~ W2 + W1 + S"),
           Zstarform = formula("Z ~  A + W1 + W2"),
           QZform = formula("Qstar_Mg ~ Z + W2 + W1+S"),
           Mstarform = formula("M ~ Z + W1 + W2"),
           Yform = formula("Y ~ M + Z + W1 + W2"))

# forms for transporting with a pooled (over S) M and Z fit
forms=list(Sform = formula("S ~  W1 + W2"),
           Aform = formula("A ~ W2 + W1 + S"),
           Zstarform = formula("Z ~  A + W1 + W2 + S"),
           QZform = formula("Qstar_Mg ~ Z + W2 + W1 + S"),
           Mstarform = formula("M ~ Z + W1 + W2 + S"),
           Yform = formula("Y ~ M + Z + W1 + W2"))


W = data[,grep("W", colnames(data))]
head(data)
data = cbind(W, A=data$A, S=data$S, Z=data$Z, M=data$M, Y=data$Y)

# define the forms for all the regressions
n=1000

Wnames = c("W1", "W2")
Wnamesalways = c("W1")

# choose no weights (all weights are 1)
data$weights = rep(1,nrow(data))
# choose randomly generated weights
data$weights = runif(nrow(data))

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

# undebug(SDE_tmle_lasso)
# undebug(get_gstarM_lasso)
# undebug(get.mediation.initdata_lasso)
# undebug(mediation.step1_lasso)
# undebug(mediation.step2_lasso)
