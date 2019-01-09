# devtools::install_github("jeremyrcoyle/sl3")
# devtools::install_github("jlstiles/SDE_transport")
library("SDEtransportsim")
library("parallel")
library("doParallel")

#Set up data generating process:
data(example_dgp)
data = data_example

# forms for non-transporting
formsNT=list(Aform = formula("A ~ W2 + W1"),
             Zstarform = formula("Z ~  A + W1 + W2"),
             QZform = formula("Qstar_Mg ~ Z + W2 + W1"),
             Mstarform = formula("M ~ Z + W1 + W2"),
             Yform = formula("Y ~ M + Z + W1 + W2"))

# forms for transporting without a pooled M
formsNP=list(Sform = formula("S ~  W1 + W2"),
           Aform = formula("A ~ W2 + W1 + S"),
           Zstarform = formula("Z ~  A + W1 + W2 + S"),
           QZform = formula("Qstar_Mg ~ Z + W2 + W1 + S"),
           Mstarform = formula("M ~ Z + W1 + W2"),
           Yform = formula("Y ~ M + Z + W1 + W2"))

# forms for transporting with a pooled M
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

Wnames = c("W1", "W2")
Wnamesalways = c("W1")
# run the tmle, no bootstrapping inference (bad for lasso anyway)
data$weights = rep(1,nrow(data))
data$weights = runif(nrow(data))

testresNT = SDE_tmle_lasso(data, formsNT, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                         B = NULL, transport = FALSE, pooledM = TRUE, gstar_S = c(1,1)) 

testres = SDE_tmle_lasso(data, forms, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                         B = NULL, transport = TRUE, pooledM = TRUE, gstar_S = c(1,1)) 

testresNP = SDE_tmle_lasso(data, formsNP, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                         B = NULL, transport = TRUE, pooledM = FALSE, gstar_S = c(1,1)) 

testresNT$CI_SDE
testres$CI_SDE
testresNP$CI_SDE

testresNT$CI_SIE
testres$CI_SIE
testresNP$CI_SIE

# debug(SDE_tmle_lasso)
# debug(SDE_tmle_lassoNT)
# debug(get_gstarM_lasso)
# undebug(mediation.step2_lasso)
