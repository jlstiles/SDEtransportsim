# devtools::install_github("jeremyrcoyle/sl3")
# devtools::install_github("jlstiles/SDE_transport")
library("SDEtransportsim")
library("parallel")
library("doParallel")

func_formsYM$func_listYMmis
func_formsYZ$func_listYZmis

forms = func_formsYM$formswell
forms = lapply(forms, formula)
# forms$Zform = formula("Z~A*W2")
# forms$Mform = formula("M~Z*W2")
truth = func_formsYM$func_listYMmis
truth

n=1e6
W = truth$f_W(n)
p_S = truth$f_S(W)
hist(p_S, breaks = 200)
S = rbinom(n, 1, p_S)
p_A = truth$f_A(S,W)
hist(p_A, breaks = 200)
A = rbinom(n, 1, p_A)
p_Z = truth$f_Z(A,S,W)
hist(p_Z, breaks = 200)
Z = rbinom(n, 1, p_Z)
p_M = truth$f_M(Z, W, S)
hist(p_M, breaks = 200)
M = rbinom(n, 1, p_M)
p_Y = truth$f_Y(M,Z,W)
hist(p_Y, breaks = 200)
Y = rbinom(n, 1, p_Y)
mean(Y)

Wnames = c("W2")
Wnamesalways = c("W2")

data = gendata.SDEtransport(1e4, truth$f_W, truth$f_S, truth$f_A, 
                   truth$f_Z, truth$f_M, truth$f_Y)
head(data)
data$weights = rep(1,nrow(data))
truth1 = truth
truth1 = NULL
res1 = try(suppressWarnings(SDE_tmle_lasso(data, forms, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                                           transport = TRUE, pooled = TRUE, gstar_S = 0, truth = truth1))) 


res2 = try(suppressWarnings(SDE_tmle_lasso_eff(data, forms, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                                           transport = TRUE, pooled = TRUE, gstar_S = 0, truth = truth1))) 

res3 = SDE_glm4(data, truth = truth1, truncate = list(lower =.0001, upper = .9999), 
                           B = NULL, forms, RCT = 0.5) 

res4 = SDE_glm_eff(data, truth = truth1, truncate = list(lower =.0001, upper = .9999), 
                B = NULL, forms, RCT = 0.5) 

result_lasso = rbind(c(res1$CI_SDE, res1$SDE0),c(res1$CI_SIE, res1$SIE0))
result_lasso_eff = rbind(c(res2$CI_SDE, res2$SDE0),c(res2$CI_SIE, res2$SIE0))
result_glm = rbind(c(res3$CI_SDE, res3$SDE_0),c(res3$CI_SIE, res3$SIE_0))
result_glm_eff = rbind(c(res4$CI_SDE, res4$SDE_0),c(res4$CI_SIE, res4$SIE_0))

res_ex = rbind(result_lasso,
result_lasso_eff,
result_glm,
result_glm_eff)
res_ex
res_ex[,3]-res_ex[,2]

debug(get_gstarM_lasso)
debug(SDE_tmle_lasso_eff)
debug(get.mediation.initdata_lasso_eff)
debug(mediation.step1_lasso_eff)
debug(mediation.step2_lasso)
debug(mediation.step1_glm)
