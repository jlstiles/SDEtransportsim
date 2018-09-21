# sl = make_learner(Lrnr_glm, family = binomial())
# n = 1e3
# truth = func_list
# data = gendata.SDEtransport(n, 
#                             f_W = truth$f_W, 
#                             f_S = truth$f_S, 
#                             f_A = truth$f_A, 
#                             f_Z = truth$f_Z, 
#                             f_M = truth$f_M, 
#                             f_Y = truth$f_Y)
# head(data)
# 
# get_gstarM_glm()
# debug(get_gstarM)
# gstar_info = get_gstarM(data = data, sl = sl, V=10, covariates = covariates, truth = func_list) 
