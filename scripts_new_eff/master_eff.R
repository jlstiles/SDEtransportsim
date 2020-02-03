
boots = 500
DGP = 1
load("func_forms1.RData")

if (DGP == 1) {
  func_list = func_formsYM$func_listYMmis
  forms = func_formsYM
}

if (DGP == 2) {
  func_list = func_formsYZ$func_listYZmis
  forms = func_formsYZ
}

if (DGP == 3) {
  func_list = func_formsYS$func_listYSmis
  forms = func_formsYS
}

form_name = paste0("forms", type)
forms = forms[[form_name]]

sim_kara = function(n, forms, truth, B = 500) {
  
  data = gendata.SDEtransport(n, 
                              f_W = truth$f_W, 
                              f_S = truth$f_S, 
                              f_A = truth$f_A, 
                              f_Z = truth$f_Z, 
                              f_M = truth$f_M, 
                              f_Y = truth$f_Y)
  res = SDE_glm_seq(data, forms, RCT = 0.5, transport = TRUE, 
                                pooled = FALSE, gstar_S = 1, truth, B = 500) 
  
  res_eff = SDE_glm_eff_seq(data, forms, RCT = 0.5, transport = TRUE, 
                            pooled = FALSE, gstar_S = 1, truth, B = 500) 
  
  return(list(res= res, res_eff = res_eff))
}
  

  B = 1000
  n=100

  res100 = mclapply(1:B, FUN = function(x) sim_kara(n=100, forms=forms, truth=func_list, B = boots),
                         mc.cores = getOption("mc.cores", cores))


  B = 1000
  n=500

  res500 = mclapply(1:B, FUN = function(x) sim_kara(n=500, forms=forms, truth=func_list, B = boots),
                         mc.cores = getOption("mc.cores", cores))

  B = 1000
  n=5000
  
  res5000 = mclapply(1:B, FUN = function(x) sim_kara(n=5000, forms=forms, truth=func_list, B = boots), 
                          mc.cores = getOption("mc.cores", cores))
  
  

