library(SDEtransportsim)

boots = 500
type = "YMmis"
load("func_forms1.RData")

if (type == "YMmis") {
  func_list = func_formsYM$func_listYMmis
  forms = func_formsYM$formswell
}

if (type == "YSmis") {
  func_list = func_formsYS$func_listYSmis
  forms = func_formsYS$formswell
}

if (type == "YZmis") {
  func_list = func_formsYZ$func_listYZmis
  forms = func_formsYZ$formswell
}

system(paste0("mkdir -p ", paste0("results", type)))
sim_kara = function(n, forms, truth, B = boots) {
  
  data = gendata.SDEtransport(n, 
                              f_W = truth$f_W, 
                              f_S = truth$f_S, 
                              f_A = truth$f_A, 
                              f_Z = truth$f_Z, 
                              f_M = truth$f_M, 
                              f_Y = truth$f_Y)
  res = SDE_glm_seq(data, forms, RCT = 0.5, transport=T, 
                    pooled=T, gstar_S = 0, truth=truth, B = boots) 
  
  res_eff = SDE_glm_eff_seq(data, forms, RCT = 0.5, transport=T,
                            pooled=T, gstar_S = 0, truth=truth, B = boots)
  return(list(res= res, res_eff = res_eff))
}

  library(parallel)

  B = 1000
  n=100

  res100_well = mclapply(1:B, FUN = function(x) sim_kara(n=100, forms=forms, truth=func_list, B = boots),
                         mc.cores = getOption("mc.cores", 20L))

  save(res100_well, func_list, forms, file = paste0("results", type, "/res100_well.RData"))

  B = 1000
  n=500

  res500_well = mclapply(1:B, FUN = function(x) sim_kara(n=500, forms=forms, truth=func_list, B = boots),
                         mc.cores = getOption("mc.cores", 20L))

  save(res500_well, func_list, forms, file =  paste0("results", type, "/res500_well.RData"))

  rm("res100_well", "res500_well")
  
  B = 1000
  n=5000
  
  res5000_well = mclapply(1:B, FUN = function(x) sim_kara(n=5000, forms=forms, truth=func_list, B = boots), 
                          mc.cores = getOption("mc.cores", 20L))
  
  save(res5000_well, func_list, forms, file = paste0("results", type, "/res5000_well.RData"))
  

