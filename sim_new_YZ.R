library(SDEtransportsim)

boots = 500
type = "YZmis"
load("func_forms1.RData")

if (type == "YMmis") {
  func_list = func_formsYM$func_listYMmis
  forms = func_formsYM$formsYMmis
}

if (type == "YSmis") {
  func_list = func_formsYS$func_listYSmis
  forms = func_formsYS$formsYMmis
}

if (type == "YZmis") {
  func_list = func_formsYZ$func_listYZmis
  forms = func_formsYZ$formsYMmis
}

# system(paste0("mkdir -p ", paste0("results", type)))

sim_kara = function(n, forms, truth, B = 500) {
  
  data = gendata.SDEtransport(n, 
                              f_W = truth$f_W, 
                              f_S = truth$f_S, 
                              f_A = truth$f_A, 
                              f_Z = truth$f_Z, 
                              f_M = truth$f_M, 
                              f_Y = truth$f_Y)
  res = SDE_glm(data, truth = truth,
           truncate = list(lower =.0001, upper = .9999),
           B=B, forms = forms, RCT = 0.5)
  res_eff = SDE_glm_eff(data, truth = truth,
                        truncate = list(lower =.0001, upper = .9999),
                        B=B, forms = forms, RCT = 0.5)
  return(list(res = res, res_eff = res_eff))
}

library(parallel)

for (simmie in 1:6) {
  if (simmie == 1) {
    suffix = "well"
    func_list = func_formsYZ$func_listYZmis
    forms = func_formsYZ$formswell
  }
  if (simmie == 2)  {
    suffix = "YAwell"
    func_list = func_formsYZ$func_listYZmis
    forms = func_formsYZ$formsYAwell
  }
  if (simmie == 3)  {
    suffix = "Ymis"
    func_list = func_formsYZ$func_listYZmis
    forms = func_formsYZ$formsYmis
  }
  if (simmie == 4)  {
    suffix = "YMmis"
    func_list = func_formsYZ$func_listYZmis
    forms = func_formsYZ$formsYMmis
  }
  if (simmie == 5)  {
    suffix = "YSmis"
    func_list = func_formsYZ$func_listYZmis
    forms = func_formsYZ$formsYSmis
  }
  if (simmie == 6)  {
    suffix = "YZmis"
    func_list = func_formsYZ$func_listYZmis
    forms = func_formsYZ$formsYZmis
  }
  
path = paste0("results_eff1",type)
system(paste0("mkdir -p ", path))
  
  B = 1000
  n=100
  
  res100 = mclapply(1:B, FUN = function(x) sim_kara(n=100, forms=forms, truth=func_list, B = boots),
                          mc.cores = getOption("mc.cores", 24L))
  
  save(res100, func_list, forms, file = paste0(path,"/res100_", suffix ,".RData"))
  
  B = 1000
  n=500
  
  res500 = mclapply(1:B, FUN = function(x) sim_kara(n=500, forms=forms, truth=func_list, B = boots),
                          mc.cores = getOption("mc.cores", 24L))
  
  save(res500, func_list, forms, file = paste0(path, "/res500_", suffix ,".RData"))
  
  rm("res100", "res500")
  
  B = 1000
  n=5000
  
  res5000 = mclapply(1:B, FUN = function(x) sim_kara(n=5000, forms=forms, truth=func_list, B = boots), 
                           mc.cores = getOption("mc.cores", 24L))
  
  save(res5000, func_list, forms, file = paste0(path,"/res5000_", suffix ,".RData"))
  rm("res5000")
}

