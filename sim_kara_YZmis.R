library(SDEtransport)

for (h in 1:2) {
  if (h == 1) {
    load("func_lists9.RData")
    func_list = func_forms9$func_listYSmis
    forms = func_forms9$formsYZmis
    forms$Zstarform = "Z~A"
    nn=10
  } else {
    load("func_lists9.RData")
    func_list = func_forms9$func_listYMmis
    forms = func_forms9$formsYZmis
    forms$Zstarform = "Z~A"
    nn=11
  }
  
  sim_kara = function(n, forms, truth, B = NULL) {
    
    data = gendata.SDEtransport(n, 
                                f_W = truth$f_W, 
                                f_S = truth$f_S, 
                                f_A = truth$f_A, 
                                f_Z = truth$f_Z, 
                                f_M = truth$f_M, 
                                f_Y = truth$f_Y)
    SDE_glm4(data, truth = truth,
             truncate = list(lower =.0001, upper = .9999),
             B=B, forms = forms, RCT = 0.5)
  }
  
  library(parallel)
  
  B = 1000
  n=100
  
  res100_YZmis = mclapply(1:B, FUN = function(x) sim_kara(n=100, forms=forms, truth=func_list, B = NULL), 
                          mc.cores = getOption("mc.cores", 20L))
  
  save(res100_YZmis, func_list, forms, file = paste0("results", nn, "/res100_YZmis.RData"))
  
  B = 1000
  n=500
  
  res500_YZmis = mclapply(1:B, FUN = function(x) sim_kara(n=500, forms=forms, truth=func_list, B = NULL), 
                          mc.cores = getOption("mc.cores", 20L))
  
  save(res500_YZmis, func_list, forms, file = paste0("results", nn, "/res500_YZmis.RData"))
  
  B = 100
  n=5000
  
  res5000_YZmis = mclapply(1:B, FUN = function(x) sim_kara(n=5000, forms=forms, truth=func_list, B = NULL), 
                           mc.cores = getOption("mc.cores", 20L))
  
  
  save(res5000_YZmis, func_list, forms, file = paste0("results", nn, "/res5000_YZmis.RData"))
}
