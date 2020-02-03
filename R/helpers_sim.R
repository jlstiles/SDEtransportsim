
#' @export
gendata.SDEtransport = function(n, f_W, f_S, f_A, f_Z, f_M, f_Y) {
  W = f_W(n)
  W = as.data.frame(W)
  P_SW = f_S(W)
  S = rbinom(n,1,P_SW)
  
  # make a pscore model
  pscores = f_A(S=S,W=W)
  A = rbinom(n, 1, pscores)
  
  # make a intermediate confounder model
  pzscores = f_Z(A=A,S=S,W=W)
  Z = rbinom(n, 1, pzscores)
  
  # make an M model according to the restrictions
  Mscores = f_M(Z=Z,W=W,S=S)
  M = rbinom(n, 1, Mscores)
  
  # make a Y model according to the restrictions
  Yscores = f_Y(M=M,Z=Z,W=W)
  Y = rbinom(n, 1, Yscores)
  
  return(as.data.frame(cbind(W, S = S, A = A, Z = Z, M = M, Y = Y)))
}


#' @export
gendata.SDE = function(n, f_W, f_A, f_Z, f_M, f_Y) {
  W = f_W(n)
  # make a pscore model
  pscores = f_A(W=W)
  A = rbinom(n, 1, pscores)
  
  # make a intermediate confounder model
  pzscores = f_Z(A=A,W=W)
  Z = rbinom(n, 1, pzscores)
  
  # make an M model according to the restrictions
  Mscores = f_M(Z=Z,W=W)
  M = rbinom(n, 1, Mscores)
  
  # make a Y model according to the restrictions
  Yscores = f_Y(M=M,Z=Z,W=W)
  Y = rbinom(n, 1, Yscores)
  
  return(as.data.frame(cbind(W, A = A, Z = Z, M = M, Y = Y)))
}

#' @export
compile_SDE = function(res_list, func_list, forms)
{  
  # res_list = res10000_YAwell
  # res_list = res500_Ymis1
  # res_list = res100
  res = lapply(res_list, FUN = function(x) {
    if (is.numeric(x[[1]][1])) unlist(x) else return(0)
  })
  
  valid_draws = which(vapply(1:length(res), FUN = function(x) res[[x]][1] != 0, FUN.VALUE = TRUE))
  
  length(valid_draws)
  res = res[valid_draws]
  res = do.call(rbind, res)
  # colnames(res)
  res = as.data.frame(res)
  coverage = c(covSDE_tmle = mean((res[,37] >= res[,2]) & (res[,37] <= res[,3])), 
               covSDE_1s = mean((res[,37] >= res[,8]) & (res[,37] <= res[,9])),
               covSDE_iptw = mean((res[,37] >= res[,14]) & (res[,37] <= res[,15])),
               covSIE_tmle = mean((res[,38] >= res[,5]) & (res[,38] <= res[,6])),
               covSIE_1s = mean((res[,38] >= res[,11]) & (res[,38] <= res[,12])),
               covSIE_iptw = mean((res[,38] >= res[,17]) & (res[,38] <= res[,18])),
               covSDE_tmleB = mean((res[,37] >= res[,20]) & (res[,37] <= res[,21])), 
               covSDE_1sB = mean((res[,37] >= res[,26]) & (res[,37] <= res[,27])),
               covSDE_iptwB = mean((res[,37] >= res[,32]) & (res[,37] <= res[,33])),
               covSIE_tmleB = mean((res[,38] >= res[,23]) & (res[,38] <= res[,24])), 
               covSIE_1sB = mean((res[,38] >= res[,29]) & (res[,38] <= res[,30])),
               covSIE_iptwB = mean((res[,38] >= res[,35]) & (res[,38] <= res[,36])))
  
  # res[, c(38,23,24,29,30,35,36)]
  res_SEboot = cbind(
    (res[,21] - res[,20])/(2*1.96),
    (res[,24] - res[,23])/(2*1.96),
    (res[,27] - res[,26])/(2*1.96),
    (res[,30] - res[,29])/(2*1.96),
    (res[,33] - res[,32])/(2*1.96),
    (res[,36] - res[,35])/(2*1.96))
  
  res_SE_SDE = res[,c(39,41,43)]
  res_SE_SIE = res[,c(40,42,44)]
  
  res_SE_SDEboot = res_SEboot[,c(1,3,5)]
  res_SE_SIEboot = res_SEboot[,c(2,4,6)]
  
  res_SE_0 = res[,45:46]
  res_SE_eff0 = res[,47:48]
  efficiency_SDE = colMeans(apply(res_SE_SDE, 2, FUN = function(x) x/res_SE_eff0[,1]))
  efficiency_SIE = colMeans(apply(res_SE_SIE, 2, FUN = function(x) x/res_SE_eff0[,2]))
  # efficiency_SDE
  # efficiency_SIE
  
  efficiency_SDEboot = colMeans(apply(res_SE_SDEboot, 2, FUN = function(x) x/res_SE_eff0[,1]))
  efficiency_SIEboot = colMeans(apply(res_SE_SIEboot, 2, FUN = function(x) x/res_SE_eff0[,2]))
  # efficiency_SDEboot
  # efficiency_SIEboot
  
  res_ests_SDE = res[,c(1,7,13,37)]
  res_ests_SIE = res[,c(4,10,16,38)]
  
  res_bias_SDE = cbind(res_ests_SDE[,1] - res_ests_SDE[,4],
                       res_ests_SDE[,2] - res_ests_SDE[,4],
                       res_ests_SDE[,3] - res_ests_SDE[,4])
  
  res_bias_SIE = cbind(res_ests_SIE[,1] - res_ests_SIE[,4],
                       res_ests_SIE[,2] - res_ests_SIE[,4],
                       res_ests_SIE[,3] - res_ests_SIE[,4])
  
  coverage = rbind(coverage[1:6], coverage[7:12])
  rownames(coverage) = c("IC variance", "Bootstrap")
  colnames(coverage) = c("SDE_tmle", "SDE_EE", "SDE_iptw", "SIE_tmle", "SIE_EE", "SIE_iptw")
  
  # precentage SE of of efficiency bound
  efficiency_percentageTrue = rbind(100*c(efficiency_SDE, efficiency_SIE),
                                    100*c(efficiency_SDEboot, efficiency_SIEboot))
  
  colnames(efficiency_percentageTrue) = c("SDE_tmle", "SDE_EE", "SDE_iptw", "SIE_tmle", "SIE_EE", "SIE_iptw")
  rownames(efficiency_percentageTrue) = c("IC variance", "Bootstrap")
  
  
  # precentage bias of SE
  bias_perc_SDE = 100*unlist(lapply(lapply(1:3, FUN = function(x) res_bias_SDE[,x]/res_SE_SDE[,x]), mean))
  bias_perc_SIE = 100*unlist(lapply(lapply(1:3, FUN = function(x) res_bias_SIE[,x]/res_SE_SIE[,x]), mean))
  bias_percentageSE = c(bias_perc_SDE, bias_perc_SIE)
  
  bias_perc_SDE = 100*unlist(lapply(lapply(1:3, FUN = function(x) res_bias_SDE[,x]/res_SE_SDE[,x]), mean))
  bias_perc_SIE = 100*unlist(lapply(lapply(1:3, FUN = function(x) res_bias_SIE[,x]/res_SE_SIE[,x]), mean))
  bias  = colMeans(cbind(res_bias_SDE, res_bias_SIE))
  names(bias) = c("SDE_tmle", "SDE_EE", "SDE_iptw", "SIE_tmle", "SIE_EE", "SIE_iptw")
  
  names(bias_percentageSE) = c("SDE_tmle", "SDE_EE", "SDE_iptw", "SIE_tmle", "SIE_EE", "SIE_iptw")
  
  RMSE = c(SDE_tmle = sqrt(mean((res$SDE_0 - res[,1])^2)),
           sqrt(SDE_1s = mean((res$SDE_0 - res[,7])^2)),
           sqrt(SDE_iptw = mean((res$SDE_0 - res[,13])^2)),
           sqrt(SIE_tmle = mean((res$SIE_0 - res[,4])^2)),
           sqrt(SIE_1s = mean((res$SIE_0 - res[,10])^2)),
           sqrt(SIE_iptw = mean((res$SIE_0 - res[,16])^2)))
  return(list(results = res,
              coverage = coverage, 
              efficiency_percentageTrue = efficiency_percentageTrue,
              bias_percentageSE = bias_percentageSE,
              bias = bias,
              RMSE = RMSE,
              valid_draws = length(valid_draws),
              func_list = func_list,
              forms = forms))
  
}

#' @export
compile_SDE_noboot = function(res_list, func_list, forms)
{  
  res_list = res100t
  res = lapply(res_list, FUN = function(x) {
    if (is.numeric(x[[1]][[1]][1]) & is.numeric(x[[2]][[1]][1])) unlist(x) else return(0)
  })
  
  valid_draws = which(vapply(1:length(res), FUN = function(x) res[[x]][1] != 0, FUN.VALUE = TRUE))
  
  length(valid_draws)
  res = res[valid_draws]
  
  res = do.call(rbind, res)
  # colnames(res)
  res = as.data.frame(res)
  colnames(res)
  
  coverage = c(SDE_tmle = mean((res$res.SDE_0 >= res[,2]) & (res$res.SDE_0 <= res[,3])),
               SDE_1s = mean((res$res.SDE_0 >= res[,8]) & (res$res.SDE_0 <= res[,9])),
               SDE_iptw = mean((res$res.SDE_0 >= res[,14]) & (res$res.SDE_0 <= res[,15])),
               SIE_tmle = mean((res$res.SIE_0 >= res[,5]) & (res$res.SIE_0 <= res[,6])),
               SIE_1s = mean((res$res.SIE_0 >= res[,11]) & (res$res.SIE_0 <= res[,12])),
               SIE_iptw = mean((res$res.SIE_0 >= res[,17]) & (res$res.SIE_0 <= res[,18])))
  
  eff = 51
  coverage_eff = c(SDE_tmle_eff = mean((res$res.SDE_0 >= res[,2+eff]) & (res$res.SDE_0 <= res[,3+eff])),
               SDE_1s_eff = mean((res$res.SDE_0 >= res[,8+eff]) & (res$res.SDE_0 <= res[,9+eff])),
               SDE_iptw_eff = mean((res$res.SDE_0 >= res[,14+eff]) & (res$res.SDE_0 <= res[,15+eff])),
               SIE_tmle_eff = mean((res$res.SIE_0 >= res[,5+eff]) & (res$res.SIE_0 <= res[,6+eff])),
               SIE_1s_eff = mean((res$res.SIE_0 >= res[,11+eff]) & (res$res.SIE_0 <= res[,12+eff])),
               SIE_iptw_eff = mean((res$res.SIE_0 >= res[,17+eff]) & (res$res.SIE_0 <= res[,18+eff])))
  
  coverage = c(coverage, coverage_eff)
  coverage
  
  coveragedf = data.frame(SDE_tmle = ((res$res.SDE_0 >= res[,2]) & (res$res.SDE_0 <= res[,3])),
               SDE_1s = ((res$res.SDE_0 >= res[,8]) & (res$res.SDE_0 <= res[,9])),
               SDE_iptw = ((res$res.SDE_0 >= res[,14]) & (res$res.SDE_0 <= res[,15])),
               SIE_tmle = ((res$res.SIE_0 >= res[,5]) & (res$res.SIE_0 <= res[,6])),
               SIE_1s = ((res$res.SIE_0 >= res[,11]) & (res$res.SIE_0 <= res[,12])),
               SIE_iptw = ((res$res.SIE_0 >= res[,17]) & (res$res.SIE_0 <= res[,18])))
  
  SDE_fucked = which(!coveragedf$SDE_tmle & coveragedf$SDE_1s)
  colMeans(res[!1:100 %in% SDE_fucked,c(46:51)])
  colMeans(res[,c(46:51)])
  res[,48]
  colMeans(abs(res[1:100 %in% SDE_fucked,c(46:51)])>1)
  colMeans(abs(res[!1:100 %in% SDE_fucked,c(46:51)])>1)
  
  colMeans(res[SDE_fucked,c(64:69)])
  forms
  func_list
  res_SE_SDE = res[,c(21,23,25,66,68)]
  res_SE_SIE = res[,c(22,24,26,67,69)]
  
  res_SE_0 = res[,70:71]
  
  efficiency_SDE = colMeans(apply(res_SE_SDE, 2, FUN = function(x) x/res_SE_0[,1]))
  efficiency_SIE = colMeans(apply(res_SE_SIE, 2, FUN = function(x) x/res_SE_0[,2]))
  
  
  res_ests_SDE = res[,c(1,7,13,19)]
  res_ests_SIE = res[,c(4,10,16,20)]
  
  res_bias_SDE = cbind(res_ests_SDE[,1] - res_ests_SDE[,4],
                       res_ests_SDE[,2] - res_ests_SDE[,4],
                       res_ests_SDE[,3] - res_ests_SDE[,4])
  
  res_bias_SIE = cbind(res_ests_SIE[,1] - res_ests_SIE[,4],
                       res_ests_SIE[,2] - res_ests_SIE[,4],
                       res_ests_SIE[,3] - res_ests_SIE[,4])
  
  coverage = rbind(coverage[1:6], coverage[7:12])
  rownames(coverage) = c("IC variance", "Bootstrap")
  colnames(coverage) = c("SDE_tmle", "SDE_EE", "SDE_iptw", "SIE_tmle", "SIE_EE", "SIE_iptw")
  
  # precentage SE of of efficiency bound
  efficiency_percentageTrue = 100*c(efficiency_SDE, efficiency_SIE)
  names(efficiency_percentageTrue) = c("SDE_tmle", "SDE_EE", "SDE_iptw", "SIE_tmle", "SIE_EE", "SIE_iptw")
  
  
  # precentage bias of SE
  bias_perc_SDE = 100*unlist(lapply(lapply(1:3, FUN = function(x) res_bias_SDE[,x]/res_SE_SDE[,x]), mean))
  bias_perc_SIE = 100*unlist(lapply(lapply(1:3, FUN = function(x) res_bias_SIE[,x]/res_SE_SIE[,x]), mean))
  bias_percentageSE = c(bias_perc_SDE, bias_perc_SIE)
  
  bias_perc_SDE = 100*unlist(lapply(lapply(1:3, FUN = function(x) res_bias_SDE[,x]/res_SE_SDE[,x]), mean))
  bias_perc_SIE = 100*unlist(lapply(lapply(1:3, FUN = function(x) res_bias_SIE[,x]/res_SE_SIE[,x]), mean))
  bias  = colMeans(cbind(res_bias_SDE, res_bias_SIE))
  names(bias) = c("SDE_tmle", "SDE_EE", "SDE_iptw", "SIE_tmle", "SIE_EE", "SIE_iptw")
  
  names(bias_percentageSE) = c("SDE_tmle", "SDE_EE", "SDE_iptw", "SIE_tmle", "SIE_EE", "SIE_iptw")
  
  RMSE = c(SDE_tmle = sqrt(mean((res$SDE_0 - res[,1])^2)),
           sqrt(SDE_1s = mean((res$SDE_0 - res[,7])^2)),
           sqrt(SDE_iptw = mean((res$SDE_0 - res[,13])^2)),
           sqrt(SIE_tmle = mean((res$SIE_0 - res[,4])^2)),
           sqrt(SIE_1s = mean((res$SIE_0 - res[,10])^2)),
           sqrt(SIE_iptw = mean((res$SIE_0 - res[,16])^2)))
  
  return(list(results = res,
              coverage = coverage, 
              efficiency_percentageTrue = efficiency_percentageTrue,
              bias_percentageSE = bias_percentageSE,
              bias = bias,
              RMSE = RMSE,
              valid_draws = length(valid_draws),
              func_list = func_list,
              forms = forms))
  
}

#' @export
compile_SDE_new = function(res_list, func_list, forms, bootcut)
{  
  # res_list = res10000_YAwell
  # res_list = res500_Ymis1
  # res_list = res100
  res = lapply(res_list, FUN = function(x) {
    if (is.numeric(x[[1]][[1]][1]) & is.numeric(x[[2]][[1]][1])) unlist(x) else return(0)
  })
  
  valid_draws = which(vapply(1:length(res), FUN = function(x) res[[x]][1] != 0, FUN.VALUE = TRUE))
  
  length(valid_draws)
  res = res[valid_draws]
  
  res = do.call(rbind, res)
  # colnames(res)
  res = as.data.frame(res)
  colnames(res)
  coverage = c(covSDE_tmle_eff = mean((res[,37] >= res[,77]) & (res[,37] <= res[,78])), 
               covSDE_1s_eff = mean((res[,37] >= res[,83]) & (res[,37] <= res[,84])),
               covSDE_tmle = mean((res[,37] >= res[,2]) & (res[,37] <= res[,3])), 
               covSDE_1s = mean((res[,37] >= res[,8]) & (res[,37] <= res[,9])),
               covSIE_tmle_eff = mean((res[,38] >= res[,80]) & (res[,38] <= res[,81])),
               covSIE_1s_eff = mean((res[,38] >= res[,86]) & (res[,38] <= res[,87])),
               covSIE_tmle = mean((res[,38] >= res[,5]) & (res[,38] <= res[,6])),
               covSIE_1s = mean((res[,38] >= res[,11]) & (res[,38] <= res[,12])),
               covSDE_tmleB_eff = mean((res[,37] >= res[,89]) & (res[,37] <= res[,90])), 
               covSDE_1sB_eff = mean((res[,37] >= res[,95]) & (res[,37] <= res[,96])),
               covSDE_tmleB = mean((res[,37] >= res[,20]) & (res[,37] <= res[,21])), 
               covSDE_1sB = mean((res[,37] >= res[,26]) & (res[,37] <= res[,27])),
               covSIE_tmleB_eff = mean((res[,38] >= res[,86]) & (res[,38] <= res[,87])),
               covSIE_1sB_eff = mean((res[,38] >= res[,92]) & (res[,38] <= res[,93])),
               covSIE_tmleB = mean((res[,38] >= res[,23]) & (res[,38] <= res[,24])), 
               covSIE_1sB = mean((res[,38] >= res[,29]) & (res[,38] <= res[,30])))
  
  # res[, c(38,23,24,29,30,35,36)]
  res_SEboot = cbind(
    (res[,90] - res[,89])/(2*1.96),
    (res[,96] - res[,95])/(2*1.96),
    (res[,21] - res[,20])/(2*1.96),
    (res[,27] - res[,26])/(2*1.96),
    (res[,87] - res[,86])/(2*1.96),
    (res[,93] - res[,92])/(2*1.96),
    (res[,24] - res[,23])/(2*1.96),
    (res[,30] - res[,29])/(2*1.96))

  res_SE_SDE = res[,c(100,102,39,41)]
  res_SE_SIE = res[,c(101,103,40,42)]
  # 
  # apply(res_SE_SIE, 2, FUN = function(x) max(x))
  res_SE_SDEboot = res_SEboot[,1:4]
  res_SE_SIEboot = res_SEboot[,5:8]
  
  # apply(res_SE_SIEboot, 2,max)
  # hist(res_SE_SIEboot[,2] - res_SE_SIEboot[,1])
  
  res_SE_0 = res[,51:52]
  res_SE_eff0 = res[,53:54]
  efficiency_SDE = colMeans(apply(res_SE_SDE, 2, FUN = function(x) x/res_SE_eff0[,1]))
  efficiency_SIE = colMeans(apply(res_SE_SIE, 2, FUN = function(x) x/res_SE_eff0[,2]))
  # efficiency_SDE
  # efficiency_SIE
  
  if (bootcut) {
  efficiency_SDEboot = unlist(lapply(lapply(1:4, FUN = function(x) {
    c = res_SE_SDEboot[,x]
    (c/res_SE_eff0[,1])[c<=20]
    }), mean))
  efficiency_SIEboot = unlist(lapply(lapply(1:4,FUN = function(x) {
    c = res_SE_SIEboot[,x]
    (c/res_SE_eff0[,2])[c<=20]
  }), mean))
  } else {
  efficiency_SDEboot = colMeans(apply(res_SE_SDEboot, 2, FUN = function(x) x/res_SE_eff0[,1]))
  efficiency_SIEboot = colMeans(apply(res_SE_SIEboot, 2, FUN = function(x) x/res_SE_eff0[,2]))
  }
  # efficiency_SDEboot
  # efficiency_SIEboot
  
  res_ests_SDE = res[,c(76, 82, 1,7, 37)]
  res_ests_SIE = res[,c(79, 85, 4,10, 38)]
  
  res_bias_SDE = cbind(res_ests_SDE[,1] - res_ests_SDE[,5],
                       res_ests_SDE[,2] - res_ests_SDE[,5],
                       res_ests_SDE[,3] - res_ests_SDE[,5],
                       res_ests_SDE[,4] - res_ests_SDE[,5],
                       res_ests_SDE[,5] - res_ests_SDE[,5])
  
  res_bias_SIE = cbind(res_ests_SIE[,1] - res_ests_SIE[,5],
                       res_ests_SIE[,2] - res_ests_SIE[,5],
                       res_ests_SIE[,3] - res_ests_SIE[,5],
                       res_ests_SIE[,4] - res_ests_SIE[,5],
                       res_ests_SIE[,5] - res_ests_SIE[,5])
  
  coverage = rbind(coverage[1:8], coverage[9:16])
  rownames(coverage) = c("IC variance", "Bootstrap")
  colnames(coverage) = c("SDE_tmle_eff", "SDE_EE_eff","SDE_tmle", "SDE_EE", 
                         "SIE_tmle_eff", "SIE_EE_eff","SIE_tmle", "SIE_EE")
  
  # precentage SE of of efficiency bound
  efficiency_percentageTrue = rbind(100*c(efficiency_SDE, efficiency_SIE),
                                    100*c(efficiency_SDEboot, efficiency_SIEboot))
  
  colnames(efficiency_percentageTrue) = colnames(coverage)
  rownames(efficiency_percentageTrue) = c("IC variance", "Bootstrap")
  
  
  # precentage bias of SE
  bias_perc_SDE = 100*unlist(lapply(lapply(1:4, FUN = function(x) res_bias_SDE[,x]/res_SE_SDE[,x]), mean))
  bias_perc_SIE = 100*unlist(lapply(lapply(1:4, FUN = function(x) res_bias_SIE[,x]/res_SE_SIE[,x]), mean))
  bias_percentageSE = c(bias_perc_SDE, bias_perc_SIE)
  
  bias_perc_SDE = 100*unlist(lapply(lapply(1:4, FUN = function(x) res_bias_SDE[,x]/res_SE_SDE[,x]), mean))
  bias_perc_SIE = 100*unlist(lapply(lapply(1:4, FUN = function(x) res_bias_SIE[,x]/res_SE_SIE[,x]), mean))
  bias  = colMeans(cbind(res_bias_SDE, res_bias_SIE))
  names(bias) = c("SDE_tmle_eff", "SDE_EE_eff", "SDE_iptw", "SIE_tmle", "SIE_EE")
  
  names(bias_percentageSE) = c("SDE_tmle_eff", "SDE_EE_eff", "SDE_iptw", "SIE_tmle", "SIE_EE")
  
  RMSE = c(SDE_tmle_eff = sqrt(mean((res$res.SDE_0 - res[,76])^2)),
           SDE_1s_eff = sqrt(mean((res$res.SDE_0 - res[,82])^2)),
           SDE_tmle = sqrt(mean((res$res.SDE_0 - res[,1])^2)),
           SDE_1s = sqrt(mean((res$res.SDE_0 - res[,7])^2)),
           SIE_tmle = sqrt(mean((res$res.SIE_0 - res[,79])^2)),
           SIE_1s = sqrt(mean((res$res.SIE_0 - res[,85])^2)),
           SIE_tmle = sqrt(mean((res$res.SIE_0 - res[,4])^2)),
           SIE_1s = sqrt(mean((res$res.SIE_0 - res[,10])^2)))
  
  Out_of_bounds = 100*c(
  mean(res[,76]<= -1 | res[,70] >= 1),
  mean(res[,82]<= -1 | res[,76] >= 1),
  mean(res[,1]<= -1 | res[,1] >= 1),
  mean(res[,7]<= -1 | res[,7] >= 1),
  mean(res[,79]<= -1 | res[,73] >= 1),
  mean(res[,85]<= -1 | res[,79] >= 1),
  mean(res[,4]<= -1 | res[,4] >= 1),
  mean(res[,10]<= -1 | res[,10] >= 1))
  names(Out_of_bounds) = colnames(coverage)
  return(list(results = res,
              coverage = coverage, 
              efficiency_percentageTrue = efficiency_percentageTrue,
              bias_percentageSE = bias_percentageSE,
              bias = bias,
              RMSE = RMSE,
              Out_of_bounds = Out_of_bounds,
              valid_draws = length(valid_draws),
              func_list = func_list,
              forms = forms))
  
}
