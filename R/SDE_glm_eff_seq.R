# Simulation function
#' @export
SDE_glm_eff_seq = function(data, forms, RCT = 0.5, transport, pooled, gstar_S = 0, truth, B = 500) 
{
  
  if (length(unique(data$Y))!=2) data$Y = (data$Y - min(data$Y))/(max(data$Y) - min(data$Y))
  if (!transport) pooled = FALSE
  # get the stochastic dist of M and true params if you want 
  gstar_info = get_gstarM_glm_eff_seqT(data = data, truth = truth, forms = forms, transport = transport, 
                                       pooled = pooled, gstar_S = gstar_S)
  # get the stochastic dist of M and true params if you want 
  gstarM_astar1 = gstar_info$gstarM_astar1
  gstarM_astar0 = gstar_info$gstarM_astar0
  gstarM_astar = list(gstarM_astar0 = gstarM_astar0, gstarM_astar1 = gstarM_astar1)
  
  # perform initial fits for the first regression
  init_info = get.mediation.initdata_glm_eff_seqT(data = data, forms = forms, RCT = RCT, pooled, transport)
  Y_preds = init_info$Y_preds
  est_info = lapply(0:1, FUN = function(astar) {
    lapply(0:1, FUN = function(a) {
      # get tmle info
      # get iptw here while I'm at it
      update = mediation.step1_glm_eff_seqT(initdata = init_info$initdata, init_info$Y_preds, data = data, 
                               gstarM_astar[[astar+1]], a = a, transport)
      
      Y_Mg = get.stochasticM(gstarM_astar[[astar+1]], Y_preds[[2]], Y_preds[[3]]) 
      A_ps = init_info$initdata$A_ps
      EE_gcomp_info = mediation.step2_glm_eff_seqT(data = data, Qstar_M = Y_preds[[1]], 
                      Qstar_Mg = Y_Mg, Hm = update$Hm, A_ps = A_ps, 
                      a = a, tmle = FALSE,
                      EE = TRUE, bootstrap = FALSE, form = forms$QZform, transport)

      Qstar_Mg = get.stochasticM(gstarM_astar[[astar+1]], update$Qstar_M1, update$Qstar_M0) 
      # compute Qstar_Mg here
      tmle_info = mediation.step2_glm_eff_seqT(data = data, Qstar_M = update$Qstar_M, 
                      Qstar_Mg = Qstar_Mg, Hm = update$Hm, A_ps = A_ps, 
                      a = a, tmle = TRUE,
                      EE = FALSE, bootstrap = FALSE, form = forms$QZform, transport)
      
      # compile all estimates
      tmle_info$eps1 = update$eps
      return(list(EE_gcomp_info = EE_gcomp_info, tmle_info = tmle_info))
    }) 
  })
  
  # run the bootstrap 500 times
  # bootstrap from the data and thus the gstarM_astar1 and gstarM_astar0
  n = nrow(data)
  if (!is.null(B)) {
    boot_ests = lapply(1:B, FUN = function(x) {
    inds = sample(1:n, replace = TRUE)
    data = data[inds,]
    init_info = get.mediation.initdata_glm_eff_seqT(data = data, forms = forms, RCT = RCT, pooled, transport)
    Y_preds = init_info$Y_preds
    gstarM_astar = list(gstarM_astar0[inds], gstarM_astar1[inds])
    
    est_info = lapply(0:1, FUN = function(astar) {
      return(lapply(0:1, FUN = function(a) {
        update = mediation.step1_glm_eff_seqT(initdata = init_info$initdata, Y_preds = Y_preds, data = data, 
                                 gstarM_astar[[astar+1]], a, transport)
        iptw_info = update$est_iptw
        
        Y_Mg = get.stochasticM(gstarM_astar[[astar+1]], Y_preds[[2]], Y_preds[[3]]) 
        A_ps = init_info$initdata$A_ps
        EE_mle_info = mediation.step2_glm_eff_seqT(data = data, Qstar_M = Y_preds[[1]], 
                                      Qstar_Mg = Y_Mg, Hm = update$Hm, A_ps = A_ps, 
                                      a = a, tmle = FALSE,
                                      EE = TRUE, bootstrap = TRUE, form = forms$QZform, transport)
        
        Qstar_Mg = get.stochasticM(gstarM_astar[[astar+1]], update$Qstar_M1, update$Qstar_M0) 
        # compute Qstar_Mg here
        tmle_info = mediation.step2_glm_eff_seqT(data = data, Qstar_M = update$Qstar_M, 
                                    Qstar_Mg = Qstar_Mg, Hm = update$Hm, A_ps = A_ps, 
                                    a = a, tmle = TRUE,
                                    EE = FALSE, bootstrap = TRUE, form = forms$QZform, transport)
        # compile all estimates
        return(list(tmle_est = tmle_info, EE_est = EE_mle_info[1],
                    mle_est = EE_mle_info[2]))
      }))
    })
  })
  
    boot_ests = lapply(boot_ests, FUN = function(boot) {
      SDE = unlist(boot[[1]][[2]]) - unlist(boot[[1]][[1]])
      SIE = unlist(boot[[2]][[2]]) - unlist(boot[[1]][[2]])
      return(list(SDE, SIE))
    })
    
    boots_SDE = do.call(rbind, lapply(boot_ests, FUN = function(boot) boot[[1]]))
    boots_SIE = do.call(rbind, lapply(boot_ests, FUN = function(boot) boot[[2]]))
    
    bootSE_SDE = apply(boots_SDE, 2, sd)
    bootSE_SIE = apply(boots_SIE, 2, sd)
  }
  
  D_SDE = est_info[[1]][[2]]$tmle_info$IC - est_info[[1]][[1]]$tmle_info$IC
  D_SIE = est_info[[2]][[2]]$tmle_info$IC- est_info[[1]][[2]]$tmle_info$IC
  
  D_SDE_1s = est_info[[1]][[2]]$EE_gcomp_info$IC - est_info[[1]][[1]]$EE_gcomp_info$IC
  D_SIE_1s = est_info[[2]][[2]]$EE_gcomp_info$IC- est_info[[1]][[2]]$EE_gcomp_info$IC

  SE_SDE = sd(D_SDE)/sqrt(n)
  SE_SIE = sd(D_SIE)/sqrt(n)
  
  SE_SDE_1s = sd(D_SDE_1s)/sqrt(n)
  SE_SIE_1s = sd(D_SIE_1s)/sqrt(n)

  if (!is.null(truth)) { 
  D_SDE_0 = gstar_info$D_astar0a1_0 - gstar_info$D_astar0a0_0
  D_SIE_0 = gstar_info$D_astar1a1_0 - gstar_info$D_astar0a1_0
  
  SE_SDE_0 = sd(D_SDE_0)/sqrt(n)
  SE_SIE_0 = sd(D_SIE_0)/sqrt(n)
  
  Psi_astar0a1_0 = gstar_info$Psi_astar0a1_0
  Psi_astar0a0_0 = gstar_info$Psi_astar0a0_0
  Psi_astar1a1_0 = gstar_info$Psi_astar1a1_0
  
  } else {
    D_SDE_0 = D_SIE_0 = SE_SDE_0 = SE_SIE_0 = NULL
    Psi_astar0a1_0 = gstar_info$Psi_astar0a1_0
    Psi_astar0a0_0 = gstar_info$Psi_astar0a0_0
    Psi_astar1a1_0 = gstar_info$Psi_astar1a1_0
  }
  ests_astar0a1 =   c(tmle = est_info[[1]][[2]]$tmle_info$est,
                      EE = est_info[[1]][[2]]$EE_gcomp_info$est,
                      iptw = est_info[[1]][[2]]$iptw_info[[2]],
                      mle = est_info[[1]][[2]]$EE_gcomp_info$init_est,
                      Psi_0 = Psi_astar0a1_0)
  
  ests_astar0a0 =   c(tmle = est_info[[1]][[1]]$tmle_info$est,
                      EE = est_info[[1]][[1]]$EE_gcomp_info$est,
                      iptw = est_info[[1]][[1]]$iptw_info[[2]],
                      mle = est_info[[1]][[1]]$EE_gcomp_info$init_est,
                      Psi_0 = Psi_astar0a0_0)
  
  ests_astar1a1 =   c(tmle = est_info[[2]][[2]]$tmle_info$est,
                      EE = est_info[[2]][[2]]$EE_gcomp_info$est,
                      iptw = est_info[[2]][[2]]$iptw_info[[2]],
                      mle = est_info[[2]][[2]]$EE_gcomp_info$init_est,
                      Psi_0 = Psi_astar1a1_0)
  
  SDE_ests = ests_astar0a1 - ests_astar0a0
  SIE_ests = ests_astar1a1 - ests_astar0a1
  
  CI_SDE = c(SDE_ests[1], SDE_ests[1] - 1.96*SE_SDE, SDE_ests[1] + 1.96*SE_SDE)
  CI_SIE = c(SIE_ests[1], SIE_ests[1] - 1.96*SE_SIE, SIE_ests[1] + 1.96*SE_SIE)
  
  CI_SDE_1s = c(SDE_ests[2], SDE_ests[2] - 1.96*SE_SDE_1s , SDE_ests[2] + 1.96*SE_SDE_1s)
  CI_SIE_1s  = c(SIE_ests[2], SIE_ests[2] - 1.96*SE_SIE_1s , SIE_ests[2] + 1.96*SE_SIE_1s)

  if (!is.null(B)) {
    CI_SDE_boot = c(SDE_ests[1], SDE_ests[1] - 1.96*bootSE_SDE[1], SDE_ests[1] + 1.96*bootSE_SDE[1])
    CI_SIE_boot = c(SIE_ests[1], SIE_ests[1] - 1.96*bootSE_SIE[1], SIE_ests[1] + 1.96*bootSE_SIE[1])
    
    CI_SDE_1s_boot = c(SDE_ests[2], SDE_ests[2] - 1.96*bootSE_SDE[2] , SDE_ests[2] + 1.96*bootSE_SDE[2])
    CI_SIE_1s_boot  = c(SIE_ests[2], SIE_ests[2] - 1.96*bootSE_SIE[2] , SIE_ests[2] + 1.96*bootSE_SIE[2])

    return(list(CI_SDE = CI_SDE, CI_SIE = CI_SIE, CI_SDE_1s = CI_SDE_1s, CI_SIE_1s = CI_SIE_1s,
                CI_SDE_boot = CI_SDE_boot, 
                CI_SIE_boot = CI_SIE_boot, CI_SDE_1s_boot = CI_SDE_1s_boot, CI_SIE_1s_boot = CI_SIE_1s_boot, 
                SDE_0 = Psi_astar0a1_0 - Psi_astar0a0_0,
                SIE_0 = Psi_astar1a1_0 - Psi_astar0a1_0, 
                SE_SDE = SE_SDE, 
                SE_SIE = SE_SIE, 
                SE_SDE_1s = SE_SDE_1s,  
                SE_SIE_1s = SE_SIE_1s, 
                SE_SDE_0 = SE_SDE_0, 
                SE_SIE_0 = SE_SIE_0,
                SE_SDE_boot=bootSE_SDE[1],
                SE_SIE_boot=bootSE_SIE[1],
                SE_SDE_1s_boot=bootSE_SDE[2],
                SE_SIE_1s_boot=bootSE_SIE[2],
                SE_SDE_iptw_boot=bootSE_SDE[3],
                SE_SIE_iptw_boot=bootSE_SIE[3],
                
                ests_astar0a1 = ests_astar0a1,
                ests_astar0a0 = ests_astar0a0,
                ests_astar1a1 = ests_astar1a1,
                
                eps1_astar0a1 = est_info[[1]][[2]]$tmle_info$eps1,
                eps2_astar0a1 = est_info[[1]][[2]]$tmle_info$eps2,
                eps1_astar0a0 = est_info[[1]][[1]]$tmle_info$eps1,
                eps2_astar0a0 = est_info[[1]][[1]]$tmle_info$eps2,
                eps1_astar1a1 = est_info[[2]][[2]]$tmle_info$eps1,
                eps2_astar1a1 = est_info[[2]][[2]]$tmle_info$eps2 
    )) 
  } else {
    
    return(list(CI_SDE = CI_SDE, CI_SIE = CI_SIE, CI_SDE_1s = CI_SDE_1s, CI_SIE_1s = CI_SIE_1s,
                SDE_0 = Psi_astar0a1_0 - Psi_astar0a0_0,
                SIE_0 = Psi_astar1a1_0 - Psi_astar0a1_0, 
                SE_SDE = SE_SDE, 
                SE_SIE = SE_SIE, 
                SE_SDE_1s = SE_SDE_1s,  
                SE_SIE_1s = SE_SIE_1s, 
                SE_SDE_0 = SE_SDE_0, 
                SE_SIE_0 = SE_SIE_0,
                
                ests_astar0a1 = ests_astar0a1,
                ests_astar0a0 = ests_astar0a0,
                ests_astar1a1 = ests_astar1a1,
                
                eps1_astar0a1 = est_info[[1]][[2]]$tmle_info$eps1,
                eps2_astar0a1 = est_info[[1]][[2]]$tmle_info$eps2,
                eps1_astar0a0 = est_info[[1]][[1]]$tmle_info$eps1,
                eps2_astar0a0 = est_info[[1]][[1]]$tmle_info$eps2,
                eps1_astar1a1 = est_info[[2]][[2]]$tmle_info$eps1,
                eps2_astar1a1 = est_info[[2]][[2]]$tmle_info$eps2 
    ))
    
  }
} 


#' @export
get_gstarM_glm_eff_seqT  = function(data, truth, forms, transport, pooled, gstar_S) 
{
  W = as.data.frame(data[,grep("W", colnames(data))])
  colnames(W) = colnames(data)[grep("W", colnames(data))]
  
  Mstarform = forms$Mstarform
  Zstarform = forms$Zstarform

  if (transport & !pooled) {
    Mstarfit = glm(Mstarform, data[data$S==gstar_S,], family = binomial)
    Zstarfit = glm(Zstarform, data[data$S==gstar_S,], family = binomial)
  } else {
    Mstarfit = glm(Mstarform, data, family = binomial)
    Zstarfit = glm(Zstarform, data, family = binomial)
  }
  
  # find the truth if simulating
  f_truth0 = function(Z, W, S) {
    # W = W[1:10]
    nn = length(Z)
    # dataM1 = data.frame(cbind)
    if (transport) {
    dataM1 = cbind(Z = rep(1,nn), W, S = rep(gstar_S,nn))
    predM1 = predict(Mstarfit, newdata = dataM1, type = 'response')
    
    dataM0 = cbind(Z = rep(0,nn), W, S = rep(gstar_S,nn))
    predM0 = predict(Mstarfit, newdata = dataM0, type = 'response')
    
    dataZ = cbind(A = rep(0,nn), W, S = rep(gstar_S, nn))
    predZ = predict(Zstarfit, newdata = dataZ, type = 'response')
    } else {
      dataM1 = cbind(Z = rep(1,nn), W)
      predM1 = predict(Mstarfit, newdata = dataM1, type = 'response')
      
      dataM0 = cbind(Z = rep(0,nn), W)
      predM0 = predict(Mstarfit, newdata = dataM0, type = 'response')
      
      dataZ = cbind(A = rep(0,nn), W)
      predZ = predict(Zstarfit, newdata = dataZ, type = 'response')
    }
    
    gM = predM1*predZ + predM0*(1 - predZ)
    return(gM)
  }
  
  f_truth1 = function(Z, W, S) {
    # W = W[1:10]
    nn = length(Z)
    if (transport) {
    dataM1 = cbind(Z = rep(1,nn), W, S = rep(gstar_S,nn))
    predM1 = predict(Mstarfit, newdata = dataM1, type = 'response')
    
    dataM0 = cbind(Z = rep(0,nn), W, S = rep(gstar_S,nn))
    predM0 = predict(Mstarfit, newdata = dataM0, type = 'response')
    
    dataZ = cbind(A = rep(1,nn), W, S = rep(gstar_S, nn))
    predZ = predict(Zstarfit, newdata = dataZ, type = 'response')
    } else {
      dataM1 = cbind(Z = rep(1,nn), W)
      predM1 = predict(Mstarfit, newdata = dataM1, type = 'response')
      
      dataM0 = cbind(Z = rep(0,nn), W)
      predM0 = predict(Mstarfit, newdata = dataM0, type = 'response')
      
      dataZ = cbind(A = rep(1,nn), W)
      predZ = predict(Zstarfit, newdata = dataZ, type = 'response')
    }
    gM = predM1*predZ + predM0*(1 - predZ)
    return(gM)
  }
  
  gstarM_astar1 = with(data, f_truth1(Z=Z, W=W, S=S))
  gstarM_astar0 = with(data, f_truth0(Z=Z, W=W, S=S))
  
  if (!is.null(truth)) {
    if (transport) {
    data_pop_astar0a1 = gendata.SDEtransport(5*1e6, f_W = truth$f_W, 
                                             f_S = truth$f_S, f_A = function(S,W) 1, 
                                             f_Z = truth$f_Z, f_M = f_truth0, f_Y = truth$f_Y) 
    data_pop_astar0a0 = gendata.SDEtransport(5*1e6, f_W = truth$f_W, 
                                             f_S = truth$f_S, f_A = function(S,W) 0, 
                                             f_Z = truth$f_Z, f_M = f_truth0, f_Y = truth$f_Y) 
    
    data_pop_astar1a1 = gendata.SDEtransport(5*1e6, f_W = truth$f_W, 
                                             f_S = truth$f_S, f_A = function(S,W) 1, 
                                             f_Z = truth$f_Z, f_M = f_truth1, f_Y = truth$f_Y) 
    
    Psi_astar0a1_0 = mean(data_pop_astar0a1$Y[data_pop_astar0a1$S==0]) 
    Psi_astar0a0_0 = mean(data_pop_astar0a0$Y[data_pop_astar0a0$S==0])
    Psi_astar1a1_0 = mean(data_pop_astar1a1$Y[data_pop_astar1a1$S==0]) 
    PS0_0 = mean(data_pop_astar1a1$S==0)   
    
    # get the true IC's
    S_ps0 = with(data, truth$f_S(W=W))
    M_ps0 = with(data, truth$f_M(Z=Z, W=W, S=S))
    Z_ps0 = with(data, truth$f_Z(A=A, W=W, S=S))
    A_ps0 = with(data, truth$f_A(W=W, S=S))
    Z_psWS0 = with(data, truth$f_Z(A=1, W=W, S=S)*A_ps0 + 
                     truth$f_Z(A=0, W=W, S=S)*(1-A_ps0))
    
    
    get_cc_0 = function(data, gstarM_astar, a) {
      df_ZS0 = data
      df_ZS0$S = 0
      df_ZS0$A = a
      ZS0_ps0 = with(df_ZS0, truth$f_Z(A=A, W=W, S=S))
      with(data, ((S == 1)*
                    ((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))*
                    ((Z == 1)*ZS0_ps0 + (Z == 0)*(1 - ZS0_ps0))*(1 - S_ps0))/
             (((M == 1)*M_ps0 + (M == 0)*(1 - M_ps0))*
                ((Z == 1)*Z_psWS0 + (Z == 0)*(1 - Z_psWS0))*S_ps0
              *PS0_0))
    }
    
    get_Hz_0 = function(a) with(data, (A == a)*(S == 0)/(A*A_ps0 + (1 - A)*(1 - A_ps0))/PS0_0)
    
    Hm_astar0a1_0 = get_cc_0(data = data, gstarM_astar = gstarM_astar0, a = 1)
    Hm_astar0a0_0 = get_cc_0(data = data, gstarM_astar = gstarM_astar0, a = 0)
    Hm_astar1a1_0 = get_cc_0(data = data, gstarM_astar = gstarM_astar1, a = 1)
    
    Hz_astar0a1_0 = get_Hz_0(1)
    Hz_astar0a0_0 = get_Hz_0(0)
    Hz_astar1a1_0 = get_Hz_0(1)
    
    get_trueIC = function(gstar_M, a, Psi_0, Hm_0, Hz_0) {
      p_0Z = with(data, truth$f_Z(A=a, W=W, S=S))
      Q_0Y = with(data, truth$f_Y(M=M, Z=Z, W=W))
      Q_0YM1 = with(data, truth$f_Y(M=1, Z=Z, W=W))
      Q_0YM0 = with(data, truth$f_Y(M=0, Z=Z, W=W))
      
      Q_0YM1Z1 = with(data, truth$f_Y(M=1, Z=1, W=W))
      Q_0YM0Z1 = with(data, truth$f_Y(M=0, Z=1, W=W))
      Q_0YM1Z0 = with(data, truth$f_Y(M=1, Z=0, W=W))
      Q_0YM0Z0 = with(data, truth$f_Y(M=0, Z=0, W=W))
      
      Q_ghat_astar = with(data, Q_0YM1*gstar_M + Q_0YM0*(1 - gstar_M)) 
      Q_ghat_astarZ1 = with(data, Q_0YM1Z1*gstar_M + Q_0YM0Z1*(1 - gstar_M)) 
      Q_ghat_astarZ0 = with(data, Q_0YM1Z0*gstar_M + Q_0YM0Z0*(1 - gstar_M)) 
      
      Q_ghat_astarW = with(data, Q_ghat_astarZ1*p_0Z + Q_ghat_astarZ0*(1-p_0Z))
      D_0Y = with(data, Hm_0*(Y - Q_0Y))
      D_0Z = Hz_0*(Q_ghat_astar - Q_ghat_astarW)
      D_0W = with(data, (Q_ghat_astarW - Psi_0)*(S ==0)/PS0_0)
      D_0 = D_0Y + D_0Z + D_0W
      return(D_0)
    }
    } else {
      data_pop_astar0a1 = gendata.SDEtransport_alt(5*1e6, f_W = truth$f_W, 
                                               f_A = function(W) 1, 
                                               f_Z = truth$f_Z, f_M = f_truth0, f_Y = truth$f_Y) 
      data_pop_astar0a0 = gendata.SDEtransport_alt(5*1e6, f_W = truth$f_W, 
                                               f_A = function(W) 0, 
                                               f_Z = truth$f_Z, f_M = f_truth0, f_Y = truth$f_Y) 
      
      data_pop_astar1a1 = gendata.SDEtransport_alt(5*1e6, f_W = truth$f_W, 
                                               f_A = function(W) 1, 
                                               f_Z = truth$f_Z, f_M = f_truth1, f_Y = truth$f_Y) 
      
      Psi_astar0a1_0 = mean(data_pop_astar0a1$Y) 
      Psi_astar0a0_0 = mean(data_pop_astar0a0$Y)
      Psi_astar1a1_0 = mean(data_pop_astar1a1$Y) 
      
      # get the true IC's
      M_ps0 = with(data, truth$f_M(Z=Z, W=W))
      Z_ps0 = with(data, truth$f_Z(A=A, W=W))
      A_ps0 = with(data, truth$f_A(W=W))
      Z_psW = with(data, truth$f_Z(A=1, W=W)*A_ps0 + 
                       truth$f_Z(A=0, W=W)*(1-A_ps0))
      
      
      get_cc_0 = function(data, gstarM_astar, a) {
        df_ZAa = data
        df_ZAa$A = a
        ZAa_ps0 = with(df_ZAa, truth$f_Z(A=A, W=W))
        with(data, (((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))*
                      ((Z == 1)*ZAa_ps0 + (Z == 0)*(1 - ZAa_ps0))/
               (((M == 1)*M_ps0 + (M == 0)*(1 - M_ps0))*
                  ((Z == 1)*Z_psW + (Z == 0)*(1 - Z_psW)))))
      }
      
      
      get_Hz_0 = function(a) with(data, (A == a)/(A*A_ps0 + (1 - A)*(1 - A_ps0)))

      Hm_astar0a1_0 = get_cc_0(data = data, gstarM_astar = gstarM_astar0, a = 1)
      Hm_astar0a0_0 = get_cc_0(data = data, gstarM_astar = gstarM_astar0, a = 0)
      Hm_astar1a1_0 = get_cc_0(data = data, gstarM_astar = gstarM_astar1, a = 1)
      
      Hz_astar0a1_0 = get_Hz_0(1)
      Hz_astar0a0_0 = get_Hz_0(0)
      Hz_astar1a1_0 = get_Hz_0(1)
      
      get_trueIC = function(gstar_M, a, Psi_0, Hm_0, Hz_0) {
        p_0Z = with(data, truth$f_Z(A=a, W=W))
        Q_0Y = with(data, truth$f_Y(M=M, Z=Z, W=W))
        Q_0YM1 = with(data, truth$f_Y(M=1, Z=Z, W=W))
        Q_0YM0 = with(data, truth$f_Y(M=0, Z=Z, W=W))
        
        Q_0YM1Z1 = with(data, truth$f_Y(M=1, Z=1, W=W))
        Q_0YM0Z1 = with(data, truth$f_Y(M=0, Z=1, W=W))
        Q_0YM1Z0 = with(data, truth$f_Y(M=1, Z=0, W=W))
        Q_0YM0Z0 = with(data, truth$f_Y(M=0, Z=0, W=W))
        
        Q_ghat_astar = with(data, Q_0YM1*gstar_M + Q_0YM0*(1 - gstar_M)) 
        Q_ghat_astarZ1 = with(data, Q_0YM1Z1*gstar_M + Q_0YM0Z1*(1 - gstar_M)) 
        Q_ghat_astarZ0 = with(data, Q_0YM1Z0*gstar_M + Q_0YM0Z0*(1 - gstar_M)) 
        
        Q_ghat_astarW = with(data, Q_ghat_astarZ1*p_0Z + Q_ghat_astarZ0*(1-p_0Z))
        D_0Y = with(data, Hm_0*(Y - Q_0Y))
        D_0Z = Hz_0*(Q_ghat_astar - Q_ghat_astarW)
        D_0W = with(data, (Q_ghat_astarW - Psi_0))
        D_0 = D_0Y + D_0Z + D_0W
        return(D_0)
      }
    }
    
    D_astar0a1_0 = get_trueIC(gstar_M = gstarM_astar0, a = 1, 
                              Psi_0 = Psi_astar0a1_0, Hm_0 = Hm_astar0a1_0,
                              Hz_0 = Hz_astar0a1_0)
    
    D_astar0a0_0 = get_trueIC(gstar_M = gstarM_astar0, a = 0, 
                              Psi_0 = Psi_astar0a0_0, Hm_0 = Hm_astar0a0_0,
                              Hz_0 = Hz_astar0a0_0)
    
    D_astar1a1_0 = get_trueIC(gstar_M = gstarM_astar1, a = 1, 
                              Psi_0 = Psi_astar1a1_0, Hm_0 = Hm_astar1a1_0,
                              Hz_0 = Hz_astar1a1_0)
    
    return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0, 
                Psi_astar0a1_0 = Psi_astar0a1_0, Psi_astar0a0_0 = Psi_astar0a0_0, 
                Psi_astar1a1_0 = Psi_astar1a1_0,  
                D_astar0a1_0 = D_astar0a1_0, D_astar0a0_0 = D_astar0a0_0,
                D_astar1a1_0 = D_astar1a1_0, 
                Hm_astar0a1_0 = Hm_astar0a1_0, 
                Hm_astar0a0_0 = Hm_astar0a0_0,
                Hm_astar1a1_0 = Hm_astar1a1_0,
                Hz_astar0a1_0 = Hz_astar0a1_0,
                Hz_astar0a0_0 = Hz_astar0a0_0,
                Hz_astar1a1_0 = Hz_astar1a1_0))
  } else return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0))
}

#' @export
get.mediation.initdata_glm_eff_seqT = function(data, forms, RCT = 0.5, pooled, transport) {
  
  if (pooled | !transport) {
    Zform = forms$Zstarform
    Mform = forms$Mstarform
  } else {
    Zform = forms$Zform
    Mform = forms$Mform
  }
  if (transport) {

  df_YM1S1 = data
  df_YM1S1$M = 1
  df_YM1S1$S = 1
  
  df_YM0S1 = df_YM1S1
  df_YM0S1$M = 0
  
  df_ZA1S0 = df_ZA0S0 = data
  df_ZA1S0$S = df_ZA0S0$S = 0
  df_ZA1S0$A = 1
  df_ZA0S0$A = 0
  
  Yform = forms$Yform
  Yfit = glm(formula = Yform, data = data[data$S==1,], family = binomial())
  
  Mfit = glm(formula = Mform, data = data, family = binomial())
  Zfit = glm(formula = Zform, data = data, family = binomial())
  
  if (is.null(RCT)) { 
    Aform = forms$Aform
    Afit = glm(formula = Aform, data = data, family = binomial())
    A_ps = predict(Afit, type = 'response')
  } else A_ps = RCT
  
  # Sform = paste0("S ~ ", paste(covariates$covariates_S, collapse = "+"))
  Sform = forms$Sform
  Sfit = glm(formula = Sform, data = data, family = binomial())
  
  # propensity scores
  S_ps = predict(Sfit, type = 'response')
  PS0 = mean(data$S == 0)
  if (is.null(RCT)) A_ps = predict(Afit, type = 'response') else A_ps = RCT
  df_ZA0 = df_ZA1 = data 
  df_ZA0$A = 0
  df_ZA1$A = 1
  Z_psWS = predict(Zfit, newdata = df_ZA0, type = 'response')*(1 - A_ps) +
    predict(Zfit, newdata = df_ZA1, type = 'response')*A_ps
  ZA1S0_ps = predict(Zfit, newdata = df_ZA1S0, type = 'response')
  ZA0S0_ps = predict(Zfit, newdata = df_ZA0S0, type = 'response')
  M_ps = predict(Mfit, newdata = data, type = 'response')
  
  Y_init = pmin(pmax(predict(Yfit, newdata = data, type = 'response'), 0.001), .999)
  Y_init_M1 = pmin(pmax(predict(Yfit, newdata = df_YM1S1, type = 'response'), 0.001), .999)
  Y_init_M0 = pmin(pmax(predict(Yfit, newdata = df_YM0S1, type = 'response'), 0.001), .999)
  return(list(initdata = list(M_ps = M_ps, ZA1S0_ps = ZA1S0_ps, ZA0S0_ps = ZA0S0_ps,
                              Z_psWS = Z_psWS, 
                              A_ps = A_ps, S_ps = S_ps, PS0 = PS0), 
              Y_preds = list(Y_init = Y_init, 
                             Y_init_M1 = Y_init_M1, 
                             Y_init_M0 = Y_init_M0)))
} else {

  df_YM1 = data
  df_YM1$M = 1
  
  df_YM0 = df_YM1
  df_YM0$M = 0
  
  df_ZA1 = df_ZA0 = data
  df_ZA1$A = 1
  df_ZA0$A = 0
  
  Yform = forms$Yform
  Yfit = glm(formula = Yform, data = data, family = binomial())
  Mfit = glm(formula = Mform, data = data, family = binomial())
  Zfit = glm(formula = Zform, data = data, family = binomial())
  
  if (is.null(RCT)) { 
    Aform = forms$Aform
    Afit = glm(formula = Aform, data = data, family = binomial())
    A_ps = predict(Afit, type = 'response')
  } else A_ps = RCT
  
  # propensity scores
  if (is.null(RCT)) A_ps = predict(Afit, type = 'response') else A_ps = RCT
  df_ZA0 = df_ZA1 = data 
  df_ZA0$A = 0
  df_ZA1$A = 1
  Z_psW = predict(Zfit, newdata = df_ZA0, type = 'response')*(1 - A_ps) +
    predict(Zfit, newdata = df_ZA1, type = 'response')*A_ps
  ZA1_ps = predict(Zfit, newdata = df_ZA1, type = 'response')
  ZA0_ps = predict(Zfit, newdata = df_ZA0, type = 'response')
  M_ps = predict(Mfit, newdata = data, type = 'response')
  
  Y_init = pmin(pmax(predict(Yfit, newdata = data, type = 'response'), 0.001), .999)
  Y_init_M1 = pmin(pmax(predict(Yfit, newdata = df_YM1, type = 'response'), 0.001), .999)
  Y_init_M0 = pmin(pmax(predict(Yfit, newdata = df_YM0, type = 'response'), 0.001), .999)
  return(list(initdata = list(M_ps = M_ps, ZA1_ps = ZA1_ps, ZA0_ps = ZA0_ps,
                              Z_psW = Z_psW, 
                              A_ps = A_ps), 
              Y_preds = list(Y_init = Y_init, 
                             Y_init_M1 = Y_init_M1, 
                             Y_init_M0 = Y_init_M0)))
  
}
}


#' @export
mediation.step1_glm_eff_seqT = function(initdata, Y_preds, data, gstarM_astar, a, transport) {
  if (transport) {
  if (a==1) ZS0_ps = initdata$ZA1S0_ps else ZS0_ps = initdata$ZA0S0_ps
  H = with(data, with(initdata, ((S == 1)*
                                   ((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))*
                                   ((Z == 1)*ZS0_ps + (Z == 0)*(1 - ZS0_ps))*(1 - S_ps))/
                        (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                           ((Z == 1)*Z_psWS + (Z == 0)*(1 - Z_psWS))*S_ps*PS0)))
  
  # updates
  Qfit = try(glm(data$Y ~ -1 + H + offset(qlogis(Y_preds$Y_init)), 
                 family = binomial), silent = TRUE)
  
  if (class(Qfit)[1]=="try-error") eps = 0 else eps = Qfit$coefficients
  
  return(list(Qstar_M  = plogis(qlogis(Y_preds$Y_init) + H*eps),
              Qstar_M1 = plogis(qlogis(Y_preds$Y_init_M1) + H*eps),
              Qstar_M0 = plogis(qlogis(Y_preds$Y_init_M0) + H*eps),
              Hm = H,
              eps = eps))
  } else {
  
    if (a==1) Z_ps = initdata$ZA1_ps else Z_ps = initdata$ZA0_ps
    H = with(data, with(initdata, (((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))*
                                     ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps)))/
                          (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                             ((Z == 1)*Z_psW + (Z == 0)*(1 - Z_psW)))))
    
    # updates
    Qfit = try(glm(data$Y ~ -1 + H + offset(qlogis(Y_preds$Y_init)), family = binomial,
                   weights = H), silent = TRUE)
    
    if (class(Qfit)[1]=="try-error") eps = 0 else eps = Qfit$coefficients
    
    return(list(Qstar_M  = plogis(qlogis(Y_preds$Y_init) + H*eps),
                Qstar_M1 = plogis(qlogis(Y_preds$Y_init_M1) + H*eps),
                Qstar_M0 = plogis(qlogis(Y_preds$Y_init_M0) + H*eps),
                Hm = H,
                eps = eps))  
  }
}


#' @export
mediation.step2_glm_eff_seqT = function(data, Qstar_M, Qstar_Mg, Hm, A_ps, a, tmle = TRUE,
                               EE = FALSE, bootstrap = FALSE, form, transport) {
  
  if (transport) PS0 = mean(data$S==0)
  df_QZ = data
  df_QZ$Qstar_Mg = Qstar_Mg
  
  # QZform = formula(paste0("Qstar_Mg ~ ", paste(covariates$covariates_QZ, collapse = "+")))
  QZform = form
  QZfit = glm(formula = QZform, data = df_QZ[df_QZ$A==a, ], family = binomial())
  # compute the clever covariate and update if tmle
  if (tmle) {
    if (transport) Hz = with(data, (A == a)*(S == 0)/(A*A_ps + (1 - A)*(1 - A_ps))/PS0) else {
      Hz = with(data, (A == a)/(A*A_ps + (1 - A)*(1 - A_ps)))
    }  
    QZ_preds_a = pmin(pmax(predict(QZfit, newdata = data, type = 'response'), .001), .999)
    # update
    QZfit_tmle = try(glm(Qstar_Mg ~ 1 + offset(qlogis(QZ_preds_a)), family = binomial,
                         weights = Hz), silent = TRUE)
    if (class(QZfit_tmle)[1]=="try-error") eps2 = 0 else eps2 = QZfit_tmle$coefficients
    
    QZstar_a = plogis(qlogis(QZ_preds_a) + eps2)
    if (transport) est = mean(QZstar_a[data$S==0]) else {
      est = mean(QZstar_a)
    }
    if(!bootstrap) { 
      D_Y = with(data, Hm*(Y - Qstar_M))
      D_Z = Hz*(Qstar_Mg - QZstar_a)
      if (transport) D_W = with(data, (QZstar_a - est)*(S ==0)/PS0) else {
        D_W = with(data, (QZstar_a - est))
      }
      D = D_Y + D_Z + D_W
    }
  } 
  
  # regress if EE or mle, EE updates the estimate, mle does not
  if (EE) {
    QZstar_a = pmin(pmax(predict(QZfit, newdata = data, type = 'response'), .001), .999) 
    if (transport) init_est = mean(QZstar_a[data$S==0]) else {
      init_est = mean(QZstar_a)
    }
    D_Y1s = with(data, Hm*(Y - Qstar_M))
    if (transport) Hz = with(data, (A == a)*(S == 0)/(A*A_ps + (1 - A)*(1 - A_ps))/PS0) else {
      Hz = with(data, (A == a)/(A*A_ps + (1 - A)*(1 - A_ps)))
    }
    D_Z1s = Hz*(Qstar_Mg - QZstar_a)
    if (transport) D_W1s = with(data, (QZstar_a - init_est)*(S ==0)/PS0) else {
      D_W1s = with(data, (QZstar_a - init_est))
    }
    D = D_Y1s + D_Z1s + D_W1s
    # update the estimate
    est = init_est + mean(D)
  }
  
  if (!bootstrap) {
    if (tmle) return(list(est = est, IC = D, eps2 = eps2)) else {
      return(list(est = est, IC = D, init_est = init_est))
    }
  } else {
    if (tmle) return(list(est = est)) else {
      return(list(est = est, init_est = init_est))
    }
  }
}

