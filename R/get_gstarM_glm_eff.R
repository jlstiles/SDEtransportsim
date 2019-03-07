#' @export
get_gstarM_glm_eff  = function(data, forms, transport,  pooled,  gstar_S, truth) 
{
  Mstarform = forms$Mstarform
  Zstarform = forms$Zstarform
  
  dataZ1 = dataZ0 = dataA1 = dataA0 = data
  # Here is where we need to define the parameter as to the mechanism used for M and Z,
  # Note, if !pooled then these lines do not matter
  if (pooled) {
    dataA1$S = dataA0$S = gstar_S
    dataZ1$S = dataZ0$S = gstar_S
  }
  
  dataZ1$Z = 1
  dataZ0$Z = 0
  dataA1$A = 1
  dataA0$A = 0
  
  if (transport & !pooled) {
    Mstarfit = glm(formula = Mstarform, data = data[data$S==gstar_S, ], family = "binomial")
    Zstarfit = glm(formula = Zstarform, data = data[data$S==gstar_S, ], family = "binomial")
  } else {
    Mstarfit = glm(formula = Mstarform, data = data, family = "binomial")
    Zstarfit = glm(formula = Zstarform, data = data, family = "binomial")
  }
  
  predMz1 = predict(Mstarfit, newdata = dataZ1, type = 'response')
  predMz0 = predict(Mstarfit, newdata = dataZ0, type = 'response')
  
  predZa0 = predict(Zstarfit, newdata = dataA0, type = 'response')
  predZa1 = predict(Zstarfit, newdata = dataA1, type = 'response')
  
  gstarM_astar0 = predMz1*predZa0 + predMz0*(1 - predZa0)
  gstarM_astar1 = predMz1*predZa1 + predMz0*(1 - predZa1)
  
  if (!is.null(truth)) {
    ftruth = function(f_W, f_S, setA, f_Z, f_Y, transport, pooled) {
      W = f_W(2e6)
      W = as.data.frame(W)
      
      if (transport) {
        P_SW = f_S(W)
        S = rbinom(2e6,1,P_SW) 
        PS0_0 = mean(S==0)
        A = setA[1]
        pzscores = f_Z(A=A,S=S,W=W)
        Z = rbinom(2e6, 1, pzscores)
        popZ1 = popZ0 = popAastar = data.frame(cbind(W, S=S,A=A,Z=Z))
      } else {
        A = setA[1]
        pzscores = f_Z(A=A,W=W)
        Z = rbinom(2e6, 1, pzscores)
        popZ1 = popZ0 = popAastar = data.frame(cbind(W,A=A,Z=Z))
      }
      
      # note, pooled is autoset to false if not transporting.  If pooled and transporting we
      # need to set S, otherwise it is autohandled by the fits being done on subset by gstar_S
      if (pooled) {
        popAastar$S = gstar_S
        popZ1$S = popZ0$S = gstar_S
      }
      
      popZ1$Z = 1
      popZ0$Z = 0
      popAastar$A = setA[2]
      
      # regardless of pooling or not the model.matrix will keep the vars according to the formula
      # so fine whether pooling or not
      pop_predZastar = predict(Zstarfit, newdata = popAastar, type = 'response')
      pop_predMz1 = predict(Mstarfit, newdata = popZ1, type = 'response')
      pop_predMz0 = predict(Mstarfit, newdata = popZ0, type = 'response')
      
      gstarM = pop_predMz1*pop_predZastar + pop_predMz0*(1 - pop_predZastar)
      
      M = rbinom(2e6, 1, gstarM)
      # make a Y model according to the restrictions
      Yscores = f_Y(M=M,Z=Z,W=W)
      Y = rbinom(2e6, 1, Yscores)
      
      if (transport) return(c(Psi = mean(Y[S==0]), PS0_0 = PS0_0)) else return(mean(Y))
    }
    
    
    info_astar0a1_0 = ftruth(truth$f_W, truth$f_S, setA = c(1,0), truth$f_Z, truth$f_Y, transport, pooled)
    info_astar0a0_0 = ftruth(truth$f_W, truth$f_S, setA = c(0,0), truth$f_Z, truth$f_Y, transport, pooled)
    info_astar1a1_0 = ftruth(truth$f_W, truth$f_S, setA = c(1,1), truth$f_Z, truth$f_Y, transport, pooled)
    Psi_astar0a1_0 = info_astar0a1_0[1]
    Psi_astar0a0_0 = info_astar0a0_0[1]
    Psi_astar1a1_0 = info_astar1a1_0[1]
    
    SDE_0 = Psi_astar0a1_0 - Psi_astar0a0_0
    SIE_0 = Psi_astar1a1_0 - Psi_astar0a1_0
    
    if (transport) PS0_0 = info_astar1a1_0[2]  
    
    # get the true IC's
    W = data[,grep("W", colnames(data))]
    if (transport) {
      S_ps0 = with(data, truth$f_S(W=W))
      M_ps0 = with(data, truth$f_M(Z=Z, W=W, S=1))
      Z_ps0 = with(data, truth$f_Z(A=A, W=W, S=S))
      A_ps0 = with(data, truth$f_A(W=W, S=S))
      Z_psWS0 = with(data, truth$f_Z(A=1, W=W, S=S)*A_ps0 + 
                       truth$f_Z(A=0, W=W, S=S)*(1-A_ps0))
    } else {
      M_ps0 = with(data, truth$f_M(Z=Z, W=W))
      Z_ps0 = with(data, truth$f_Z(A=A, W=W))
      A_ps0 = with(data, truth$f_A(W=W))
      Z_psW = with(data, truth$f_Z(A=1, W=W)*A_ps0 + 
                     truth$f_Z(A=0, W=W)*(1-A_ps0))
    }
    if (transport) {
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
    } else {
      get_cc_0 = function(data, gstarM_astar, a) {
        df_ZA0 = data
        df_ZA0$A = a
        ZA0_ps0 = with(df_ZA0, truth$f_Z(A=A, W=W))
        with(data, (((M == 1)*gstarM_astar0 + (M == 0)*(1 - gstarM_astar0))
                    *(ZA0_ps0*Z+ (1 - Z)*(1 - ZA0_ps0))/
                      (((M == 1)*M_ps0 + (M == 0)*(1 - M_ps0))*
                         (Z_psW*Z+ (1 - Z)*(1 - Z_psW)))))
      }
    }
    
    if (transport) {
      get_Hz_0 = function(a) with(data, (A == a)*(S == 0)/(A*A_ps0 + (1 - A)*(1 - A_ps0))/PS0_0)
    } else {
      get_Hz_0 = function(a) with(data, (A == a)/(A*A_ps0 + (1 - A)*(1 - A_ps0)))
    }
    
    Hm_astar0a1_0 = get_cc_0(data = data, gstarM_astar = gstarM_astar0, a = 1)
    Hm_astar0a0_0 = get_cc_0(data = data, gstarM_astar = gstarM_astar0, a = 0)
    Hm_astar1a1_0 = get_cc_0(data = data, gstarM_astar = gstarM_astar1, a = 1)
    
    Hz_astar0a1_0 = get_Hz_0(1)
    Hz_astar0a0_0 = get_Hz_0(0)
    Hz_astar1a1_0 = get_Hz_0(1)
    
    get_trueIC = function(gstar_M, a, Psi_0, Hm_0, Hz_0) {
      if (transport) p_0Z = with(data, truth$f_Z(A=a, W=W, S=S)) else {
        p_0Z = with(data, truth$f_Z(A=a, W=W))
      }
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
      if (transport) D_0W = with(data, (Q_ghat_astarW - Psi_0)*(S ==0)/PS0_0) else {
        D_0W = with(data, (Q_ghat_astarW - Psi_0))
      }
      D_0 = D_0Y + D_0Z + D_0W
      return(D_0)
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
    
    SE_SDE_0 = sd(D_astar0a1_0 - D_astar0a0_0)/sqrt(nrow(data))
    SE_SIE_0 = sd(D_astar1a1_0 - D_astar0a1_0)/sqrt(nrow(data))
    
    return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0,
                SDE_0 = SDE_0, SIE_0 = SIE_0, SE_SDE_0 = SE_SDE_0, SE_SIE_0 = SE_SIE_0))
  } else {
    return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0))
  }
}

