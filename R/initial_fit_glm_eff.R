#' @export
get.mediation.initdata_glm_eff = function(data, forms, RCT = 0.5, transport, 
                                           pooled, gstarM_astar) {
  
  # data = cbind(data$W,A=data$A, Z=data$Z, M=data$M, Y=data$Y)
  if (pooled | !transport) {
    Zform = forms$Zstarform
    Mform = forms$Mstarform
  } else {
    Zform = forms$Zform
    Mform = forms$Mform
  }
  if (transport) { 
    # This is to avoid NA's being a pain, we won't use these outcomes
    data$Y[data$S!=1] = 2
    
    df_ZA1S0 = df_ZA0S0 = data
    df_ZA1S0$S = df_ZA0S0$S = 0
    df_ZA0S0$A = 0
    df_ZA1S0$A = 1
    df_ZA1S1 = df_ZA1S0
    df_ZA0S1 = df_ZA0S0
    df_ZA0S1$S = 1
    df_ZA1S1$S = 1
    df_S1 = data
    df_S1$S = 1
    Sform = forms$Sform
  }
  
  df_YM1 = data
  df_YM1$M = 1
  
  df_YM0 = df_YM1
  df_YM0$M = 0
  
  df_YM0Z0 = df_YM0Z1 = df_YM0
  df_YM0Z0$Z = 0 
  df_YM0Z1$Z = 1
  df_YM1Z0 = df_YM1Z1 = df_YM1
  df_YM1Z0$Z = 0 
  df_YM1Z1$Z = 1
  
  df_ZA0 = df_ZA1 = data
  df_ZA0$A = 0
  df_ZA1$A = 1
  
  Yform = forms$Yform
  
  if (transport) {
    Yfit = glm(formula = Yform, data = data[data$S==1,], family = "binomial")
  } else {
    Yfit = glm(formula = Yform, data = data, family = "binomial")
  }
  
  Mfit = glm(formula = Mform, data = data, family = "binomial")
  
  if (!transport) {
    Zfit = glm(formula = Zform, data = data, family = "binomial")
    ZA1_ps = predict(Zfit, newdata = df_ZA1, type = 'response')
    ZA0_ps = predict(Zfit, newdata = df_ZA0, type = 'response')
  }
  if (transport) {
    Zfit = glm(formula = Zform, data = data, family = "binomial")
    Z_ps = predict(Zfit, type = 'response')
    ZA1_ps = predict(Zfit, newdata = df_ZA1, type = 'response')
    ZA0_ps = predict(Zfit, newdata = df_ZA0, type = 'response')
    ZA1S0_ps = predict(Zfit, newdata = df_ZA1S0, type = 'response')
    ZA0S0_ps = predict(Zfit, newdata = df_ZA0S0, type = 'response')
    ZA1S1_ps = predict(Zfit, newdata = df_ZA1S1, type = 'response')
    ZA0S1_ps = predict(Zfit, newdata = df_ZA0S1, type = 'response')
    MS1_ps = predict(Mfit, newdata = df_S1, type = 'response')
    Sfit = glm(formula = Sform, data = data, family = "binomial")
    S_ps = predict(Sfit, type = 'response')
    PS0 = mean(data$S==0)
  } else {
    Zfit = glm(formula = Zform, data = data, family = "binomial")
    Z_ps = predict(Zfit, type = 'response')
    ZA1_ps = predict(Zfit, newdata = df_ZA1, type = 'response')
    ZA0_ps = predict(Zfit, newdata = df_ZA0, type = 'response')
  }
  # propensity scores
  if (is.null(RCT)) { 
    Aform = forms$Aform
    Afit = glm(formula = Aform, data = data, family = "binomial")
    A_ps = predict(Afit, type = 'response')
    AS1_ps = predict(Afit, df_S1, type = 'response')
  } else A_ps = AS1_ps = RCT
  # as clev cov is 0 otherwise 
  M_ps = predict(Mfit, type = 'response')
  
  # Predict Y for whole data, also with M = 1 and 0
  Q = pmin(pmax(predict(Yfit, newdata = data, type = 'response'), 0.001), .999)
  QM1Z0 = pmin(pmax(predict(Yfit, newdata = df_YM1Z0, type = 'response'), 0.001), .999)
  QM1Z1 = pmin(pmax(predict(Yfit, newdata = df_YM1Z1, type = 'response'), 0.001), .999)
  QM0Z0 = pmin(pmax(predict(Yfit, newdata = df_YM0Z0, type = 'response'), 0.001), .999)
  QM0Z1 = pmin(pmax(predict(Yfit, newdata = df_YM0Z1, type = 'response'), 0.001), .999)
  
  gstarM_astar1 = gstarM_astar$gstarM_astar1
  gstarM_astar0 = gstarM_astar$gstarM_astar0
  
  Z_Wps = ZA1_ps*A_ps + ZA0_ps*(1 - A_ps)
  Q1Wa1 = QM1Z1*gstarM_astar1 + QM0Z1*(1-gstarM_astar1)
  Q0Wa1 = QM1Z0*gstarM_astar1 + QM0Z0*(1-gstarM_astar1)
  Q1Wa0 = QM1Z1*gstarM_astar0 + QM0Z1*(1-gstarM_astar0)
  Q0Wa0 = QM1Z0*gstarM_astar0 + QM0Z0*(1-gstarM_astar0)
  
  if (transport) {
    return(list(MS1_ps = MS1_ps, ZA1S0_ps = ZA1S0_ps, ZA0S0_ps = ZA0S0_ps, 
                ZA1S1_ps = ZA1S1_ps, ZA0S1_ps = ZA0S1_ps,
                ZA1_ps = ZA1_ps, ZA0_ps = ZA0_ps, A_ps = A_ps, AS1_ps = AS1_ps, 
                S_ps = S_ps, PS0 = PS0, Z_ps = Z_ps, Q = Q, 
                gstarM_astar1 = gstarM_astar1,
                gstarM_astar0 = gstarM_astar0, 
                Z_Wps = Z_Wps,
                QM1Z1 = QM1Z1, 
                QM1Z0 = QM1Z0,
                QM0Z1 = QM0Z1,
                QM0Z0 = QM0Z0,
                Q1Wa1 = Q1Wa1, 
                Q0Wa1 = Q0Wa1,
                Q1Wa0 = Q1Wa0,
                Q0Wa0 = Q0Wa0,
                Y = data$Y,
                M = data$M,
                Z = data$Z,
                A = data$A,
                S = data$S))
  } else {
    return(list(M_ps = M_ps, ZA1_ps = ZA1_ps, ZA0_ps = ZA0_ps, 
                A_ps = A_ps, Z_ps = Z_ps,Q = Q, Z_Wps = Z_Wps,
                gstarM_astar1 = gstarM_astar1,
                gstarM_astar0 = gstarM_astar0,
                QM1Z1 = QM1Z1, 
                QM1Z0 = QM1Z0,
                QM0Z1 = QM0Z1,
                QM0Z0 = QM0Z0,
                Q1Wa1 = Q1Wa1, 
                Q0Wa1 = Q0Wa1,
                Q1Wa0 = Q1Wa0,
                Q0Wa0 = Q0Wa0,
                Y = data$Y,
                M = data$M,
                Z = data$Z,
                A = data$A))  
  }
  
}

