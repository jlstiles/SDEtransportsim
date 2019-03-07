#' @export
mediation.steps_glm_eff1 = function(initdata, transport, bootstrap, max_iter) {
  tmledata = initdata
  calculateHQ = function(tmledata, setM, setZ, setS) {
    if (is.null(setM)) setM = tmledata$M
    if (is.null(setZ)) setZ = tmledata$Z
    if (is.null(setS)) setS = tmledata$S
    H_SDE = with(tmledata, ((setS == 1)*
                              ((setM == 1)*gstarM_astar0 + (setM == 0)*(1 - gstarM_astar0))*
                              ((ZA1S0_ps*setZ + (1 - setZ)*(1 - ZA1S0_ps))-
                                 (ZA0S0_ps*setZ + (1 - setZ)*(1 - ZA0S0_ps)))*(1 - S_ps))/
                   (((setM == 1)*M_ps + (setM == 0)*(1 - M_ps))*
                      ((setZ == 1)*Z_Wps + (setZ == 0)*(1 - Z_Wps))*S_ps*PS0))
    
    H_SIE = with(tmledata, ((S == 1)*
                              (((setM == 1)*gstarM_astar1 + (setM == 0)*(1 - gstarM_astar1))-
                                 ((setM == 1)*gstarM_astar0 + (setM == 0)*(1 - gstarM_astar0)))*
                              (ZA1S0_ps*setZ+ (1 - setZ)*(1 - ZA1S0_ps))*(1 - S_ps))/
                   (((setM == 1)*M_ps + (setM == 0)*(1 - M_ps))*
                      ((setZ == 1)*Z_Wps + (setZ == 0)*(1 - Z_Wps))*S_ps*PS0))
    
    
    return(matrix(c(H_SDE,H_SIE), ncol = 2))
  }
  
  calculateHZ = function(tmledata, setA, setS) {
    if (is.null(setA)) setA = tmledata$A
    if (is.null(setS)) setS = tmledata$S
    HZ_SDE = with(tmledata, (setS==0)*((setA == 1)/(A_ps*setA + (1 - setA)*(1 - A_ps))/PS0-
                                         (setA == 0)/(A_ps*setA + (1 - setA)*(1 - A_ps))/PS0)*(Q1Wa0 - Q0Wa0))
    
    HZ_SIE = with(tmledata, (setS==0)*(setA == 1)/A_ps/PS0*(Q1Wa1 - Q0Wa1 - Q1Wa0 + Q0Wa0))
    
    return(matrix(c(HZ_SDE, HZ_SIE), ncol = 2))
  }
  
  n = length(tmledata$Q)
  for (i in 1:max_iter){
    its = i
    HQ = calculateHQ(tmledata, setM = NULL, setZ = NULL, setS = NULL)[,1]
    HZ = calculateHZ(tmledata, setA = NULL, setS = NULL)[,1]
    
    Qa1WS = with(tmledata, Q1Wa0*ZA1_ps + Q0Wa0*(1-ZA1_ps))
    Qa0WS = with(tmledata, Q1Wa0*ZA0_ps + Q0Wa0*(1-ZA0_ps))
    
    Psi_SDE =  with(tmledata, mean((Qa1WS - Qa0WS)[S==0]))
    
    Qastar1WS = with(tmledata, Q1Wa1*ZA1_ps + Q0Wa1*(1-ZA1_ps))
    # Qastar0WS = with(tmledata, Q1Wa0*ZA1_ps + Q0Wa0*(1-ZA1_ps))
    
    Psi_SIE =  with(tmledata, mean((Qastar1WS - Qa1WS)[S==0]))
    
    DY = with(tmledata, HQ*as.vector(Y - Q))
    DZ = with(tmledata, HZ*as.vector(Z - Z_ps))
    DW = with(tmledata,
              matrix(c((S==0)/PS0*(Qa1WS - Qa0WS - Psi_SDE), 
                       (S==0)/PS0*(Qastar1WS - Qa1WS - Psi_SIE)), 
                     ncol = 2))[,1]
    D = DY + DZ + DW 
    
    # sigma = apply(D, 2, sd)/(n)
    sigma = sd(D)/n
    
    if (i == 1) {
      EE_info = list(
        Psi_SDE_EE = Psi_SDE + mean(D),
        # Psi_SIE_EE = Psi_SIE + mean(D[,2]),
        D = D,
        SE_EE = sigma*sqrt(n)
      )
    }
    ED = mean(D)
    tol = abs(ED) <= sigma
    if (tol) break
    
    Hdot = c(HQ, HZ)
    
    
    OS = with(tmledata, c(Q[S==1], Z_ps))
    OC = with(tmledata, c(Y[S==1], Z))
    cc = c(with(tmledata, Hdot[1:n][S==1]),Hdot[1:n+n])
    flucfit = try(suppressWarnings(glm(OC ~ -1 + cc + offset(qlogis(OS)), family = binomial)), 
                  silent = TRUE)
    
    if (class(flucfit)[1]=="try-error") eps = 0 else eps = flucfit$coefficients
    
    HM1Z1 = calculateHQ(tmledata, setM = 1, setZ = 1, setS = 1)[,1]
    HM1Z0 = calculateHQ(tmledata, setM = 1, setZ = 0, setS = 1)[,1]
    HM0Z1 = calculateHQ(tmledata, setM = 0, setZ = 1, setS = 1)[,1]
    HM0Z0 = calculateHQ(tmledata, setM = 0, setZ = 0, setS = 1)[,1]
    
    HA1S0 = calculateHZ(tmledata, setA = 1, setS = 0)[,1]
    HA0S0 = calculateHZ(tmledata, setA = 0, setS = 0)[,1]
    
    HA1 = calculateHZ(tmledata, setA = 1, setS = NULL)[,1]
    HA0 = calculateHZ(tmledata, setA = 0, setS = NULL)[,1]
    
    QandZ = with(tmledata, plogis(qlogis(c(Q,Z_ps)) + Hdot*eps))
    tmledata$Q = QandZ[1:n]
    tmledata$Z_ps = QandZ[1:n+n]
    
    tmledata$QM1Z1 = with(tmledata, plogis(qlogis(QM1Z1) + HM1Z1*eps))
    tmledata$QM1Z0 = with(tmledata, plogis(qlogis(QM1Z0) + HM1Z0*eps))
    tmledata$QM0Z1 = with(tmledata, plogis(qlogis(QM0Z1) + HM0Z1*eps))
    tmledata$QM0Z0 = with(tmledata, plogis(qlogis(QM0Z0) + HM0Z0*eps))
    
    
    tmledata$ZA1S0_ps = with(tmledata, plogis(qlogis(ZA1S0_ps) + HA1S0*eps))
    tmledata$ZA0S0_ps = with(tmledata, plogis(qlogis(ZA0S0_ps) + HA0S0*eps))
    
    tmledata$ZA1_ps = with(tmledata, plogis(qlogis(ZA1_ps) + HA1*eps))
    tmledata$ZA0_ps = with(tmledata, plogis(qlogis(ZA0_ps) + HA0*eps))
    tmledata$Z_Wps = with(tmledata, ZA1_ps*A_ps + ZA0_ps*(1 - A_ps))
    
    tmledata$Q1Wa1 = with(tmledata, QM1Z1*gstarM_astar1 + QM0Z1*(1-gstarM_astar1))
    tmledata$Q0Wa1 = with(tmledata, QM1Z0*gstarM_astar1 + QM0Z0*(1-gstarM_astar1))
    tmledata$Q1Wa0 = with(tmledata, QM1Z1*gstarM_astar0 + QM0Z1*(1-gstarM_astar0))
    tmledata$Q0Wa0 = with(tmledata, QM1Z0*gstarM_astar0 + QM0Z0*(1-gstarM_astar0))
  }
  if (bootstrap) {return(c(Psi_SDE = Psi_SDE,  
                           Psi_SDE_EE = EE_info$Psi_SDE, Psi_SIE_EE = EE_info$Psi_SIE, its = its))
  } else {
    return(list(tmledata = tmledata, initdata = initdata, Psi_SDE = Psi_SDE, 
                ED = ED, SE = sigma*sqrt(n), D = D, EE_info = EE_info, its = its))
  }
}

  
