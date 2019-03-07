#' @export
mediation.steps_glm = function(initdata, transport, bootstrap, max_iter) {
  tmledata = initdata
  calculateHQ = function(tmledata, setM, setZ, setA, setS) {
    if (is.null(setM)) setM = tmledata$M
    if (is.null(setZ)) setZ = tmledata$Z
    if (is.null(setA)) setA = tmledata$A

    if (transport) {
      if (is.null(setS)) setS = tmledata$S
      H_SDE = with(tmledata, ((setS == 1)*((setA ==1)*(setZ*ZA1S0_ps + (1 - setZ)*(1 - ZA1S0_ps))-
                                             (setA==0)*(setZ*ZA0S0_ps + (1 - setZ)*(1 - ZA0S0_ps)))*
                                ((setM == 1)*gstarM_astar0 + (setM == 0)*(1 - gstarM_astar0))*
                                (1 - S_ps))/((AS1_ps*setA + (1 - setA)*(1 - AS1_ps))*
                                               ((setM == 1)*MS1_ps + (setM == 0)*(1 - MS1_ps))*
                                               ((setZ == 1)*ZS1_ps + (setZ == 0)*(1 - ZS1_ps))*S_ps*PS0))
    } else {
      H_SDE = with(tmledata, (((setA ==1)-(setA==0))/(A_ps*setA + (1 - setA)*(1 - A_ps))*
                                ((setM == 1)*gstarM_astar0 + (setM == 0)*(1 - gstarM_astar0))/
                                ((setM == 1)*M_ps + (setM == 0)*(1 - M_ps))))
    }
    
    if (transport) {
      H_SIE = with(tmledata, ((setS == 1)*(setA==1)*
                                (((setM == 1)*gstarM_astar1 + (setM == 0)*(1 - gstarM_astar1))-
                                   ((setM == 1)*gstarM_astar0 + (setM == 0)*(1 - gstarM_astar0)))*
                                (setZ*ZA1S0_ps + (1 - setZ)*(1 - ZA1S0_ps))*
                                (1 - S_ps))/
                     (((setM == 1)*MS1_ps + (setM == 0)*(1 - MS1_ps))*
                        ((setZ == 1)*ZS1_ps + (setZ == 0)*(1 - ZS1_ps))*
                        (AS1_ps*setA + (1 - setA)*(1 - AS1_ps))*S_ps*PS0))
    } else {
      H_SIE = with(tmledata, ((setA ==1)*
                                (((setM == 1)*gstarM_astar1 + (setM == 0)*(1 - gstarM_astar1))-
                                   ((setM == 1)*gstarM_astar0 + (setM == 0)*(1 - gstarM_astar0)))/
                                (((setM == 1)*M_ps + (setM == 0)*(1 - M_ps))*
                                   (A_ps*setA + (1 - setA)*(1 - A_ps)))))
    }
    
    return(matrix(c(H_SDE,H_SIE), ncol = 2))
  }
  
  calculateHZ = function(tmledata, setA, setS) {
    if (is.null(setA)) setA = tmledata$A
      if (transport) { 
        if (is.null(setS)) setS = tmledata$S
        HZ_SDE = with(tmledata, (setS==0)*((setA == 1)/A_ps/PS0-
                                             (setA == 0)/(1 - A_ps)/PS0)*(Q1Wa0 - Q0Wa0))
        
      } else {
        HZ_SDE = with(tmledata, ((setA == 1)/A_ps-
                                   (setA == 0)/(1 - A_ps))*(Q1Wa0 - Q0Wa0))
      }
      if (transport) {  
        HZ_SIE = with(tmledata, (setS==0)*(setA == 1)/A_ps/PS0*(Q1Wa1 - Q0Wa1 - Q1Wa0 + Q0Wa0))
        
      } else {
        HZ_SIE = with(tmledata, (setA == 1)/A_ps*(Q1Wa1 - Q0Wa1 - Q1Wa0 + Q0Wa0))
      }
    return(matrix(c(HZ_SDE, HZ_SIE), ncol = 2))
  }
  
  n = length(tmledata$Q)
  for (i in 1:max_iter){
    its = i
    HQ = calculateHQ(tmledata, setM = NULL, setZ = NULL, setA = NULL, setS = NULL)
    HZ = calculateHZ(tmledata, setA = NULL, setS = NULL)
  
    clfm = function(cc, dir) {
      H = cc %*% dir
    }
    
    if (transport) {
      Qa1WS = with(tmledata, Q1Wa0*ZA1S0_ps + Q0Wa0*(1-ZA1S0_ps))
      Qa0WS = with(tmledata, Q1Wa0*ZA0S0_ps + Q0Wa0*(1-ZA0S0_ps))
      
      Psi_SDE =  with(tmledata, mean((Qa1WS - Qa0WS)[S==0]))
      
      Qastar1WS = with(tmledata, Q1Wa1*ZA1S0_ps + Q0Wa1*(1-ZA1S0_ps))
      # Qastar0WS = with(tmledata, Q1Wa0*ZA1_ps + Q0Wa0*(1-ZA1_ps))
      
      Psi_SIE =  with(tmledata, mean((Qastar1WS - Qa1WS)[S==0]))
      
      DY = with(tmledata, HQ*as.vector(Y - Q))
      DZ = with(tmledata, HZ*as.vector(Z - Z_ps))
      DW = with(tmledata,
                matrix(c((S==0)/PS0*(Qa1WS - Qa0WS - Psi_SDE), 
                         (S==0)/PS0*(Qastar1WS - Qa1WS - Psi_SIE)), 
                       ncol = 2))
      D = DY + DZ + DW 
 
      sigma = apply(D, 2, sd)/n 
      
      if (i == 1) {
        EE_info = list(
          Psi_SDE_EE = Psi_SDE + mean(D[,1]),
          Psi_SIE_EE = Psi_SIE + mean(D[,2]),
          D = D,
          SE_EE = sigma*sqrt(n)
          )
      }
      ED = colMeans(D)
      tol = all(abs(ED) <= sigma)
      if (tol) break
      
      dir = ED/sqrt(ED[1]^2 + ED[2]^2)
      Hdot = clfm(rbind(HQ,HZ), dir)
        
    } else {
      Qa1W = with(tmledata, Q1Wa0*ZA1_ps + Q0Wa0*(1-ZA1_ps))
      Qa0W = with(tmledata, Q1Wa0*ZA0_ps + Q0Wa0*(1-ZA0_ps))
      
      Psi_SDE =  mean(Qa1W - Qa0W)
      
      Qastar1W = with(tmledata, Q1Wa1*ZA1_ps + Q0Wa1*(1-ZA1_ps))
      # Qastar0W = with(tmledata, Q1Wa0*ZA1_ps + Q0Wa0*(1-ZA1_ps))
      
      Psi_SIE =  with(tmledata, mean(Qastar1W - Qa1W))
      
      DY = with(tmledata, HQ*as.vector(Y - Q))
      DZ = with(tmledata, HZ*as.vector(Z - Z_ps))
      DW = with(tmledata,
                matrix(c((Qa1W - Qa0W - Psi_SDE),(Qastar1W - Qa1W - Psi_SIE)), 
                       ncol = 2))
      D = DY + DZ + DW 
      sigma = apply(D, 2, sd)/(n) 
      if (i == 1) {
        EE_info = list(
          Psi_SDE_EE = Psi_SDE + mean(D[,1]),
          Psi_SIE_EE = Psi_SIE + mean(D[,2]),
          D = D,
          SE_EE = sigma*sqrt(n)
        )
      }
      
      
      ED = colMeans(D)
      tol = all(abs(ED) <= sigma)
      if (tol) break
      
      dir = ED/sqrt(ED[1]^2 + ED[2]^2)
      Hdot = clfm(rbind(HQ,HZ), dir)
    }
    
    OS = with(tmledata, c(Q[S==1], Z_ps))
    OC = with(tmledata, c(Y[S==1], Z))
    cc = c(with(tmledata, c(Hdot[1:n][S==1])),Hdot[1:n+n])
    flucfit = try(suppressWarnings(glm(OC ~ -1 + cc + offset(qlogis(OS)), family = binomial)), 
                  silent = TRUE)
    
    if (class(flucfit)[1]=="try-error") eps = 0 else eps = flucfit$coefficients
    
    HM1Z1 = clfm(calculateHQ(tmledata, setM = 1, setZ = 1, setA = NULL, setS = 1),dir)
    HM1Z0 = clfm(calculateHQ(tmledata, setM = 1, setZ = 0, setA = NULL, setS = 1),dir)
    HM0Z1 = clfm(calculateHQ(tmledata, setM = 0, setZ = 1, setA = NULL, setS = 1),dir)
    HM0Z0 = clfm(calculateHQ(tmledata, setM = 0, setZ = 0, setA = NULL, setS = 1),dir)

    tmledata$QM1Z1 = with(tmledata, plogis(qlogis(QM1Z1) + HM1Z1*eps))
    tmledata$QM1Z0 = with(tmledata, plogis(qlogis(QM1Z0) + HM1Z0*eps))
    tmledata$QM0Z1 = with(tmledata, plogis(qlogis(QM0Z1) + HM0Z1*eps))
    tmledata$QM0Z0 = with(tmledata, plogis(qlogis(QM0Z0) + HM0Z0*eps))
    
    if (transport) {
      HA1S0 = clfm(calculateHZ(tmledata, setA = 1, setS = 0),dir)
      tmledata$ZA1S0_ps = with(tmledata, plogis(qlogis(ZA1S0_ps) + HA1S0*eps))
      HA0S0 = clfm(calculateHZ(tmledata, setA = 0, setS = 0),dir)
      tmledata$ZA0S0_ps = with(tmledata, plogis(qlogis(ZA0S0_ps) + HA0S0*eps))
    } else {
      HA1 = clfm(calculateHZ(tmledata, setA = 1, setS = NULL),dir)
      HA0 = clfm(calculateHZ(tmledata, setA = 0, setS = NULL),dir)
      tmledata$ZA1_ps = with(tmledata, plogis(qlogis(ZA1_ps) + HA1*eps))
      tmledata$ZA0_ps = with(tmledata, plogis(qlogis(ZA0_ps) + HA0*eps))
    }
    
    tmledata$Q1Wa1 = with(tmledata, QM1Z1*gstarM_astar1 + QM0Z1*(1-gstarM_astar1))
    tmledata$Q0Wa1 = with(tmledata, QM1Z0*gstarM_astar1 + QM0Z0*(1-gstarM_astar1))
    tmledata$Q1Wa0 = with(tmledata, QM1Z1*gstarM_astar0 + QM0Z1*(1-gstarM_astar0))
    tmledata$Q0Wa0 = with(tmledata, QM1Z0*gstarM_astar0 + QM0Z0*(1-gstarM_astar0))
    
    QandZ = with(tmledata, plogis(qlogis(c(Q,Z_ps)) + Hdot*eps))
    tmledata$Q = QandZ[1:n]
    tmledata$Z_ps = QandZ[1:n+n]
  }
  if (bootstrap) {return(c(Psi_SDE = Psi_SDE, Psi_SIE = Psi_SIE, 
                         Psi_SDE_EE = EE_info$Psi_SDE, Psi_SIE_EE = EE_info$Psi_SIE, its = its))
  } else {
    return(list(tmledata = tmledata, initdata = initdata, Psi_SDE = Psi_SDE, Psi_SIE = Psi_SIE,
                ED = ED, SE = sigma*sqrt(n), D = D, EE_info = EE_info,its = its))
  }
}
  
