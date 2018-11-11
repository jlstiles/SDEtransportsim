#' @title SDE_tmle_lasso
#' @description computes the sequential regression, targeted maximum likelihood estimate
#' for the stochastic direct effect or stochastic indirect effect using lasso. Note, this is
#' a non-transport parameter.
#' @param data, data.frame where confounders have any names but the following, which must be 
#' reserved as follows: A = treatment, Z = intermediate confounder, M = mediator and Y is the outcome.
#' @param forms, list of formulas. Include for each necessary model for outcome, 
#' called Yform for outcome Y, QZform for outcome Qstar_Mg, Mform, Zform, Aform (can be NULL if RCT)
#' is selected as TRUE. 
#' @param RCT either NULL or a value, if null, then the Aform is used to fit the propensity score, 
#' otherwise propensity scores are set to RCT.
#' @param B, the number of bootstraps, default is NULL and should not be changed since this will
#' be invalid for use with lasso anyway.
#' @return  a list with a CI's for SDE and SIE for the means under (a*,a) combos (0,0), (0,1), (1,1) 
#' and the epsilons for both sequential regressions for those three parameters
#' @example /inst/example_SDE_lassoNOtransport.R 
#' @export
SDE_tmle_lasso = function(data, forms, RCT = 0.5, B = NULL, truth = NULL) 
{
  # get the stochastic dist of M and true params if you want 
  gstar_info = get_gstarM_lasso(data = data,forms = forms)
  gstarM_astar1 = gstar_info$gstarM_astar1
  gstarM_astar0 = gstar_info$gstarM_astar0
  gstarM_astar = list(gstarM_astar0 = gstarM_astar0, gstarM_astar1 = gstarM_astar1)
  
  # perform initial fits for the first regression
  init_info = get.mediation.initdata_lasso(data = data, forms = forms, RCT = RCT)
  Y_preds = init_info$Y_preds
  
  est_info = lapply(0:1, FUN = function(astar) {
    lapply(0:1, FUN = function(a) {
      # get tmle info
      # get iptw here while I'm at it
      update = mediation.step1_lasso(initdata = init_info$initdata, init_info$Y_preds, data = data, 
                               gstarM_astar[[astar+1]], a, iptw = FALSE)
      
      Y_Mg = get.stochasticM(gstarM_astar[[astar+1]], Y_preds[[2]], Y_preds[[3]]) 
      A_ps = init_info$initdata$A_ps

      Qstar_Mg = get.stochasticM(gstarM_astar[[astar+1]], update$Qstar_M1, update$Qstar_M0) 
      # compute Qstar_Mg here
      tmle_info = mediation.step2_lasso(data = data, Qstar_M = update$Qstar_M, 
                      Qstar_Mg = Qstar_Mg, Hm = update$Hm, A_ps = A_ps, 
                      a = a, tmle = TRUE,
                      EE = FALSE, bootstrap = FALSE, form = forms$QZform)
      
      # compile all estimates
      tmle_info$eps1 = update$eps
      return(tmle_info)
    }) 
  })
  
  # run the bootstrap 500 times
  # bootstrap from the data and thus the gstarM_astar1 and gstarM_astar0
  n = nrow(data)
  if (!is.null(B)) {
    boot_ests = lapply(1:B, FUN = function(x) {
    inds = sample(1:n, replace = TRUE)
    data = data[inds, ]
    init_info = get.mediation.initdata_lasso(data = data, forms = forms, RCT = RCT)
    Y_preds = init_info$Y_preds
    gstarM_astar = list(gstarM_astar0[inds], gstarM_astar1[inds])
    
    est_info = lapply(0:1, FUN = function(astar) {
      return(lapply(0:1, FUN = function(a) {
        update = mediation.step1_lasso(initdata = init_info$initdata, Y_preds = Y_preds, data = data, 
                                 gstarM_astar[[astar+1]], a, iptw = FALSE)
        
        Y_Mg = get.stochasticM(gstarM_astar[[astar+1]], Y_preds[[2]], Y_preds[[3]]) 
        A_ps = init_info$initdata$A_ps
        
        Qstar_Mg = get.stochasticM(gstarM_astar[[astar+1]], update$Qstar_M1, update$Qstar_M0) 
        # compute Qstar_Mg here
        tmle_info = mediation.step2_lasso(data = data, Qstar_M = update$Qstar_M, 
                                    Qstar_Mg = Qstar_Mg, Hm = update$Hm, A_ps = A_ps, 
                                    a = a, tmle = TRUE,
                                    EE = FALSE, bootstrap = TRUE, form = forms$QZform)
        # compile all estimates
        return(tmle_info)
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
  
  D_SDE = est_info[[1]][[2]]$IC - est_info[[1]][[1]]$IC
  D_SIE = est_info[[2]][[2]]$IC- est_info[[1]][[2]]$IC

  SE_SDE = sd(D_SDE)/sqrt(n)
  SE_SIE = sd(D_SIE)/sqrt(n)

  if (!is.null(truth)) { 
  D_SDE_0 = gstar_info$D_astar0a1_0 - gstar_info$D_astar0a0_0
  D_SIE_0 = gstar_info$D_astar1a1_0 - gstar_info$D_astar0a1_0
  
  SE_SDE_0 = sd(D_SDE_0)/sqrt(n)
  SE_SIE_0 = sd(D_SIE_0)/sqrt(n)
  
  Psi_astar0a1_0 = gstar_info$Psi_astar0a1_0
  Psi_astar0a0_0 = gstar_info$Psi_astar0a0_0
  Psi_astar1a1_0 = gstar_info$Psi_astar1a1_0
  
  SDE_0 = Psi_astar0a1_0 - Psi_astar0a0_0
  SIE_0 = Psi_astar1a1_0 - Psi_astar0a1_0 
  
  } else {
    D_SDE_0 = D_SIE_0 = SE_SDE_0 = SE_SIE_0 = Psi_astar0a1_0 = 
      Psi_astar0a0_0 = Psi_astar1a1_0 = SDE_0 = SIE_0 = NULL
  }
  ests_astar0a1 =   c(tmle = est_info[[1]][[2]]$est,
                      Psi_0 = Psi_astar0a1_0)
  
  ests_astar0a0 =   c(tmle = est_info[[1]][[1]]$est,
                      Psi_0 = Psi_astar0a0_0)
  
  ests_astar1a1 =   c(tmle = est_info[[2]][[2]]$est,
                      Psi_0 = Psi_astar1a1_0)
  
  SDE_ests = ests_astar0a1 - ests_astar0a0
  SIE_ests = ests_astar1a1 - ests_astar0a1
  
  CI_SDE = c(SDE_ests[1], SDE_ests[1] - 1.96*SE_SDE, SDE_ests[1] + 1.96*SE_SDE)
  CI_SIE = c(SIE_ests[1], SIE_ests[1] - 1.96*SE_SIE, SIE_ests[1] + 1.96*SE_SIE)
  
  if (!is.null(B)) {
    CI_SDE_boot = c(SDE_ests[1], SDE_ests[1] - 1.96*bootSE_SDE[1], SDE_ests[1] + 1.96*bootSE_SDE[1])
    CI_SIE_boot = c(SIE_ests[1], SIE_ests[1] - 1.96*bootSE_SIE[1], SIE_ests[1] + 1.96*bootSE_SIE[1])

    return(list(CI_SDE = CI_SDE, CI_SIE = CI_SIE,
                CI_SDE_boot = CI_SDE_boot,
                CI_SIE_boot = CI_SIE_boot,
                SE_SDE = SE_SDE, 
                SE_SIE = SE_SIE, 
                ests_astar0a1 = ests_astar0a1,
                ests_astar0a0 = ests_astar0a0,
                ests_astar1a1 = ests_astar1a1,
                
                eps1_astar0a1 = est_info[[1]][[2]]$eps1,
                eps2_astar0a1 = est_info[[1]][[2]]$eps2,
                eps1_astar0a0 = est_info[[1]][[1]]$eps1,
                eps2_astar0a0 = est_info[[1]][[1]]$eps2,
                eps1_astar1a1 = est_info[[2]][[2]]$eps1,
                eps2_astar1a1 = est_info[[2]][[2]]$eps2 
    )) 
  } else {
    
    return(list(CI_SDE = CI_SDE, 
                CI_SIE = CI_SIE, 
                SE_SDE = SE_SDE, 
                SE_SIE = SE_SIE, 
                ests_astar0a1 = ests_astar0a1,
                ests_astar0a0 = ests_astar0a0,
                ests_astar1a1 = ests_astar1a1,
                
                eps1_astar0a1 = est_info[[1]][[2]]$eps1,
                eps2_astar0a1 = est_info[[1]][[2]]$eps2,
                eps1_astar0a0 = est_info[[1]][[1]]$eps1,
                eps2_astar0a0 = est_info[[1]][[1]]$eps2,
                eps1_astar1a1 = est_info[[2]][[2]]$eps1,
                eps2_astar1a1 = est_info[[2]][[2]]$eps2 
    ))
    
  }
} 

#' @export
get_gstarM_lasso  = function(data, truth, forms) 
{
  # W = data$W
  # nn = nrow(W)
  nn = nrow(data)
  Mstarform = forms$Mstarform
  Zstarform = forms$Zstarform
  # fit M for S = 1
  # data = cbind(W,A=data$A, Z=data$Z, M=data$M, Y=data$Y)
  dataZ1 = dataZ0 = dataA1 = dataA0 = data
  dataZ1$Z = 1
  dataZ0$Z = 0
  dataA1$A = 1
  dataA0$A = 0
  
  dataMstar = model.matrix(Mstarform, data)[,-1]
  Mstarfit = cv.glmnet(dataMstar, data$M, family = "binomial")
  
  dataZstar = model.matrix(Zstarform, data)[,-1]
  Zstarfit = cv.glmnet(dataZstar, data$Z, family = "binomial")

  dataMz1 = model.matrix(Mstarform, dataZ1)[,-1]
  predMz1 = predict(Mstarfit, newx = dataMz1, type = 'response')
  
  dataMz0 = model.matrix(Mstarform, dataZ0)[,-1]
  predMz0 = predict(Mstarfit, newx = dataMz0, type = 'response')
  
  dataZa0 = model.matrix(Zstarform, dataA0)[,-1]
  predZa0 = predict(Zstarfit, newx = dataZa0, type = 'response')
  
  dataZa1 = model.matrix(Zstarform, dataA1)[,-1]
  predZa1 = predict(Zstarfit, newx = dataZa1, type = 'response')
  
  gstarM_astar0 = predMz1*predZa0 + predMz0*(1 - predZa0)
  gstarM_astar1 = predMz1*predZa1 + predMz0*(1 - predZa1)
  
  return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0))
}


#' @export
get.mediation.initdata_lasso = function(data, forms, RCT = 0.5) {

  # data = cbind(data$W,A=data$A, Z=data$Z, M=data$M, Y=data$Y)
  df_YM1 = data
  df_YM1$M = 1
  
  df_YM0 = df_YM1
  df_YM0$M = 0
  Yform = forms$Yform
  Mform = forms$Mstarform
  Zform = forms$Zstarform
  Aform = forms$Aform
  # convert to model matrices to match for prediction
  df_YM1 = model.matrix(Yform, df_YM1)[,-1]
  df_YM0 = model.matrix(Yform, df_YM0)[,-1]
  
  df_Z = data

  Yform = forms$Yform
  dataY = model.matrix(Yform, data)[,-1]
  Yfit = cv.glmnet(dataY, data$Y, family = "binomial")
  
  # Mform = paste0("M ~ ", paste(covariates$covariates_M, collapse = "+"))
  dataM = model.matrix(Mform, data)[,-1]
  Mfit = cv.glmnet(dataM, data$M, family = "binomial")
  
  # Zform = paste0("Z ~ ", paste(covariates$covariates_Z, collapse = "+"))
  dataZ = model.matrix(Zform, data)[,-1]
  Zfit = cv.glmnet(dataZ, data$Z, family = "binomial")
  
  # propensity scores
  if (is.null(RCT)) { 
    dataA = model.matrix(Aform, data)[,-1]
    Afit = cv.glmnet(dataA, data$A, family = "binomial")
    A_ps = predict(Afit, newx = dataA, type = 'response')
  } else A_ps = RCT

  # as clev cov is 0 otherwise 
  M_ps = predict(Mfit, newx = dataM, type = 'response')

  # Predict Y for whole data, also with M = 1 and 0
  Y_init = pmin(pmax(predict(Yfit, newx = dataY, type = 'response'), 0.001), .999)
  Y_init_M1 = pmin(pmax(predict(Yfit, newx = df_YM1, type = 'response'), 0.001), .999)
  Y_init_M0 = pmin(pmax(predict(Yfit, newx = df_YM0, type = 'response'), 0.001), .999)
  return(list(initdata = list(M_ps = M_ps, A_ps = A_ps), 
              Y_preds = list(Y_init = Y_init, 
                             Y_init_M1 = Y_init_M1, 
                             Y_init_M0 = Y_init_M0)))
}



#' @export
mediation.step1_lasso = function(initdata, Y_preds, data, gstarM_astar, a, iptw = TRUE) {
  H = with(data, with(initdata, ((A == a)*((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))/
                        (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*(A_ps*A + (1 - A)*(1 - A_ps))))))
  
  # updates
  Qfit = try(glm(data$Y ~ 1 + offset(qlogis(Y_preds$Y_init)), family = binomial,
                 weights = H), silent = TRUE)
  
  if (class(Qfit)=="try-error") eps = 0 else eps = Qfit$coefficients
  
  est_iptw = mean(data$Y*H)
  IC_iptw = data$Y*H - est_iptw
  if (iptw) {
    return(list(Qstar_M  = plogis(qlogis(Y_preds$Y_init) + eps),
              Qstar_M1 = plogis(qlogis(Y_preds$Y_init_M1) + eps),
              Qstar_M0 = plogis(qlogis(Y_preds$Y_init_M0) + eps),
              IC_iptw = IC_iptw, est_iptw = est_iptw,
              Hm = H,
              eps = eps))
  } else {
    return(list(Qstar_M  = plogis(qlogis(Y_preds$Y_init) + eps),
                Qstar_M1 = plogis(qlogis(Y_preds$Y_init_M1) + eps),
                Qstar_M0 = plogis(qlogis(Y_preds$Y_init_M0) + eps),
                Hm = H,
                eps = eps))
  }
}


#' @export
mediation.step2_lasso = function(data, Qstar_M, Qstar_Mg, Hm, A_ps, a, tmle = TRUE,
                               EE = FALSE, bootstrap = FALSE, form) {
  QZform = form
  Qstar_Mg = as.vector(Qstar_Mg)
  df = cbind(data, Qstar_Mg = Qstar_Mg)
  
  # QZform = formula(paste0("Qstar_Mg ~ ", paste(covariates$covariates_QZ, collapse = "+")))
  df_QZ = model.matrix(QZform, df)[,-1]
  QZfit = cv.glmnet(df_QZ[data$A==a, ], Qstar_Mg[data$A==a], family = "gaussian")
  # compute the clever covariate and update if tmle
  if (tmle) {
    Hz = with(data, (A == a)/(A*A_ps + (1 - A)*(1 - A_ps)))  
    QZ_preds_a = pmin(pmax(predict(QZfit, newx = df_QZ,  s="lambda.1se"), .001), .999)
    # update
    QZfit_tmle = try(glm(Qstar_Mg ~ 1 + offset(qlogis(QZ_preds_a)), family = binomial,
                         weights = Hz), silent = TRUE)
    if (class(QZfit_tmle)=="try-error") eps2 = 0 else eps2 = QZfit_tmle$coefficients
    
    QZstar_a = plogis(qlogis(QZ_preds_a) + eps2)
    est = mean(QZstar_a)
    if(!bootstrap) { 
      D_Y = with(data, Hm*(Y - Qstar_M))
      D_Z = Hz*(Qstar_Mg - QZstar_a)
      D_W = with(data, (QZstar_a - est))
      D = D_Y + D_Z + D_W
    }
  } 
  
  # regress if EE or mle, EE updates the estimate, mle does not
  if (EE) {
    QZstar_a = pmin(pmax(predict(QZfit, newx = df_QZ, s="lambda.1se"), .001), .999) 
    init_est = mean(QZstar_a)
    D_Y1s = with(data, Hm*(Y - Qstar_M))
    Hz = with(data, (A == a)/(A*A_ps + (1 - A)*(1 - A_ps)))
    D_Z1s = Hz*(Qstar_Mg - QZstar_a)
    D_W1s = with(data, QZstar_a - init_est)
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

