
#' @title SDE_tmle3
#' @description computes the sequential regression, targeted maximum likelihood estimate
#' for the stochastic direct effect or stochastic indirect effect when the outcome and 
#' mediator model are only available on site 1 (S = 1).  This is a data adaptive parameter
#' as the stochastic direct effect has a model for the mediator is determined on the data for 
#' site 1. 
#' @param data, data.frame of variables in time ordering from left to right
#' @param a, the treatment intervention of interest
#' @param a_star, the treatment intervention of the stochastic model for Mediator, M
#' @param sl the sl3 superlearner defined, see sl3 documentation for defining a superlearner
#' and the example below
#' @param V number of folds for cross-validation (fixed to 10 for now)
#' @param covariates, list of covariates for each necessary model, going backwards from the 
#' outcome, Y, M, Z, A, W, S where S is the binary site, W are confounders, A is the treatment
#' Z is the intermediary confounder (binary) and M is the mediator, Y is the outcome. 
#' @param truth for testing purposes input a list with names f_W, f_S, f_Z and f_Y models representing
#' and the corresponding elements are functions of appropriate variables so as to be able to 
#' generate the truth and check the estimator's performance.  Default is NULL
#' @param truncate, a list with elements lower and upper to truncate the various p-scores  
#' not functional at present
#' @return  a list with a CI for the estimate, and estimate using linear main terms MLE 
#' gcomp formula (est_mle), the influence curve (IC), the superlearner coefficients for
#' the Y model and the QZ model (SL_coef)
#' @export
#' @example /inst/tester_SDE_tmle.R
SDE_tmle3 = function(data, sl, V=10, covariates, truth = NULL, 
                     truncate = list(lower =.0001, upper = .9999), glm_only = TRUE,
                     iptw = TRUE, onestep = TRUE, B = 500) 
{
  
  # n = 1e3
  # data = gendata.SDEtransport(n, f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
  # 
  if (glm_only) sl = make_learner(Lrnr_glm_fast, family = binomial())
  
  L = truncate$lower
  U = truncate$upper
  
  # get the stochastic dist of M and true params
  gstar_info = get_gstarM(data = data, sl, V=10, covariates = covariates, truth = truth)
  gstarM_astar1 = gstar_info$gstarM_astar1 
  gstarM_astar0 = gstar_info$gstarM_astar0 
  Psi_astar0a1_0 = gstar_info$Psi_astar0a1
  Psi_astar0a0_0 = gstar_info$Psi_astar0a0
  Psi_astar1a1_0 = gstar_info$Psi_astar1a1
  PS0_0 = gstar_info$PS0_0
  
  # get the true IC's
  df_ZS0 = data
  df_ZS0$S = 0
  
  S_ps0 = with(data, truth$f_S(W=W))
  ZS0_ps0 = with(df_ZS0, truth$f_Z(A=A, W=W, S=S))
  M_ps0 = with(data, truth$f_M(Z=Z, W=W, S=1))
  Z_ps0 = with(data, truth$f_Z(A=A, W=W, S=S))
  A_ps0 = with(data, truth$f_A(W=W, S=S))
  
  get_cc_0 = function(data, gstarM_astar, a) {
    with(data, ((S == 1)*(A == a)*
                  ((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))*
                  ((Z == 1)*ZS0_ps0 + (Z == 0)*(1 - ZS0_ps0))*(1 - S_ps0))/
           (((M == 1)*M_ps0 + (M == 0)*(1 - M_ps0))*
              ((Z == 1)*Z_ps0 + (Z == 0)*(1 - Z_ps0))*
              ((A == 1)*A_ps0 + (A == 0)*(1 - A_ps0))*S_ps0
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
  
  D_astar0a1_0 = get_trueIC(gstar_M = gstarM_astar0, a = 1, 
                            Psi_0 = Psi_astar0a1_0, Hm_0 = Hm_astar0a1_0,
                            Hz_0 = Hz_astar0a1_0)
  
  D_astar0a0_0 = get_trueIC(gstar_M = gstarM_astar0, a = 0, 
                            Psi_0 = Psi_astar0a0_0, Hm_0 = Hm_astar0a0_0,
                            Hz_0 = Hz_astar0a0_0)
  
  D_astar1a1_0 = get_trueIC(gstar_M = gstarM_astar1, a = 1, 
                            Psi_0 = Psi_astar1a1_0, Hm_0 = Hm_astar1a1_0,
                            Hz_0 = Hz_astar1a1_0)
  
  
  ####
  ####
  # Here we call a function to compute the estimates based on stochastic M
  ####
  ####
  get_estimates = function(data, gstarM_astar1, gstarM_astar0, bootstrap) {
    
    W = data[,grep("W", colnames(data))]
    df_YM1S1 = data
    df_YM1S1$M = 1
    df_YM1S1$S = 1
    
    df_YM0S1 = df_YM1S1
    df_YM0S1$M = 0
    
    df_ZS0 = data
    df_ZS0$S = 0
    
    # to be used for fitting Y
    task_YsubS1 <- sl3_Task$new(
      data = data.table::copy(data[data$S==1,]),
      covariates = covariates$covariates_Y,
      outcome = "Y"
    )
    
    # Used for predicting Y
    task_Y <- sl3_Task$new(
      data = data.table::copy(data),
      covariates = covariates$covariates_Y,
      outcome = "Y"
    )
    
    task_YM1S1 <- sl3_Task$new(
      data = data.table::copy(df_YM1S1),
      covariates = covariates$covariates_Y,
      outcome = "Y"
    )
    
    task_YM0S1 <- sl3_Task$new(
      data = data.table::copy(df_YM0S1),
      covariates = covariates$covariates_Y,
      outcome = "Y"
    )
    
    # used for fitting M
    task_MsubS1 <- sl3_Task$new(
      data = data.table::copy(data[data$S == 1,]),
      covariates = covariates$covariates_M,
      outcome = "M"
    )
    
    # used for predicting M
    task_MS1 <- sl3_Task$new(
      data = data.table::copy(df_YM1S1),
      covariates = covariates$covariates_M,
      outcome = "M"
    )
    
    
    task_Z <- sl3_Task$new(
      data = data.table::copy(data),
      covariates = covariates$covariates_Z,
      outcome = "Z"
    )
    
    task_ZS0 <- sl3_Task$new(
      data = data.table::copy(df_ZS0),
      covariates = covariates$covariates_Z,
      outcome = "Z"
    )
    
    task_A <- sl3_Task$new(
      data = data.table::copy(data),
      covariates = covariates$covariates_A,
      outcome = "A"
    )
    
    task_ZS1 <- sl3_Task$new(
      data = data.table::copy(df_YM1S1),
      covariates = covariates$covariates_Z,
      outcome = "Z"
    )
    
    task_S <- sl3_Task$new(
      data = data.table::copy(data),
      covariates = covariates$covariates_S,
      outcome = "S"
    )
    
    # Y, M, Z, A, S fits
    Sfit = sl$train(task_S)
    Afit = sl$train(task_A)
    Zfit = sl$train(task_Z)
    
    # fitting M and Y on subsets where S is 1.  This is all that is needed for M and Y
    Mfit = sl$train(task_MsubS1)
    Yfit = sl$train(task_YsubS1)
    
    # propensity scores
    S_ps = Sfit$predict()
    PS0 = mean(data$S == 0)
    A_ps = Afit$predict()
    Z_ps = Zfit$predict()
    ZS0_ps = Zfit$predict(task_ZS0)
    # might as well predict for S = 1 on whole data because only S = 1 subset is relevant
    # as clev cov is 0 otherwise 
    M_ps = Mfit$predict(task_MS1)
    
    # 1st clever cov FOR SDE, ALSO NEED THE ONE FOR SIE
    get_cc = function(data, gstarM_astar, a) {
      with(data, ((S == 1)*(A == a)*
                    ((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))*
                    ((Z == 1)*ZS0_ps + (Z == 0)*(1 - ZS0_ps))*(1 - S_ps))/
             (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                (A_ps*A + (1 - A)*(1 - A_ps))*S_ps*PS0))
    }
    
    
    Hm_astar0a1 = get_cc(data = data, gstarM_astar = gstarM_astar0, a = 1)
    Hm_astar0a0 = get_cc(data = data, gstarM_astar = gstarM_astar0, a = 0)
    Hm_astar1a1 = get_cc(data = data, gstarM_astar = gstarM_astar1, a = 1)
    
    # Predict Y for whole data, also with M = 1 and 0, S does not matter since
    # Y is not a function of that variable so won't be used for prediction
    Y_preds = Yfit$predict(task_Y)
    Y_preds_M1 = Yfit$predict(task_YM1S1)
    Y_preds_M0 = Yfit$predict(task_YM0S1)
    
    # updates
    update_fcn = function(weights) {
      Qfit = try(glm(data$Y ~ 1 + offset(qlogis(Y_preds)), family = binomial,
                 weights = weights), silent = TRUE)

      if (class(Qfit)=="try-error") eps = 0 else eps = Qfit$coefficients
      
      Qstar_M  = plogis(qlogis(Y_preds) + eps)
      Qstar_M1 = plogis(qlogis(Y_preds_M1) + eps)
      Qstar_M0 = plogis(qlogis(Y_preds_M0) + eps)
      df = data.frame(Qstar_M = Qstar_M, Qstar_M0 = Qstar_M0, Qstar_M1 = Qstar_M1)
      return(list(df = df, eps = eps))
    }
    
    QstarY_astar0a1_info = update_fcn(weights = Hm_astar0a1)
    QstarY_astar0a0_info = update_fcn(weights = Hm_astar0a0)
    QstarY_astar1a1_info = update_fcn(weights = Hm_astar1a1)
    
    QstarY_astar0a1_df = QstarY_astar0a1_info$df
    QstarY_astar0a0_df = QstarY_astar0a0_info$df 
    QstarY_astar1a1_df = QstarY_astar1a1_info$df
    
    QstarMg_astar0a1 = QstarY_astar0a1_df[,3]*gstarM_astar0 + 
        QstarY_astar0a1_df[,2]*(1 - gstarM_astar0)
    
    QstarMg_astar0a0 = QstarY_astar0a0_df[,3]*gstarM_astar0 + 
       QstarY_astar0a0_df[,2]*(1 - gstarM_astar0)
    
    QstarMg_astar1a1 = QstarY_astar1a1_df[,3]*gstarM_astar1 + 
       QstarY_astar1a1_df[,2]*(1 - gstarM_astar1)
    
    YMg_astar0a1 = YMg_astar0a0 = Y_preds_M1*gstarM_astar0 + 
       Y_preds_M0*(1 - gstarM_astar0)
    
    YMg_astar1a1 = Y_preds_M1*gstarM_astar1 + 
        Y_preds_M0*(1 - gstarM_astar1)
    
    # NEXT REGRESSION
    regress_step2 = function(Qstar_M, Qstar_Mg, Y_Mg, Hm, a) {
      df_QZ = data
      df_QZ$Qstar_Mg = Qstar_Mg
      df_QZ$Y_Mg = Y_Mg  
      
      # for tmle 
      task_QZsubAa <- sl3_Task$new(
        data = data.table::copy(df_QZ[df_QZ$A == a, ]),
        covariates = covariates$covariates_QZ,
        outcome = "Qstar_Mg"
      )
      
      # for standard gcomp
      task_QZ_YsubAa <- sl3_Task$new(
        data = data.table::copy(df_QZ[df_QZ$A == a, ]),
        covariates = covariates$covariates_QZ,
        outcome = "Y_Mg"
      )
      
      task_data <- sl3_Task$new(
        data = data.table::copy(data),
        covariates = covariates$covariates_QZ,
        outcome = "Y"
      )
      
      QZfit = sl$train(task_QZsubAa)
      YZfit = sl$train(task_QZ_YsubAa)
      
      # compute the clever covariate 2
      Hz = with(data, (A == a)*(S == 0)/(A*A_ps + (1 - A)*(1 - A_ps))/PS0)

      # estimates predicted on whole data.  Since they were fit on A = a, A is not
      # a covariate in these predictions as A = a is therefore implicit
      QZ_preds_a = pmin(pmax(QZfit$predict(task_data), .001), .999)
      YZ_preds_a = pmin(pmax(YZfit$predict(task_data), .001), .999)
      
      # update
      QZfit_tmle = try(glm(Qstar_Mg ~ 1 + offset(qlogis(QZ_preds_a)), family = binomial,
                       weights = Hz), silent = TRUE)
      if (class(QZfit_tmle)=="try-error") eps2 = 0 else eps2 = QZfit_tmle$coefficients

      QZstar_a = plogis(qlogis(QZ_preds_a) + eps2)
      
      # compute the parameter estimate
      est = mean(QZstar_a[data$S==0])
      est_1s_init = mean(YZ_preds_a[data$S==0])
      est_iptw = mean(with(data, Y*Hm))
      est_mle = mean(YZ_preds_a[data$S==0])
      
      D_Y = with(data, Hm*(Y - Qstar_M))
      D_Y1s = with(data, Hm*(Y - Y_preds))
      
      D_Z = Hz*(Qstar_Mg - QZstar_a)
      D_Z1s = Hz*(Y_Mg - YZ_preds_a)
      
      D_W = with(data, (QZstar_a - est)*(S ==0)/PS0)
      D_W1s = with(data, (YZ_preds_a - est_1s_init)*(S ==0)/PS0)
      
      D = D_Y + D_Z + D_W
      D_1s = D_Y1s + D_Z1s + D_W1s
      D_iptw = with(data, Y*Hm) - est_iptw
      est_1s = est_1s_init + mean(D_1s)
      
      if (!bootstrap) {
        n = nrow(data)
        SE = sd(D)/sqrt(n)
        SE_1s = sd(D_1s)/sqrt(n)
        SE_iptw = sd(D_iptw)/sqrt(n)
        
        CI = c(est = est, left = est - 1.96*SE, right = est + 1.96*SE)
        CI_1s = c(est_1s, left = est_1s - 1.96*SE_1s, right = est_1s + 1.96*SE_1s)
        CI_iptw = c(est_iptw, left = est_iptw - 1.96*SE_iptw, right = est_iptw + 1.96*SE_iptw)
        
        return(list(est = est, est_1s = est_1s, est_iptw = est_iptw,
                    est_mle = est_mle, IC = D, IC_1s = D_1s, IC_iptw = D_iptw, eps2 = eps2))
      } else {
        return(list(est = est, est_1s = est_1s, est_iptw = est_iptw, est_mle = est_mle))
      }
    }
    
    # undebug(regress_step2)
    info_astar0a1 = regress_step2(Qstar_M = QstarY_astar0a1_df$Qstar_M,
                                  Qstar_Mg = QstarMg_astar0a1, 
                                  Y_Mg = YMg_astar0a1, 
                                  Hm = Hm_astar0a1, 
                                  a = 1)
    
    info_astar0a0 = regress_step2(Qstar_M = QstarY_astar0a0_df$Qstar_M,
                                  Qstar_Mg = QstarMg_astar0a0, 
                                  Y_Mg = YMg_astar0a0, 
                                  Hm = Hm_astar0a0, 
                                  a = 0)
    
    info_astar1a1 = regress_step2(Qstar_M = QstarY_astar1a1_df$Qstar_M,
                                  Qstar_Mg = QstarMg_astar1a1, 
                                  Y_Mg = YMg_astar1a1, 
                                  Hm = Hm_astar1a1, 
                                  a = 1)
    
    return(list(info_astar0a1 = info_astar0a1, 
                info_astar0a0 = info_astar0a0, 
                info_astar1a1 = info_astar1a1,
                eps_astar0a1 = QstarY_astar0a1_info$eps, 
                eps_astar0a0 = QstarY_astar0a0_info$eps, 
                eps_astar1a1 = QstarY_astar1a1_info$eps
                ))
  }
  
  # get the estimates and IC's
  info = get_estimates(data, gstarM_astar1, gstarM_astar0, bootstrap = FALSE)
  
  # run the bootstrap 500 times
  # bootstrap from the data and thus the gstarM_astar1 and gstarM_astar0
  n = nrow(data)
  boots = lapply(1:B, FUN = function(x) {
    # print (i)
    # x = sample(1:1000000,1)
    # set.seed(x)
    inds = sample(1:n, replace = TRUE)
    data = data[inds,]
    gstarM_astar1 = gstarM_astar1[inds]
    gstarM_astar0 = gstarM_astar0[inds]
    get_estimates(data, gstarM_astar1, gstarM_astar0, bootstrap = TRUE)
  })

  boot_ests = lapply(boots, FUN = function(boot) {
    SDE = unlist(boot[[1]]) - unlist(boot[[2]])
    SIE = unlist(boot[[3]]) - unlist(boot[[1]])
    return(list(SDE, SIE))
  })

  boots_SDE = do.call(rbind, lapply(boot_ests, FUN = function(boot) boot[[1]]))
  boots_SIE = do.call(rbind, lapply(boot_ests, FUN = function(boot) boot[[2]]))

  bootSE_SDE = apply(boots_SDE, 2, sd)
  bootSE_SIE = apply(boots_SIE, 2, sd)

  D_SDE = info$info_astar0a1$IC - info$info_astar0a0$IC
  D_SIE = info$info_astar1a1$IC - info$info_astar0a1$IC
  
  D_SDE_1s = info$info_astar0a1$IC_1s - info$info_astar0a0$IC_1s
  D_SIE_1s = info$info_astar1a1$IC_1s - info$info_astar0a1$IC_1s
  
  D_SDE_iptw = info$info_astar0a1$IC_iptw - info$info_astar0a0$IC_iptw
  D_SIE_iptw = info$info_astar1a1$IC_iptw - info$info_astar0a1$IC_iptw
  
  D_SDE_0 = D_astar0a1_0 - D_astar0a0_0
  D_SIE_0 = D_astar1a1_0 - D_astar0a1_0
  
  SE_SDE = sd(D_SDE)/sqrt(n)
  SE_SIE = sd(D_SIE)/sqrt(n)
  
  SE_SDE_1s = sd(D_SDE_1s)/sqrt(n)
  SE_SIE_1s = sd(D_SIE_1s)/sqrt(n)
  
  SE_SDE_iptw = sd(D_SDE_iptw)/sqrt(n)
  SE_SIE_iptw = sd(D_SIE_iptw)/sqrt(n)
  
  SE_SDE_0 = sd(D_SDE_0)/sqrt(n)
  SE_SIE_0 = sd(D_SIE_0)/sqrt(n)
  
  ests_astar0a1 =   c(info$info_astar0a1$est,
                      info$info_astar0a1$est_1s,
                      info$info_astar0a1$est_iptw,
                      info$info_astar0a1$est_mle,
                      Psi_astar0a1_0)
  
  ests_astar0a0 =   c(info$info_astar0a0$est,
                      info$info_astar0a0$est_1s,
                      info$info_astar0a0$est_iptw,
                      info$info_astar0a0$est_mle,
                      Psi_astar0a0_0)
  
  ests_astar1a1 =   c(info$info_astar1a1$est,
                      info$info_astar1a1$est_1s,
                      info$info_astar1a1$est_iptw,
                      info$info_astar1a1$est_mle,
                      Psi_astar1a1_0)
  
  SDE_ests = ests_astar0a1 - ests_astar0a0
  SIE_ests = ests_astar1a1 - ests_astar0a1
  
  CI_SDE = c(SDE_ests[1], SDE_ests[1] - 1.96*SE_SDE, SDE_ests[1] + 1.96*SE_SDE)
  CI_SIE = c(SIE_ests[1], SIE_ests[1] - 1.96*SE_SIE, SIE_ests[1] + 1.96*SE_SIE)
  
  CI_SDE_1s = c(SDE_ests[2], SDE_ests[2] - 1.96*SE_SDE_1s , SDE_ests[2] + 1.96*SE_SDE_1s)
  CI_SIE_1s  = c(SIE_ests[2], SIE_ests[2] - 1.96*SE_SIE_1s , SIE_ests[2] + 1.96*SE_SIE_1s)
  
  CI_SDE_iptw = c(SDE_ests[3], SDE_ests[3] - 1.96*SE_SDE_iptw, SDE_ests[3] + 1.96*SE_SDE_iptw)
  CI_SIE_iptw = c(SIE_ests[3], SIE_ests[3] - 1.96*SE_SIE_iptw, SIE_ests[3] + 1.96*SE_SIE_iptw)
  
  CI_SDE_boot = c(SDE_ests[1], SDE_ests[1] - 1.96*bootSE_SDE[1], SDE_ests[1] + 1.96*bootSE_SDE[1])
  CI_SIE_boot = c(SIE_ests[1], SIE_ests[1] - 1.96*bootSE_SIE[1], SIE_ests[1] + 1.96*bootSE_SIE[1])
  
  CI_SDE_1s_boot = c(SDE_ests[2], SDE_ests[2] - 1.96*bootSE_SDE[2] , SDE_ests[2] + 1.96*bootSE_SDE[2])
  CI_SIE_1s_boot  = c(SIE_ests[2], SIE_ests[2] - 1.96*bootSE_SIE[2] , SIE_ests[2] + 1.96*bootSE_SIE[2])
  
  CI_SDE_iptw_boot = c(SDE_ests[3], SDE_ests[3] - 1.96*bootSE_SDE[3], SDE_ests[3] + 1.96*bootSE_SDE[3])
  CI_SIE_iptw_boot = c(SIE_ests[3], SIE_ests[3] - 1.96*bootSE_SIE[3], SIE_ests[3] + 1.96*bootSE_SIE[3])
  
  return(list(CI_SDE = CI_SDE, CI_SIE = CI_SIE, CI_SDE_1s = CI_SDE_1s, CI_SIE_1s = CI_SIE_1s,
  CI_SDE_iptw = CI_SDE_iptw, CI_SIE_iptw = CI_SIE_iptw, CI_SDE_boot = CI_SDE_boot, 
  CI_SIE_boot = CI_SIE_boot, CI_SDE_1s_boot = CI_SDE_1s_boot, CI_SIE_1s_boot = CI_SIE_1s_boot, 
  CI_SDE_iptw_boot = CI_SDE_iptw_boot, CI_SIE_iptw_boot = CI_SIE_iptw_boot,
  SDE_0 = Psi_astar0a1_0 - Psi_astar0a0_0,
  SDE_0 = Psi_astar1a1_0 - Psi_astar0a1_0, 
  SE_SDE = SE_SDE, 
  SE_SIE = SE_SIE, 
  SE_SDE_1s = SE_SDE_1s,  
  SE_SIE_1s = SE_SIE_1s, 
  SE_SDE_iptw = SE_SDE_iptw,   
  SE_SIE_iptw = SE_SIE_iptw, 
  SE_SDE_0 = SE_SDE_0, 
  SE_SIE_0 = SE_SIE_0,
  eps2_astar0a1 = info$info_astar0a1$eps2,
  eps2_astar0a0 = info$info_astar0a0$eps2,
  eps2_astar1a1 = info$info_astar1a1$eps2,
  eps_astar0a1 = info$eps_astar0a1, 
  eps_astar0a0 = info$eps_astar0a0, 
  eps_astar1a1 = info$eps_astar1a1)) 
} 

#' @title get_gstarM
#' @description computes the sequential regression, targeted maximum likelihood estimate
#' for the stochastic direct effect or stochastic indirect effect when the outcome and 
#' mediator model are only available on site 1 (S = 1).  This is a data adaptive parameter
#' as the stochastic direct effect has a model for the mediator is determined on the data for 
#' site 1. 
#' @param data, data.frame of variables in time ordering from left to right
#' @param a, the treatment intervention of interest
#' @param a_star, the treatment intervention of the stochastic model for Mediator, M
#' @param sl the sl3 superlearner defined, see sl3 documentation for defining a superlearner
#' and the example below
#' @param V number of folds for cross-validation (fixed to 10 for now)
#' @param covariates, list of covariates for each necessary model, going backwards from the 
#' outcome, Y, M, Z, A, W, S where S is the binary site, W are confounders, A is the treatment
#' Z is the intermediary confounder (binary) and M is the mediator, Y is the outcome. 
#' @param truth for testing purposes input a list with names f_W, f_S, f_Z and f_Y models representing
#' and the corresponding elements are functions of appropriate variables so as to be able to 
#' generate the truth and check the estimator's performance.  Default is NULL
#' @param truncate, a list with elements lower and upper to truncate the various p-scores  
#' not functional at present
#' @return  a list with a CI for the estimate, and estimate using linear main terms MLE 
#' gcomp formula (est_mle), the influence curve (IC), the superlearner coefficients for
#' the Y model and the QZ model (SL_coef)
#' @export
#' @example /inst/tester_SDE_tmle.R
get_gstarM  = function(data, sl, V=10, covariates, truth) 
{
  W = data[,grep("W", colnames(data))]
  task_Mstar <- sl3_Task$new(
    data = data.table::copy(data[data$S == 1,]),
    covariates = covariates$covariates_M,
    outcome = "M"
  )
  
  df_MZ1 = df_MZ0 = data
  df_MZ1$Z = 1
  df_MZ0$Z = 0
  
  task_MZ1 <- sl3_Task$new(
    data = data.table::copy(df_MZ1),
    covariates = covariates$covariates_M,
    outcome = "M"
  )
  
  task_MZ0 <- sl3_Task$new(
    data = data.table::copy(df_MZ0),
    covariates = covariates$covariates_M,
    outcome = "M"
  )
  
  task_Zstar <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  df_ZA1 = df_ZA0 = data
  df_ZA1$A = 1
  df_ZA0$A = 0
  
  task_ZA1 <- sl3_Task$new(
    data = data.table::copy(df_ZA1),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  task_ZA0 <- sl3_Task$new(
    data = data.table::copy(df_ZA0),
    covariates = covariates$covariates_Z,
    outcome = "Z"
  )
  
  Mstarfit = sl$train(task_Mstar)
  Zstarfit = sl$train(task_Zstar)
  
  # find the truth if simulating
  f_truth0 = function(Z, W, S) {
    # W = W[1:10]
    nn = nrow(W)
    # dataM1 = data.frame(cbind)
    dataM1 = cbind(Z = rep(1,nn), W, M = rep(1,nn))
    taskM1 =   sl3_Task$new(
      data = data.table::copy(dataM1),
      covariates = covariates$covariates_M,
      outcome = "M"
    )
    predM1 = Mstarfit$predict(taskM1)
    
    dataM0 = cbind(Z = rep(0,nn), W, M = rep(1,nn))
    taskM0 =   sl3_Task$new(
      data = data.table::copy(dataM0),
      covariates = covariates$covariates_M,
      outcome = "M"
    )
    
    predM0 = Mstarfit$predict(taskM0)
    
    dataZ = cbind(A = rep(0,nn), W, S = rep(1, nn), Z = rep(1,nn))
    
    taskZ <- sl3_Task$new(
      data = data.table::copy(dataZ),
      covariates = covariates$covariates_Z,
      outcome = "Z"
    )
    
    predZ = Zstarfit$predict(taskZ)
    gM = predM1*predZ + predM0*(1 - predZ)
    return(gM)
  }
  
  f_truth1 = function(Z, W, S) {
    # W = W[1:10]
    nn = nrow(W)
    # dataM1 = data.frame(cbind)
    dataM1 = cbind(Z = rep(1,nn), W, M = rep(1,nn))
    taskM1 =   sl3_Task$new(
      data = data.table::copy(dataM1),
      covariates = covariates$covariates_M,
      outcome = "M"
    )
    predM1 = Mstarfit$predict(taskM1)
    
    dataM0 = cbind(Z = rep(0,nn), W, M = rep(1,nn))
    taskM0 =   sl3_Task$new(
      data = data.table::copy(dataM0),
      covariates = covariates$covariates_M,
      outcome = "M"
    )
    
    predM0 = Mstarfit$predict(taskM0)
    
    dataZ = cbind(A = rep(1,nn), W, S = rep(1, nn), Z = rep(1,nn))
    
    taskZ <- sl3_Task$new(
      data = data.table::copy(dataZ),
      covariates = covariates$covariates_Z,
      outcome = "Z"
    )
    
    predZ = Zstarfit$predict(taskZ)
    gM = predM1*predZ + predM0*(1 - predZ)
    return(gM)
  }
  
  data_pop_astar0a1 = gendata.SDEtransport(5*1e6, f_W = truth$f_W, f_S = truth$f_S, f_A = function(S,W) 1, 
                                           f_Z = truth$f_Z, f_M = f_truth0, f_Y = truth$f_Y) 
  data_pop_astar0a0 = gendata.SDEtransport(5*1e6, f_W = truth$f_W, f_S = truth$f_S, f_A = function(S,W) 0, 
                                           f_Z = truth$f_Z, f_M = f_truth0, f_Y = truth$f_Y) 
  
  data_pop_astar1a1 = gendata.SDEtransport(5*1e6, f_W = truth$f_W, f_S = truth$f_S, f_A = function(S,W) 1, 
                                           f_Z = truth$f_Z, f_M = f_truth1, f_Y = truth$f_Y) 
  
  assign("gstarM_astar1", with(data, f_truth1(Z=Z, W=W, S=S)))
  assign("gstarM_astar0", with(data, f_truth0(Z=Z, W=W, S=S)))
  
  
  Psi_astar0a1 = mean(data_pop_astar0a1$Y[data_pop_astar0a1$S==0]) 
  Psi_astar0a0 = mean(data_pop_astar0a0$Y[data_pop_astar0a0$S==0])
  Psi_astar1a1 = mean(data_pop_astar1a1$Y[data_pop_astar1a1$S==0]) 
  PS0_0 = mean(data_pop_astar1a1$S==0)   
  
  
  return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0, 
              Psi_astar0a1 = Psi_astar0a1, Psi_astar0a0 = Psi_astar0a0, 
              Psi_astar1a1= Psi_astar1a1, PS0_0 = PS0_0)) 
}

