
#' @title SDE_tmle2
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
SDE_tmle2 = function(data, a, a_star, sl, V=10, covariates, truth = NULL, 
                    truncate = list(lower =.0001, upper = .9999), glm_only = TRUE,
                    iptw = TRUE, onestep = TRUE) 
  {
  
  n = 1e4
  data = gendata.SDEtransport(n, f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
  
  W = data[,grep("W", colnames(data))]
  
  if (glm_only) sl = make_learner(Lrnr_glm_fast, family = binomial())
  
  L = truncate$lower
  U = truncate$upper
  
  task_Mstar <- sl3_Task$new(
    data = data.table::copy(data[data$S == 1,]),
    covariates = covariates$covariates_M,
    outcome = "M",
    outcome_type = "binomial"
  )
  
  df_MZ1 = df_MZ0 = data
  df_MZ1$Z = 1
  df_MZ0$Z = 0
  
  task_MZ1 <- sl3_Task$new(
    data = data.table::copy(df_MZ1),
    covariates = covariates$covariates_M,
    outcome = "M",
    outcome_type = "binomial"
  )
  
  task_MZ0 <- sl3_Task$new(
    data = data.table::copy(df_MZ0),
    covariates = covariates$covariates_M,
    outcome = "M",
    outcome_type = "binomial"
  )
  
  
  
  task_Zstar <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Z,
    outcome = "Z",
    outcome_type = "binomial"
  )
  
  df_ZA1 = df_ZA0 = data
  df_ZA1$A = 1
  df_ZA0$A = 0
  
  task_ZA1 <- sl3_Task$new(
    data = data.table::copy(df_ZA1),
    covariates = covariates$covariates_Z,
    outcome = "Z",
    outcome_type = "binomial"
  )
  
  task_ZA0 <- sl3_Task$new(
    data = data.table::copy(df_ZA0),
    covariates = covariates$covariates_Z,
    outcome = "Z",
    outcome_type = "binomial"
  )
  
  Mstarfit = sl$train(task_Mstar)
  Zstarfit = sl$train(task_Zstar)
  
  # find the truth if simulating
  if (!is.null(truth)) {
    for (a_star in 0:1) {
    f_truth = function(Z, W, S) {
      # W = W[1:10]
      nn = nrow(W)
      # dataM1 = data.frame(cbind)
      dataM1 = cbind(Z = rep(1,nn), W, M = rep(1,nn))
      taskM1 =   sl3_Task$new(
        data = data.table::copy(dataM1),
        covariates = covariates$covariates_M,
        outcome = "M",
        outcome_type = "binomial"
      )
      predM1 = Mstarfit$predict(taskM1)
      
      dataM0 = cbind(Z = rep(0,nn), W, M = rep(1,nn))
      taskM0 =   sl3_Task$new(
        data = data.table::copy(dataM0),
        covariates = covariates$covariates_M,
        outcome = "M",
        outcome_type = "binomial"
      )
      
      predM0 = Mstarfit$predict(taskM0)
      
      dataZ = cbind(A = rep(a_star,nn), W, S = rep(1, nn), Z = rep(1,nn))
      
      taskZ <- sl3_Task$new(
        data = data.table::copy(dataZ),
        covariates = covariates$covariates_Z,
        outcome = "Z",
        outcome_type = "binomial"
      )
      
      predZ = Zstarfit$predict(taskZ)
      gM = predM1*predZ + predM0*(1 - predZ)
      return(gM)
    }
    assign(paste0("f_truth",a_star), f_truth)
    }
    
    data_pop1 = gendata.SDEtransport(5*1e6, f_W = truth$f_W, f_S = truth$f_S, f_A = function(S,W) 1, 
                                    f_Z = truth$f_Z, f_M = f_truth0, f_Y = truth$f_Y) 
    data_pop0 = gendata.SDEtransport(5*1e6, f_W = truth$f_W, f_S = truth$f_S, f_A = function(S,W) 0, 
                                     f_Z = truth$f_Z, f_M = f_truth0, f_Y = truth$f_Y) 
    data_pop11 = gendata.SDEtransport(5*1e6, f_W = truth$f_W, f_S = truth$f_S, f_A = function(S,W) 1, 
                                     f_Z = truth$f_Z, f_M = f_truth1, f_Y = truth$f_Y) 
      
    SDE_0 = mean(data_pop1$Y[data_pop1$S==0]) - mean(data_pop0$Y[data_pop0$S==0])
    SIE_0 = mean(data_pop11$Y[data_pop11$S==0]) - mean(data_pop1$Y[data_pop1$S==0])
    PS0_0 = mean(data_pop1$S==0)   
    }
  
  
  pred_MZ1 = Mstarfit$predict(task_MZ1)
  pred_MZ0 = Mstarfit$predict(task_MZ0)
  pred_ZA1 =  Zstarfit$predict(task_ZA1)
  gstarMastar1_ps = pred_MZ1*pred_ZA1 + pred_MZ0*(1 - pred_ZA1)
  
  pred_ZA0 = Zstarfit$predict(task_ZA0)
  gstarMastar0_ps = pred_MZ1*pred_ZA0 + pred_MZ0*(1 - pred_ZA0)
  
  df_YM1S1 = data
  df_YM1S1$M = 1
  df_YM1S1$S = 1
  
  df_YM0S1 = df_YM1S1
  df_YM0S1$M = 0
  
  df_ZS0 = data
  df_ZS0$S = 0

  task_YS1 <- sl3_Task$new(
    data = data.table::copy(data[data$S==1,]),
    covariates = covariates$covariates_Y,
    outcome = "Y",
    outcome_type = "binomial"
  )
  
  task_Y <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Y,
    outcome = "Y",
    outcome_type = "binomial"
  )
  
  task_YM1S1 <- sl3_Task$new(
    data = data.table::copy(df_YM1S1),
    covariates = covariates$covariates_Y,
    outcome = "Y",
    outcome_type = "binomial"
  )
  
  task_YM0S1 <- sl3_Task$new(
    data = data.table::copy(df_YM0S1),
    covariates = covariates$covariates_Y,
    outcome = "Y",
    outcome_type = "binomial"
  )
  
  
  task_ZS0 <- sl3_Task$new(
    data = data.table::copy(df_ZS0),
    covariates = covariates$covariates_Z,
    outcome = "Z",
    outcome_type = "binomial"
  )
  
  task_M <- sl3_Task$new(
    data = data.table::copy(data[data$S == 1,]),
    covariates = covariates$covariates_M,
    outcome = "M",
    outcome_type = "binomial"
  )
  
  task_MS1 <- sl3_Task$new(
    data = data.table::copy(df_YM1S1),
    covariates = covariates$covariates_M,
    outcome = "M",
    outcome_type = "binomial"
  )
  
  
  task_Z <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Z,
    outcome = "Z",
    outcome_type = "binomial"
  )
  
  task_ZS1 <- sl3_Task$new(
    data = data.table::copy(df_YM1S1),
    covariates = covariates$covariates_Z,
    outcome = "Z",
    outcome_type = "binomial"
  )
  
  Mfit = sl$train(task_M)
  Zfit = sl$train(task_Z)
  
  Z_S0_ps = Zfit$predict(task_ZS0)
  Z_S0_ps0 = with(df_ZS0, truth$f_Z(A=A, W=W, S=S))
  
  task_S <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_S,
    outcome = "S",
    outcome_type = "binomial"
  )
  
  # compute probs S = 0 given W
  Sfit = sl$train(task_S)
  S0_preds = 1 - Sfit$predict()
  S0_preds0 = 1 - with(data, truth$f_S(W=W))
  
  M_ps = rep(.5, nrow(data))
  M_ps[data$S==1] = Mfit$predict(task_M)
  
  M_S1_ps = Mfit$predict(task_MS1)
  M_ps0 = with(data, truth$f_M(Z=Z, W=W, S=1))
  
  # compute Z preds 
  Z_ps = Zfit$predict()
  Z_ps0 = with(data, f_Z(A=A, W=W, S=S))
  
  Z_S1_ps = Zfit$predict(task_ZS1)
  Z_S1_ps0 = with(df_YM1S1, truth$f_Z(A=A, W=W, S=S))
  
  # compute A=1 preds for S = 1
  
  task_A <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_A,
    outcome = "A",
    outcome_type = "binomial"
  )
  df_AS1 = data
  df_AS1$S = 1

  task_AS1 <- sl3_Task$new(
    data = data.table::copy(df_AS1),
    covariates = covariates$covariates_A,
    outcome = "A",
    outcome_type = "binomial"
  )

  Afit = sl$train(task_A)
  A_ps = Afit$predict()
  A_ps0 = with(data, truth$f_A(W=W, S=S))
  
  A_S1_ps = Afit$predict(task_AS1)
  A_S1_ps0 = with(df_AS1, truth$f_A(W=W, S=S))
  
  #compute prob S = 1 given W
  S1_preds = 1 - S0_preds
  S1_preds0 = 1 - S0_preds0
  
  # compute prob S = 0
  PS0 = mean(data$S == 0)
  
  # 1st clever cov FOR SDE, ALSO NEED THE ONE FOR SIE
  
  Hm_SDE = with(data, ((S == 1)*(A == 1)*
                     ((M == 1)*gstarMastar0_ps + (M == 0)*(1 - gstarMastar0_ps))*
                     ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
                     S0_preds)/
              (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                 ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                  A_ps*S1_preds*PS0) -
              ((S == 1)*(A == 0)*
                 ((M == 1)*gstarMastar0_ps + (M == 0)*(1 - gstarMastar0_ps))*
                 ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
                 S0_preds)/
              (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                 ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                 (1 - A_ps)*S1_preds*PS0))
  
  Hm_SIE = with(data, ((S == 1)*(A == 1)*
                         ((M == 1)*gstarMastar1_ps + (M == 0)*(1 - gstarMastar1_ps))*
                         ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
                         S0_preds)/
                  (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                     ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                     A_ps*S1_preds*PS0) -
                  ((S == 1)*(A == 1)*
                     ((M == 1)*gstarMastar0_ps + (M == 0)*(1 - gstarMastar0_ps))*
                     ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
                     S0_preds)/
                  (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                     ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                     A_ps*S1_preds*PS0))
  
                
  
  Hm_0_SDE = with(data, ((S == 1)*(A == 1)*
                       ((M == 1)*gstarMastar0_ps + (M == 0)*(1 - gstarMastar0_ps))*
                       ((Z == 1)*Z_S0_ps0 + (Z == 0)*(1 - Z_S0_ps0))*
                       S0_preds0)/
                (((M == 1)*M_ps0 + (M == 0)*(1 - M_ps0))*
                   ((Z == 1)*Z_ps0 + (Z == 0)*(1 - Z_ps0))*
                   ((A == 1)*A_ps0 + (A == 0)*(1 - A_ps0))*S1_preds0
                 *PS0_0) -
                ((S == 1)*(A == 0)*
                   ((M == 1)*gstarMastar0_ps + (M == 0)*(1 - gstarMastar0_ps))*
                   ((Z == 1)*Z_S0_ps0 + (Z == 0)*(1 - Z_S0_ps0))*
                   S0_preds0)/
                (((M == 1)*M_ps0 + (M == 0)*(1 - M_ps0))*
                   ((Z == 1)*Z_ps0 + (Z == 0)*(1 - Z_ps0))*
                    A_ps0*S1_preds0
                 *PS0_0))
  
  Hm_0_SIE = with(data, ((S == 1)*(A == 1)*
                         ((M == 1)*gstarMastar1_ps + (M == 0)*(1 - gstarMastar1_ps))*
                         ((Z == 1)*Z_S0_ps0 + (Z == 0)*(1 - Z_S0_ps0))*
                         S0_preds0)/
                  (((M == 1)*M_ps0 + (M == 0)*(1 - M_ps0))*
                     ((Z == 1)*Z_ps0 + (Z == 0)*(1 - Z_ps0))*
                     A_ps0*S1_preds0*PS0_0) -
                  ((S == 1)*(A == 1)*
                     ((M == 1)*gstarMastar0_ps + (M == 0)*(1 - gstarMastar0_ps))*
                     ((Z == 1)*Z_S0_ps0 + (Z == 0)*(1 - Z_S0_ps0))*
                     S0_preds0)/
                  (((M == 1)*M_ps0 + (M == 0)*(1 - M_ps0))*
                     ((Z == 1)*Z_ps0 + (Z == 0)*(1 - Z_ps0))*
                     A_ps0*S1_preds0*PS0_0))

  # for clever covariate we set S = 1 because it is not affecting the outcome, then we 
  # intervene on A and M via stoc intervention
  Hm1_SDE = with(df_YM1S1, ((S == 1)*(A == 1)*
       gstarMastar0_ps*
      ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
      S0_preds)/
      (M_S1_ps*
         ((Z == 1)*Z_S1_ps + (Z == 0)*(1 - Z_S1_ps))*
          A_S1_ps*S1_preds*PS0) - 
        ((S == 1)*(A == 0)*
           gstarMastar0_ps*
           ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
           S0_preds)/
        (M_S1_ps*
           ((Z == 1)*Z_S1_ps + (Z == 0)*(1 - Z_S1_ps))*
            (1 - A_S1_ps)*S1_preds*PS0))
  
  Hm1_SIE = with(df_YM1S1, ((S == 1)*(A == 1)*
                         gstarMastar1_ps*
                         ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
                         S0_preds)/
                  (M_S1_ps*
                     ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                     A_ps*S1_preds*PS0) -
                  ((S == 1)*(A == 1)*
                     gstarMastar0_ps*
                     ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
                     S0_preds)/
                  (M_S1_ps*
                     ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                     A_ps*S1_preds*PS0))
  
  Hm0_SDE = with(df_YM0S1, ((S == 1)*(A == 1)*
       (1 - gstarMastar0_ps)*
      ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
      S0_preds)/
      (((M == 1)*M_S1_ps + (M == 0)*(1 - M_S1_ps))*
         ((Z == 1)*Z_S1_ps + (Z == 0)*(1 - Z_S1_ps))*
         A_S1_ps*S1_preds*PS0) - 
        ((S == 1)*(A == 0)*
           (1 - gstarMastar0_ps)*
           ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
           S0_preds)/
        (((M == 1)*M_S1_ps + (M == 0)*(1 - M_S1_ps))*
           ((Z == 1)*Z_S1_ps + (Z == 0)*(1 - Z_S1_ps))*
           (1 - A_S1_ps)*S1_preds*PS0))
  
  Hm0_SIE = with(df_YM0S1, ((S == 1)*(A == 1)*
                          (1 - gstarMastar1_ps)*
                          ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
                          S0_preds)/
                   ((1 - M_S1_ps)*
                      ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                      A_ps*S1_preds*PS0) -
                   ((S == 1)*(A == 1)*
                      (1 - gstarMastar0_ps)*
                      ((Z == 1)*Z_S0_ps + (Z == 0)*(1 - Z_S0_ps))*
                      S0_preds)/
                   ((1 - M_S1_ps)*
                      ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                      A_ps*S1_preds*PS0))
  
  
  # H_clfm  
  # estimates
  Yfit = sl$train(task_YS1)
  Y_preds = Yfit$predict(task_Y)
  
  for (a in 1:100) {
  Y_preds = Qstar_M
  Y_preds_M1 = Yfit$predict(task_YM1S1)
  Y_preds_M0 = Yfit$predict(task_YM0S1)
  
  # compute the clfm clever covariate
  H = matrix(c(Hm_SDE, Hm_SIE), ncol = 2)
  D_Yinit = apply(H, 2, FUN = function(col) col*(data$Y - Y_preds))
  PnD_Yinit = colMeans(D_Yinit)
  colMeans(D_Yinit)
  norm_PnD = sqrt(sum(PnD_Yinit^2))
  H_clfm = H %*% PnD_Yinit/norm_PnD
  
  lrnr_glm = make_learner(Lrnr_glm_fast, family = binomial())
  Yfit_mle = lrnr_glm$train(task_YS1)
  Y_mle_preds = Yfit_mle$predict(task_Y)
  Y_mle_preds_M1 = Yfit_mle$predict(task_YM1S1)
  Y_mle_preds_M0 = Yfit_mle$predict(task_YM0S1)
  
  # updates
  Qfit = glm(data$Y ~ H_clfm - 1 + offset(qlogis(Y_preds)), family = binomial)
  eps = Qfit$coefficients
  
  cc_sign = sign(H_clfm)
  Qfit_wt = glm(data$Y ~ cc_sign - 1 + offset(qlogis(Y_preds)), family = binomial, 
                weights = abs(H_clfm))
  eps_wt = Qfit_wt$coefficients
  
  Qstar_M  = plogis(qlogis(Y_preds) + eps*H_clfm)
  Qstar_M_wt  = plogis(qlogis(Y_preds) + eps_wt*cc_sign)
  
  D_Ynew = apply(H, 2, FUN = function(col) col*(data$Y - Qstar_M))
  D_Ynew_wt = apply(H, 2, FUN = function(col) col*(data$Y - Qstar_M_wt))
  
  mean((data$Y - Qstar_M)*H_clfm)
  mean((data$Y - Qstar_M_wt)*H_clfm)
  
  # print(colMeans(D_Ynew_wt))
  print(colMeans(D_Ynew))
}
  
  apply(D_Ynew, 2,sd)/n
  colMeans(D_Ynew_wt)
  apply(D_Ynew_wt, 2,sd)/n
  
  Qstar_M1 = plogis(qlogis(Y_preds_M1) + eps*Hm1_SDE)
  Qstar_M0 = plogis(qlogis(Y_preds_M0) + eps*Hm0_SDE)
  
  
  # THIS IS WHERE WE ONLY HAVE COMMON LIK for SDE but not SIE
  # perform the stochastic intervention on Qstar, fcn of Z, W, S.  These are inits
  Qstar_Mg = Qstar_M1*gstarM_ps + Qstar_M0*(1 - gstarM_ps)
  Y_Mg1s = Y_preds_M1*gstarM_ps + Y_preds_M0*(1 - gstarM_ps)
  Y_Mg = Y_mle_preds_M1*gstarM_ps + Y_mle_preds_M0*(1 - gstarM_ps)
  # NEXT REGRESSION
  
  # regress on Z,W,S
  
  df_QZ = data
  df_QZ$Qstar_Mg = Qstar_Mg
  df_QZ$Y_Mg = Y_Mg  
  df_QZS0a = df_QZ
  df_QZS0a$S = 0 
  df_QZS0a$A = a
  
  task_QZ <- sl3_Task$new(
    data = data.table::copy(df_QZ),
    covariates = covariates$covariates_Z,
    outcome = "Qstar_Mg",
    outcome_type = "gaussian"
  )
  
  task_QZS0a <- sl3_Task$new(
    data = data.table::copy(df_QZS0a),
    covariates = covariates$covariates_Z,
    outcome = "Qstar_Mg",
    outcome_type = "gaussian"
  )
  
  task_QZ_Y <- sl3_Task$new(
    data = data.table::copy(df_QZ),
    covariates = covariates$covariates_Z,
    outcome = "Y_Mg",
    outcome_type = "guassian"
  )
  
  QZfit = sl$train(task_QZ)
  # QZfit = lglm$train(task_QZ)
  YZfit = lrnr_glm$train(task_QZ_Y)

  # compute the clever covariate 2
  A_preds = data$A*A_ps + (1 - data$A)*(1 - A_ps)
  A_preds0 = with(data, A*A_ps0 + (1 - A)*(1 - A_ps0))
  
  Hz = (data$A == a)*(data$S == 0)/A_preds/PS0
  if (a == 1) Hza = (data$S == 0)/A_ps/PS0 else {
    Hza = (data$S == 0)/(1-A_ps)/PS0
  }
  
  Hz_0 = (data$A == a)*(data$S == 0)/A_preds0/PS0_0
  if (a == 1) Hza_0 = (data$S == 0)/A_ps0/PS0_0 else {
    Hza_0 = (data$S == 0)/(1-A_ps0)/PS0_0
  }
  
  
  # estimates
  QZ_preds = pmin(pmax(QZfit$predict(), .001), .999)
  
  QZ_preds_a = pmin(pmax(QZfit$predict(task_QZS0a), .001), .999)
  YZ_preds_a = pmin(pmax(YZfit$predict(task_QZS0a), .001), .999)
  
  gstar_M = with(data, f_truth(Z=Z, W=W, S=S))
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
  
  # update
  QZfit_tmle = glm(Qstar_Mg ~ Hz - 1 + offset(qlogis(QZ_preds)), family = binomial)
  eps2 = QZfit_tmle$coefficients
  QZstar_a = plogis(qlogis(QZ_preds_a) + Hza*eps2)
  
  # compute the parameter estimate
  est = mean(QZstar_a[data$S==0])
  est_1s_init = mean(YZ_preds_a[data$S==0])
  est_iptw = mean(with(data, Y*Hm))
  est_mle = mean(YZ_preds_a[data$S==0])
  
  D_Y = with(data, Hm*(Y - Qstar_M))
  D_Y1s = with(data, Hm*(Y - Y_preds))
  D_0Y = with(data, Hm_0*(Y - Q_0Y))
  
  D_Z = Hz*(Qstar_Mg - QZstar_a)
  D_Z1s = Hz*(Y_Mg1s - YZ_preds_a)
  D_0Z = Hz_0*(Q_ghat_astar - Q_ghat_astarW)
  
  D_W = with(data, (QZstar_a - est)*(S ==0)/PS0)
  D_W1s = with(data, (YZ_preds_a - est_1s_init))
  D_0W = with(data, (Q_ghat_astarW - Psi_0)*(S ==0)/PS0_0)
  
  D = D_Y + D_Z + D_W
  D_1s = D_Y1s + D_Z1s + D_W1s
  D_0 = D_0Y + D_0Z + D_0W
  D_iptw = with(data, Y*Hm) - est_iptw
  est_1s = est_1s_init + mean(D_1s)
  
  n = nrow(data)
  SE = sd(D)/sqrt(n)
  SE_1s = sd(D_1s)/sqrt(n)
  SE_iptw = sd(D_iptw)/sqrt(n)
  SE_0 = sd(D_0)/sqrt(n)
  
  CI = c(est = est, left = est - 1.96*SE, right = est + 1.96*SE)
  CI_1s = c(est_1s, left = est_1s - 1.96*SE_1s, right = est_1s + 1.96*SE_1s)
  CI_iptw = c(est_iptw, left = est_iptw - 1.96*SE_iptw, right = est_iptw + 1.96*SE_iptw)
  return(list(CI = CI, CI_1s, CI_iptw = CI_iptw,
              est_mle = est_mle, IC = D, IC1s = D_1s, IC_0 = D_0, SE = SE, SE_1s = SE_1s, 
              SE_iptw = SE_iptw, SE_0 = SE_0, 
              SL_coef = list(Y = Yfit$coefficients, QZ = QZfit$coefficients),
              Psi_0 = ifelse(is.null(truth), NULL, Psi_0)))
}

