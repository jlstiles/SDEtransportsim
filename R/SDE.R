
#' @export
fcn.SDE = function(data, 
                   f_M0, 
                   f_M1, 
                   Mstarfit, 
                   Zstarfit, 
                   Yfit,
                   Mfit,
                   Zfit,
                   Afit,
                   Sfit,
                   a, a_star) {
  # compute the first clever covariate
  
  if (a_star == 1) gstarM_ps = f_M1(Z=data$Z, W = data$W, S = data$S) else
    gstarM_ps = f_M0(Z=data$Z, W = data$W, S = data$S)
  
  # compute Z preds for S = 0
  data_M1S1a = data_M0S1a = data
  data_M1S1a$M = 1
  data_M1S1a$A = a
  data_M1S1a$S = 1
  
  data_M0S1a$M = 0
  data_M0S1a$A = a
  data_M0S1a$S = 1
  
  data_S0a = data
  data_S0a$S = 0
  data_S0a$A = a
  
  Z0a_ps = predict(Zfit,newdata = data_S0a, type = 'response')
  
  # compute probs S = 0 given W
  S0_preds = 1 - predict(Sfit, type = 'response')
  
  # compute M preds 
  M_ps = predict(Mfit, newdata = data, type = 'response')
  M1a_ps = predict(Mfit, newdata = data_M1S1a, type = 'response')
  
  # compute Z preds 
  Z_ps = predict(Zfit, type = 'response')
  Z1a_ps = predict(Zfit, newdata = data_M1S1a, type = 'response')
  
  # compute A=1 preds for S = 1
  A_ps = predict(Afit, type = 'response')
  A1_ps = predict(Afit, newdata = data_M1S1a, type = 'response')
  
  #compute prob S = 1 given W
  S1_preds = 1 - S0_preds
  
  # compute prob S = 0
  PS0 = mean(data$S == 0)
  
  # 1st clever cov
  
  Hm = with(data, ((S == 1)*(A == a)*
                     ((M == 1)*gstarM_ps + (M == 0)*(1 - gstarM_ps))*
                     ((Z == 1)*Z0a_ps + (Z == 0)*(1 - Z0a_ps))*
                     S0_preds)/
              (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                 ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                 ((A == 1)*A_ps + (A == 0)*(1 - A_ps))*S1_preds*PS0))
  
  # for clever covariate we set S = 1 because it is not affecting the outcome, then we 
  # intervene on A and M via stoc intervention
  Hm1 = with(data_M1S1a, (
    (S == 1)*(A == a)*
      ((M == 1)*gstarM_ps + (M == 0)*(1 - gstarM_ps))*
      ((Z == 1)*Z0a_ps + (Z == 0)*(1 - Z0a_ps))*
      S0_preds)/
      (((M == 1)*M1a_ps + (M == 0)*(1 - M1a_ps))*
         ((Z == 1)*Z1a_ps + (Z == 0)*(1 - Z1a_ps))*
         ((A == 1)*A1_ps + (A == 0)*(1 - A1_ps))*S1_preds*PS0))
  
  Hm0 = with(data_M0S1a, (
    (S == 1)*(A == a)*
      ((M == 1)*gstarM_ps + (M == 0)*(1 - gstarM_ps))*
      ((Z == 1)*Z0a_ps + (Z == 0)*(1 - Z0a_ps))*
      S0_preds)/
      (((M == 1)*M1a_ps + (M == 0)*(1 - M1a_ps))*
         ((Z == 1)*Z1a_ps + (Z == 0)*(1 - Z1a_ps))*
         ((A == 1)*A1_ps + (A == 0)*(1 - A1_ps))*S1_preds*PS0))
  
  # estimates
  Y_preds = predict(Yfit, newdata = data, type = 'response')
  Y_preds_M1 = predict(Yfit, newdata = data_M1S1a, type = 'response')
  Y_preds_M0 = predict(Yfit, newdata = data_M0S1a, type = 'response')
  
  # updates
  Qfit = glm(data$Y ~ Hm - 1 + offset(qlogis(Y_preds)), family = binomial)
  eps = Qfit$coefficients
  
  Qstar_M  = plogis(qlogis(Y_preds) + eps*Hm)
  Qstar_M1 = plogis(qlogis(Y_preds_M1) + eps*Hm1)
  Qstar_M0 = plogis(qlogis(Y_preds_M0) + eps*Hm0)
  
  # perform the stochastic intervention on Qstar, fcn of Z, W, S.  These are inits
  Qstar_Mg = Qstar_M1*gstarM_ps + Qstar_M0*(1 - gstarM_ps)
  Y_Mg = Y_preds_M1*gstarM_ps + Y_preds_M0*(1 - gstarM_ps)
  # NEXT REGRESSION
  
  # regress on Z,W,S
  QZfit = glm(Qstar_Mg ~ A + W + S, data = data, family = 'binomial')
  YZfit = glm(Y_Mg ~ A + W + S, data = data, family = 'binomial')
  
  # compute the clever covariate 2
  A_preds = data$A*A_ps + (1 - data$A)*(1 - A_ps)
  
  Hz = (data$A == a)*(data$S == 0)/A_preds/PS0 
  if (a == 1) Hza = (data$S == 0)/A_ps/PS0 else Hza = (data$S == 0)/(1-A_ps)/PS0
  
  # estimates
  QZ_preds = predict(QZfit, newdata = data, type = 'response')
  QZ_preds_a = predict(QZfit, newdata = data_S0a, type = 'response')
  
  YZ_preds_a = predict(YZfit, newdata = data_S0a, type = 'response')
  
  # update
  QZfit_tmle = glm(Qstar_Mg ~ Hz - 1 + offset(qlogis(QZ_preds)), family = binomial)
  eps2 = QZfit_tmle$coefficients
  QZstar_a = plogis(qlogis(QZ_preds_a) + Hza*eps2)
  
  # compute the parameter estimate
  est = mean(QZstar_a[data$S==0])
  est_mle = mean(YZ_preds_a[data$S==0])
  
  D_Y = with(data, Hm*(Y - Qstar_M))
  D_Z = Hz*(Qstar_Mg - QZstar_a)
  D_W = with(data, (QZstar_a - est)*(S ==0)/PS0)
  
  D = D_Y + D_Z + D_W
  return(list(est = est, est_mle = est_mle, IC = D))
}


#' @title SDE_tmle
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
SDE_tmle = function(data, a, a_star, sl, V=10, covariates, truth = NULL, 
                    truncate = list(lower =.0001, upper = .9999), glm_only = TRUE) {
  # compute the first clever covariate
  # data_ipcw
  # data_ipcw[, censored_at_t_obs := as.numeric(!Delta & t_disc == t_obs)]
  # V=2
  # dd = data
  # folds = make_folds(n, V=2)
  if (glm_only) sl = make_learner(Lrnr_glm_fast, family = binomial())
  
  L = truncate$lower
  U = truncate$upper
  
  task_Mstar <- sl3_Task$new(
    data = data.table::copy(data[data$S == 1,]),
    covariates = covariates$covariates_M,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "M",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  df_M1 = df_M0 = data
  df_M1$Z = 1
  df_M0$Z = 0
  
  task_M1 <- sl3_Task$new(
    data = data.table::copy(df_M1),
    covariates = covariates$covariates_M,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "M",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  task_M0 <- sl3_Task$new(
    data = data.table::copy(df_M0),
    covariates = covariates$covariates_M,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "M",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  
  
  task_Zstar <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Z,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Z",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  df_Z1 = df_Z0 = data
  df_Z1$A = 1
  df_Z0$A = 0
  
  task_Z1 <- sl3_Task$new(
    data = data.table::copy(df_Z1),
    covariates = covariates$covariates_Z,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Z",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  task_Z0 <- sl3_Task$new(
    data = data.table::copy(df_Z0),
    covariates = covariates$covariates_Z,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Z",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  Mstarfit = sl$train(task_Mstar)
  Zstarfit = sl$train(task_Zstar)
  
  # find the truth if simulating
  if (!is.null(truth)) {
  f_truth = function(Z, W, S) {
    # W = W[1:10]
    nn = length(W)
    dataM1 = data.frame(Z = rep(1,nn), W = W, M = rep(1,nn))
    taskM1 =   sl3_Task$new(
      data = data.table::copy(dataM1),
      covariates = covariates$covariates_M,
      outcome = "M",
      outcome_type = "binomial"
    )
    predM1 = Mstarfit$predict(taskM1)
    
    dataM0 = data.frame(Z = rep(0,nn), W = W, M = rep(1,nn))
    taskM0 =   sl3_Task$new(
      data = data.table::copy(dataM0),
      covariates = covariates$covariates_M,
      outcome = "M",
      outcome_type = "binomial"
    )
    
    predM0 = Mstarfit$predict(taskM0)
    
    dataZ = data.frame(A = rep(a_star,nn), W = W, S = rep(1, nn), Z = rep(1,nn))
    
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
  
  data_pop = gendata.SDEtransport(5*1e6, f_W = truth$f_W, f_S = truth$f_S, f_A = function(S,W) a, 
                                  f_Z = truth$f_Z, f_M = f_truth, f_Y = truth$f_Y)
  Psi_0 = mean(data_pop$Y[data_pop$S==0])
  PS0_0 = mean(data_pop$S)   
  }

  if (a_star == 1) {
    predM1 = Mstarfit$predict(task_M1)
    predM0 = Mstarfit$predict(task_M0)
    predZ =  Zstarfit$predict(task_Z1)
    gstarM_ps = predM1*predZ + predM0*(1 - predZ)
  } else {
    predM1 = Mstarfit$predict(task_M1)
    predM0 = Mstarfit$predict(task_M0)
    predZ = Zstarfit$predict(task_Z0)
    gstarM_ps = predM1*predZ + predM0*(1 - predZ)
  }
  
  
  
  # compute Z preds for S = 0
  df_M1S1a = df_M0S1a = data
  df_M1S1a$M = 1
  df_M1S1a$A = a
  df_M1S1a$S = 1
  
  df_M0S1a$M = 0
  df_M0S1a$A = a
  df_M0S1a$S = 1
  
  df_S0a = data
  df_S0a$S = 0
  df_S0a$A = a
  
  
  task_YS1 <- sl3_Task$new(
    data = data.table::copy(data[data$S==1,]),
    covariates = covariates$covariates_Y,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Y",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  task_Y <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Y,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Y",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  task_M1S1a <- sl3_Task$new(
    data = data.table::copy(df_M1S1a),
    covariates = covariates$covariates_Y,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Y",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  task_M0S1a <- sl3_Task$new(
    data = data.table::copy(df_M0S1a),
    covariates = covariates$covariates_Y,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Y",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  task_S0a <- sl3_Task$new(
    data = data.table::copy(df_S0a),
    covariates = covariates$covariates_Z,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Z",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  task_M <- sl3_Task$new(
    data = data.table::copy(data[data$S == 1,]),
    covariates = covariates$covariates_M,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "M",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  taskM_M1S1a <- sl3_Task$new(
    data = data.table::copy(df_M1S1a),
    covariates = covariates$covariates_M,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "M",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  taskM_M0S1a <- sl3_Task$new(
    data = data.table::copy(df_M0S1a),
    covariates = covariates$covariates_M,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "M",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  task_Z <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_Z,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Z",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  taskZ_M1S1a <- sl3_Task$new(
    data = data.table::copy(df_M1S1a),
    covariates = covariates$covariates_Z,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Z",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  taskZ_M0S1a <- sl3_Task$new(
    data = data.table::copy(df_M0S1a),
    covariates = covariates$covariates_Z,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Z",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  Mfit = sl$train(task_M)
  Zfit = sl$train(task_Z)
  
  Z0a_ps = Zfit$predict(task_S0a)
  Z0a_ps0 = with(df_S0a, truth$f_Z(A=A, W=W, S=S))
  
  task_S <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_S,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "S",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  # compute probs S = 0 given W
  Sfit = sl$train(task_S)
  S0_preds = 1 - Sfit$predict()
  S0_preds0 = 1 - with(data, truth$f_S(W=W))
  
  # compute M preds
  M_ps = rep(10, nrow(data))
  M_ps[data$S==1] = Mfit$predict(task_M)
  M1a_ps = Mfit$predict(taskM_M1S1a)
  M_ps0 = with(data, truth$f_M(Z=Z, W=W, S=1))
    
  # compute Z preds 
  Z_ps = Zfit$predict()
  Z_ps0 = with(data, f_Z(A=A, W=W, S=S))
  Z1a_ps = Zfit$predict(taskZ_M1S1a)
  Z1a_ps0 = with(df_M1S1a, truth$f_Z(A=A, W=W, S=S))
  
  # compute A=1 preds for S = 1
  
  task_A <- sl3_Task$new(
    data = data.table::copy(data),
    covariates = covariates$covariates_A,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "A",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  df_S1a = data
  df_S1a$S = 1
  df_S1a$A = a
  
  task_AS1a <- sl3_Task$new(
    data = data.table::copy(df_S1a),
    covariates = covariates$covariates_A,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "A",
    # id = "ID",
    # folds = folds,
    outcome_type = "binomial"
  )
  
  Afit = sl$train(task_A)
  A_ps = Afit$predict()
  A_ps0 = with(data, truth$f_A(W=W, S=S))
  
  A1_ps = Afit$predict(task_AS1a)
  A1_ps0 = with(df_S1a, truth$f_A(W=W, S=S))
  
  #compute prob S = 1 given W
  S1_preds = 1 - S0_preds
  S1_preds0 = 1 - S0_preds0
  
  # compute prob S = 0
  PS0 = mean(data$S == 0)
  
  # 1st clever cov
  
  Hm = with(data, ((S == 1)*(A == a)*
                     ((M == 1)*gstarM_ps + (M == 0)*(1 - gstarM_ps))*
                     ((Z == 1)*Z0a_ps + (Z == 0)*(1 - Z0a_ps))*
                     S0_preds)/
              (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                 ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                 ((A == 1)*A_ps + (A == 0)*(1 - A_ps))*S1_preds*PS0))
  
  Hm_0 = with(data, ((S == 1)*(A == a)*
                     ((M == 1)*gstarM_ps + (M == 0)*(1 - gstarM_ps))*
                     ((Z == 1)*Z0a_ps0 + (Z == 0)*(1 - Z0a_ps0))*
                     S0_preds0)/
              (((M == 1)*M_ps0 + (M == 0)*(1 - M_ps0))*
                 ((Z == 1)*Z_ps0 + (Z == 0)*(1 - Z_ps0))*
                 ((A == 1)*A_ps0 + (A == 0)*(1 - A_ps0))*S1_preds*PS0_0))
  
  # computing true clever covariate

  
  # for clever covariate we set S = 1 because it is not affecting the outcome, then we 
  # intervene on A and M via stoc intervention
  Hm1 = with(df_M1S1a, (
    (S == 1)*(A == a)*
      ((M == 1)*gstarM_ps + (M == 0)*(1 - gstarM_ps))*
      ((Z == 1)*Z0a_ps + (Z == 0)*(1 - Z0a_ps))*
      S0_preds)/
         (((M == 1)*M1a_ps + (M == 0)*(1 - M1a_ps))*
         ((Z == 1)*Z1a_ps + (Z == 0)*(1 - Z1a_ps))*
         ((A == 1)*A1_ps + (A == 0)*(1 - A1_ps))*S1_preds*PS0))
  
  Hm0 = with(df_M0S1a, (
    (S == 1)*(A == a)*
      ((M == 1)*gstarM_ps + (M == 0)*(1 - gstarM_ps))*
      ((Z == 1)*Z0a_ps + (Z == 0)*(1 - Z0a_ps))*
      S0_preds)/
      (((M == 1)*M1a_ps + (M == 0)*(1 - M1a_ps))*
         ((Z == 1)*Z1a_ps + (Z == 0)*(1 - Z1a_ps))*
         ((A == 1)*A1_ps + (A == 0)*(1 - A1_ps))*S1_preds*PS0))
  
  # estimates
  Yfit = sl$train(task_YS1)
  Y_preds = Yfit$predict(task_Y)
  Y_preds_M1 = Yfit$predict(task_M1S1a)
  Y_preds_M0 = Yfit$predict(task_M0S1a)
  
  lrnr_glm = make_learner(Lrnr_glm_fast, family = binomial())
  Yfit_mle = lrnr_glm$train(task_YS1)
  Y_mle_preds = Yfit_mle$predict(task_Y)
  Y_mle_preds_M1 = Yfit_mle$predict(task_M1S1a)
  Y_mle_preds_M0 = Yfit_mle$predict(task_M0S1a)
  
  # updates
  Qfit = glm(data$Y ~ Hm - 1 + offset(qlogis(Y_preds)), family = binomial)
  eps = Qfit$coefficients
  
  Qstar_M  = plogis(qlogis(Y_preds) + eps*Hm)
  Qstar_M1 = plogis(qlogis(Y_preds_M1) + eps*Hm1)
  Qstar_M0 = plogis(qlogis(Y_preds_M0) + eps*Hm0)
  
  # perform the stochastic intervention on Qstar, fcn of Z, W, S.  These are inits
  Qstar_Mg = Qstar_M1*gstarM_ps + Qstar_M0*(1 - gstarM_ps)
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
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Qstar_Mg",
    # id = "ID",
    # folds = folds,
    outcome_type = "gaussian"
  )
  
  task_QZS0a <- sl3_Task$new(
    data = data.table::copy(df_QZS0a),
    covariates = covariates$covariates_Z,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Qstar_Mg",
    # id = "ID",
    # folds = folds,
    outcome_type = "gaussian"
  )
  
  task_QZ_Y <- sl3_Task$new(
    data = data.table::copy(df_QZ),
    covariates = covariates$covariates_Z,
    # covariates = colnames(dt_in)[str_detect(colnames(dt_in), "W")],
    outcome = "Y_Mg",
    # id = "ID",
    # folds = folds,
    outcome_type = "guassian"
  )
  
  QZfit = sl$train(task_QZ)
  # QZfit = lglm$train(task_QZ)
  YZfit = lrnr_glm$train(task_QZ_Y)
  # pp = glm(Qstar_Mg~A+W+S, data = df_QZ, family = binomial)
  # QZ_preds = predict(pp, type = 'response')
  # 
  # QZ_preds_a = predict(pp, newdata = df_QZS0a, type = 'response')
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
  YZ_preds_a = pmin(pmax(YZfit$predict(task_QZ_Y), .001), .999)
  
  gstar_M = with(data, f_truth(Z=Z, W=W, S=S))
  p_0Z = with(data, truth$f_Z(A=a, W=W, S=S))
  Q_0Y = with(data, truth$f_Y(M=M, Z=Z, W=W))
  
  Q_ghat_astar = Q_0Y*gstar_M
  Q_ghat_astarW = Q_ghat_astar*p_0Z
  
  # update
  QZfit_tmle = glm(Qstar_Mg ~ Hz - 1 + offset(qlogis(QZ_preds)), family = binomial)
  eps2 = QZfit_tmle$coefficients
  QZstar_a = plogis(qlogis(QZ_preds_a) + Hza*eps2)
  
  # compute the parameter estimate
  est = mean(QZstar_a[data$S==0])
  est_mle = mean(YZ_preds_a[data$S==0])
  
  D_Y = with(data, Hm*(Y - Qstar_M))
  D_0Y = with(data, Hm_0*(Y - Q_0Y))
  
  D_Z = Hz*(Qstar_Mg - QZstar_a)
  D_0Z = Hz_0*(Q_ghat_astar - Q_ghat_astarW)
  
  D_W = with(data, (QZstar_a - est)*(S ==0)/PS0)
  D_0W = with(data, (Q_ghat_astarW - Psi_0)*(S ==0)/PS0_0)
  
  D = D_Y + D_Z + D_W
  D_0 = D_0Y + D_0Z + D_0W
  
  n = nrow(data)
  SE = sd(D)/sqrt(n)
  SE_0 = sd(D_0)/sqrt(n)
  
  CI = c(est = est, left = est - 1.96*SE, right = est + 1.96*SE)
  
  return(list(CI = CI, est_mle = est_mle, IC = D, IC_0 = D_0, SE = SE, SE_0 = SE_0, 
              SL_coef = list(Y = Yfit$coefficients, QZ = QZfit$coefficients),
              Psi_0 = ifelse(is.null(truth), NULL, Psi_0)))
}

