
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

#' @export
SDE_tmle = function(data, a, a_star, sl, V=10, covariates) {
  # compute the first clever covariate
  # data_ipcw
  # data_ipcw[, censored_at_t_obs := as.numeric(!Delta & t_disc == t_obs)]
  # V=2
  # dd = data
  # folds = make_folds(n, V=2)
  
  task_Mstar <- sl3_Task$new(
    data = data.table::copy(data[data$S == 1,]),
    covariates = covariates$covariates_Mstar,
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
    data = data.table::copy(data[data$S == 1,]),
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
  
  a_star = 0
  if (a_star == 1) {
    predM1 = Mstarfit$predict(task_M1)
    predM0 = Mstarfit$predict(task_M0)
    predZ = Zstarfit$predict(task_Z1)
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
    data = data.table::copy(data),
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
  
  # compute M preds 
  M_ps = Mfit$predict(task_M)
  M1a_ps = Mfit$predict(taskM_M1S1a)
  
  # compute Z preds 
  Z_ps = Zfit$predict()
  Z1a_ps = Zfit$predict(taskZ_M1S1a)
  
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
  A1_ps = Afit$predict(task_AS1a)
  
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
  Yfit = sl$train(task_Y)
  Y_preds = Yfit$predict()
  Y_preds_M1 = Yfit$predict(task_M1S1a)
  Y_preds_M0 = Yfit$predict(task_M0S1a)
  
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
  YZfit = sl$train(task_QZ_Y)
  # pp = glm(Qstar_Mg~A+W+S, data = df_QZ, family = binomial)
  # QZ_preds = predict(pp, type = 'response')
  # 
  # QZ_preds_a = predict(pp, newdata = df_QZS0a, type = 'response')
  # compute the clever covariate 2
  A_preds = data$A*A_ps + (1 - data$A)*(1 - A_ps)
  
  Hz = (data$A == a)*(data$S == 0)/A_preds/PS0 
  if (a == 1) Hza = (data$S == 0)/A_ps/PS0 else Hza = (data$S == 0)/(1-A_ps)/PS0
  
  # estimates
  QZ_preds = pmin(pmax(QZfit$predict(), .001), .999)
  
  QZ_preds_a = pmin(pmax(QZfit$predict(task_QZS0a), .001), .999)
  YZ_preds_a = pmin(pmax(YZfit$predict(task_QZ_Y), .001), .999)
  
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
  return(list(est = est, est_mle = est_mle, IC = D, 
              SL_coef = list(Y = Yfit$coefficients, QZ = QZfit$coefficients)))
}

