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
  
  # df_Z = data
  
  Yform = forms$Yform
  dataY = model.matrix(Yform, data)[,-1]
  Yfit = cv.glmnet(dataY, data$Y, family = "binomial")
  
  # Mform = paste0("M ~ ", paste(covariates$covariates_M, collapse = "+"))
  dataM = model.matrix(Mform, data)[,-1]
  Mfit = cv.glmnet(dataM, data$M, family = "binomial")
  
  # Zform = paste0("Z ~ ", paste(covariates$covariates_Z, collapse = "+"))
  # dataZ = model.matrix(Zform, data)[,-1]
  # Zfit = cv.glmnet(dataZ, data$Z, family = "binomial")
  
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

#' @export
bound = function(x, lower, upper) {
  pmin(pmax(lower, x), upper)
}

#' @export
get.stochasticM = function(gstarM_astar, Y_preds1, Y_preds0) {
  Y_preds1*gstarM_astar + Y_preds0*(1 - gstarM_astar)
}

