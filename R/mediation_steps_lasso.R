#' @export
mediation.step1_lasso = function(initdata, Y_preds, data, gstarM_astar, a, transport) {
  
  if (!transport) {
    H = with(data, with(initdata, ((A == a)*data$weights/sum(data$weights)*((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))/
                                     (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*(A_ps*A + (1 - A)*(1 - A_ps))))))
  } else {
    H = with(data, with(initdata, ((S == 1)*(A == a)*data$weights/sum(data$weights[data$S==1])*
                                     ((M == 1)*gstarM_astar + (M == 0)*(1 - gstarM_astar))*
                                     ((Z == 1)*ZS0_ps + (Z == 0)*(1 - ZS0_ps))*(1 - S_ps))/
                          (((M == 1)*M_ps + (M == 0)*(1 - M_ps))*
                             ((Z == 1)*Z_ps + (Z == 0)*(1 - Z_ps))*
                             (A_ps*A + (1 - A)*(1 - A_ps))*S_ps*PS0)))
  }
  
  # updates
  Qfit = try(suppressWarnings(glm(data$Y ~ 1 + offset(qlogis(Y_preds$Y_init)), family = binomial,
                                  weights = H)), silent = TRUE)
  
  if (class(Qfit)[1]=="try-error") eps = 0 else eps = Qfit$coefficients
  
  est_iptw = mean(data$Y*H)
  IC_iptw = data$Y*H - est_iptw
  return(list(Qstar_M  = plogis(qlogis(Y_preds$Y_init) + eps),
              Qstar_M1 = plogis(qlogis(Y_preds$Y_init_M1) + eps),
              Qstar_M0 = plogis(qlogis(Y_preds$Y_init_M0) + eps),
              IC_iptw = IC_iptw, est_iptw = est_iptw,
              Hm = H,
              eps = eps))
}


#' @export
mediation.step2_lasso = function(data, Qstar_M, Qstar_Mg, Hm, A_ps, a, tmle = TRUE,
                                 EE = FALSE, bootstrap = FALSE, form, transport) {
  
  if (transport) PS0 = mean(data$S==0)
  QZform = form
  Qstar_Mg = as.vector(Qstar_Mg)
  df = cbind(data, Qstar_Mg = Qstar_Mg)
  
  cl<-makePSOCKcluster(4)
  registerDoParallel(cl)
  # QZform = formula(paste0("Qstar_Mg ~ ", paste(covariates$covariates_QZ, collapse = "+")))
  df_QZ = model.matrix(QZform, df)[,-1]
  QZfit = cv.glmnet(df_QZ[data$A==a, ], Qstar_Mg[data$A==a], family = "gaussian", parallel=TRUE)
  stopCluster(cl)
  # compute the clever covariate and update if tmle
  if (!transport) {
    norm_wts = sum(data$weights)
    Hz = data$weights/norm_wts*with(data, (A == a)/(A*A_ps + (1 - A)*(1 - A_ps)))
  } else {
    norm_wts = sum(data$weights[data$S==0])
    Hz = data$weights/norm_wts*with(data, (A == a)*(S == 0)/(A*A_ps + (1 - A)*(1 - A_ps))/PS0)  
  }
  
  if (tmle) {
    QZ_preds_a = pmin(pmax(predict(QZfit, newx = df_QZ,  s="lambda.1se"), .001), .999)
    # update
    QZfit_tmle = try(glm(Qstar_Mg ~ 1 + offset(qlogis(QZ_preds_a)), family = binomial,
                         weights = Hz), silent = TRUE)
    if (class(QZfit_tmle)=="try-error") eps2 = 0 else eps2 = QZfit_tmle$coefficients
    
    QZstar_a = plogis(qlogis(QZ_preds_a) + eps2)
    if (transport) {
      est = mean(sum(data$S==0)*(QZstar_a*data$weights/norm_wts)*(data$S==0)/PS0)
    } else {
      est = mean(nrow(data)*QZstar_a*data$weights/norm_wts)
    }
    
    if(!bootstrap) { 
      D_Y = with(data, Hm*(Y - Qstar_M))
      D_Z = Hz*(Qstar_Mg - QZstar_a)
      if (transport) {
        wts = data$weights/sum(data$weights[data$S==0])
        D_W = with(data, (sum(data$S==0)*QZstar_a*data$weights/norm_wts - est)*(S ==0)/PS0)
      } else {
        wts = data$weights/sum(data$weights)
        D_W = with(data, (nrow(data)*QZstar_a*data$weights/norm_wts - est))
      }
      D = D_Y + D_Z + D_W
    }
  } 
  
  # regress if EE or mle, EE updates the estimate, mle does not
  if (EE) {
    QZstar_a = pmin(pmax(predict(QZfit, newx = df_QZ, s="lambda.1se"), .001), .999) 
    if (transport) {
      init_est = mean(sum(data$S==0)*(QZstar_a*data$weights/norm_wts)*(data$S==0)/PS0)
    } else {
      init_est = mean(nrow(data)*QZstar_a*data$weights/norm_wts)
    }
    D_Y1s = with(data, Hm*(Y - Qstar_M))
    D_Z1s = Hz*(Qstar_Mg - QZstar_a)
    if (transport) {
      D_W1s = with(data, (sum(data$S==0)*QZstar_a*data$weights/norm_wts - init_est)*(S ==0)/PS0)
    } else {
      D_W1s = with(data, (QZstar_a*data$weights - init_est*norm_wts)/norm_wts)
    }
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
