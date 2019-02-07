#' @export
get.mediation.initdata_lasso_eff = function(data, forms, RCT = 0.5, Wnames, Wnamesalways, 
                                        transport, pooled) {
  
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
    Sform = forms$Sform
  }
  
  df_YM1 = data
  df_YM1$M = 1
  
  df_YM0 = df_YM1
  df_YM0$M = 0

  df_ZA0 = df_ZA1 = data
  df_ZA0$A = 0
  df_ZA1$A = 1
  
  Yform = forms$Yform
  
  # convert to model matrices to match for prediction
  data_YM1 = model.matrix(Yform, df_YM1)[,-1]
  data_YM0 = model.matrix(Yform, df_YM0)[,-1]

  data_Y = model.matrix(Yform, data)[,-1]
  pfac<-rep(1, ncol(data_Y))
  pfac[which( colnames(data_Y[,c("Z", "M", Wnames)]) %in% c("Z", "M", Wnamesalways) )]<-0
  
  cl<-makePSOCKcluster(4)
  registerDoParallel(cl)
 
   if (transport) {
    x = data_Y[data$S==1, ]
    y = data$Y[data$S==1]
    weights = data$weights[data$S==1]
  } else {
    x = data_Y
    y = data$Y
    weights = data$weights
  }

  Yfit = cv.glmnet(x, y, family = "binomial", penalty.factor=pfac, parallel=TRUE)
  
  data_M = model.matrix(Mform, data)[,-1]
  pfac<-rep(1, ncol(data_M))
  pfac[which( colnames(data_M[,c("Z", Wnames)]) %in% c("Z", Wnamesalways) )]<-0
  Mfit = cv.glmnet(data_M, data$M, family = "binomial", penalty.factor=pfac, parallel=TRUE)
  
  if (!transport) {
  data_Z = model.matrix(Zform, data)[,-1]
  pfac<-rep(1, ncol(data_Z))
  pfac[which( colnames(data_Z[,c("A", Wnames)]) %in% c("A", Wnamesalways) )]<-0
  Zfit = cv.glmnet(data_Z, data$Z, family = "binomial", penalty.factor=pfac, parallel=TRUE)
  data_ZA1 = model.matrix(Zform, df_ZA1)[,-1]
  data_ZA0 = model.matrix(Zform, df_ZA0)[,-1]
  ZA1_ps = predict(Zfit, newx = data_ZA1, type = 'response', s="lambda.min")
  ZA0_ps = predict(Zfit, newx = data_ZA0, type = 'response', s="lambda.min")
  }
  if (transport) {
    data_Z = model.matrix(Zform, data)[,-1]
    pfac<-rep(1, ncol(data_Z))
    pfac[which( colnames(data_Z[,c("A", Wnames)]) %in% c("A", Wnamesalways) )]<-0
    Zfit = cv.glmnet(data_Z, data$Z, family = "binomial", penalty.factor=pfac, parallel=TRUE)
    data_ZA1 = model.matrix(Zform, df_ZA1)[,-1]
    data_ZA0 = model.matrix(Zform, df_ZA0)[,-1]
    ZA1_ps = predict(Zfit, newx = data_ZA1, type = 'response', s="lambda.min")
    ZA0_ps = predict(Zfit, newx = data_ZA0, type = 'response', s="lambda.min")
    data_ZA1S0 = model.matrix(Zform, df_ZA1S0)[,-1]
    data_ZA0S0 = model.matrix(Zform, df_ZA0S0)[,-1]
    ZA1S0_ps = predict(Zfit, newx = data_ZA1S0, type = 'response', s="lambda.min")
    ZA0S0_ps = predict(Zfit, newx = data_ZA0S0, type = 'response', s="lambda.min")
    data_S = model.matrix(forms$Sform, data)[,-1]
    
    if (is.vector(data_S)) {
      data_S = data.frame(S = data$S, W = data_S)
      Sfit = glm(formula = S~., data = data_S, family = binomial())
      S_ps = predict(Sfit, type = 'response')
    } else {
    Sfit = cv.glmnet(data_S, data$S, family = "binomial", parallel=TRUE)
    S_ps = predict(Sfit, newx = data_S, type = 'response', s="lambda.min")
    }
    PS0 = mean(data$S==0)
   }
  stopCluster(cl)
  # propensity scores
  if (is.null(RCT)) { 
    Aform = forms$Aform
    data_A = model.matrix(Aform, data)[,-1]
    Afit = cv.glmnet(data_A, data$A, family = "binomial")
    A_ps = predict(Afit, newx = data_A, type = 'response', s="lambda.min")
  } else A_ps = RCT
  Z_ps = ZA1_ps*A_ps + ZA0_ps*(1 - A_ps)
  # as clev cov is 0 otherwise 
  M_ps = predict(Mfit, newx = data_M, type = 'response', s="lambda.min")
  
  # Predict Y for whole data, also with M = 1 and 0
  Y_init = pmin(pmax(predict(Yfit, newx = data_Y, type = 'response', s="lambda.min"), 0.001), .999)
  Y_init_M1 = pmin(pmax(predict(Yfit, newx = data_YM1, type = 'response', s="lambda.min"), 0.001), .999)
  Y_init_M0 = pmin(pmax(predict(Yfit, newx = data_YM0, type = 'response', s="lambda.min"), 0.001), .999)
  
  if (transport) {
    return(list(initdata = list(M_ps = M_ps, ZA1S0_ps = ZA1S0_ps, ZA0S0_ps = ZA0S0_ps, 
                                Z_ps = Z_ps, A_ps = A_ps, 
                                S_ps = S_ps, PS0 = PS0), 
                Y_preds = list(Y_init = Y_init, 
                               Y_init_M1 = Y_init_M1, 
                               Y_init_M0 = Y_init_M0)))
  } else {
    return(list(initdata = list(M_ps = M_ps, ZA1_ps = ZA1_ps, ZA0_ps = ZA0_ps, 
                                Z_ps = Z_ps, A_ps = A_ps), 
                Y_preds = list(Y_init = Y_init, 
                               Y_init_M1 = Y_init_M1, 
                               Y_init_M0 = Y_init_M0)))    
  }

}

