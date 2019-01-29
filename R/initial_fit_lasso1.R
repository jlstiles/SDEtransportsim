#' @export
get.mediation.initdata_lasso1 = function(data, forms, RCT = 0.5, Wnames, Wnamesalways, 
                                        transport, pooled) {
  
  # data = cbind(data$W,A=data$A, Z=data$Z, M=data$M, Y=data$Y)
  if (transport) { 
    # This is to avoid NA's being a pain, we won't use these outcomes
    data$Y[data$S!=1] = 2
    
    df_ZS0 = data
    df_ZS0$S = 0
    if (pooled) Zform = forms$Zstarform else Zform = forms$Zform
    Sform = forms$Sform
  }
  
  df_YM1 = data
  df_YM1$M = 1
  
  df_YM0 = df_YM1
  df_YM0$M = 0

  Yform = forms$Yform
  if (pooled | !transport) Mform = forms$Mstarform else Mform = forms$Mform
  Aform = forms$Aform
  
  # convert to model matrices to match for prediction
  df_YM1 = model.matrix(Yform, df_YM1)[,-1]
  df_YM0 = model.matrix(Yform, df_YM0)[,-1]

  dataY = model.matrix(Yform, data)[,-1]
  
  cl<-makePSOCKcluster(4)
  registerDoParallel(cl)
 
   if (transport) {
    x = dataY[data$S==1, ]
    y = data$Y[data$S==1]
    weights = data$weights[data$S==1]
  } else {
    x = dataY
    y = data$Y
    weights = data$weights
  }

  Yfit = glmnet(x, y, family = "binomial")
  
  dataM = model.matrix(Mform, data)[,-1]

  x = dataM
  y = data$M
  Mfit = glmnet(x, y, family = "binomial")
  
  if (transport) {
    dataZ = model.matrix(Zform, data)[,-1]
    Zfit = glmnet(dataZ, data$Z, family = "binomial")
    Z_ps = predict(Zfit, newx = dataZ, type = 'response')
    dataZS0 = model.matrix(Zform, df_ZS0)[,-1]
    ZS0_ps = predict(Zfit, newx = dataZS0, type = 'response')
    dataS = model.matrix(Sform, data)[,-1]
    if (is.vector(dataS)) {
      dataS = data.frame(S = data$S, W = dataS)
      Sfit = glm(formula = S~., data = dataS, family = binomial())
      S_ps = predict(Sfit, type = 'response')
    } else {
    Sfit = glmnet(dataS, data$S, family = "binomial")
    S_ps = predict(Sfit, newx = dataS, type = 'response')
    }
    PS0 = mean(data$S==0)
   }
  stopCluster(cl)
  # propensity scores
  if (is.null(RCT)) { 
    dataA = model.matrix(Aform, data)[,-1]
    Afit = glmnet(dataA, data$A, family = "binomial")
    A_ps = predict(Afit, newx = dataA, type = 'response')
  } else A_ps = RCT
  
  # as clev cov is 0 otherwise 
  M_ps = predict(Mfit, newx = dataM, type = 'response')
  
  # Predict Y for whole data, also with M = 1 and 0
  Y_init = pmin(pmax(predict(Yfit, newx = dataY, type = 'response'), 0.001), .999)
  Y_init_M1 = pmin(pmax(predict(Yfit, newx = df_YM1, type = 'response'), 0.001), .999)
  Y_init_M0 = pmin(pmax(predict(Yfit, newx = df_YM0, type = 'response'), 0.001), .999)
  
  if (transport) {
    return(list(initdata = list(M_ps = M_ps, ZS0_ps = ZS0_ps, Z_ps = Z_ps, A_ps = A_ps, 
                                S_ps = S_ps, PS0 = PS0), 
                Y_preds = list(Y_init = Y_init, 
                               Y_init_M1 = Y_init_M1, 
                               Y_init_M0 = Y_init_M0)))
  } else {
    return(list(initdata = list(M_ps = M_ps, A_ps = A_ps), 
                Y_preds = list(Y_init = Y_init, 
                               Y_init_M1 = Y_init_M1, 
                               Y_init_M0 = Y_init_M0)))    
  }

}

