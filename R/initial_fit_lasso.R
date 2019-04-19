#' @export
get.mediation.initdata_lasso = function(data, forms, RCT = 0.5, Wnames, Wnamesalways, 
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
  pfac<-rep(1, ncol(dataY))
  pfac[which( colnames(dataY[,c("Z", "M", Wnames)]) %in% c("Z", "M", Wnamesalways) )]<-0
  
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
  if (length(unique(data$Y))>2) Yfamily = "gaussian" else Yfamily = "binomial"
  Yfit = cv.glmnet(x, y, family = "binomial", penalty.factor=pfac, parallel=TRUE)
  
  dataM = model.matrix(Mform, data)[,-1]
  pfac<-rep(1, ncol(dataM))
  pfac[which( colnames(dataM[,c("Z", Wnames)]) %in% c("Z", Wnamesalways) )]<-0

  x = dataM
  y = data$M
  Mfit = cv.glmnet(x, y, family = "binomial", penalty.factor=pfac, parallel=TRUE)
  
  if (transport) {
    dataZ = model.matrix(Zform, data)[,-1]
    pfac<-rep(1, ncol(dataZ))
    pfac[which( colnames(dataZ[,c("A", Wnames)]) %in% c("A", Wnamesalways) )]<-0
    Zfit = cv.glmnet(dataZ, data$Z, family = "binomial", penalty.factor=pfac, parallel=TRUE)
    Z_ps = predict(Zfit, newx = dataZ, type = 'response', s="lambda.min")
    dataZS0 = model.matrix(Zform, df_ZS0)[,-1]
    ZS0_ps = predict(Zfit, newx = dataZS0, type = 'response', s="lambda.min")
    dataS = model.matrix(forms$Sform, data)[,-1]
    if (is.vector(dataS)) {
      dataS = data.frame(S = data$S, W = dataS)
      Sfit = glm(formula = S~., data = dataS, family = binomial())
      S_ps = predict(Sfit, type = 'response')
    } else {
    Sfit = cv.glmnet(dataS, data$S, family = "binomial", parallel=TRUE)
    S_ps = predict(Sfit, newx = dataS, type = 'response', s="lambda.min")
    }
    PS0 = mean(data$S==0)
   }
  stopCluster(cl)
  # propensity scores
  if (is.null(RCT)) { 
    dataA = model.matrix(Aform, data)[,-1]
    Afit = cv.glmnet(dataA, data$A, family = "binomial")
    A_ps = predict(Afit, newx = dataA, type = 'response', s="lambda.min")
  } else A_ps = RCT
  
  # as clev cov is 0 otherwise 
  M_ps = predict(Mfit, newx = dataM, type = 'response', s="lambda.min")
  
  # Predict Y for whole data, also with M = 1 and 0
  Y_init = pmin(pmax(predict(Yfit, newx = dataY, type = 'response', s="lambda.min"), 0.001), .999)
  Y_init_M1 = pmin(pmax(predict(Yfit, newx = df_YM1, type = 'response', s="lambda.min"), 0.001), .999)
  Y_init_M0 = pmin(pmax(predict(Yfit, newx = df_YM0, type = 'response', s="lambda.min"), 0.001), .999)
  
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

