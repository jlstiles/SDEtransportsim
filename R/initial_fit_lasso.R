#' @export
get.mediation.initdata_lasso = function(data, forms, RCT = 0.5, Wnames, Wnamesalways) {
  
  # Need to set a stock outcome for Y just to make the model matrix, in case Y is NA for S!=1.
  data$Y[data$S!=1] = data$M[data$S!=1] = 2
  df_YM1S1 = data
  df_YM1S1$M = 1
  df_YM1S1$S = 1
  
  df_YM0S1 = df_YM1S1
  df_YM0S1$M = 0
  
  df_ZS0 = data
  df_ZS0$S = 0

  
  Yform = forms$Yform
  Mform = forms$Mstarform
  Zform = forms$Zstarform
  Aform = forms$Aform
  Sform = forms$Sform
  
  # convert to model matrices to match for prediction
  df_YM1S1 = model.matrix(Yform, df_YM1S1)[,-1]
  df_YM0S1 = model.matrix(Yform, df_YM0S1)[,-1]
  
  # df_Z = data
  
  dataY = model.matrix(Yform, data)[,-1]
  #KER:
  pfac<-rep(1, ncol(dataY))
  pfac[which( colnames(dataY[,c("Z", "M", Wnames)]) %in% c("Z", "M", Wnamesalways) )]<-0
  cl<-makePSOCKcluster(4)
  registerDoParallel(cl)
  Yfit = cv.glmnet(dataY[data$S==1, ], data$Y[data$S==1], family = "binomial", 
                   weights = data$weights[data$S==1], 
                   penalty.factor=pfac, parallel=TRUE)
  
  # Mform = paste0("M ~ ", paste(covariates$covariates_M, collapse = "+"))
  dataM = model.matrix(Mform, data)[,-1]
  pfac<-rep(1, ncol(dataM))
  pfac[which( colnames(dataM[,c("Z", Wnames)]) %in% c("Z", Wnamesalways) )]<-0
  Mfit = cv.glmnet(dataM[data$S==1, ], data$M[data$S==1], family = "binomial", 
                   weights = data$weights[data$S==1], 
                   penalty.factor=pfac, parallel=TRUE)
  
  dataZ = model.matrix(Zform, data)[,-1]
  pfac<-rep(1, ncol(dataZ))
  pfac[which( colnames(dataZ[,c("A", Wnames)]) %in% c("A", Wnamesalways) )]<-0
  Zfit = cv.glmnet(dataZ, data$Z, family = "binomial", weights=data$weights, penalty.factor=pfac, 
                   parallel=TRUE)
  
  dataS = model.matrix(Sform, data)[,-1]
  Sfit = cv.glmnet(dataS, data$S, family = "binomial", weights=data$weights, parallel=TRUE)
  
  stopCluster(cl)
  # propensity scores
  if (is.null(RCT)) { 
    dataA = model.matrix(Aform, data)[,-1]
    Afit = cv.glmnet(dataA, data$A, family = "binomial")
    A_ps = predict(Afit, newx = dataA, type = 'response', s="lambda.1se")
  } else A_ps = RCT
  
  # set to predicted values and over full data set for cc is 0 when S==0 anyway.  
  M_ps = predict(Mfit, newx = dataM, type = 'response', s="lambda.1se")
  Z_ps = predict(Zfit, newx = dataZ, type = 'response', s="lambda.1se")
  
  dataZS0 = model.matrix(Zform, df_ZS0)[,-1]
  ZS0_ps = predict(Zfit, newx = dataZS0, type = 'response', s="lambda.1se")
  S_ps = predict(Sfit, newx = dataS, type = 'response', s="lambda.1se")
  PS0 = mean(data$S==0)
  
  # Predict Y for whole data, also with M = 1 and 0.  
  Y_init = pmin(pmax(predict(Yfit, newx = dataY, type = 'response', s="lambda.1se"), 0.001), .999)
  Y_init_M1 = pmin(pmax(predict(Yfit, newx = df_YM1S1, type = 'response', s="lambda.1se"), 0.001), .999)
  Y_init_M0 = pmin(pmax(predict(Yfit, newx = df_YM0S1, type = 'response', s="lambda.1se"), 0.001), .999)
  return(list(initdata = list(M_ps = M_ps, ZS0_ps = ZS0_ps, Z_ps = Z_ps, A_ps = A_ps, 
                              S_ps = S_ps, PS0 = PS0), 
              Y_preds = list(Y_init = Y_init, 
                             Y_init_M1 = Y_init_M1, 
                             Y_init_M0 = Y_init_M0)))
}


