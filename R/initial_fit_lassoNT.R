#' @export
get.mediation.initdata_lassoNT = function(data, forms, RCT = 0.5, Wnames, Wnamesalways) {
  
  # data = cbind(data$W,A=data$A, Z=data$Z, M=data$M, Y=data$Y)
  df_YM1 = data
  df_YM1$M = 1
  
  df_YM0 = df_YM1
  df_YM0$M = 0
  Yform = forms$Yform
  Mform = forms$Mstarform
  # Zform = forms$Zstarform
  Aform = forms$Aform
  # convert to model matrices to match for prediction
  df_YM1 = model.matrix(Yform, df_YM1)[,-1]
  df_YM0 = model.matrix(Yform, df_YM0)[,-1]
  
  # df_Z = data
  
  dataY = model.matrix(Yform, data)[,-1]
  #KER:
  pfac<-rep(1, ncol(dataY))
  pfac[which( colnames(dataY[,c("Z", "M", Wnames)]) %in% c("Z", "M", Wnamesalways) )]<-0
  cl<-makePSOCKcluster(4)
  registerDoParallel(cl)
  Yfit = cv.glmnet(dataY, data$Y, family = "binomial", weights=data$weights, 
                   penalty.factor=pfac, parallel=TRUE)
  
  # Mform = paste0("M ~ ", paste(covariates$covariates_M, collapse = "+"))
  dataM = model.matrix(Mform, data)[,-1]
  pfac<-rep(1, ncol(dataM))
  pfac[which( colnames(dataM[,c("Z", Wnames)]) %in% c("Z", Wnamesalways) )]<-0
  Mfit = cv.glmnet(dataM, data$M, family = "binomial", weights=data$weights, 
                   penalty.factor=pfac, parallel=TRUE)
  
  # Zform = paste0("Z ~ ", paste(covariates$covariates_Z, collapse = "+"))
  # dataZ = model.matrix(Zform, data)[,-1]
  # pfac<-rep(1, ncol(dataZ))
  # pfac[which( colnames(dataZ[,c("A", Wnames)]) %in% c("A", Wnamesalways) )]<-0
  # Zfit = cv.glmnet(dataZ, data$Z, family = "binomial", weights=data$weights, penalty.factor=pfac, parallel=TRUE)
  
  stopCluster(cl)
  # propensity scores
  if (is.null(RCT)) { 
    dataA = model.matrix(Aform, data)[,-1]
    Afit = cv.glmnet(dataA, data$A, family = "binomial")
    A_ps = predict(Afit, newx = dataA, type = 'response', s="lambda.1se")
  } else A_ps = RCT
  
  # as clev cov is 0 otherwise 
  M_ps = predict(Mfit, newx = dataM, type = 'response', s="lambda.1se")
  
  # Predict Y for whole data, also with M = 1 and 0
  Y_init = pmin(pmax(predict(Yfit, newx = dataY, type = 'response', s="lambda.1se"), 0.001), .999)
  Y_init_M1 = pmin(pmax(predict(Yfit, newx = df_YM1, type = 'response', s="lambda.1se"), 0.001), .999)
  Y_init_M0 = pmin(pmax(predict(Yfit, newx = df_YM0, type = 'response', s="lambda.1se"), 0.001), .999)
  return(list(initdata = list(M_ps = M_ps, A_ps = A_ps), 
              Y_preds = list(Y_init = Y_init, 
                             Y_init_M1 = Y_init_M1, 
                             Y_init_M0 = Y_init_M0)))
}

