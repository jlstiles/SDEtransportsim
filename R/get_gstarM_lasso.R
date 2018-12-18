#' @export
get_gstarM_lasso  = function(data, forms, Wnames, Wnamesalways, transport = TRUE, pooledM) 
{
  # W = data$W
  # nn = nrow(W)
  nn = nrow(data)
  Mstarform = forms$Mstarform
  Zstarform = forms$Zstarform
  # fit M for S = 1
  # data = cbind(W,A=data$A, Z=data$Z, M=data$M, Y=data$Y)
  # in case data$M = NA for S=1, we make sure model.matrix works on full data set
  if (transport & !pooledM) data$M[data$S==0] = 2 
  dataZ1 = dataZ0 = dataA1 = dataA0 = data
  dataZ1$Z = 1
  dataZ0$Z = 0
  dataA1$A = 1
  dataA0$A = 0
  
  #KER
  #update to add any more covariates
  #KER I actually think that Wnames should be Wnames should change to be unique per model
  cl<-makePSOCKcluster(4)
  registerDoParallel(cl)
  #dataMstar=data.matrix(data[,c("Z", Wnames)])
  #
  dataMstar = model.matrix(Mstarform, data)[,-1]

  pfac<-rep(1, ncol(dataMstar))
  pfac[which(colnames(dataMstar) %in% c("Z", Wnamesalways) )]<-0
  if (transport & !pooledM) {
    Mstarfit = cv.glmnet(dataMstar[data$S==1, ], data$M[data$S==1], family = "binomial", 
                         penalty.factor=pfac, parallel=TRUE)
  } else {
    Mstarfit = cv.glmnet(dataMstar, data$M, family = "binomial", 
                         penalty.factor=pfac, parallel=TRUE)
  }

  
  #KER note to self that data$M is updated to be the correct M based on line 87 of illustrativeanalysisfunction.R
  #Mstarfit = cv.glmnet(dataMstar, data$M, family = "binomial")
  
  #dataZstar=data.matrix(data[,c("A", Wnames)])
  dataZstar = model.matrix(Zstarform, data)[,-1]
  #Ker
  pfac<-rep(1, ncol(dataZstar))
  pfac[which( colnames(dataZstar) %in% c("A", Wnamesalways) )]<-0
  Zstarfit = cv.glmnet(dataZstar, data$Z, family = "binomial", 
                       penalty.factor=pfac, parallel=TRUE)
  #Zstarfit = cv.glmnet(dataZstar, data$Z, family = "binomial")
  
  dataMz1 = model.matrix(Mstarform, dataZ1)[,-1]
  predMz1 = predict(Mstarfit, newx = dataMz1, type = 'response', s="lambda.1se")
  
  dataMz0 = model.matrix(Mstarform, dataZ0)[,-1]
  predMz0 = predict(Mstarfit, newx = dataMz0, type = 'response', s="lambda.1se")
  
  dataZa0 = model.matrix(Zstarform, dataA0)[,-1]
  predZa0 = predict(Zstarfit, newx = dataZa0, type = 'response', s="lambda.1se")
  
  dataZa1 = model.matrix(Zstarform, dataA1)[,-1]
  predZa1 = predict(Zstarfit, newx = dataZa1, type = 'response', s="lambda.1se")
  
  gstarM_astar0 = predMz1*predZa0 + predMz0*(1 - predZa0)
  gstarM_astar1 = predMz1*predZa1 + predMz0*(1 - predZa1)
  
  stopCluster(cl)
  
  return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0))
}

