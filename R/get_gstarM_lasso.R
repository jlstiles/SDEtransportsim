#' @export
get_gstarM_lasso  = function(data, forms, Wnames, Wnamesalways, transport = TRUE, 
                             pooled, 
                             gstar_S = 1) 
{
  nn = nrow(data)
  Mstarform = forms$Mstarform
  Zstarform = forms$Zstarform
  
  dataZ1 = dataZ0 = dataA1 = dataA0 = data
  # Here is where we need to define the parameter as to the mechanism used for M and Z,
  # Note, if !pooled then these lines do not matter
  if (pooled) {
  dataA1$S = dataA0$S = gstar_S
  dataZ1$S = dataZ0$S = gstar_S
  }
  
  dataZ1$Z = 1
  dataZ0$Z = 0
  dataA1$A = 1
  dataA0$A = 0

  cl<-makePSOCKcluster(4)
  registerDoParallel(cl)

  dataMstar = model.matrix(Mstarform, data)[,-1]
  dataZstar = model.matrix(Zstarform, data)[,-1]
  
 # make sure the following gets forced into Mstar fit
  pfacM<-rep(1, ncol(dataMstar))
  pfacM[which(colnames(dataMstar) %in% c("Z", Wnamesalways) )]<-0
  
 # make sure the following gets forced into Zstar fit
  pfacZ<-rep(1, ncol(dataZstar))
  pfacZ[which( colnames(dataZstar) %in% c("A", Wnamesalways) )]<-0
  
  if (transport & !pooled) {
    Mstarfit = cv.glmnet(dataMstar[data$S==gstar_S, ], data$M[data$S==gstar_S], family = "binomial", 
                         penalty.factor=pfacM, parallel=TRUE)
    Zstarfit = cv.glmnet(dataZstar[data$S==gstar_S, ], data$Z[data$S==gstar_S], family = "binomial", 
                         penalty.factor=pfacZ, parallel=TRUE)
  } else {
    Mstarfit = cv.glmnet(dataMstar, data$M, family = "binomial", 
                         penalty.factor=pfacM, parallel=TRUE)
    Zstarfit = cv.glmnet(dataZstar, data$Z, family = "binomial", 
                         penalty.factor=pfacZ, parallel=TRUE)
  }

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

