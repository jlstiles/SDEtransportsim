#' @export
get_gstarM_lasso  = function(data, forms, Wnames, Wnamesalways, transport, 
                             pooled, 
                             gstar_S,
                             truth) 
{
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
  pfacZ[which(colnames(dataZstar) %in% c("A", Wnamesalways))]<-0
  
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
  
  stopCluster(cl)
  
  dataMz1 = model.matrix(Mstarform, dataZ1)[,-1]
  predMz1 = predict(Mstarfit, newx = dataMz1, type = 'response', s="lambda.min")
  
  dataMz0 = model.matrix(Mstarform, dataZ0)[,-1]
  predMz0 = predict(Mstarfit, newx = dataMz0, type = 'response', s="lambda.min")
  
  dataZa0 = model.matrix(Zstarform, dataA0)[,-1]
  predZa0 = predict(Zstarfit, newx = dataZa0, type = 'response', s="lambda.min")
  
  dataZa1 = model.matrix(Zstarform, dataA1)[,-1]
  predZa1 = predict(Zstarfit, newx = dataZa1, type = 'response', s="lambda.min")
  
  gstarM_astar0 = predMz1*predZa0 + predMz0*(1 - predZa0)
  gstarM_astar1 = predMz1*predZa1 + predMz0*(1 - predZa1)
  
  if (!is.null(truth)) {
    ftruth = function(f_W, f_S, setA, f_Z, f_Y, transport, pooled) {
      W = f_W(2e6)
      W = as.data.frame(W)
      
      if (transport) {
        P_SW = f_S(W)
        S = rbinom(2e6,1,P_SW) 
        A = setA[1]
        pzscores = f_Z(A=A,S=S,W=W)
        Z = rbinom(2e6, 1, pzscores)
        popZ1 = popZ0 = popAastar = data.frame(cbind(W, S=S,A=A,Z=Z))
      } else {
        A = setA[1]
        pzscores = f_Z(A=A,W=W)
        Z = rbinom(2e6, 1, pzscores)
        popZ1 = popZ0 = popAastar = data.frame(cbind(W,A=A,Z=Z))
      }
    
      popZ1$M = popZ0$M = 2
      
      # note, pooled is autoset to false if not transporting.  If pooled and transporting we
      # need to set S, otherwise it is autohandled by the fits being done on subset by gstar_S
      if (pooled) {
        popAastar$S = gstar_S
        popZ1$S = popZ0$S = gstar_S
      }
      
      popZ1$Z = 1
      popZ0$Z = 0
      popAastar$A = setA[2]
      
      # regardless of pooling or not the model.matrix will keep the vars according to the formula
      # so fine whether pooling or not
      popZastar = model.matrix(Zstarform, popAastar)[,-1]
      pop_predZastar = predict(Zstarfit, newx = popZastar, type = 'response', s="lambda.min")
      
      popMz1 = model.matrix(Mstarform, popZ1)[,-1]
      pop_predMz1 = predict(Mstarfit, newx = popMz1, type = 'response', s="lambda.min")
      
      popMz0 = model.matrix(Mstarform, popZ0)[,-1]
      pop_predMz0 = predict(Mstarfit, newx = popMz0, type = 'response', s="lambda.min")
      
      gstarM = pop_predMz1*pop_predZastar + pop_predMz0*(1 - pop_predZastar)
      
      M = rbinom(2e6, 1, gstarM)
      # make a Y model according to the restrictions
      Yscores = f_Y(M=M,Z=Z,W=W)
      Y = rbinom(2e6, 1, Yscores)
      if (transport) return(mean(Y[S==0])) else return(mean(Y))
    }
    SDE0 = ftruth(truth$f_W, truth$f_S, setA = c(1,0), truth$f_Z, truth$f_Y, transport, pooled) - 
      ftruth(truth$f_W, truth$f_S, setA = c(0,0), truth$f_Z, truth$f_Y, transport, pooled)
    
    SIE0 = ftruth(truth$f_W, truth$f_S, setA = c(1,1), truth$f_Z, truth$f_Y, transport, pooled) - 
      ftruth(truth$f_W, truth$f_S, setA = c(1,0), truth$f_Z, truth$f_Y, transport, pooled)
  }
  
  if (!is.null(truth)) {
    return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0,
                SDE0 = SDE0, SIE0 = SIE0))
  } else {
    return(list(gstarM_astar1 = gstarM_astar1, gstarM_astar0 = gstarM_astar0))
  }
}

