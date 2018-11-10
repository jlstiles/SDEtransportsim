
#' @export
gendata.SDEtransport = function(n, f_W, f_S, f_A, f_Z, f_M, f_Y) {
  W = f_W(n)
  W = as.data.frame(W)
  P_SW = f_S(W)
  S = rbinom(n,1,P_SW)
  
  # make a pscore model
  pscores = f_A(S=S,W=W)
  A = rbinom(n, 1, pscores)
  
  # make a intermediate confounder model
  pzscores = f_Z(A=A,S=S,W=W)
  Z = rbinom(n, 1, pzscores)
  
  # make an M model according to the restrictions
  Mscores = f_M(Z=Z,W=W,S=S)
  M = rbinom(n, 1, Mscores)
  
  # make a Y model according to the restrictions
  Yscores = f_Y(M=M,Z=Z,W=W)
  Y = rbinom(n, 1, Yscores)
  
  return(as.data.frame(cbind(W, S = S, A = A, Z = Z, M = M, Y = Y)))
}

#' @export
bound = function(x, lower, upper) {
  pmin(pmax(lower, x), upper)
}

#' @export
get.stochasticM = function(gstarM_astar, Y_preds1, Y_preds0) {
    Y_preds1*gstarM_astar + Y_preds0*(1 - gstarM_astar)
}


#' @export
gendata.SDEtransport_alt = function(n, f_W, f_A, f_Z, f_M, f_Y) {
  W = f_W(n)
  # make a pscore model
  pscores = f_A(W=W)
  A = rbinom(n, 1, pscores)
  
  # make a intermediate confounder model
  pzscores = f_Z(A=A,W=W)
  Z = rbinom(n, 1, pzscores)
  
  # make an M model according to the restrictions
  Mscores = f_M(Z=Z,W=W)
  M = rbinom(n, 1, Mscores)
  
  # make a Y model according to the restrictions
  Yscores = f_Y(M=M,Z=Z,W=W)
  Y = rbinom(n, 1, Yscores)
  
  return(as.data.frame(cbind(W, A = A, Z = Z, M = M, Y = Y)))
}
