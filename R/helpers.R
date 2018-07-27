#' @export
gendata.SDEtransport = function(n, f_W, f_S, f_A, f_Z, f_M, f_Y) {
  W = f_W(n)
  P_SW = f_S(W)
  S = rbinom(n,1,P_SW)
  
  # make a pscore model
  pscores = f_A(S,W)
  A = rbinom(n, 1, pscores)
  
  # make a intermediate confounder model
  pzscores = f_Z(A,S,W)
  Z = rbinom(n, 1, pzscores)
  
  # make an M model according to the restrictions
  Mscores = f_M(Z,W,S)
  M = rbinom(n, 1, Mscores)
  
  # make a Y model according to the restrictions
  Yscores = f_Y(M,Z,W)
  Y = rbinom(n, 1, Yscores)
  
  return(data.frame(W = W, S = S, A = A, Z = Z, M = M, Y = Y))
}
