library(cateSurvival)

# Functions to generate data for transport:

# make W different for the two sites:
n=1e6
W = rnorm(n)
f_S = function(W) plogis(W +.7)
P_SW = f_S(W)
S = rbinom(n,1,P_SW)
mean(S)

# make a pscore model
f_A = function(S,W) plogis(-.2*S-.2*W+.7)
pscores = f_A(S,W)
max(pscores)
min(pscores)
hist(pscores, 200)
A = rbinom(n, 1, pscores)

# make a intermediate confounder model
f_Z = function(A,S,W) plogis(.1*S-.4*W+2*A-.8)
pzscores = f_Z(A,S,W)
hist(pzscores,200)
Z = rbinom(n, 1, pzscores)

# make an M model according to the restrictions
f_M = function(Z,W,S) plogis(-.14*S+.25*W+2*Z+.6)
Mscores = f_M(Z,W,S)
hist(Mscores, 200)
M = rbinom(n, 1, Mscores)

# make a Y model according to the restrictions
f_Y = function(M,Z,W) plogis(1*M+.25*W+1*Z-1)
Yscores = f_Y(M,Z,W)
Y = rbinom(n, 1, Yscores)
hist(Yscores, 200)

# pack these functions into a DGP
gendata = function(n, f_SW, f_A, f_Z, f_M, f_Y) {
  W = rnorm(n)
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

data_pop = gendata(5*1e6, f_SW = f_SW, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
mean(data_pop$Y)

data = gendata(1000, f_SW = f_SW, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
# Let's generate the truth!!!
# figure out the chosen formula for M

Mfit = glm(M ~ Z + W, data = data[data$S==1, ], family = binomial())
Zfit = glm(Z ~ A + W + S, data = data, family = binomial())

f_M1 = function(Z, W, S) {
  # W = W[1:10]
  nn = length(W)
  dataM1 = data.frame(Z = rep(1,nn), W = W)
  predM1 = predict(Mfit, newdata = dataM1, type = 'response')
  dataM0 = data.frame(Z = rep(0,nn), W = W)
  predM0 = predict(Mfit, newdata = dataM0, type = 'response')
  dataZ = data.frame(A = rep(1,nn), W = W, S = rep(1, nn))
  predZ = predict(Zfit, newdata = dataZ, type = 'response')
  gM = predM1*predZ + predM0*(1 - predZ)
  return(gM)
}  

f_M0 = function(Z, W, S) {
  # W = W[1:10]
  nn = length(W)
  dataM1 = data.frame(Z = rep(1,nn), W = W)
  predM1 = predict(Mfit, newdata = dataM1, type = 'response')
  dataM0 = data.frame(Z = rep(0,nn), W = W)
  predM0 = predict(Mfit, newdata = dataM0, type = 'response')
  dataZ = data.frame(A = rep(0,nn), W = W, S = rep(1, nn))
  predZ = predict(Zfit, newdata = dataZ, type = 'response')
  gM = predM1*predZ + predM0*(1 - predZ)
  return(gM)
}  
# We can now set A to 1, and stochastically intervene for A = 0 to generate M and our outcome for pop

data_popS11 = gendata(2*1e6, f_SW = f_SW, f_A = function(S,W) 1, f_Z = f_Z, f_M = f_M1, f_Y = f_Y)
head(data_popS11)
mean(data_popS11$Y)

data_popS10 = gendata(2*1e6, f_SW = f_SW, f_A = function(S,W) 1, f_Z = f_Z, f_M = f_M0, f_Y = f_Y)
head(data_popS10)
mean(data_popS10$Y)

data_popS01 = gendata(2*1e6, f_SW = f_SW, f_A = function(S,W) 0, f_Z = f_Z, f_M = f_M1, f_Y = f_Y)
head(data_popS01)
mean(data_popS01$Y)

data_popS00 = gendata(2*1e6, f_SW = f_SW, f_A = function(S,W) 0, f_Z = f_Z, f_M = f_M0, f_Y = f_Y)
head(data_popS00)
mean(data_popS00$Y)

