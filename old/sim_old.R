# library(cateSurvival)
library(Simulations)
library(SDEtransport)
# Functions to generate data for transport:

# make W different for the two sites:
n=1e6

f_W = function(n) rnorm(n)
W = f_W(n)
f_S = function(W) plogis(W +.7)
P_SW = f_S(W)
S = rbinom(n,1,P_SW)
mean(S)
hist(P_SW, breaks = 200)
# make a pscore model

f_A = function(S,W) plogis(-.7*S-.9*W +.7)

f_A = function(S,W) plogis(-.6*S-.7*W +.17)
pscores = f_A(S,W)
max(pscores)
min(pscores)
hist(pscores, 200)
A = rbinom(n, 1, pscores)

# make a intermediate confounder model

f_Z = function(A,S,W) plogis(.1*S-1*W +.72*A -.8)

f_Z = function(A,S,W) plogis(.1*S-.4*W+1*A-.3)
pzscores = f_Z(A,S,W)
hist(pzscores,200)
Z = rbinom(n, 1, pzscores)

# make an M model according to the restrictions

# make an M model according to the restrictions
f_M = function(Z,W,S) plogis(-.14*S + 1*W + 1.2*Z +.1)
Mscores = f_M(Z,W,S)
hist(Mscores, 200)
M = rbinom(n, 1, Mscores)

# make a Y model according to the restrictions
f_Y = function(M,Z,W) {
  plogis(1*M + 1.5*W*Z + 1*Z - W - .7*M*W - 1)
  # plogis(1*M  + Z - .68*W - 1)
  # plogis(1*M  - Z*M - .68*W^2*Z - .3*Z - (W > .7)*M*Z + (W < -.4)*M - 1)
  # plogis(1*M + 1.5*W + 1*Z - 1)
  # plogis(1*M -.5*cos(W) +.68*cos(W)*Z - .38*Z*M + (W < -.4)*M - 1)
  # plogis(1*M + 1.5*W + 1*Z - 1)
  # plogis(W * M + 3 * cos(W) * Z - 0.4)
}

Yscores = f_Y(M,Z,W)
Y = rbinom(n, 1, Yscores)
hist(Yscores, 200)
min(Y*log(Yscores)+(1-Y)*log(1-Yscores))

# pack these functions into a DGP

n = 1e5
data = gendata.SDEtransport(n, f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)

# Let's generate the truth!!!
# figure out the chosen formula for M

Mstarfit = glm(M ~ Z + W , data = data[data$S==1, ], family = binomial())
Zstarfit = glm(Z ~ A + W + S, data = data, family = binomial())

f_M1 = function(Z, W, S) {
  # W = W[1:10]
  nn = length(W)
  dataM1 = data.frame(Z = rep(1,nn), W = W)
  predM1 = predict(Mstarfit, newdata = dataM1, type = 'response')
  dataM0 = data.frame(Z = rep(0,nn), W = W)
  predM0 = predict(Mstarfit, newdata = dataM0, type = 'response')
  dataZ = data.frame(A = rep(1,nn), W = W, S = rep(1, nn))
  predZ = predict(Zstarfit, newdata = dataZ, type = 'response')
  gM = predM1*predZ + predM0*(1 - predZ)
  return(gM)
}  

f_M0 = function(Z, W, S) {
  # W = W[1:10]
  nn = length(W)
  dataM1 = data.frame(Z = rep(1,nn), W = W)
  predM1 = predict(Mstarfit, newdata = dataM1, type = 'response')
  dataM0 = data.frame(Z = rep(0,nn), W = W)
  predM0 = predict(Mstarfit, newdata = dataM0, type = 'response')
  dataZ = data.frame(A = rep(0,nn), W = W, S = rep(1, nn))
  predZ = predict(Zstarfit, newdata = dataZ, type = 'response')
  gM = predM1*predZ + predM0*(1 - predZ)
  return(gM)
}  
# We can now set A to 1, and stochastically intervene for A = 0 to generate M and our outcome for pop

data_popS11 = gendata.SDEtransport(2*1e6,f_W, f_S = f_S, f_A = function(S,W) 1, f_Z = f_Z, f_M = f_M1, f_Y = f_Y)
true11 = mean(data_popS11$Y[data_popS11$S==0])

data_popS10 = gendata.SDEtransport(2*1e6, f_W, f_S = f_S, f_A = function(S,W) 1, f_Z = f_Z, f_M = f_M0, f_Y = f_Y)
true10 = mean(data_popS10$Y[data_popS10$S==0])

data_popS01 = gendata.SDEtransport(2*1e6, f_W, f_S = f_S, f_A = function(S,W) 0, f_Z = f_Z, f_M = f_M1, f_Y = f_Y)
true01 = mean(data_popS01$Y[data_popS01$S==0])

data_popS00 = gendata.SDEtransport(2*1e6, f_W, f_S = f_S, f_A = function(S,W) 0, f_Z = f_Z, f_M = f_M0, f_Y = f_Y)
true00 = mean(data_popS00$Y[data_popS00$S==0])

rm(data_popS11, data_popS10, data_popS01, data_popS00)
# computing a TMLE for each of these parameters
# start with a regression of Y on the past given S = 1
Yfit = glm(Y ~ M + Z + W, data = data[data$S==1,], family = binomial)
# Yfit = glm(Y ~ 1, data = data[data$S==1,], family = binomial)
Mfit = Mstarfit
# Mfit = glm(M ~ Z + W + S, data = data, family = binomial)
# Mfit = glm(M ~ 1, data = data, family = binomial)
Zfit = Zstarfit
# Zfit = glm(Z ~ 1, data = data, family = binomial)
Afit = glm(A ~ W + S, data = data, family = binomial)
# Afit = glm(A ~ M, data = data, family = binomial)
Sfit = glm(S ~ W, data = data, family = binomial)
# Sfit = glm(S ~ 1, data = data, family = binomial)
# undebug(fcn.SDE)

res = append(lapply(0:1, FUN = function(a_star) SDE_tmle1(data = data, 
                                                          a = 0, 
                                                          a_star = a_star, 
                                                          sl = sl, 
                                                          covariates = covariates,
                                                          truth = list(f_W = f_W, f_A = f_A, f_S = f_S, f_Z = f_Z, f_M = f_M, f_Y = f_Y),
                                                          truncate = list(lower =.0001, upper = .9999))),
lapply(0:1, FUN = function(a_star) SDE_tmle1(data = data, 
                                             a = 1, 
                                             a_star = a_star, 
                                             sl = sl, 
                                             covariates = covariates,
                                             truth = list(f_W = f_W, f_A = f_A, f_S = f_S, f_Z = f_Z, f_M = f_M, f_Y = f_Y),
                                             truncate = list(lower =.0001, upper = .9999))))

unlist(lapply(res, FUN = function(x) x[[1]][1]))
unlist(lapply(res, FUN = function(x) x$Psi_0))
c(true00, true01, true10, true11)
unlist(lapply(res, FUN = function(x) x$est_mle))
unlist(lapply(res, FUN = function(x) mean(x$IC)))


func_list = list(f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
covariates = list(covariates_S = c("W"),
                  covariates_A = c("S","W"),
                  covariates_Z = c("S","A","W"),
                  covariates_M = c("Z","W"),
                  covariates_Y = c("M","Z","W"))

sl = make_learner(Lrnr_glm_fast, family = binomial())

undebug(SDE_tmle)
sim_SDE = function(n, func_list, a, a_star) {
  f_W = func_list$f_W
  f_S = func_list$f_S
  f_A = func_list$f_A
  f_Z = func_list$f_Z
  f_M = func_list$f_M
  f_Y = func_list$f_Y
  data = gendata.SDEtransport(n, f_W = f_W, f_S = f_S, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)
  res = SDE_tmle1(data = data, a = a, a_star = a_star, sl = sl, covariates = covariates,
                 truth = list(f_W = f_W, f_A = f_A, f_S = f_S, f_Z = f_Z, f_M = f_M, f_Y = f_Y),
                 truncate = list(lower =.0001, upper = .9999))
  
  # res1 = SDE_tmle1(data = data, a = a, a_star = a_star, sl = sl, covariates = covariates,
  #                truth = list(f_W = f_W, f_A = f_A, f_S = f_S, f_Z = f_Z, f_M = f_M, f_Y = f_Y),
  #                truncate = list(lower =.0001, upper = .9999))
  
  # tmle est
  return(list(CI = res$CI, LR_est = res$est_mle, IC = res$IC, IC0 = res$IC_0, SL_coef = res$SL_coef,
              Psi_0 = res$Psi_0, SE = res$SE, SE0 = res$SE_0))
  
}

library("parallel")
library("ggplot2")
library(cowplot)
B = 1000
n=1000
a=0
a_star = 1

res = mclapply(1:B, FUN = function(x) sim_SDE(n, func_list, a = a, a_star = a_star), 
               mc.cores = getOption("mc.cores", 2L))


se = lapply(res, FUN = function(x) c(x$SE, x$SE0))
se = do.call(rbind, se)
colMeans(se)

# save(res, func_list, a, a_star, sl, file = "resSDE100.RData")
# load("resSDE1.RData")
res_est = do.call(rbind, lapply(res, FUN = function(x) c(x[[1]],x[[2]], x$Psi_0)))
cover_tmle = mean((res_est[,5] >= res_est[,2]) & (res_est[,5] <= res_est[,3]))
cover_tmle

inds = c(1,4)
mseDA_tmle = mean((res_est[,1] - res_est[,5])^2)
mseDA_lr = mean((res_est[,4] - res_est[,5])^2)
biasDA_tmle = mean((res_est[,1] - res_est[,5]))
biasDA_lr = mean((res_est[,4] - res_est[,5]))

psi0_avg = mean(res_est[,5])
mse = t(apply(res_est[,inds], 2, perf, psi0_avg))
rownames(mse) = c("TMLE", "LR")
mse
mseDA_tmle 
mseDA_lr 
biasDA_tmle
biasDA_lr
# real SE
mean(res_est[,3] - res_est[,2])/4
sd(res_est[,1])



types = c("TMLE", "LR")[order(inds)]
type = c(rep("TMLE",B), rep("LR", B))
ests = unlist(lapply(inds, FUN = function(x) res_est[,x]))

inds = inds[order(types)]
colors = c("green","blue")
plot_df = data.frame(ests = ests, type = type)

plot = ggplot(plot_df, aes(x = ests, fill = type)) + geom_density(alpha=.5)+
  ggtitle("Sampling Dists for SDE transport estimates")+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  geom_vline(xintercept = mean(res_est[,5]),color="black")+
  geom_vline(xintercept=mean(res_est[,inds[1]]),color = colors[1])+
  geom_vline(xintercept=mean(res_est[,inds[2]]),color = colors[2])

  capt = paste0("avg Truth is at black vline. ",B, " simulations were performed, n = ", n,".\n",
                 "Coverage for TMLE is ", round(cover_tmle,4), ". W is 1-d, a = ",a,
                " and a_star = ", a_star
               )
  plot = ggdraw(add_sub(plot ,capt, x= 0, y = 0.5, hjust = 0, vjust = 0.5,
                        vpadding = grid::unit(1, "lines"), fontfamily = "", 
                        fontface = "plain",colour = "black", size = 10, angle = 0, 
                        lineheight = 0.9))
  plot

hist(res_est[,1])
debug(fcn.SDE)
sd(res_est[,1])
mean(res_est[3]-res_est[,2])
library(latex2exp)
latex2exp(X<Y)

plot(x, xlim=c(0, 4), ylim=c(0, 10), 
     xlab='x', ylab=TeX('$\\alpha  x^\\alpha$, where $\\alpha \\in 1\\ldots 5$'), 
     type='n', main=TeX('Using $\\LaTeX$ for plotting in base graphics!'))

x <- seq(0, 4, length.out=100)
alpha <- 1:5

plot(x, xlim=c(0, 4), ylim=c(0, 10), 
     xlab='x', ylab=TeX('$\\alpha  x^\\alpha$, where $\\alpha \\in 1\\ldots 5$'), 
     type='n', main=TeX('Using $\\LaTeX$ for plotting in base graphics!'))

invisible(sapply(alpha, function(a) lines(x, a*x^a, col=a)))

legend('topleft', legend=TeX(sprintf("$\\alpha = %d$", alpha)), 
       lwd=1, col=alpha)


