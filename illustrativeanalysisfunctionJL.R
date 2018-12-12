# devtools::install_github("jeremyrcoyle/sl3")
devtools::install_github("jlstiles/SDEtransportsim")
library("SDEtransportsim")
setwd("/Users/kararudolph/Documents/K99R00/R00/transportSDE")
source("transportSDEfunctions.R")
library("glmnet")
#Set up data generating process:
# data generating process for 2-d W
f_W = function(n) {
  W1 = rnorm(n)
  gender = rbinom(n, 1, 0.5)
  W2 = rnorm(n)
  data.frame(W1=W1, W2=W2, Wgender = gender)
}

# make a pscore model
f_A = function(W) {
  with(W, plogis(-.7*W1 + 0.1*W2 + .17))
}
# make a intermediate confounder model
f_Z = function(A,W) {
  df = as.data.frame(cbind(W, A = A))
  with(df, plogis(.4*W1 - W2 + .2*Wgender + 1*A*Wgender -.3))
}
# make an M model according to the restrictions
f_M = function(Z,W) {
  df = as.data.frame(cbind(Z=Z, W))
  with(df, plogis(1*W1 - W2 + .2*Wgender- 1.2*Z +.1))
}
# make a Y model according to the restrictions, main terms linear logistic reg.
# plug-in is biased and not robust like tmle
f_Y = function(M,Z,W) {
  df = as.data.frame(cbind(M = M, Z = Z, W))
  with(df, plogis(W2*M + 3*cos(W1)*Z-.4*Wgender - .4))
}

# generate n random samples
n = 1e4
# set.seed(1)
df = gendata.SDEtransport_alt(n, f_W = f_W, f_A = f_A, f_Z = f_Z, f_M = f_M, f_Y = f_Y)

# setting up a multinomial site variable
df$S = apply(rmultinom(n, 1, c(.25,.25,.25,.25)), 2, FUN = function(x) which(x==1))

# adding two more outcome cols so three total
df$Y1 = df$Y
df$Y2 = df$Y

# adding four more mediator cols so 5 total
df$M1 = df$M
df$M2 = df$M
df$M3 = df$M
df$M4 = df$M

#KER weights
df$weights<-rnorm(nrow(df), 1,1)
# These can be easily gotten from your data.frame 
# KERthe wcols column names need to include  all possible interaction variables too
# Wnames<-colnames(df[1:2])
Wnames = c("W1", "W2")
Wnamesalways =c("W1", "W2") 

library(doParallel)

#Ker 
# Wnamesalways = c("svy_ethrace", "svy_age_imp")

#KER note to rename Z variable with Z and A variable with A, weights variable with weights

# make formulas with main terms and interactions
forms = list(
  Aform = NULL, 
  Zstarform = formula(paste0("Z ~ (", paste(c(Wnames, "A"), "", collapse = "+"), ")^2")), 
  Mstarform = formula(paste0("M~ (", paste(c(Wnames, "Z"), "", collapse = "+"), ")^2")),
  QZform = formula(paste0("Qstar_Mg ~ (", paste(c(Wnames), "", collapse = "+"), ")^2")), 
  Yform = formula(paste0("Y ~ (", paste(c(Wnames, "Z", "M"), "", collapse = "+"), ")^2")) 
)

#KER

Wtmp = c("W1", "W2", "W3")
subsetM1 = Wtmp[c(1,3)]
formM = formula(paste0("Y~", paste(paste0(subsetM1,"*M"), "", collapse = "+")))

forms = list(
  Aform = NULL, 
  Zstarform = formula(paste0("Z ~ ", paste(c(Wnames, "A",paste0(Wnames,":A")), "", collapse = "+") )), 
  Mstarform = formula(paste0("M~ ", paste(c(Wnames, "Z",paste0(Wnames,":Z")), "", collapse = "+") )),
  QZform = formula(paste0("Qstar_Mg ~ ", paste(c(Wnames), "", collapse = "+") )), 
  Yform = formula(paste0("Y ~ ", paste(c(Wnames, "Z", "M",paste0(Wnames,":Z"),paste0(Wnames,":M")), "", collapse = "+") )) 
)

# The variables for subsetting and for different mediator-oc combos
gender = c(0,1)
site = c(1,2,3,4)
mediator = c("M", "M1", "M2", "M3", "M4")
outcome = c("Y", "Y1", "Y2")

results = lapply(gender, FUN = function(g) {
  lapply(site, FUN = function(s) {
    lapply(mediator, FUN = function(med) {
      lapply(outcome, FUN = function(oc) {
         # KR: I'm going to put mediator and outcome first because their spots may vary
         # KER edit for Wgender to be colname of gender and S to be colname of site
        data = subset(df, Wgender==g & S==s, select = c(med, oc, Wnames, "A", "Z", "weights"))
        #data = subset(df, Wgender==g & S==s, select = c("W1", "W2", "A", "Z", med, oc))
        # replace colnames so formulas all work

        colnames(data)[1:2] = c("M", "Y")
        # the main function here
        res = SDE_tmle_lasso(data, forms, RCT = 0.5, B = NULL, Wnames = Wnames, Wnamesalways = Wnamesalways,
                             transport = FALSE)
        return(res)
      })
    })
  })
})

results_df = data.frame(matrix(rep(NA, 120*8), nrow = 120))
# store results in a dataframe
j=1
for (g in 1:length(gender)) {
  for (s in 1:length(site)) {
    for (med in 1:length(mediator)) {
      for (oc in 1:length(outcome)) {
        
        results_df[(2*j-1),] = c(results[[g]][[s]][[med]][[oc]]$CI_SDE, 
                                 gender[g], site[s], mediator[med], outcome[oc], "SDE") 
        results_df[(2*j),] = c(results[[g]][[s]][[med]][[oc]]$CI_SIE, 
                               gender[g], site[s], mediator[med], outcome[oc], "SIE") 
        j = j + 1
      }
    }
  }
}

results_df[,1:3] = apply(results_df[,1:3], 2, FUN = function(x) round(as.numeric(x), 5))    
colnames(results_df) = c("est", "left", "right", "gender", "site", "mediator", "outcome", "param")
subset(results_df, param=="SDE"&gender=="female")

# name the list easily as you like
gender = c("female", "male")
site = c("s1", "s2", "s3", "s4")
mediator = c("M1", "M2", "M3", "M4", "M5")
outcome = c("Y1", "Y2", "Y3")
for (a in 1:2) {
  names(results)[a] = gender[a]
  for (b in 1:4) {
    names(results[[a]])[b] = site[b]
    for (c in 1:5) {
      names(results[[a]][[b]])[c] = mediator[c]
      for (d in 1:3)
        names(results[[a]][[b]][[c]])[d] = outcome[d]
    }
  }
}

results$female$s1$M4$Y2
