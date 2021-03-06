% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SDE_glm.R
\name{SDE_glm}
\alias{SDE_glm}
\title{SDE_tmle_glm}
\usage{
SDE_glm(data, forms, RCT = 0.5, transport, pooled, gstar_S = 1, truth,
  B = NULL, max_iter = 100)
}
\arguments{
\item{data, }{data.frame where confounders have any names but the following, which must be
reserved as follows: A = treatment, Z = intermediate confounder, M = mediator and Y is the outcome.}

\item{forms, }{list of formulas. Include for each necessary model for outcome,
called Yform for outcome Y, QZform for outcome Qstar_Mg, Mform, Zform, Aform (can be NULL if RCT)
is selected as TRUE.}

\item{RCT}{either NULL or a value, if null, then the Aform is used to fit the propensity score,
otherwise propensity scores are set to RCT.}

\item{transport}{if true you are transporting to site S=0}

\item{pooled}{set to TRUE if you wish to define the stochastic intervention by the mechanism for
the mediator defined by pooling the regression across both sites.  Otherwise the stochastic intervention
will only be defined by the subset for S = gstar_S for both M and Z.}

\item{gstar_S}{set to 0 or 1 depending on which site you want to use to define the stochastic
intervention}
}
\value{
a list with a CI's for SDE and SIE for the means under (a*,a) combos (0,0), (0,1), (1,1)
and the epsilons for both sequential regressions for those three parameters
}
\description{
computes the sequential regression, targeted maximum likelihood estimate
for the stochastic direct effect or stochastic indirect effect using lasso. Note, this is
a non-transport parameter.
}
\examples{
# uncomment these if testing the lasso functions without loading the package
# library("parallel")
# library("doParallel")
# library("glmnet")

#Set up data generating process:
data(example_dgp)
data = data_example

# forms for non-transporting contains no S variable
formsNT=list(Aform = formula("A ~ W2 + W1"),
             Zstarform = formula("Z ~  A + W1 + W2"),
             QZform = formula("Qstar_Mg ~ Z + W2 + W1"),
             Mstarform = formula("M ~ Z + W1 + W2"),
             Yform = formula("Y ~ M + Z + W1 + W2"))

# forms for transporting without a pooled (over S) M or Z fit
formsNP=list(Sform = formula("S ~  W1 + W2"),
           Aform = formula("A ~ W2 + W1 + S"),
           Zstarform = formula("Z ~  A + W1 + W2"),
           QZform = formula("Qstar_Mg ~ Z + W2 + W1+S"),
           Mstarform = formula("M ~ Z + W1 + W2"),
           Yform = formula("Y ~ M + Z + W1 + W2"))

# forms for transporting with a pooled (over S) M and Z fit
forms=list(Sform = formula("S ~  W1 + W2"),
           Aform = formula("A ~ W2 + W1 + S"),
           Zstarform = formula("Z ~  A + W1 + W2 + S"),
           QZform = formula("Qstar_Mg ~ Z + W2 + W1 + S"),
           Mstarform = formula("M ~ Z + W1 + W2 + S"),
           Yform = formula("Y ~ M + Z + W1 + W2"))


W = data[,grep("W", colnames(data))]
head(data)
data = cbind(W, A=data$A, S=data$S, Z=data$Z, M=data$M, Y=data$Y)

# define the forms for all the regressions
n=1000

Wnames = c("W1", "W2")
Wnamesalways = c("W1")

# choose no weights (all weights are 1)
data$weights = rep(1,nrow(data))
# choose randomly generated weights
data$weights = runif(nrow(data))

i=1
res = data.frame(matrix(rep(NA,36),nrow = 6))
colnames(res) = c("SDE", "left", "right","SIE", "left", "right")
rownames(res) = c("Trans, pooled, gstar_S=0", 
                  "Trans, pooled, gstar_S=1", 
                  "Trans, not pooled, gstar_S=0", 
                  "Trans, not pooled, gstar_S=1",
                  "Not trans, gstar_S=0", 
                  "Not trans, gstar_S=1")

for (pooled in c(FALSE, TRUE)) {
  for (gstar_S in c(0,1)) {
    if (pooled) forms = forms else forms = formsNP
    ans = SDE_tmle_lasso(data, forms, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                         B = NULL, transport = TRUE, pooled = pooled, gstar_S = gstar_S) 
    res[i,]=c(ans$CI_SDE,ans$CI_SIE)
    i = i + 1
  }
}

for (gstar_S in c(0,1)) {
  ans = SDE_tmle_lasso(data, formsNT, RCT = 0.5, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                       B = NULL, transport = FALSE, pooled = pooled, gstar_S = gstar_S) 
  res[i,]=c(ans$CI_SDE,ans$CI_SIE)
  i = i + 1
}

res

# undebug(SDE_tmle_lasso)
# undebug(get_gstarM_lasso)
# undebug(get.mediation.initdata_lasso)
# undebug(mediation.step1_lasso)
# undebug(mediation.step2_lasso)
}
