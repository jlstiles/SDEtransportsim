#' @title SDE_tmle_lasso
#' @description computes the sequential regression, targeted maximum likelihood estimate
#' for the stochastic direct effect or stochastic indirect effect using lasso. Note, this is
#' a non-transport parameter.
#' @param data, data.frame where confounders have any names but the following, which must be 
#' reserved as follows: A = treatment, Z = intermediate confounder, M = mediator and Y is the outcome.
#' @param forms, list of formulas. Include for each necessary model for outcome, 
#' called Yform for outcome Y, QZform for outcome Qstar_Mg, Mform, Zform, Aform (can be NULL if RCT)
#' is selected as TRUE. 
#' @param RCT either NULL or a value, if null, then the Aform is used to fit the propensity score, 
#' otherwise propensity scores are set to RCT.
#' @param Wnames
#' @param Wnamesalways
#' @param transport if true you are transporting to site S=0
#' @param pooled set to TRUE if you wish to define the stochastic intervention by the mechanism for 
#' the mediator defined by pooling the regression across both sites.  Otherwise the stochastic intervention
#' will only be defined by the subset for S = gstar_S for both M and Z.   
#' @param gstar_S set to 0 or 1 depending on which site you want to use to define the stochastic
#' intervention
#' @param truth can use a truth generator to check matching the truth
#' @return  a list with a CI's for SDE and SIE for the means under (a*,a) combos (0,0), (0,1), (1,1) 
#' and the epsilons for both sequential regressions for those three parameters
#' 
#'
SDE_tmle_lasso1 = function(data, forms, RCT = 0.5,Wnames, Wnamesalways, transport,
                          pooled, gstar_S = 1, truth = NULL) 
{
    if (!transport) pooled = FALSE
    # get the stochastic dist of M and true params if you want 
    gstar_info = get_gstarM_lasso1(data = data, forms = forms, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                                  transport = transport, pooled = pooled, gstar_S = gstar_S, truth= truth)
    gstarM_astar1 = gstar_info$gstarM_astar1 
    gstarM_astar0 = gstar_info$gstarM_astar0
    gstarM_astar = list(gstarM_astar0 = gstarM_astar0, gstarM_astar1 = gstarM_astar1)
    
    # perform initial fits for the first regression
    
    init_info = get.mediation.initdata_lasso1(data = data, forms = forms, RCT = RCT, Wnames = Wnames, 
                                             Wnamesalways = Wnamesalways, transport = transport, 
                                             pooled = pooled)
    
    Y_preds = init_info$Y_preds
    
    est_info = lapply(0:1, FUN = function(astar) {
      lapply(0:1, FUN = function(a) {
        # get tmle info
        # get iptw here while I'm at it
        update = mediation.step1_lasso1(initdata = init_info$initdata, init_info$Y_preds, data = data, 
                                       gstarM_astar[[astar+1]], a, transport = transport)
        
        Y_Mg = get.stochasticM(gstarM_astar[[astar+1]], Y_preds[[2]], Y_preds[[3]]) 
        A_ps = init_info$initdata$A_ps
        
        Qstar_Mg = get.stochasticM(gstarM_astar[[astar+1]], update$Qstar_M1, update$Qstar_M0) 
        # compute Qstar_Mg here
        tmle_info = mediation.step2_lasso1(data = data, Qstar_M = update$Qstar_M, 
                                          Qstar_Mg = Qstar_Mg, Hm = update$Hm, A_ps = A_ps, 
                                          a = a, tmle = TRUE,
                                          EE = FALSE, bootstrap = FALSE, form = forms$QZform,
                                          transport = transport)
        
        # compile all estimates
        tmle_info$eps1 = update$eps
        return(tmle_info)
      }) 
    })
    
    # run the bootstrap 500 times
    # bootstrap from the data and thus the gstarM_astar1 and gstarM_astar0
    n = nrow(data)
    
    D_SDE = est_info[[1]][[2]]$IC - est_info[[1]][[1]]$IC
    D_SIE = est_info[[2]][[2]]$IC- est_info[[1]][[2]]$IC
    
    # weighted variance as per weights option in data matrix
    
    SE_SDE = sd(D_SDE)/sqrt(n)
    SE_SIE = sd(D_SIE)/sqrt(n)

    ests_astar0a1 =   est_info[[1]][[2]]$est
    ests_astar0a0 =   est_info[[1]][[1]]$est
    ests_astar1a1 =   est_info[[2]][[2]]$est
    
    SDE_ests = ests_astar0a1 - ests_astar0a0
    SIE_ests = ests_astar1a1 - ests_astar0a1
    
    CI_SDE = c(SDE_ests, SDE_ests - 1.96*SE_SDE, SDE_ests + 1.96*SE_SDE)
    CI_SIE = c(SIE_ests, SIE_ests - 1.96*SE_SIE, SIE_ests + 1.96*SE_SIE)
    
    if (!is.null(truth)) {   
      return(list(CI_SDE = CI_SDE, 
                  CI_SIE = CI_SIE,
                  SDE0 = gstar_info$SDE0,
                  SIE0 = gstar_info$SIE0,
                  SE_SDE = SE_SDE, 
                  SE_SIE = SE_SIE, 
                  ests_astar0a1 = ests_astar0a1,
                  ests_astar0a0 = ests_astar0a0,
                  ests_astar1a1 = ests_astar1a1,
                  eps1_astar0a1 = est_info[[1]][[2]]$eps1,
                  eps2_astar0a1 = est_info[[1]][[2]]$eps2,
                  eps1_astar0a0 = est_info[[1]][[1]]$eps1,
                  eps2_astar0a0 = est_info[[1]][[1]]$eps2,
                  eps1_astar1a1 = est_info[[2]][[2]]$eps1,
                  eps2_astar1a1 = est_info[[2]][[2]]$eps2 
      )) 
    } else {
      
      return(list(CI_SDE = CI_SDE, 
                  CI_SIE = CI_SIE, 
                  SE_SDE = SE_SDE, 
                  SE_SIE = SE_SIE, 
                  ests_astar0a1 = ests_astar0a1,
                  ests_astar0a0 = ests_astar0a0,
                  ests_astar1a1 = ests_astar1a1,
                  eps1_astar0a1 = est_info[[1]][[2]]$eps1,
                  eps2_astar0a1 = est_info[[1]][[2]]$eps2,
                  eps1_astar0a0 = est_info[[1]][[1]]$eps1,
                  eps2_astar0a0 = est_info[[1]][[1]]$eps2,
                  eps1_astar1a1 = est_info[[2]][[2]]$eps1,
                  eps2_astar1a1 = est_info[[2]][[2]]$eps2 
      ))
      
    }
  } 


