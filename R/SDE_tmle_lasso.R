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
#' @param B, the number of bootstraps, default is NULL and should not be changed since this will
#' be invalid for use with lasso anyway.
#' @param Wnames
#' @param Wnamesalways
#' @param transport if true you are transporting to site S=0
#' @param pooledM set to TRUE if you wish to define the stochastic intervention by the mechanism for 
#' the mediator defined by pooling the regression across both sites.  Otherwise the stochastic intervention
#' will only be defined by the fit for S = 1.
#' @param truth set permanently to NULL, not used
#' @return  a list with a CI's for SDE and SIE for the means under (a*,a) combos (0,0), (0,1), (1,1) 
#' and the epsilons for both sequential regressions for those three parameters
#' @example /inst/example_SDE_lasso.R 
#' @export
SDE_tmle_lasso = function(data, forms, RCT = 0.5, B = NULL, Wnames, Wnamesalways, transport = TRUE,
                          pooledM = TRUE, truth = NULL) 
{
    # get the stochastic dist of M and true params if you want 
    gstar_info = get_gstarM_lasso(data = data, forms = forms, Wnames = Wnames, Wnamesalways = Wnamesalways, 
                                  transport = transport, pooledM = pooledM)
    gstarM_astar1 = gstar_info$gstarM_astar1
    gstarM_astar0 = gstar_info$gstarM_astar0
    gstarM_astar = list(gstarM_astar0 = gstarM_astar0, gstarM_astar1 = gstarM_astar1)
    
    # perform initial fits for the first regression
    
    init_info = get.mediation.initdata_lasso(data = data, forms = forms, RCT = RCT, Wnames = Wnames, 
                                             Wnamesalways = Wnamesalways, transport = transport,
                                             pooledM = pooledM)
    
    Y_preds = init_info$Y_preds
    
    est_info = lapply(0:1, FUN = function(astar) {
      lapply(0:1, FUN = function(a) {
        # get tmle info
        # get iptw here while I'm at it
        update = mediation.step1_lasso(initdata = init_info$initdata, init_info$Y_preds, data = data, 
                                       gstarM_astar[[astar+1]], a, transport = transport)
        
        Y_Mg = get.stochasticM(gstarM_astar[[astar+1]], Y_preds[[2]], Y_preds[[3]]) 
        A_ps = init_info$initdata$A_ps
        
        Qstar_Mg = get.stochasticM(gstarM_astar[[astar+1]], update$Qstar_M1, update$Qstar_M0) 
        # compute Qstar_Mg here
        tmle_info = mediation.step2_lasso(data = data, Qstar_M = update$Qstar_M, 
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
    if (!is.null(B)) {
      boot_ests = lapply(1:B, FUN = function(x) {
        inds = sample(1:n, replace = TRUE)
        data = data[inds, ]
        
        init_info = get.mediation.initdata_lasso(data = data, forms = forms, RCT = RCT, Wnames = Wnames, 
                                                 Wnamesalways = Wnamesalways, transport = transport, 
                                                 pooledM = pooledM)
        Y_preds = init_info$Y_preds
        gstarM_astar = list(gstarM_astar0[inds], gstarM_astar1[inds],
                            Wnames = Wnames, Wnamesalways = Wnamesalways, 
                            transport = transport, pooledM = pooledM)
        
        est_info = lapply(0:1, FUN = function(astar) {
          return(lapply(0:1, FUN = function(a) {
            update = mediation.step1_lasso(initdata = init_info$initdata, Y_preds = Y_preds, data = data, 
                                           gstarM_astar[[astar+1]], a, transport = transport)
            
            Y_Mg = get.stochasticM(gstarM_astar[[astar+1]], Y_preds[[2]], Y_preds[[3]]) 
            A_ps = init_info$initdata$A_ps
            
            Qstar_Mg = get.stochasticM(gstarM_astar[[astar+1]], update$Qstar_M1, update$Qstar_M0) 
            # compute Qstar_Mg here
            tmle_info = mediation.step2_lasso(data = data, Qstar_M = update$Qstar_M, 
                                              Qstar_Mg = Qstar_Mg, Hm = update$Hm, A_ps = A_ps, 
                                              a = a, tmle = TRUE,
                                              EE = FALSE, bootstrap = TRUE, form = forms$QZform,
                                              transport = transport)
            # compile all estimates
            return(tmle_info)
          }))
        })
      })
      
      boot_ests = lapply(boot_ests, FUN = function(boot) {
        SDE = unlist(boot[[1]][[2]]) - unlist(boot[[1]][[1]])
        SIE = unlist(boot[[2]][[2]]) - unlist(boot[[1]][[2]])
        return(list(SDE, SIE))
      })
      
      boots_SDE = do.call(rbind, lapply(boot_ests, FUN = function(boot) boot[[1]]))
      boots_SIE = do.call(rbind, lapply(boot_ests, FUN = function(boot) boot[[2]]))
      
      bootSE_SDE = apply(boots_SDE, 2, sd)
      bootSE_SIE = apply(boots_SIE, 2, sd)
    }
    
    D_SDE = est_info[[1]][[2]]$IC - est_info[[1]][[1]]$IC
    D_SIE = est_info[[2]][[2]]$IC- est_info[[1]][[2]]$IC
    
    SE_SDE = sd(D_SDE)/sqrt(n)
    SE_SIE = sd(D_SIE)/sqrt(n)
    
    if (!is.null(truth)) { 
      D_SDE_0 = gstar_info$D_astar0a1_0 - gstar_info$D_astar0a0_0
      D_SIE_0 = gstar_info$D_astar1a1_0 - gstar_info$D_astar0a1_0
      
      SE_SDE_0 = sd(D_SDE_0)/sqrt(n)
      SE_SIE_0 = sd(D_SIE_0)/sqrt(n)
      
      Psi_astar0a1_0 = gstar_info$Psi_astar0a1_0
      Psi_astar0a0_0 = gstar_info$Psi_astar0a0_0
      Psi_astar1a1_0 = gstar_info$Psi_astar1a1_0
      
      SDE_0 = Psi_astar0a1_0 - Psi_astar0a0_0
      SIE_0 = Psi_astar1a1_0 - Psi_astar0a1_0 
      
    } else {
      D_SDE_0 = D_SIE_0 = SE_SDE_0 = SE_SIE_0 = Psi_astar0a1_0 = 
        Psi_astar0a0_0 = Psi_astar1a1_0 = SDE_0 = SIE_0 = NULL
    }
    ests_astar0a1 =   c(tmle = est_info[[1]][[2]]$est,
                        Psi_0 = Psi_astar0a1_0)
    
    ests_astar0a0 =   c(tmle = est_info[[1]][[1]]$est,
                        Psi_0 = Psi_astar0a0_0)
    
    ests_astar1a1 =   c(tmle = est_info[[2]][[2]]$est,
                        Psi_0 = Psi_astar1a1_0)
    
    SDE_ests = ests_astar0a1 - ests_astar0a0
    SIE_ests = ests_astar1a1 - ests_astar0a1
    
    CI_SDE = c(SDE_ests[1], SDE_ests[1] - 1.96*SE_SDE, SDE_ests[1] + 1.96*SE_SDE)
    CI_SIE = c(SIE_ests[1], SIE_ests[1] - 1.96*SE_SIE, SIE_ests[1] + 1.96*SE_SIE)
    
    if (!is.null(B)) {
      CI_SDE_boot = c(SDE_ests[1], SDE_ests[1] - 1.96*bootSE_SDE[1], SDE_ests[1] + 1.96*bootSE_SDE[1])
      CI_SIE_boot = c(SIE_ests[1], SIE_ests[1] - 1.96*bootSE_SIE[1], SIE_ests[1] + 1.96*bootSE_SIE[1])
      
      return(list(CI_SDE = CI_SDE, CI_SIE = CI_SIE,
                  CI_SDE_boot = CI_SDE_boot,
                  CI_SIE_boot = CI_SIE_boot,
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


