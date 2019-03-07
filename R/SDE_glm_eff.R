#' @title SDE_tmle_glm_effT
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
#' @example /inst/example_SDE_lasso.R 
#' @export
SDE_glm_eff = function(data, forms, RCT = 0.5,transport, pooled, gstar_S = 1, truth = NULL, B = NULL, 
                       max_iter = 100) 
{
    if (!transport) pooled = FALSE
    # get the stochastic dist of M and true params if you want 
    gstar_info = get_gstarM_glm_eff(data = data, forms = forms, 
                                  transport = transport, pooled = pooled, gstar_S = gstar_S, truth= truth)
    gstarM_astar1 = gstar_info$gstarM_astar1 
    gstarM_astar0 = gstar_info$gstarM_astar0
    gstarM_astar = list(gstarM_astar0 = gstarM_astar0, gstarM_astar1 = gstarM_astar1)
    
    # perform initial fits for the first regression
    
    initdata = get.mediation.initdata_glm_eff(data = data, forms = forms, RCT = RCT, transport = transport, 
                                             pooled = pooled, gstarM_astar)

    est_info = mediation.steps_glm_eff(initdata, transport, bootstrap = FALSE, max_iter)
    Psi_SDE = est_info$Psi_SDE
    Psi_SIE = est_info$Psi_SIE
    Psi_SDE_EE = est_info$EE_info$Psi_SDE
    Psi_SIE_EE = est_info$EE_info$Psi_SIE    
    # bootstrap from the data and thus the gstarM_astar1 and gstarM_astar0
    n = nrow(data)
    if (!is.null(B)) {
      boot_ests = lapply(1:B, FUN = function(x) {
        inds = sample(1:n, replace = TRUE)
        data = data[inds,]
        gstarM_astar = list(gstarM_astar0 = gstarM_astar0[inds], gstarM_astar1 = gstarM_astar1[inds])
        
        initdata = get.mediation.initdata_glm_eff(data = data, forms = forms, RCT = RCT, transport = transport, 
                                                   pooled = pooled, gstarM_astar)
        
        est_info = mediation.steps_glm_eff(initdata, transport, bootstrap = TRUE, max_iter)
            # compile all estimates
        return(est_info)
      }
      )
      its_boot = mean(unlist(lapply(boot_ests, FUN = function(x) x[5])))
      boot_ests = do.call(rbind, boot_ests)
      boot_SE = apply(boot_ests, 2, sd)
      SE_SDE_boot = boot_SE[1]
      SE_SIE_boot = boot_SE[2]
      SE_SDE_EE_boot = boot_SE[3]
      SE_SIE_EE_boot = boot_SE[4]
      CI_SDE_boot = c(Psi_SDE, Psi_SDE - 1.96*SE_SDE_boot, Psi_SDE + 1.96*SE_SDE_boot)
      CI_SIE_boot = c(Psi_SIE, Psi_SIE - 1.96*SE_SIE_boot, Psi_SIE + 1.96*SE_SIE_boot)
      CI_SDE_EE_boot = c(Psi_SDE_EE, Psi_SDE_EE - 1.96*SE_SDE_EE_boot, Psi_SDE_EE + 1.96*SE_SDE_EE_boot)
      CI_SIE_EE_boot = c(Psi_SIE_EE, Psi_SIE_EE - 1.96*SE_SIE_EE_boot, Psi_SIE_EE + 1.96*SE_SIE_EE_boot)
    } else {
      SE_SDE_boot = NULL
      SE_SIE_boot = NULL
      SE_SDE_EE_boot =NULL
      SE_SIE_EE_boot = NULL
      CI_SDE_boot = NULL
      CI_SIE_boot = NULL
      CI_SDE_EE_boot = NULL
      CI_SIE_EE_boot = NULL
      its_boot = NULL
    }
    
    D_SDE = est_info$D[,1]
    D_SIE = est_info$D[,2]
    D_SDE_EE = est_info$EE_info$D[,1]
    D_SIE_EE = est_info$EE_info$D[,2]
    
    # weighted variance as per weights option in data matrix
    
    SE_SDE = est_info$SE[1]
    SE_SIE = est_info$SE[2]
    SE_SDE_EE = est_info$EE_info$SE[1]
    SE_SIE_EE = est_info$EE_info$SE[2]
    
    CI_SDE = c(Psi_SDE, Psi_SDE - 1.96*SE_SDE, Psi_SDE + 1.96*SE_SDE)
    CI_SIE = c(Psi_SIE, Psi_SIE - 1.96*SE_SIE, Psi_SIE + 1.96*SE_SIE)
    
    CI_SDE_EE = c(Psi_SDE_EE, Psi_SDE_EE - 1.96*SE_SDE_EE, Psi_SDE_EE + 1.96*SE_SDE_EE)
    CI_SIE_EE = c(Psi_SIE_EE, Psi_SIE_EE - 1.96*SE_SIE_EE, Psi_SIE_EE + 1.96*SE_SIE_EE)
    
    if (is.null(truth)) {   
      return(list(CI_SDE = CI_SDE, 
                  CI_SIE = CI_SIE,
                  CI_SDE_EE = CI_SDE_EE, 
                  CI_SIE_EE = CI_SIE_EE,
                  SE_SDE = SE_SDE, 
                  SE_SIE = SE_SIE,
                  SE_SDE_EE = SE_SDE_EE, 
                  SE_SIE_EE = SE_SIE_EE, 
                  CI_SDE_boot = CI_SDE_boot, 
                  CI_SIE_boot = CI_SIE_boot,
                  CI_SDE_EE_boot = CI_SDE_EE_boot, 
                  CI_SIE_EE_boot = CI_SIE_EE_boot,
                  SE_SDE_boot = SE_SDE_boot, 
                  SE_SIE_boot = SE_SIE_boot,
                  SE_SDE_EE_boot = SE_SDE_EE_boot, 
                  SE_SIE_EE_boot = SE_SIE_EE_boot,
                  its = est_info$its,
                  its_boot = its_boot
      )) 
    } else {
      
      return(list(CI_SDE = CI_SDE, 
                  CI_SIE = CI_SIE,
                  CI_SDE_EE = CI_SDE_EE, 
                  CI_SIE_EE = CI_SIE_EE,
                  SE_SDE = SE_SDE, 
                  SE_SIE = SE_SIE,
                  SE_SDE_EE = SE_SDE_EE, 
                  SE_SIE_EE = SE_SIE_EE, 
                  SE_SDE_boot = SE_SDE_boot, 
                  SE_SIE_boot = SE_SIE_boot,
                  SE_SDE_EE_boot = SE_SDE_EE_boot, 
                  SE_SIE_EE_boot = SE_SIE_EE_boot,
                  CI_SDE_boot = CI_SDE_boot, 
                  CI_SIE_boot = CI_SIE_boot,
                  CI_SDE_EE_boot = CI_SDE_EE_boot, 
                  CI_SIE_EE_boot = CI_SIE_EE_boot,
                  its = est_info$its,
                  its_boot = its_boot,
                  SDE_0 = gstar_info$SDE_0,
                  SIE_0 = gstar_info$SIE_0,
                  SE_SDE_0 = gstar_info$SE_SDE_0,
                  SE_SIE_0 = gstar_info$SE_SIE_0
      ))
      
    }
  } 


