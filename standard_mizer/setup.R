# In this file we are setting up the functions for
# - the dynamics of the resources
# - the contribution of the resources to the encounter rate
# - the senescence mortality
# and provide a function for setting up Asta's model with
# benthos and algae.
#
# This file is sourced by the file run.R that runs the 
# climate change scenarios

library(tidyverse)
library(mizerExperimental)
library(mizerStarvation)

background_semichemostat <- function(params, n_other, rates, dt, component,
                                     ...) {
    c <- params@other_params[[component]]
    # name of interaction parameter for this component in species_params
    interaction_component <- paste0("interaction_", component)
    interaction <- params@species_params[[interaction_component]]
    mort <- as.vector(interaction  %*% rates$pred_rate)
    tmp <- c$rate * c$capacity / (c$rate + mort)
    return(tmp - (tmp - n_other[[component]]) * exp(-(c$rate + mort) * dt))
}

background_encounter <- function(params, n, n_pp, n_other, ...) {
    idx_sp <- (length(params@w_full) - length(params@w) + 1):length(params@w_full)
    prey <- outer(params@species_params$interaction_resource, n_pp) +
        outer(params@species_params$interaction_aa, n_other$aa) +
        outer(params@species_params$interaction_bb, n_other$bb)
    prey[, idx_sp] <- prey[, idx_sp] + params@interaction %*% n
    prey <- sweep(prey, 2, params@w_full * params@dw_full, "*")
    avail_energy <- Re(base::t(mvfft(base::t(params@ft_pred_kernel_e) * 
                                         mvfft(base::t(prey)),
                                     inverse = TRUE))) / length(params@w_full)
    avail_energy <- avail_energy[, idx_sp, drop = FALSE]
    avail_energy[avail_energy < 1e-18] <- 0
    
    params@search_vol * avail_energy
}



#parameters for senescence mortality as used in Law et al. 2009
k.sm <- 0.1 # mortality per year at the threshold size (should be 0.5 originally)
xsw <- 0.95 # proportion of w_inf at which mortality is at k.sm (should be 0.9)
sen.e <- 3  # exponent of the senescence mortality (larger value will give 
# steeper increase in the last few sizes) (should be 3)

sen_mort <- function(sppParams, params, k.sm, xsw, sen.e) {
    sen.mort.m = array(0, dim = c(length(sppParams$species),length(params@dw)))
    for (i in 1:length(sppParams$species)) {
        mu_Sen = k.sm * 10^(sen.e*(log(params@w) - log(xsw*sppParams$w_inf[i])))
        sen.mort.m[i,] <- mu_Sen    
    }
    # For really small species, like Trachinops, Pictilabrus and urchins 
    # predation mortality will be so high that senescence mortality is unlikely 
    # to be the case and perhaps should not even be applied
    sen.mort.m[which(params@species_params$w_inf < 400),] <- 0
    return(sen.mort.m)
}

newAstaParams <- function(sp, interaction, 
                          temperature = 12,
                          Ea = 0.63,
                          t_ref = 12,
                          no_w = 200, 
                          kappa = 2,
                          lambda = 2.15,
                          kappa_ben = 6, 
                          lambda_ben = 1.9,
                          stable_abund,
                          stable_pl,
                          stable_ben,
                          stable_alg) {
    # temperature factor
    temperature <- temperature + 273 # converting to Kelvin from Celcius
    t_ref <- t_ref + 273
    temperatureScalar <- exp(-Ea / 8.617332e-5 * (1/temperature - 1/t_ref))
    
    # We rescale rates. Mortality is rescaled futher down
    sp$ks <- sp$ks * temperatureScalar
    sp$h <- sp$h * temperatureScalar
    sp$gamma <- sp$gamma * temperatureScalar
    
    # We are choosing the smallest plankton size to agree with what Asta's
    # setup function produces
    min_w_pp <- 0.95e-10
    
    params <- newMultispeciesParams(
        sp, interaction = interaction, no_w = no_w, min_w_pp = min_w_pp,
        n = 2/3, r_pp = 1, kappa = kappa, lambda = lambda, w_pp_cutoff = 1)
    
    # Add starvation mortality with default parameter
    params <- setStarvation(params)
    
    # Add senescence mortality and rescale with temperature factor 
    z0 <- getExtMort(params) + 
        sen_mort(mariaParams, params, k.sm, xsw, sen.e)
    params <- setExtMort(params, z0 = z0 * temperatureScalar)
    
    # Add macroalgae 
    kappa <- 16
    lambda <- 1.6
    r <- 2
    n <- 2/3
    max <- 50
    min <- 1e-3
    rate <- r * params@w_full^(n - 1)
    capacity <- kappa * params@w_full^(-lambda)
    capacity[params@w_full > max] <- 0
    capacity[params@w_full < min] <- 0
    params <- setComponent(params = params, component = "aa",
                           initial_value = stable_alg,
                           dynamics_fun =  "background_semichemostat",
                           component_params = list(rate = rate,
                                                   capacity = capacity))
    
    # Add benthic resource
    kappa <- kappa_ben
    lambda <- lambda_ben
    r <- 1
    n <- 2/3
    max <- 5
    min <- 1e-3
    rate <- r * params@w_full^(n - 1)
    capacity <- kappa * params@w_full^(-lambda)
    capacity[params@w_full > max] <- 0
    capacity[params@w_full < min] <- 0
    params <- setComponent(params = params, component = "bb",
                           initial_value = stable_ben,
                           dynamics_fun =  "background_semichemostat",
                           component_params = list(rate = rate,
                                                   capacity = capacity))
    
    # Include these extra resources in the encounter rate
    params <- setRateFunction(params, "Encounter", "background_encounter")
    
    # Update initial abundances 
    params@initial_n <- stable_abund
    params@initial_n_pp <- stable_pl
    
    params
}


