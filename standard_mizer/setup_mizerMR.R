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
library(mizerMR)

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
        n = 2/3)
    
    # Add starvation mortality with default parameter
    params <- setStarvation(params)
    
    # Add senescence mortality and rescale with temperature factor 
    z0 <- getExtMort(params) + 
        sen_mort(mariaParams, params, k.sm, xsw, sen.e)
    params <- setExtMort(params, z0 = z0 * temperatureScalar)
    
    # Set up resources
    resource_params(params) <- data.frame(
        resource = c("pl", "aa", "bb"),
        lambda = c(lambda, 1.6, lambda_ben),
        kappa = c(kappa, 16, kappa_ben),
        r_pp = c(1, 2, 1),
        w_min = c(NA, 1e-3, 1e-3),
        w_max = c(1, 50, 5)
    )
    resource_interaction(params)[, 1] <- sp$avail_PP
    resource_interaction(params)[, 2] <- sp$avail_AA
    resource_interaction(params)[, 3] <- sp$avail_BB
    
    # Update initial abundances 
    initialN(params) <- stable_abund
    initialNResource(params)[1, ] <- stable_pl
    initialNResource(params)[2, ] <- stable_alg
    initialNResource(params)[3, ] <- stable_ben
    
    params
}
