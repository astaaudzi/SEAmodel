# In this file we are running the climate change scenarios
# 
# The code is mostly copied from the "Actual scenario runs" section of
# Asta's "BenPel_SEAmodel_dbV2.Rmd" file, just slightly adjusted to run
# with standard mizer. Places where I made importantchanges are marked with 
# comments starting with GWD

# We source the functions needed to set up Asta's model in standard mizer
source("standard_mizer/setup.R")


# Species params ----
# GWD: We have to make some changes to Asta's species parameter file:
# - Change in names of resource interaction parameters
# - Absorb the growth efficiency into the existing parameters
load(file = "../modelParams/mariaParamsMs.RData") #species parameter file 
mariaParams <- mariaParams %>% 
    rename(interaction_resource = avail_PP,
           interaction_aa = avail_AA,
           interaction_bb = avail_BB) %>% 
    # get rid of separate growth efficiency
    mutate(alpha = alpha * alpha_g,
           ks = ks * alpha_g,
           erepro = erepro / alpha_g)

# GWD: We change the maturity curve to correspond to an exponent u = 5 
mariaParams$w_mat25 <- mariaParams$w_mat/3^(1/5)

# We load interaction matrix and add rownames
load(file = "../modelParams/inter_N19.RData")
dimnames(inter)[[1]] <- dimnames(inter)[[2]]

# Load initial values
load(file = "../modelParams/naa_N19.RData")
load(file = "../modelParams/nbb_N19.RData")
load(file = "../modelParams/abund_N19.RData")
load(file = "../modelParams/npp_N19.RData")

# Make a rund with the base scenario just to check that things are working
# 
# params <- newAstaParams(mariaParams, interaction = inter,
#                         stable_abund = stable_abund,
#                         stable_pl = stable_pl,
#                         stable_ben = stable_ben,
#                         stable_alg = stable_alg)
# 
# tasm1 <- project(params, t_max = 10, effort = 0, dt = 0.2)
# plot(tasm1)


### Load scenarios ----

load(file = "../modelParams/params28ms.RData")
accepted1 <- as.data.frame(params28)

# List of seven productivity scenarios. Each of them will be run 4 times: with
# and without fishing, with and without heating. Note, I used very approximate
# array of fishing mortality (shown below). This can be updated

#if no fishing is applied 
fish.mort.0 <- rep(0, times = 17)
#a vector of low fishing mortalities
fish.mort.v <-  c(0.15, 0.15, 0, 0.1, 0.15, 0.1, 0.15, 0.1, 0.1, 0.15, 0.15, 
                  0.15, 0.15, 0.15, 0.1, 0.15, 0.15)


prod_scen<-list( 
    "baseline_12_f0"=       c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.0), 
    "more_plankt_12_f0"=    c("kappa" =2.6, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.0), 
    "less_plankt_12_f0"=    c("kappa" =1.5, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.0), 
    "small_plankt_12_f0"=   c("kappa" =2, "kappa_ben"=6, "lambda"=2.18, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.0), 
    "large_plankt_12_f0"=   c("kappa" =2, "kappa_ben"=6, "lambda"=2.12, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.0), 
    "more_benth_12_f0"=     c("kappa" =2, "kappa_ben"=9, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.0), 
    "less_benth_12_f0"=     c("kappa" =2, "kappa_ben"=4, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.0), 
    "small_benth_12_f0"=    c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=2.0, "temp" = 12, "fmort" = fish.mort.0),
    "large_benth_12_f0"=    c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.8, "temp" = 12, "fmort" = fish.mort.0),
    
    "baseline_14_f0"=       c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.0), 
    "more_plankt_14_f0"=    c("kappa" =2.6, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.0), 
    "less_plankt_14_f0"=    c("kappa" =1.5, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.0), 
    "small_plankt_14_f0"=   c("kappa" =2, "kappa_ben"=6, "lambda"=2.18, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.0), 
    "large_plankt_14_f0"=   c("kappa" =2, "kappa_ben"=6, "lambda"=2.12, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.0), 
    "more_benth_14_f0"=     c("kappa" =2, "kappa_ben"=9, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.0), 
    "less_benth_14_f0"=     c("kappa" =2, "kappa_ben"=4, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.0), 
    "small_benth_14_f0"=    c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=2.0, "temp" = 14.5, "fmort" = fish.mort.0), 
    "large_benth_14_f0"=    c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.8, "temp" = 14.5, "fmort" = fish.mort.0), 
    
    "baseline_12_f2"=       c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.v), 
    "more_plankt_12_f2"=    c("kappa" =2.6, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.v), 
    "less_plankt_12_f2"=    c("kappa" =1.5, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.v), 
    "small_plankt_12_f2"=   c("kappa" =2, "kappa_ben"=6, "lambda"=2.18, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.v), 
    "large_plankt_12_f2"=   c("kappa" =2, "kappa_ben"=6, "lambda"=2.12, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.v), 
    "more_benth_12_f2"=     c("kappa" =2, "kappa_ben"=9, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.v), 
    "less_benth_12_f2"=     c("kappa" =2, "kappa_ben"=4, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 12, "fmort" = fish.mort.v), 
    "small_benth_12_f2"=    c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=2.0, "temp" = 12, "fmort" = fish.mort.v),
    "large_benth_12_f2"=    c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.8, "temp" = 12, "fmort" = fish.mort.v),
    
    "baseline_14_f2"=       c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.v), 
    "more_plankt_14_f2"=    c("kappa" =2.6, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.v), 
    "less_plankt_14_f2"=    c("kappa" =1.5, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.v), 
    "small_plankt_14_f2"=   c("kappa" =2, "kappa_ben"=6, "lambda"=2.18, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.v), 
    "large_plankt_14_f2"=   c("kappa" =2, "kappa_ben"=6, "lambda"=2.12, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.v), 
    "more_benth_14_f2"=     c("kappa" =2, "kappa_ben"=9, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.v), 
    "less_benth_14_f2"=     c("kappa" =2, "kappa_ben"=4, "lambda"=2.15, "lambda_ben"=1.9, "temp" = 14.5, "fmort" = fish.mort.v), 
    "small_benth_14_f2"=    c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=2.0, "temp" = 14.5, "fmort" = fish.mort.v), 
    "large_benth_14_f2"=    c("kappa" =2, "kappa_ben"=6, "lambda"=2.15, "lambda_ben"=1.8, "temp" = 14.5, "fmort" = fish.mort.v)) 

### Setup run parameters ----

## do you want a plot for each run, while the simulations are going? This takes
## more time but is sometimes fun to watch and keep an eye if there are
## oscillating dynamics
showplot = F

## setup run time 
tmax = 100
#number of size groups 
no_size_groups = 200
#timestep used in the integration 
dt = 0.2
#our reference temperature 
temp_ref = 12

## fixed background params- THESE ARE NOT CHANGING IN SCENARIOS
w_pp_cutoff = 1 #g
r_pp = 1 # rate of regeneration
min_w_pp = 1e-10 #g

w_bb_cutoff = 5 #
r_bb = 1 # 
min_w_bb = 0.001 # 0.01

#algae are not really size structured food resources, but are treated as such here for now 
kappa_alg = 16 #intercept assuming g/m2
lambda_alg = 1.6 #slope is much more shallow, to allow for lots of large kelp, but size structure does not really make much sense here
w_aa_cutoff = 50 
r_aa = 2 
min_w_aa = 0.001

### Setup arrays to store data ----

# We will store data on: 
# 1) Numbers at size to look at species spectra and so on
# 2) Numbers at size for all background spectra, just in case we need it
# 3) Proportion of Large fish indicator 
# 4) Mean weight above some thresholds: two thresholds used are MeanWeight of 
#    mature fish (maturation size is a threshold) and MeanWeight of fish above 
#    5cm in length
# 5) Calculate community slopes: for all species, for fish only 
#    (exclude lobsters and urchins), and for four functional groups 
#    separately - benthivores, planktivores, herbivores (includes urchins) 
#    and predators (includes lobsters)

#number of scenarios in a list
scen_num <- length(lengths(prod_scen))
n_spp <- length(params@species_params$species)

## 200, 326 or 120 in these dataframes refers to the numbers of size groups 
# in different arrays
numb_all <-array(data=NA,c(scen_num,n_spp,200,dim(params28)[1])) 
biom_all <- array(data=NA, c(scen_num,n_spp, dim(params28)[1]))
plankt_all <- array(data=NA, c(scen_num,326, dim(params28)[1]))
bent_all <- array(data=NA, c(scen_num,120, dim(params28)[1]))
alg_all <- array(data=NA, c(scen_num,153, dim(params28)[1]))
PropLarFish_all <- array(data= NA, c(scen_num,n_spp, dim(params28)[1]))
Yield_all <- array(data= NA, c(scen_num,n_spp, dim(params28)[1]))
#Yield_all[c(1:14),,] <- 0 #no fishing scenarios will have zero yield naturally
MeanWgtMat_all <- array(data= NA, c(scen_num,n_spp, dim(params28)[1]))
#MeanWgtHalfMat_all <- array(data= NA, c(scen_num,n_spp, dim(params28)[1]))
MeanWgtAbove5_all <- array(data= NA, c(scen_num,n_spp, dim(params28)[1]))
MeanWgtAbove10_all <- array(data= NA, c(scen_num,n_spp, dim(params28)[1]))
CommSlope_all <- array(data= NA, c(scen_num,3,dim(params28)[1]))
CommSlopeFish_all <- array(data= NA, c(scen_num,3,dim(params28)[1]))
CommSlopeBent_all <- array(data= NA, c(scen_num,3,dim(params28)[1]))
CommSlopePred_all <- array(data= NA, c(scen_num,3,dim(params28)[1]))
CommSlopePlan_all <- array(data= NA, c(scen_num,3,dim(params28)[1]))
CommSlopeHerb_all <- array(data= NA, c(scen_num,3,dim(params28)[1]))

## setup different species groups for community slope calcualtions
fishonly <- as.character(mariaParams$species[c(1:15)])
bentiv <- as.character(mariaParams$species[which(mariaParams$funcgr == "omni")])
plankt <- as.character(mariaParams$species[which(mariaParams$funcgr == "plankt")])
predat <- as.character(mariaParams$species[which(mariaParams$funcgr == "predat")])
herbiv <- as.character(mariaParams$species[which(mariaParams$funcgr == "herbi")])


### All simulations ----

# This will run 36 scenarios with 29 parameter combinations each and save
# outputs into data arrays

## loop through 36 scenarios

for (scen in 1:scen_num) {
    
    ## resource params
    kappa = as.numeric(prod_scen[[scen]]["kappa"])
    kappa_ben = as.numeric(prod_scen[[scen]]["kappa_ben"])  
    lambda = as.numeric(prod_scen[[scen]]["lambda"])
    lambda_ben = as.numeric(prod_scen[[scen]]["lambda_ben"])
    
    ## get temperature value for the scenario 
    temp_run <- as.numeric(prod_scen[[scen]]["temp"])
    
    ## for each scenario loop through 29 parameter combinations   
    for (iter in 1:dim(params28)[1]) {
        
        print("scenario = ")
        print (scen)
        print(prod_scen[scen])
        print("iteration out of 29")
        print(iter)  
        ## update parameters
        
        mariaParams$r_max <- as.numeric(accepted1[iter,c(1:17)])
        mariaParams$gamma <- as.numeric(accepted1[iter,c(18:34)])
        availUr <- as.numeric(accepted1[iter,35])
        availUrLob <- as.numeric(accepted1[iter,36])
        availSchooling <- as.numeric(accepted1[iter,37])
        
        # overwrite default interaction values with this (only for predators that feed on fish, i.e. have values > 0)
        inter[c(which(inter[,which(mariaParams$species == "T_caudimaculatus")] >0)),which(mariaParams$species == "T_caudimaculatus")] <- availSchooling
        inter[c(which(inter[,which(mariaParams$species == "P_laticlavius")] >0)),which(mariaParams$species == "P_laticlavius")] <- availSchooling
        inter[c(which(inter[,which(mariaParams$species == "C_rasor")] >0)),which(mariaParams$species == "C_rasor")] <- availSchooling
        inter[c(which(inter[,which(mariaParams$species == "urchins")] >0)),which(mariaParams$species == "urchins")] <- availUr
        inter[which(mariaParams$species == "lobsters"), which(mariaParams$species == "urchins")] <- availUrLob
        
        
        ### setup again with new parameter values
        # GWD: We use our own setup function here, which already
        # implements the temperature rescalings as well as setting up
        # the resource parameters, adding the senescence mortality,
        # and setting initial abundances.
        params <- newAstaParams(mariaParams, 
                                interaction = inter,
                                temperature = temp_run,
                                no_w = no_size_groups, 
                                kappa = kappa, 
                                lambda = lambda, 
                                kappa_ben = kappa_ben, 
                                lambda_ben = lambda_ben,
                                stable_abund = stable_abund,
                                stable_pl = stable_pl,
                                stable_ben = stable_ben,
                                stable_alg = stable_alg)
        # GWD: Rather than having one gear per species and setting a different
        # effort for each, we set the fishing mortality via the catchability
        # and then use effort 1.
        # set catchability
        fish.mort <- array(as.numeric(prod_scen[[scen]][6:22]), dim = c(1, 17))
        # fish.mort <- c(0.3, 0.2, 0, 0.1, 0.2, 0.1, 0.3, 0.1, 0.1, 0.2, 0.2, 0.1, 0.1, 0.2, 0.1, 0, 0.3)
        params@catchability[] <- fish.mort
        
        #run the model
        tmax <- 100
        dt <- 0.2
        # GWD: use standard project function. We don't need to supply
        # temperature here because we have already rescaled the rates
        # when we set up the params object.
        tasm1 <- project(params, t_max = tmax, effort = 1, dt = dt)
        
        ### Calculate various statistics inside the run to avoid saving massive model objects
        if (showplot == T) {
            plot(tasm1)
        }
        #get relative biomasses 
        #need to take an average of the last 30 years because it is oscilating
        biomass <- apply((getBiomass(tasm1)[c((tmax-29):tmax),]),2,mean)
        #biomass <- getBiomass(tasm1)[tmax,]
        #numbers <- tasm1@n[tmax,,] 
        numbers <- apply((tasm1@n[c((tmax-29):tmax),,]),c(2,3),mean)
        
        PropLarFish <- rep(NA,17)
        MeanWgtMat <- rep(NA,17)
        MeanWgtAbove5 <- rep(NA,17)
        MeanWgtAbove10 <- rep(NA,17)
        
        for (xx in 1:length(mariaParams$species)) {
            PropLarFish[xx] <- mean(getProportionOfLargeFish(tasm1,species = xx, threshold_w = mariaParams$w_mat[xx])[(tmax-29):tmax])
            MeanWgtMat[xx]  <- mean(getMeanWeight(tasm1, min_w = mariaParams$w_mat[xx], max_w = mariaParams$w_inf[xx])[(tmax-29):tmax])
            MeanWgtAbove5[xx]  <- mean(getMeanWeight(tasm1, min_w = mariaParams$cm5[xx], max_w = mariaParams$w_inf[xx])[(tmax-29):tmax])
            MeanWgtAbove10[xx]  <- mean(getMeanWeight(tasm1, min_w = mariaParams$cm10[xx], max_w = mariaParams$w_inf[xx])[(tmax-29):tmax])
            
        }
        
        #get mean values of various statistics over the last 30 years and save them in arrays 
        ttemp <- tasm1@n_pp[c((tmax-29):tmax),which(params@w_full > min_w_pp & params@w_full < w_pp_cutoff)]
        plankt_all[scen,,iter] <- apply(ttemp, 2, mean)
        
        # GWD changes to next 4 lines because extra resources are saved differently
        ttemp <- Reduce("+", tasm1@n_other[(tmax-29):tmax, "bb"])[which(params@w_full > min_w_bb & params@w_full < w_bb_cutoff)]
        bent_all[scen,,iter] <- ttemp
        ttemp <- Reduce("+", tasm1@n_other[(tmax-29):tmax, "aa"])[which(params@w_full > min_w_aa & params@w_full < w_aa_cutoff)]
        alg_all[scen,,iter] <- ttemp
        
        numb_all[scen,,,iter] <- numbers
        biom_all[scen,,iter] <- biomass
        PropLarFish_all[scen,,iter] <- PropLarFish
        MeanWgtMat_all[scen,,iter] <- MeanWgtMat
        MeanWgtAbove5_all[scen,,iter] <- MeanWgtAbove5
        MeanWgtAbove10_all[scen,,iter] <- MeanWgtAbove10
        Yield_all[scen,,iter] <- as.numeric(apply((getYield(tasm1)[c((tmax-29):tmax),]),2,mean))
        CommSlope_all[scen,,iter] <- as.numeric(apply((getCommunitySlope(tasm1)[c((tmax-29):tmax),]),2,mean))
        CommSlopeFish_all[scen,,iter] <- as.numeric(apply((getCommunitySlope(tasm1, species = fishonly)[c((tmax-29):tmax),]),2,mean))
        CommSlopeBent_all[scen,,iter] <- as.numeric(apply((getCommunitySlope(tasm1, species = bentiv)[c((tmax-29):tmax),]),2,mean))
        CommSlopePlan_all[scen,,iter] <- as.numeric(apply((getCommunitySlope(tasm1, species = plankt)[c((tmax-29):tmax),]),2,mean))
        CommSlopeHerb_all[scen,,iter] <- as.numeric(apply((getCommunitySlope(tasm1, species = herbiv)[c((tmax-29):tmax),]),2,mean))
        CommSlopePred_all[scen,,iter] <- as.numeric(apply((getCommunitySlope(tasm1, species = predat)[c((tmax-29):tmax),]),2,mean))
        
    }
    
}

### Save results ----

save(numb_all, file = "outputs/numb.RData")
save(biom_all, file = "outputs/biom.RData")
save(plankt_all, file = "outputs/plankt.RData")
save(bent_all, file = "outputs/bent.RData")
save(alg_all, file = "outputs/alg.RData")
save(PropLarFish_all, file = "outputs/PropLargFish.RData")
save(MeanWgtMat_all, file = "outputs/MeanWgtMat.RData")
save(MeanWgtAbove5_all, file = "outputs/MeanWgtAbove5.RData")
save(MeanWgtAbove10_all, file = "outputs/MeanWgtAbove10.RData")
save(Yield_all, file = "outputs/Yield.RData")
save(CommSlope_all, file = "outputs/CommSlope.RData")
save(CommSlopeFish_all, file = "outputs/CommSlopeFish.RData")
save(CommSlopeBent_all, file = "outputs/CommSlopeBent.RData")
save(CommSlopePlan_all, file = "outputs/CommSlopePlan.RData")
save(CommSlopeHerb_all, file = "outputs/CommSlopeHerb.RData")
save(CommSlopePred_all, file = "outputs/CommSlopePred.RData")

#saveRDS(foodWebStats, file = "../altRuns28/foodwebstats.rds")
