# ==============================================================================
# ===== SARS-CoV-X Emergence: ABC Wuhan ========================================
# ==============================================================================

# ------------------------------------------------------------------------------
# ----- 0. Initialisation ------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 0.1. Description -------------------------------------------------------

# The following uses an ABC process to approximate the posterior distributions
# of the parameters R0_wuhan and recovery_wuhan.

# Priors:
#   - R0 parameters: lognormal with meanlog = 1, sdlog = 0.5
#   - recovery parameters: logit-normal with underlying normal mean = -1 and sd = 0.6
#   - waning parameters: logit-normal with underlying normal mean = -5 and sd = 2

# Run conditions:
#   - 2.5 million iterations recommended. It is advisible to parallelise this
#     process by duplicating this file and changing the file_prefix object. For
#     example, 100 duplicate files can be run in parallel with file_prefix of
#     1 to 100, each with 25000 iterations.

# File paths:
# This script uses the here() library to set paths dynamically. Path errors may occur
# if this script is opened in an already running Rstudio session. To resolve this,
# close Rstudio and reopen it by double-clicking on this file.

# ----- 0.2. Dependencies ------------------------------------------------------

file_prefix <- 1

Sys.sleep((file_prefix - 1) * 2.5)

library(odin.dust); library(tidyverse); library(socialmixr); library(abind); library(here)

# ----- 0.3 Convenience Functions ----------------------------------------------

add <- function(list, name, object){
  list[[name]] <- object
  list[order(names(list))]
}

# ----- 0.4. Load Odin Model ---------------------------------------------------

model_path <- here("models", "Odin_Model.R")
model <- odin_dust(model_path)

# ------------------------------------------------------------------------------
# ----- 1. Model Parameterisation ----------------------------------------------
# ------------------------------------------------------------------------------

# ----- 1.1. Iterations and Duration -------------------------------------------

parameters <- list()

parameters <- add(parameters, "model_iterations", 10)

parameters <- add(parameters, "model_duration", 450)

parameters <- add(parameters, "model_times", seq(1, parameters$model_duration, by = 1))

# ----- 1.2. Scotland Population Structure -------------------------------------

scotland <- list()

with(scotland, {
  groups <- read_csv(here("data", "pop_structure.csv"))$age_group
  values <- read_csv(here("data", "pop_structure.csv"))[,2:6]
  for (i in c(1:ncol(values))){
    new <- values[[i]]
    names(new) <- groups
    scotland <<- add(scotland, colnames(values)[i], new)
  }
  
  ages <- read_csv(here("data", "pop_individual_ages.csv"))$age
  n <- read_csv(here("data", "pop_individual_ages.csv"))$n_by_indiv_years
  names(n) <- ages
  
  scotland <<- add(scotland, "n_by_year", n)
  
})

scotland <<- add(scotland, "groups", list(
  "0-5" = c(1:5),
  "5-10" = c(6:10),
  "10-15" = c(11:15),
  "15-20" = c(16:20),
  "20-25" = c(21:25),
  "25-30" = c(26:30),
  "30-35" = c(31:35),
  "35-40" = c(36:40),
  "40-45" = c(41:45),
  "45-50" = c(46:50),
  "50-55" = c(51:55),
  "55-60" = c(56:60),
  "60-65" = c(61:65),
  "65-70" = c(66:70),
  "70-75" = c(71:75),
  "75+" = c(76:86)))

scotland <<- add(scotland, "n_by_year_cutoff", c(scotland$n_by_year[1:85], "85+" = sum(scotland$n_by_year[86:91])))

# ----- 1.3. Scotland Vacination Rates -----------------------------------------

scotland <- add(scotland, "vaccination_dose1", as.matrix(read_csv(here("data", "vaccination_dose1.csv"))[,2:17])[1:parameters$model_duration, ])
scotland <- add(scotland, "vaccination_dose2",
                as.matrix(read_csv(here("data", "vaccination_dose2.csv"))[,2:17])[1:parameters$model_duration, ])

# ----- 1.4. Scotland Social Structure -----------------------------------------

scotland <- add(scotland, "contact_matrices", readRDS(here("data", "contact_matrices.rds")))

# ----- 1.5 SARS-CoV-2 Variant Prevalences -------------------------------------

viruses <- list()

viruses <- add(viruses, "prevalence_wuhan", read_csv(here("data", "SARS2_variant_prevalences.csv"))$wuhan[1:parameters$model_duration])
viruses <- add(viruses, "prevalence_alpha", read_csv(here("data", "SARS2_variant_prevalences.csv"))$alpha[1:parameters$model_duration])
viruses <- add(viruses, "prevalence_delta", read_csv(here("data", "SARS2_variant_prevalences.csv"))$delta[1:parameters$model_duration])
viruses <- add(viruses, "prevalence_omicron", read_csv(here("data", "SARS2_variant_prevalences.csv"))$omicron[1:parameters$model_duration])

# ----- 1.6 SARS-CoV-2 Variant Phenotypes --------------------------------------

with(viruses, {
  phenotypes <- read_csv(here("data", "SARS2_infection_phenotypes_by_age.csv"))
  
  for (v in unique(phenotypes$virus)){
    current <- filter(phenotypes, virus == v)
    viruses <<- add(viruses, sprintf("u_death_%s", v), current$u_death)
    viruses <<- add(viruses, sprintf("u_incubation_%s", v), current$u_incubation)
    viruses <<- add(viruses, sprintf("disability_weighting_%s", v), current$disability_weighting)
    
  }
  
})

# ----- 1.7 SARS-CoV-2 Immunity Phenotypes -------------------------------------

immunity <- list()

with(immunity, {
  
  reinfection <- read_csv(here("data", "SARS2_immunity_reinfection_by_age.csv"))
  for(v in reinfection$virus){
    immunity <<- add(immunity, sprintf("immunity_reinfection_%s", v), filter(reinfection, virus == v)[,3:21] %>% as.matrix())
  }
  
  virulence <- read_csv(here("data", "SARS2_immunity_virulence_by_age.csv"))
  for(v in virulence$virus){
    immunity <<- add(immunity, sprintf("immunity_virulence_%s", v), filter(virulence, virus == v)[,3:21] %>% as.matrix())
  }
  
  waning <- read_csv(here("data", "SARS2_waning_rates.csv"))
  
  for (i in c(1:nrow(waning))){
    immunity <<- add(immunity, sprintf("u_waning_%s", waning$immunity[[i]]), rep(waning$u_waning[[i]], 16))
  }
  
})

# ----- 1.8 Immunity by Age Functions ------------------------------------------

immune_effectiveness <- function(z) {
  function(x){
    y <- 58.100686 / (1 + exp(-(0.031678) * (x - 7.635491))) - 
      48.95726 / (1 + exp(-0.11386 * (x - 83.1))) + 
      -42.46706 + 1.05007 * z
    pmax(pmin(y, 100), 0)
  }
}

find_z_immune_effectiveness <- function(target){
  if (target == 1){target <- target - 1e-8}
  effectiveness_auc <- function(z) {
    (integrate(immune_effectiveness(z), lower = 0, upper = 80)$value / 80) / 100
  }
  effectiveness_delta <- function(z) {
    current_auc <- effectiveness_auc(z)
    current_auc - target
  }
  effectiveness_z <- uniroot(effectiveness_delta, interval = c(-25, 125), tol = 0.0001)$root
  return(effectiveness_z)
}

# ----- 1.9 Infection-related Mortality Rate Function -------------------------

virus_death_rate <- function(x, z) {
  
  z <- z * 100
  
  rate <- 0.0224086 / (1 + exp(-0.0975006 * (x - (95.0328598 - 3.4256239 * z))))
  rate <- pmax(pmin(rate, 1), 0)
  return(rate)
}

# ----- 1.10 Transmission Rate Function ----------------------------------------

transmission_rate <- function(R0, incubation_rate, recovery_rate, mortality_natural, mortality_infection, contact_rate) {
  rate <-
    ((R0 * (incubation_rate + mortality_natural) *
       (recovery_rate + mortality_natural + mortality_infection)) / incubation_rate) / contact_rate
  return(rate)
}

# ------------------------------------------------------------------------------
# ----- 2. Model Run -----------------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1 Priors -------------------------------------------------------------

prior_R0_wuhan <- c(mean = 1, sd = 0.5)
prior_R0_alpha <- c(mean = 1, sd = 0.5)
prior_R0_delta <- c(mean = 1, sd = 0.5)
prior_R0_omicron <- c(mean = 1, sd = 0.5)

prior_recovery_wuhan <- c(mean = -1, sd = 0.6)
prior_recovery_alpha <- c(mean = -1, sd = 0.6)
prior_recovery_delta <- c(mean = -1, sd = 0.6)
prior_recovery_omicron <- c(mean = -1, sd = 0.6)

prior_waning_wuhan <- c(meanlog = -5, sdlog = 2)
prior_waning_alpha <- c(meanlog = -5, sdlog = 2)
prior_waning_delta <- c(meanlog = -5, sdlog = 2)
prior_waning_omicron <- c(meanlog = -5, sdlog = 2)

iterations <- 2500000

# ----- 2.2 Model Execution ----------------------------------------------------

if(exists(here("models", "ABC", "ABC_samples_Wuhan"))){
  rm(here("models", "ABC", "ABC_samples_Wuhan"))
}

dir.create(here("models", "ABC", "ABC_samples_Wuhan"))
  
for (s in (1:iterations)){
  
  if(exists("prevalences")){
    rm(prevalences)
  }
  
  # ----- Prior Sampling -----
  
  R0_wuhan <- rlnorm(1, meanlog = prior_R0_wuhan[1], sdlog = prior_R0_wuhan[2])
  R0_alpha <- rlnorm(1, meanlog = prior_R0_alpha[1], sdlog = prior_R0_alpha[2])
  R0_delta <- rlnorm(1, meanlog = prior_R0_delta[1], sdlog = prior_R0_delta[2])
  R0_omicron <- rlnorm(1, meanlog = prior_R0_omicron[1], sdlog = prior_R0_omicron[2])
  
  recovery_wuhan <- plogis(rnorm(1, mean = prior_recovery_wuhan[1], sd = prior_recovery_wuhan[2]))
  recovery_alpha <- plogis(rnorm(1, mean = prior_recovery_alpha[1], sd = prior_recovery_alpha[2]))
  recovery_delta <- plogis(rnorm(1, mean = prior_recovery_delta[1], sd = prior_recovery_delta[2]))
  recovery_omicron <- plogis(rnorm(1, mean = prior_recovery_omicron[1], sd = prior_recovery_omicron[2]))
  
  waning_wuhan <- plogis(rnorm(1, mean = prior_waning_wuhan[1], sd = prior_waning_wuhan[2]))
  waning_alpha <- plogis(rnorm(1, mean = prior_waning_alpha[1], sd = prior_waning_alpha[2]))
  waning_delta <- plogis(rnorm(1, mean = prior_waning_delta[1], sd = prior_waning_delta[2]))
  waning_omicron <- plogis(rnorm(1, mean = prior_waning_omicron[1], sd = prior_waning_omicron[2]))

    # ----- Transmission Rate Calculation -----
  
  death_natural <- 10.84157 / (1000 * 365)
  
  u_transmission_wuhan <- transmission_rate(R0_wuhan, viruses$u_incubation_wuhan, recovery_wuhan, death_natural, viruses$u_death_wuhan, mean(scotland$contact_matrices[,,1:1200]) * sum(scotland$n))
  u_transmission_alpha <- transmission_rate(R0_alpha, viruses$u_incubation_alpha, recovery_alpha, death_natural, viruses$u_death_alpha, mean(scotland$contact_matrices[,,1:1200]) * sum(scotland$n))
  u_transmission_delta <- transmission_rate(R0_delta, viruses$u_incubation_delta, recovery_delta, death_natural, viruses$u_death_delta, mean(scotland$contact_matrices[,,1:1200]) * sum(scotland$n))
  u_transmission_omicron <- transmission_rate(R0_omicron, viruses$u_incubation_omicron, recovery_omicron, death_natural, viruses$u_death_omicron, mean(scotland$contact_matrices[,,1:1200]) * sum(scotland$n))

print(mean(scotland$contact_matrices))
  
  # ----- Odin Model Call -----
  
  cat("\n")
  print(sprintf("Iteration: %s", s))
  
  pb <- txtProgressBar(min = 0, max = parameters$model_iterations, style = 3, width = 10)
  
  for (trial in c(1:parameters$model_iterations)){
    
    model_run <- model$new(pars = list(blank = rep(0, 16),
                                       duration = parameters$model_duration,
                                       N_age = 16,
                                       contact_matrices = scotland$contact_matrices[,,1:parameters$model_duration],
                                       S_ini = scotland$n,
                                       u_birth = scotland$u_birth,
                                       u_death = scotland$u_death,
                                       u_aging = 1/(365*5),
                                       n_migration = scotland$n_migration,
                                       
                                       u_V1 = scotland$vaccination_dose1,
                                       u_V2 = scotland$vaccination_dose2,
                                       
                                       u_waning_V1 = immunity$u_waning_V1,
                                       u_waning_V2 = immunity$u_waning_V2,
                                       
                                       prevalence_W = viruses$prevalence_wuhan,
                                       prevalence_A = viruses$prevalence_alpha,
                                       prevalence_D = viruses$prevalence_delta,
                                       prevalence_O = viruses$prevalence_omicron,
                                       
                                       u_trans_W = u_transmission_wuhan,
                                       u_trans_A = u_transmission_alpha,
                                       u_trans_D = u_transmission_delta,
                                       u_trans_O = u_transmission_omicron,
                                       u_trans_X = rep(0,16),
                                       
                                       immunity_inf_W = immunity$immunity_reinfection_wuhan,
                                       immunity_inf_A = immunity$immunity_reinfection_alpha,
                                       immunity_inf_D = immunity$immunity_reinfection_delta,
                                       immunity_inf_O = immunity$immunity_reinfection_omicron,
                                       immunity_inf_X = immunity$immunity_reinfection_wuhan,
                                       
                                       u_death_W = viruses$u_death_wuhan,
                                       u_death_A = viruses$u_death_alpha,
                                       u_death_D = viruses$u_death_delta,
                                       u_death_O = viruses$u_death_omicron,
                                       u_death_X = rep(0, 16),
                                       
                                       immunity_vir_W = immunity$immunity_virulence_wuhan,
                                       immunity_vir_A = immunity$immunity_virulence_alpha,
                                       immunity_vir_D = immunity$immunity_virulence_delta,
                                       immunity_vir_O = immunity$immunity_virulence_omicron,
                                       immunity_vir_X = immunity$immunity_virulence_wuhan,
                                       
                                       u_incub_W = viruses$u_incubation_wuhan,
                                       u_incub_A = viruses$u_incubation_alpha,
                                       u_incub_D = viruses$u_incubation_delta,
                                       u_incub_O = viruses$u_incubation_omicron,
                                       u_incub_X = rep(0, 16),
                                       
                                       u_recov_W = rep(recovery_wuhan,16),
                                       u_recov_A = rep(recovery_alpha,16),
                                       u_recov_D = rep(recovery_delta,16),
                                       u_recov_O = rep(recovery_omicron,16),
                                       u_recov_X = rep(0, 16),
                                       
                                       u_waning_W = rep(waning_wuhan,16),
                                       u_waning_A = rep(waning_alpha,16),
                                       u_waning_D = rep(waning_delta,16),
                                       u_waning_O = rep(waning_omicron,16),
                                       u_waning_X = rep(0,16),
                                       time_intro_W = 1,
                                       amount_intro_W = rep(5,16),
                                       time_intro_X = 0,
                                       amount_intro_X = c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)),
                           time = 1L,
                           n_particles = 1,
                           n_threads = 1L)
    
    output <- array(NA, dim = c(model_run$info()$len, 1, parameters$model_duration))
    
    for (i in c(1:(parameters$model_duration))){
      output[ , , i] <- model_run$run(i)
    }
    
    prevalence <- output[2,,]/output[4,,]
    prevalence <- ifelse(is.finite(prevalence), prevalence, 0)
    
    prev <- data.frame(trial = trial,
                       prevalence = paste(prevalence, collapse = ","))
    
    if(!exists("prevalences")){
      prevalences <- prev
    } else {
      prevalences <- rbind(prevalences, prev)
    }
    
    setTxtProgressBar(pb, trial)
    
  }
  
  virus <- data.frame(R0_wuhan = R0_wuhan,
                      R0_alpha = R0_alpha,
                      R0_delta = R0_delta,
                      R0_omicron = R0_omicron,
                      recovery_wuhan = recovery_wuhan,
                      recovery_alpha = recovery_alpha,
                      recovery_delta = recovery_delta,
                      recovery_omicron = recovery_omicron,
                      waning_wuhan = waning_wuhan,
                      waning_alpha = waning_alpha,
                      waning_delta = waning_delta,
                      waning_omicron = waning_omicron,
                      prevalence = paste(prevalences$prevalence, collapse = ";"))
  
  cat("\n")
  print(virus[,c(1:12)])
  
  outfile <- here("models", "ABC", "ABC_samples_Wuhan", sprintf("/%s_%s.csv", file_prefix, s))
  
  write_csv(virus, outfile)
  
}
