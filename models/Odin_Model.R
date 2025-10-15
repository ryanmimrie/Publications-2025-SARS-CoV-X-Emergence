# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~ SARS-X Emergence Odin Model ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# ------------------------------------------------------------------------------
# ----- 0. Initialization ------------------------------------------------------
# ------------------------------------------------------------------------------
# ----- 0.1. Description -------------------------------------------------------
# ----- 0.2. Initial Compartment Values ----------------------------------------

initial(tracker_IX)    <- 0
initial(tracker_I2)    <- 0
initial(tracker_V)     <- 0
initial(tracker_N)     <- 0

initial(deaths_IX[])    <- blank[i]
initial(days_IXS[])     <- blank[i]
initial(days_IXV1[])    <- blank[i]
initial(days_IXV2[])    <- blank[i]
initial(days_IXROS[])   <- blank[i]
initial(days_IXROV1[])  <- blank[i]
initial(days_IXROV2[])  <- blank[i]
initial(N[])            <- S_ini[i]
initial(S[])            <- S_ini[i]
initial(V1[])           <- blank[i]
initial(V2[])           <- blank[i]
initial(EWS[])          <- blank[i]
initial(IWS[])          <- blank[i]
initial(RWS[])          <- blank[i]
initial(EWV1[])         <- blank[i]
initial(IWV1[])         <- blank[i]
initial(RWV1[])         <- blank[i]
initial(EWV2[])         <- blank[i]
initial(IWV2[])         <- blank[i]
initial(RWV2[])         <- blank[i]
initial(EAS[])          <- blank[i]
initial(IAS[])          <- blank[i]
initial(RAS[])          <- blank[i]
initial(EARWS[])        <- blank[i]
initial(IARWS[])        <- blank[i]
initial(EAV1[])         <- blank[i]
initial(IAV1[])         <- blank[i]
initial(RAV1[])         <- blank[i]
initial(EARWV1[])       <- blank[i]
initial(IARWV1[])       <- blank[i]
initial(EAV2[])         <- blank[i]
initial(IAV2[])         <- blank[i]
initial(RAV2[])         <- blank[i]
initial(EARWV2[])       <- blank[i]
initial(IARWV2[])       <- blank[i]
initial(EDS[])          <- blank[i]
initial(IDS[])          <- blank[i]
initial(RDS[])          <- blank[i]
initial(EDRWS[])        <- blank[i]
initial(IDRWS[])        <- blank[i]
initial(EDRAS[])        <- blank[i]
initial(IDRAS[])        <- blank[i]
initial(EDV1[])         <- blank[i]
initial(IDV1[])         <- blank[i]
initial(RDV1[])         <- blank[i]
initial(EDRWV1[])       <- blank[i]
initial(IDRWV1[])       <- blank[i]
initial(EDRAV1[])       <- blank[i]
initial(IDRAV1[])       <- blank[i]
initial(EDV2[])         <- blank[i]
initial(IDV2[])         <- blank[i]
initial(RDV2[])         <- blank[i]
initial(EDRWV2[])       <- blank[i]
initial(IDRWV2[])       <- blank[i]
initial(EDRAV2[])       <- blank[i]
initial(IDRAV2[])       <- blank[i]
initial(EOS[])          <- blank[i]
initial(IOS[])          <- blank[i]
initial(ROS[])          <- blank[i]
initial(EORWS[])        <- blank[i]
initial(IORWS[])        <- blank[i]
initial(EORAS[])        <- blank[i]
initial(IORAS[])        <- blank[i]
initial(EORDS[])        <- blank[i]
initial(IORDS[])        <- blank[i]
initial(EOV1[])         <- blank[i]
initial(IOV1[])         <- blank[i]
initial(ROV1[])         <- blank[i]
initial(EORWV1[])       <- blank[i]
initial(IORWV1[])       <- blank[i]
initial(EORAV1[])       <- blank[i]
initial(IORAV1[])       <- blank[i]
initial(EORDV1[])       <- blank[i]
initial(IORDV1[])       <- blank[i]
initial(EOV2[])         <- blank[i]
initial(IOV2[])         <- blank[i]
initial(ROV2[])         <- blank[i]
initial(EORWV2[])       <- blank[i]
initial(IORWV2[])       <- blank[i]
initial(EORAV2[])       <- blank[i]
initial(IORAV2[])       <- blank[i]
initial(EORDV2[])       <- blank[i]
initial(IORDV2[])       <- blank[i]
initial(EXS[])          <- blank[i]
initial(IXS[])          <- blank[i]
initial(RXS[])          <- blank[i]
initial(EXV1[])         <- blank[i]
initial(IXV1[])         <- blank[i]
initial(RXV1[])         <- blank[i]
initial(EXV2[])         <- blank[i]
initial(IXV2[])         <- blank[i]
initial(RXV2[])         <- blank[i]
initial(EORXS[])        <- blank[i]
initial(IORXS[])        <- blank[i]
initial(EORXV1[])       <- blank[i]
initial(IORXV1[])       <- blank[i]
initial(EORXV2[])       <- blank[i]
initial(IORXV2[])       <- blank[i]
initial(EXROS[])        <- blank[i]
initial(IXROS[])        <- blank[i]
initial(EXROV1[])       <- blank[i]
initial(IXROV1[])       <- blank[i]
initial(EXROV2[])       <- blank[i]
initial(IXROV2[])       <- blank[i]
initial(ROXS[])         <- blank[i]
initial(ROXV1[])        <- blank[i]
initial(ROXV2[])        <- blank[i]

# ----- 0.3. User Provided Parameters ------------------------------------------

blank[]                <- user()
duration               <- user()
N_age                  <- user()
contact_matrices[,,]   <- user()
S_ini[]                <- user()
u_birth[]              <- user()
u_death[]              <- user()
u_aging                <- user()
n_migration[]          <- user()
u_V1[,]                <- user()
u_V2[,]                <- user()
u_waning_V1[]          <- user()
u_waning_V2[]          <- user()
prevalence_W[]         <- user()
prevalence_A[]         <- user()
prevalence_D[]         <- user()
prevalence_O[]         <- user()
u_trans_W[]            <- user()
u_trans_A[]            <- user()
u_trans_D[]            <- user()
u_trans_O[]            <- user()
u_trans_X[]            <- user()
immunity_inf_W[,]      <- user()
immunity_inf_A[,]      <- user()
immunity_inf_D[,]      <- user()
immunity_inf_O[,]      <- user()
immunity_inf_X[,]      <- user()
u_death_W[]            <- user()
u_death_A[]            <- user()
u_death_D[]            <- user()
u_death_O[]            <- user()
u_death_X[]            <- user()
immunity_vir_W[,]      <- user()
immunity_vir_A[,]      <- user()
immunity_vir_D[,]      <- user()
immunity_vir_O[,]      <- user()
immunity_vir_X[,]      <- user()
u_incub_W[]            <- user()
u_incub_A[]            <- user()
u_incub_D[]            <- user()
u_incub_O[]            <- user()
u_incub_X[]            <- user()
u_recov_W[]            <- user()
u_recov_A[]            <- user()
u_recov_D[]            <- user()
u_recov_O[]            <- user()
u_recov_X[]            <- user()
u_waning_W[]           <- user()
u_waning_A[]           <- user()
u_waning_D[]           <- user()
u_waning_O[]           <- user()
u_waning_X[]           <- user()
time_intro_W           <- user()
amount_intro_W[]       <- user()
time_intro_X           <- user()
amount_intro_X[]       <- user()

# ----- 0.4. Object Dimensions -------------------------------------------------
# ----- 0.4.1. User Parameters -----

dim(blank)            <- N_age
dim(contact_matrices) <- c(N_age,N_age,duration)
dim(S_ini)            <- N_age
dim(u_birth)          <- N_age
dim(u_death)          <- N_age
dim(n_migration)      <- N_age
dim(u_V1)             <- c(duration,N_age)
dim(u_V2)             <- c(duration,N_age)
dim(u_waning_V1)      <- N_age
dim(u_waning_V2)      <- N_age
dim(prevalence_W)     <- duration
dim(prevalence_A)     <- duration
dim(prevalence_D)     <- duration
dim(prevalence_O)     <- duration
dim(u_trans_W)      <- N_age
dim(u_trans_A)      <- N_age
dim(u_trans_D)      <- N_age
dim(u_trans_O)      <- N_age
dim(u_trans_X)      <- N_age
dim(immunity_inf_W)   <- c(N_age,19)
dim(immunity_inf_A)   <- c(N_age,19)
dim(immunity_inf_D)   <- c(N_age,19)
dim(immunity_inf_O)   <- c(N_age,19)
dim(immunity_inf_X)   <- c(N_age,19)
dim(u_death_W)      <- N_age
dim(u_death_A)      <- N_age
dim(u_death_D)      <- N_age
dim(u_death_O)      <- N_age
dim(u_death_X)      <- N_age
dim(immunity_vir_W)   <- c(N_age,19)
dim(immunity_vir_A)   <- c(N_age,19)
dim(immunity_vir_D)   <- c(N_age,19)
dim(immunity_vir_O)   <- c(N_age,19)
dim(immunity_vir_X)   <- c(N_age,19)
dim(u_incub_W)        <- N_age
dim(u_incub_A)        <- N_age
dim(u_incub_D)        <- N_age
dim(u_incub_O)        <- N_age
dim(u_incub_X)        <- N_age
dim(u_recov_W)        <- N_age
dim(u_recov_A)        <- N_age
dim(u_recov_D)        <- N_age
dim(u_recov_O)        <- N_age
dim(u_recov_X)        <- N_age
dim(u_waning_W)       <- N_age
dim(u_waning_A)       <- N_age
dim(u_waning_D)       <- N_age
dim(u_waning_O)       <- N_age
dim(u_waning_X)       <- N_age
dim(amount_intro_W)   <- N_age
dim(amount_intro_X)   <- N_age

# ----- 0.4.2. Trackers & Compartments -----

dim(deaths_IX)   <- N_age
dim(days_IXS)    <- N_age
dim(days_IXV1)   <- N_age
dim(days_IXV2)   <- N_age
dim(days_IXROS)  <- N_age
dim(days_IXROV1) <- N_age
dim(days_IXROV2) <- N_age
dim(N)           <- N_age
dim(S)           <- N_age
dim(V1)          <- N_age
dim(V2)          <- N_age
dim(EWS)         <- N_age
dim(IWS)         <- N_age
dim(RWS)         <- N_age
dim(EWV1)        <- N_age
dim(IWV1)        <- N_age
dim(RWV1)        <- N_age
dim(EWV2)        <- N_age
dim(IWV2)        <- N_age
dim(RWV2)        <- N_age
dim(EAS)         <- N_age
dim(IAS)         <- N_age
dim(RAS)         <- N_age
dim(EARWS)       <- N_age
dim(IARWS)       <- N_age
dim(EAV1)        <- N_age
dim(IAV1)        <- N_age
dim(RAV1)        <- N_age
dim(EARWV1)      <- N_age
dim(IARWV1)      <- N_age
dim(EAV2)        <- N_age
dim(IAV2)        <- N_age
dim(RAV2)        <- N_age
dim(EARWV2)      <- N_age
dim(IARWV2)      <- N_age
dim(EDS)         <- N_age
dim(IDS)         <- N_age
dim(RDS)         <- N_age
dim(EDRWS)       <- N_age
dim(IDRWS)       <- N_age
dim(EDRAS)       <- N_age
dim(IDRAS)       <- N_age
dim(EDV1)        <- N_age
dim(IDV1)        <- N_age
dim(RDV1)        <- N_age
dim(EDRWV1)      <- N_age
dim(IDRWV1)      <- N_age
dim(EDRAV1)      <- N_age
dim(IDRAV1)      <- N_age
dim(EDV2)        <- N_age
dim(IDV2)        <- N_age
dim(RDV2)        <- N_age
dim(EDRWV2)      <- N_age
dim(IDRWV2)      <- N_age
dim(EDRAV2)      <- N_age
dim(IDRAV2)      <- N_age
dim(EOS)         <- N_age
dim(IOS)         <- N_age
dim(ROS)         <- N_age
dim(EORWS)       <- N_age
dim(IORWS)       <- N_age
dim(EORAS)       <- N_age
dim(IORAS)       <- N_age
dim(EORDS)       <- N_age
dim(IORDS)       <- N_age
dim(EOV1)        <- N_age
dim(IOV1)        <- N_age
dim(ROV1)        <- N_age
dim(EORWV1)      <- N_age
dim(IORWV1)      <- N_age
dim(EORAV1)      <- N_age
dim(IORAV1)      <- N_age
dim(EORDV1)      <- N_age
dim(IORDV1)      <- N_age
dim(EOV2)        <- N_age
dim(IOV2)        <- N_age
dim(ROV2)        <- N_age
dim(EORWV2)      <- N_age
dim(IORWV2)      <- N_age
dim(EORAV2)      <- N_age
dim(IORAV2)      <- N_age
dim(EORDV2)      <- N_age
dim(IORDV2)      <- N_age
dim(EXS)         <- N_age
dim(IXS)         <- N_age
dim(RXS)         <- N_age
dim(EXV1)        <- N_age
dim(IXV1)        <- N_age
dim(RXV1)        <- N_age
dim(EXV2)        <- N_age
dim(IXV2)        <- N_age
dim(RXV2)        <- N_age
dim(EORXS)       <- N_age
dim(IORXS)       <- N_age
dim(EORXV1)      <- N_age
dim(IORXV1)      <- N_age
dim(EORXV2)      <- N_age
dim(IORXV2)      <- N_age
dim(EXROS)       <- N_age
dim(IXROS)       <- N_age
dim(EXROV1)      <- N_age
dim(IXROV1)      <- N_age
dim(EXROV2)      <- N_age
dim(IXROV2)      <- N_age
dim(ROXS)        <- N_age
dim(ROXV1)       <- N_age
dim(ROXV2)       <- N_age

# ----- 0.4.3. Mean Rates -----

dim(contacts_2) <- c(N_age,N_age)
dim(contacts_X) <- c(N_age,N_age)

dim(u_S_exposed_W) <- N_age
dim(u_S_exposed_A) <- N_age
dim(u_S_exposed_D) <- N_age
dim(u_S_exposed_O) <- N_age
dim(u_S_exposed_X)  <- N_age

dim(u_V1_exposed_W) <- N_age
dim(u_V1_exposed_A) <- N_age
dim(u_V1_exposed_D) <- N_age
dim(u_V1_exposed_O) <- N_age
dim(u_V1_exposed_X)  <- N_age

dim(u_V2_exposed_W) <- N_age
dim(u_V2_exposed_A) <- N_age
dim(u_V2_exposed_D) <- N_age
dim(u_V2_exposed_O) <- N_age
dim(u_V2_exposed_X)  <- N_age

dim(u_RWS_exposed_A) <- N_age
dim(u_RWS_exposed_D) <- N_age
dim(u_RAS_exposed_D) <- N_age
dim(u_RWS_exposed_O) <- N_age
dim(u_RAS_exposed_O) <- N_age
dim(u_RDS_exposed_O) <- N_age

dim(u_RWV1_exposed_A) <- N_age
dim(u_RWV1_exposed_D) <- N_age
dim(u_RAV1_exposed_D) <- N_age
dim(u_RWV1_exposed_O) <- N_age
dim(u_RAV1_exposed_O) <- N_age
dim(u_RDV1_exposed_O) <- N_age

dim(u_RWV2_exposed_A) <- N_age
dim(u_RWV2_exposed_D) <- N_age
dim(u_RAV2_exposed_D) <- N_age
dim(u_RWV2_exposed_O) <- N_age
dim(u_RAV2_exposed_O) <- N_age
dim(u_RDV2_exposed_O) <- N_age

dim(u_RXS_exposed_O)  <- N_age
dim(u_RXV1_exposed_O) <- N_age
dim(u_RXV2_exposed_O) <- N_age

dim(u_ROS_exposed_X)  <- N_age
dim(u_ROV1_exposed_X) <- N_age
dim(u_ROV2_exposed_X) <- N_age

dim(u_death_W_S) <- N_age
dim(u_death_A_S) <- N_age
dim(u_death_D_S) <- N_age
dim(u_death_O_S) <- N_age
dim(u_death_X_S) <- N_age

dim(u_death_W_V1) <- N_age
dim(u_death_A_V1) <- N_age
dim(u_death_D_V1) <- N_age
dim(u_death_O_V1) <- N_age
dim(u_death_X_V1) <- N_age

dim(u_death_W_V2) <- N_age
dim(u_death_A_V2) <- N_age
dim(u_death_D_V2) <- N_age
dim(u_death_O_V2) <- N_age
dim(u_death_X_V2) <- N_age

dim(u_death_A_RWS) <- N_age
dim(u_death_D_RWS) <- N_age
dim(u_death_D_RAS) <- N_age
dim(u_death_O_RWS) <- N_age
dim(u_death_O_RAS) <- N_age
dim(u_death_O_RDS) <- N_age

dim(u_death_A_RWV1) <- N_age
dim(u_death_D_RWV1) <- N_age
dim(u_death_D_RAV1) <- N_age
dim(u_death_O_RWV1) <- N_age
dim(u_death_O_RAV1) <- N_age
dim(u_death_O_RDV1) <- N_age

dim(u_death_A_RWV2) <- N_age
dim(u_death_D_RWV2) <- N_age
dim(u_death_D_RAV2) <- N_age
dim(u_death_O_RWV2) <- N_age
dim(u_death_O_RAV2) <- N_age
dim(u_death_O_RDV2) <- N_age

dim(u_death_X_ROS) <- N_age
dim(u_death_X_ROV1) <- N_age
dim(u_death_X_ROV2) <- N_age

dim(u_death_O_RXS) <- N_age
dim(u_death_O_RXV1) <- N_age
dim(u_death_O_RXV2) <- N_age

# ----- 0.4.4. Transition Probabilities -----
# ----- 0.4.4.1. Naive to Exposed -----

dim(p_S_to_EWS) <- N_age
dim(p_S_to_EAS) <- N_age
dim(p_S_to_EDS) <- N_age
dim(p_S_to_EOS) <- N_age
dim(p_S_to_EXS)  <- N_age

dim(p_V1_to_EWV1) <- N_age
dim(p_V1_to_EAV1) <- N_age
dim(p_V1_to_EDV1) <- N_age
dim(p_V1_to_EOV1) <- N_age
dim(p_V1_to_EXV1)  <- N_age

dim(p_V2_to_EWV2) <- N_age
dim(p_V2_to_EAV2) <- N_age
dim(p_V2_to_EDV2) <- N_age
dim(p_V2_to_EOV2) <- N_age
dim(p_V2_to_EXV2)  <- N_age

# ----- 0.4.4.2. Recovered to Exposed -----

dim(p_RWS_to_EARWS) <- N_age
dim(p_RWS_to_EDRWS) <- N_age
dim(p_RAS_to_EDRAS) <- N_age
dim(p_RWS_to_EORWS) <- N_age
dim(p_RAS_to_EORAS) <- N_age
dim(p_RDS_to_EORDS) <- N_age

dim(p_RWV1_to_EARWV1) <- N_age
dim(p_RWV1_to_EDRWV1) <- N_age
dim(p_RAV1_to_EDRAV1) <- N_age
dim(p_RWV1_to_EORWV1) <- N_age
dim(p_RAV1_to_EORAV1) <- N_age
dim(p_RDV1_to_EORDV1) <- N_age

dim(p_RWV2_to_EARWV2) <- N_age
dim(p_RWV2_to_EDRWV2) <- N_age
dim(p_RAV2_to_EDRAV2) <- N_age
dim(p_RWV2_to_EORWV2) <- N_age
dim(p_RAV2_to_EORAV2) <- N_age
dim(p_RDV2_to_EORDV2) <- N_age

dim(p_ROS_to_EXROS)   <- N_age
dim(p_ROV1_to_EXROV1) <- N_age
dim(p_ROV2_to_EXROV2) <- N_age

dim(p_RXS_to_EORXS)   <- N_age
dim(p_RXV1_to_EORXV1) <- N_age
dim(p_RXV2_to_EORXV2) <- N_age

# ----- 0.4.4.3. Exposed to Infected -----

dim(p_incub_W) <- N_age
dim(p_incub_A) <- N_age
dim(p_incub_D) <- N_age
dim(p_incub_O) <- N_age
dim(p_incub_X)  <- N_age

# ----- 0.4.4.4. Infected to Recovered -----

dim(p_recov_W) <- N_age
dim(p_recov_A) <- N_age
dim(p_recov_D) <- N_age
dim(p_recov_O) <- N_age
dim(p_recov_X)  <- N_age

# ----- 0.4.4.5. Recovered to Naive -----

dim(p_waning_W)  <- N_age
dim(p_waning_A)  <- N_age
dim(p_waning_D)  <- N_age
dim(p_waning_O)  <- N_age
dim(p_waning_X)   <- N_age

# ----- 0.4.4.6. Vaccination -----

dim(p_S_to_V1)  <- N_age
dim(p_V1_to_V2) <- N_age

# ----- 0.4.4.7. Vaccination Waning -----

dim(p_V1_to_S) <- N_age
dim(p_V2_to_S) <- N_age

# ----- 0.4.5.8. Births,Migrations -----

dim(p_birth) <- N_age

dim(p_migration_S)  <- N_age
dim(p_migration_V1) <- N_age
dim(p_migration_V2) <- N_age

dim(p_migration_EWS)  <- N_age
dim(p_migration_IWS)  <- N_age
dim(p_migration_RWS)  <- N_age
dim(p_migration_EWV1) <- N_age
dim(p_migration_IWV1) <- N_age
dim(p_migration_RWV1) <- N_age
dim(p_migration_EWV2) <- N_age
dim(p_migration_IWV2) <- N_age
dim(p_migration_RWV2) <- N_age

dim(p_migration_EAS)     <- N_age
dim(p_migration_IAS)     <- N_age
dim(p_migration_RAS)     <- N_age
dim(p_migration_EARWS)  <- N_age
dim(p_migration_IARWS)  <- N_age
dim(p_migration_EAV1)    <- N_age
dim(p_migration_IAV1)    <- N_age
dim(p_migration_RAV1)    <- N_age
dim(p_migration_EARWV1) <- N_age
dim(p_migration_IARWV1) <- N_age
dim(p_migration_EAV2)    <- N_age
dim(p_migration_IAV2)    <- N_age
dim(p_migration_RAV2)    <- N_age
dim(p_migration_EARWV2) <- N_age
dim(p_migration_IARWV2) <- N_age

dim(p_migration_EDS)     <- N_age
dim(p_migration_IDS)     <- N_age
dim(p_migration_RDS)     <- N_age
dim(p_migration_EDRWS)  <- N_age
dim(p_migration_IDRWS)  <- N_age
dim(p_migration_EDRAS)  <- N_age
dim(p_migration_IDRAS)  <- N_age
dim(p_migration_EDV1)    <- N_age
dim(p_migration_IDV1)    <- N_age
dim(p_migration_RDV1)    <- N_age
dim(p_migration_EDRWV1) <- N_age
dim(p_migration_IDRWV1) <- N_age
dim(p_migration_EDRAV1) <- N_age
dim(p_migration_IDRAV1) <- N_age
dim(p_migration_EDV2)    <- N_age
dim(p_migration_IDV2)    <- N_age
dim(p_migration_RDV2)    <- N_age
dim(p_migration_EDRWV2) <- N_age
dim(p_migration_IDRWV2) <- N_age
dim(p_migration_EDRAV2) <- N_age
dim(p_migration_IDRAV2) <- N_age

dim(p_migration_EOS)     <- N_age
dim(p_migration_IOS)     <- N_age
dim(p_migration_ROS)     <- N_age
dim(p_migration_EORWS)  <- N_age
dim(p_migration_IORWS)  <- N_age
dim(p_migration_EORAS)  <- N_age
dim(p_migration_IORAS)  <- N_age
dim(p_migration_EORDS)  <- N_age
dim(p_migration_IORDS)  <- N_age
dim(p_migration_EOV1)    <- N_age
dim(p_migration_IOV1)    <- N_age
dim(p_migration_ROV1)    <- N_age
dim(p_migration_EORWV1) <- N_age
dim(p_migration_IORWV1) <- N_age
dim(p_migration_EORAV1) <- N_age
dim(p_migration_IORAV1) <- N_age
dim(p_migration_EORDV1) <- N_age
dim(p_migration_IORDV1) <- N_age
dim(p_migration_EOV2)    <- N_age
dim(p_migration_IOV2)    <- N_age
dim(p_migration_ROV2)    <- N_age
dim(p_migration_EORWV2) <- N_age
dim(p_migration_IORWV2) <- N_age
dim(p_migration_EORAV2) <- N_age
dim(p_migration_IORAV2) <- N_age
dim(p_migration_EORDV2) <- N_age
dim(p_migration_IORDV2) <- N_age

dim(p_migration_EXS)  <- N_age
dim(p_migration_IXS)  <- N_age
dim(p_migration_RXS)  <- N_age
dim(p_migration_EXV1) <- N_age
dim(p_migration_IXV1) <- N_age
dim(p_migration_RXV1) <- N_age
dim(p_migration_EXV2) <- N_age
dim(p_migration_IXV2) <- N_age
dim(p_migration_RXV2) <- N_age

dim(p_migration_EORXS)  <- N_age
dim(p_migration_IORXS)  <- N_age
dim(p_migration_EORXV1) <- N_age
dim(p_migration_IORXV1) <- N_age
dim(p_migration_EORXV2) <- N_age
dim(p_migration_IORXV2) <- N_age

dim(p_migration_EXROS)  <- N_age
dim(p_migration_IXROS)  <- N_age
dim(p_migration_EXROV1) <- N_age
dim(p_migration_IXROV1) <- N_age
dim(p_migration_EXROV2) <- N_age
dim(p_migration_IXROV2) <- N_age

dim(p_migration_ROXS)   <- N_age
dim(p_migration_ROXV1)  <- N_age
dim(p_migration_ROXV2)  <- N_age

# ----- 0.4.5.9. Deaths -----

dim(p_death) <- N_age

dim(p_death_IWS) <- N_age
dim(p_death_IAS) <- N_age
dim(p_death_IDS) <- N_age
dim(p_death_IOS) <- N_age
dim(p_death_IXS)  <- N_age

dim(p_death_IWV1) <- N_age
dim(p_death_IAV1) <- N_age
dim(p_death_IDV1) <- N_age
dim(p_death_IOV1) <- N_age
dim(p_death_IXV1)  <- N_age

dim(p_death_IWV2) <- N_age
dim(p_death_IAV2) <- N_age
dim(p_death_IDV2) <- N_age
dim(p_death_IOV2) <- N_age
dim(p_death_IXV2)  <- N_age

dim(p_death_IARWS) <- N_age
dim(p_death_IDRWS) <- N_age
dim(p_death_IDRAS) <- N_age
dim(p_death_IORWS) <- N_age
dim(p_death_IORAS) <- N_age
dim(p_death_IORDS) <- N_age

dim(p_death_IARWV1) <- N_age
dim(p_death_IDRWV1) <- N_age
dim(p_death_IDRAV1) <- N_age
dim(p_death_IORWV1) <- N_age
dim(p_death_IORAV1) <- N_age
dim(p_death_IORDV1) <- N_age

dim(p_death_IARWV2) <- N_age
dim(p_death_IDRWV2) <- N_age
dim(p_death_IDRAV2) <- N_age
dim(p_death_IORWV2) <- N_age
dim(p_death_IORAV2) <- N_age
dim(p_death_IORDV2) <- N_age

dim(p_death_IXROS)  <- N_age
dim(p_death_IXROV1) <- N_age
dim(p_death_IXROV2) <- N_age

dim(p_death_IORXS)  <- N_age
dim(p_death_IORXV1) <- N_age
dim(p_death_IORXV2) <- N_age

# ----- 0.4.5. Transitions -----
# ----- 0.4.5.1. Naive to Exposed -----

dim(c_S_to_EWS) <- N_age
dim(c_S_to_EAS) <- N_age
dim(c_S_to_EDS) <- N_age
dim(c_S_to_EOS) <- N_age
dim(c_S_to_EXS)  <- N_age

dim(c_V1_to_EWV1) <- N_age
dim(c_V1_to_EAV1) <- N_age
dim(c_V1_to_EDV1) <- N_age
dim(c_V1_to_EOV1) <- N_age
dim(c_V1_to_EXV1)  <- N_age

dim(c_V2_to_EWV2) <- N_age
dim(c_V2_to_EAV2) <- N_age
dim(c_V2_to_EDV2) <- N_age
dim(c_V2_to_EOV2) <- N_age
dim(c_V2_to_EXV2)  <- N_age

# ----- 0.4.5.2. Recovered to Exposed -----

dim(c_RWS_to_EARWS) <- N_age
dim(c_RWS_to_EDRWS) <- N_age
dim(c_RAS_to_EDRAS) <- N_age
dim(c_RWS_to_EORWS) <- N_age
dim(c_RAS_to_EORAS) <- N_age
dim(c_RDS_to_EORDS) <- N_age

dim(c_RWV1_to_EARWV1) <- N_age
dim(c_RWV1_to_EDRWV1) <- N_age
dim(c_RAV1_to_EDRAV1) <- N_age
dim(c_RWV1_to_EORWV1) <- N_age
dim(c_RAV1_to_EORAV1) <- N_age
dim(c_RDV1_to_EORDV1) <- N_age

dim(c_RWV2_to_EARWV2) <- N_age
dim(c_RWV2_to_EDRWV2) <- N_age
dim(c_RAV2_to_EDRAV2) <- N_age
dim(c_RWV2_to_EORWV2) <- N_age
dim(c_RAV2_to_EORAV2) <- N_age
dim(c_RDV2_to_EORDV2) <- N_age

dim(c_ROS_to_EXROS)   <- N_age
dim(c_ROV1_to_EXROV1) <- N_age
dim(c_ROV2_to_EXROV2) <- N_age

dim(c_RXS_to_EORXS)   <- N_age
dim(c_RXV1_to_EORXV1) <- N_age
dim(c_RXV2_to_EORXV2) <- N_age

# ----- 0.4.5.3. Exposed to Infected -----

dim(c_EWS_to_IWS) <- N_age
dim(c_EAS_to_IAS) <- N_age
dim(c_EDS_to_IDS) <- N_age
dim(c_EOS_to_IOS) <- N_age
dim(c_EXS_to_IXS)   <- N_age

dim(c_EWV1_to_IWV1) <- N_age
dim(c_EAV1_to_IAV1) <- N_age
dim(c_EDV1_to_IDV1) <- N_age
dim(c_EOV1_to_IOV1) <- N_age
dim(c_EXV1_to_IXV1)   <- N_age

dim(c_EWV2_to_IWV2) <- N_age
dim(c_EAV2_to_IAV2) <- N_age
dim(c_EDV2_to_IDV2) <- N_age
dim(c_EOV2_to_IOV2) <- N_age
dim(c_EXV2_to_IXV2)   <- N_age

dim(c_EARWS_to_IARWS) <- N_age
dim(c_EDRWS_to_IDRWS) <- N_age
dim(c_EDRAS_to_IDRAS) <- N_age
dim(c_EORWS_to_IORWS) <- N_age
dim(c_EORAS_to_IORAS) <- N_age
dim(c_EORDS_to_IORDS) <- N_age

dim(c_EARWV1_to_IARWV1) <- N_age
dim(c_EDRWV1_to_IDRWV1) <- N_age
dim(c_EDRAV1_to_IDRAV1) <- N_age
dim(c_EORWV1_to_IORWV1) <- N_age
dim(c_EORAV1_to_IORAV1) <- N_age
dim(c_EORDV1_to_IORDV1) <- N_age

dim(c_EARWV2_to_IARWV2) <- N_age
dim(c_EDRWV2_to_IDRWV2) <- N_age
dim(c_EDRAV2_to_IDRAV2) <- N_age
dim(c_EORWV2_to_IORWV2) <- N_age
dim(c_EORAV2_to_IORAV2) <- N_age
dim(c_EORDV2_to_IORDV2) <- N_age

dim(c_EXROS_to_IXROS)   <- N_age
dim(c_EXROV1_to_IXROV1) <- N_age
dim(c_EXROV2_to_IXROV2) <- N_age

dim(c_EORXS_to_IORXS)   <- N_age
dim(c_EORXV1_to_IORXV1) <- N_age
dim(c_EORXV2_to_IORXV2) <- N_age

# ----- 0.4.5.4. Infected to Recovered -----

dim(c_IWS_to_RWS) <- N_age
dim(c_IAS_to_RAS) <- N_age
dim(c_IDS_to_RDS) <- N_age
dim(c_IOS_to_ROS) <- N_age
dim(c_IXS_to_RXS)   <- N_age

dim(c_IWV1_to_RWV1) <- N_age
dim(c_IAV1_to_RAV1) <- N_age
dim(c_IDV1_to_RDV1) <- N_age
dim(c_IOV1_to_ROV1) <- N_age
dim(c_IXV1_to_RXV1)   <- N_age

dim(c_IWV2_to_RWV2) <- N_age
dim(c_IAV2_to_RAV2) <- N_age
dim(c_IDV2_to_RDV2) <- N_age
dim(c_IOV2_to_ROV2) <- N_age
dim(c_IXV2_to_RXV2)   <- N_age

dim(c_IARWS_to_RAS) <- N_age
dim(c_IDRWS_to_RDS) <- N_age
dim(c_IDRAS_to_RDS) <- N_age
dim(c_IORWS_to_ROS) <- N_age
dim(c_IORAS_to_ROS) <- N_age
dim(c_IORDS_to_ROS) <- N_age

dim(c_IARWV1_to_RAV1) <- N_age
dim(c_IDRWV1_to_RDV1) <- N_age
dim(c_IDRAV1_to_RDV1) <- N_age
dim(c_IORWV1_to_ROV1) <- N_age
dim(c_IORAV1_to_ROV1) <- N_age
dim(c_IORDV1_to_ROV1) <- N_age

dim(c_IARWV2_to_RAV2) <- N_age
dim(c_IDRWV2_to_RDV2) <- N_age
dim(c_IDRAV2_to_RDV2) <- N_age
dim(c_IORWV2_to_ROV2) <- N_age
dim(c_IORAV2_to_ROV2) <- N_age
dim(c_IORDV2_to_ROV2) <- N_age

dim(c_IXROS_to_ROXS)   <- N_age
dim(c_IXROV1_to_ROXV1) <- N_age
dim(c_IXROV2_to_ROXV2) <- N_age

dim(c_IORXS_to_ROXS)   <- N_age
dim(c_IORXV1_to_ROXV1) <- N_age
dim(c_IORXV2_to_ROXV2) <- N_age

# ----- 0.4.5.5. Recovered to Naive -----

dim(c_RWS_to_S) <- N_age
dim(c_RAS_to_S) <- N_age
dim(c_RDS_to_S) <- N_age
dim(c_ROS_to_S) <- N_age
dim(c_RXS_to_S)  <- N_age

dim(c_RWV1_to_V1) <- N_age
dim(c_RAV1_to_V1) <- N_age
dim(c_RDV1_to_V1) <- N_age
dim(c_ROV1_to_V1) <- N_age
dim(c_RXV1_to_V1)  <- N_age

dim(c_RWV2_to_V2) <- N_age
dim(c_RAV2_to_V2) <- N_age
dim(c_RDV2_to_V2) <- N_age
dim(c_ROV2_to_V2) <- N_age
dim(c_RXV2_to_V2)  <- N_age

dim(c_ROXS_to_S)   <- N_age
dim(c_ROXV1_to_V1) <- N_age
dim(c_ROXV2_to_V2) <- N_age

# ----- 0.4.5.6. Vaccination -----

dim(c_S_to_V1)  <- N_age
dim(c_V1_to_V2) <- N_age

dim(c_EWS_to_EWV1)  <- N_age
dim(c_EWV1_to_EWV2) <- N_age
dim(c_EAS_to_EAV1)  <- N_age
dim(c_EAV1_to_EAV2) <- N_age
dim(c_EDS_to_EDV1)  <- N_age
dim(c_EDV1_to_EDV2) <- N_age
dim(c_EOS_to_EOV1)  <- N_age
dim(c_EOV1_to_EOV2) <- N_age

dim(c_EXS_to_EXV1)  <- N_age
dim(c_EXV1_to_EXV2) <- N_age

dim(c_IWS_to_IWV1)  <- N_age
dim(c_IWV1_to_IWV2) <- N_age
dim(c_IAS_to_IAV1)  <- N_age
dim(c_IAV1_to_IAV2) <- N_age
dim(c_IDS_to_IDV1)  <- N_age
dim(c_IDV1_to_IDV2) <- N_age
dim(c_IOS_to_IOV1)  <- N_age
dim(c_IOV1_to_IOV2) <- N_age

dim(c_IXS_to_IXV1)  <- N_age
dim(c_IXV1_to_IXV2) <- N_age

dim(c_RWS_to_RWV1)  <- N_age
dim(c_RWV1_to_RWV2) <- N_age
dim(c_RAS_to_RAV1)  <- N_age
dim(c_RAV1_to_RAV2) <- N_age
dim(c_RDS_to_RDV1)  <- N_age
dim(c_RDV1_to_RDV2) <- N_age
dim(c_ROS_to_ROV1)  <- N_age
dim(c_ROV1_to_ROV2) <- N_age

dim(c_RXS_to_RXV1)  <- N_age
dim(c_RXV1_to_RXV2) <- N_age

dim(c_EARWS_to_EARWV1) <- N_age
dim(c_EDRWS_to_EDRWV1) <- N_age
dim(c_EDRAS_to_EDRAV1) <- N_age
dim(c_EORWS_to_EORWV1) <- N_age
dim(c_EORAS_to_EORAV1) <- N_age
dim(c_EORDS_to_EORDV1) <- N_age

dim(c_EARWV1_to_EARWV2) <- N_age
dim(c_EDRWV1_to_EDRWV2) <- N_age
dim(c_EDRAV1_to_EDRAV2) <- N_age
dim(c_EORWV1_to_EORWV2) <- N_age
dim(c_EORAV1_to_EORAV2) <- N_age
dim(c_EORDV1_to_EORDV2) <- N_age

dim(c_IARWS_to_IARWV1) <- N_age
dim(c_IDRWS_to_IDRWV1) <- N_age
dim(c_IDRAS_to_IDRAV1) <- N_age
dim(c_IORWS_to_IORWV1) <- N_age
dim(c_IORAS_to_IORAV1) <- N_age
dim(c_IORDS_to_IORDV1) <- N_age

dim(c_IARWV1_to_IARWV2) <- N_age
dim(c_IDRWV1_to_IDRWV2) <- N_age
dim(c_IDRAV1_to_IDRAV2) <- N_age
dim(c_IORWV1_to_IORWV2) <- N_age
dim(c_IORAV1_to_IORAV2) <- N_age
dim(c_IORDV1_to_IORDV2) <- N_age

dim(c_EORXS_to_EORXV1)  <- N_age
dim(c_EORXV1_to_EORXV2) <- N_age

dim(c_EXROS_to_EXROV1)  <- N_age
dim(c_EXROV1_to_EXROV2) <- N_age

dim(c_IORXS_to_IORXV1)  <- N_age
dim(c_IORXV1_to_IORXV2) <- N_age

dim(c_IXROS_to_IXROV1)  <- N_age
dim(c_IXROV1_to_IXROV2) <- N_age

dim(c_ROXS_to_ROXV1)  <- N_age
dim(c_ROXV1_to_ROXV2) <- N_age

# ----- 0.4.5.7. Vaccination Waning -----

dim(c_V1_to_S) <- N_age
dim(c_V2_to_S) <- N_age

dim(c_EWV1_to_EWS) <- N_age
dim(c_EWV2_to_EWS) <- N_age
dim(c_EAV1_to_EAS) <- N_age
dim(c_EAV2_to_EAS) <- N_age
dim(c_EDV1_to_EDS) <- N_age
dim(c_EDV2_to_EDS) <- N_age
dim(c_EOV1_to_EOS) <- N_age
dim(c_EOV2_to_EOS) <- N_age

dim(c_EXV1_to_EXS) <- N_age
dim(c_EXV2_to_EXS) <- N_age

dim(c_IWV1_to_IWS) <- N_age
dim(c_IWV2_to_IWS) <- N_age
dim(c_IAV1_to_IAS) <- N_age
dim(c_IAV2_to_IAS) <- N_age
dim(c_IDV1_to_IDS) <- N_age
dim(c_IDV2_to_IDS) <- N_age
dim(c_IOV1_to_IOS) <- N_age
dim(c_IOV2_to_IOS) <- N_age

dim(c_IXV1_to_IXS) <- N_age
dim(c_IXV2_to_IXS) <- N_age

dim(c_RWV1_to_RWS) <- N_age
dim(c_RWV2_to_RWS) <- N_age
dim(c_RAV1_to_RAS) <- N_age
dim(c_RAV2_to_RAS) <- N_age
dim(c_RDV1_to_RDS) <- N_age
dim(c_RDV2_to_RDS) <- N_age
dim(c_ROV1_to_ROS) <- N_age
dim(c_ROV2_to_ROS) <- N_age

dim(c_RXV1_to_RXS) <- N_age
dim(c_RXV2_to_RXS) <- N_age

dim(c_EARWV1_to_EARWS) <- N_age
dim(c_EDRWV1_to_EDRWS) <- N_age
dim(c_EDRAV1_to_EDRAS) <- N_age
dim(c_EORWV1_to_EORWS) <- N_age
dim(c_EORAV1_to_EORAS) <- N_age
dim(c_EORDV1_to_EORDS) <- N_age

dim(c_EARWV2_to_EARWS) <- N_age
dim(c_EDRWV2_to_EDRWS) <- N_age
dim(c_EDRAV2_to_EDRAS) <- N_age
dim(c_EORWV2_to_EORWS) <- N_age
dim(c_EORAV2_to_EORAS) <- N_age
dim(c_EORDV2_to_EORDS) <- N_age

dim(c_IARWV1_to_IARWS) <- N_age
dim(c_IDRWV1_to_IDRWS) <- N_age
dim(c_IDRAV1_to_IDRAS) <- N_age
dim(c_IORWV1_to_IORWS) <- N_age
dim(c_IORAV1_to_IORAS) <- N_age
dim(c_IORDV1_to_IORDS) <- N_age

dim(c_IARWV2_to_IARWS) <- N_age
dim(c_IDRWV2_to_IDRWS) <- N_age
dim(c_IDRAV2_to_IDRAS) <- N_age
dim(c_IORWV2_to_IORWS) <- N_age
dim(c_IORAV2_to_IORAS) <- N_age
dim(c_IORDV2_to_IORDS) <- N_age

dim(c_EORXV1_to_EORXS) <- N_age
dim(c_EORXV2_to_EORXS) <- N_age

dim(c_EXROV1_to_EXROS) <- N_age
dim(c_EXROV2_to_EXROS) <- N_age

dim(c_IORXV1_to_IORXS) <- N_age
dim(c_IORXV2_to_IORXS) <- N_age

dim(c_IXROV1_to_IXROS) <- N_age
dim(c_IXROV2_to_IXROS) <- N_age

dim(c_ROXV1_to_ROXS) <- N_age
dim(c_ROXV2_to_ROXS) <- N_age

# ----- 0.4.5.8. Births,Migrations -----

dim(c_birth) <- N_age

dim(c_migration_S)  <- N_age
dim(c_migration_V1) <- N_age
dim(c_migration_V2) <- N_age

dim(c_migration_EWS)  <- N_age
dim(c_migration_IWS)  <- N_age
dim(c_migration_RWS)  <- N_age
dim(c_migration_EWV1) <- N_age
dim(c_migration_IWV1) <- N_age
dim(c_migration_RWV1) <- N_age
dim(c_migration_EWV2) <- N_age
dim(c_migration_IWV2) <- N_age
dim(c_migration_RWV2) <- N_age

dim(c_migration_EAS)     <- N_age
dim(c_migration_IAS)     <- N_age
dim(c_migration_RAS)     <- N_age
dim(c_migration_EARWS)  <- N_age
dim(c_migration_IARWS)  <- N_age
dim(c_migration_EAV1)    <- N_age
dim(c_migration_IAV1)    <- N_age
dim(c_migration_RAV1)    <- N_age
dim(c_migration_EARWV1) <- N_age
dim(c_migration_IARWV1) <- N_age
dim(c_migration_EAV2)    <- N_age
dim(c_migration_IAV2)    <- N_age
dim(c_migration_RAV2)    <- N_age
dim(c_migration_EARWV2) <- N_age
dim(c_migration_IARWV2) <- N_age

dim(c_migration_EDS)     <- N_age
dim(c_migration_IDS)     <- N_age
dim(c_migration_RDS)     <- N_age
dim(c_migration_EDRWS)  <- N_age
dim(c_migration_IDRWS)  <- N_age
dim(c_migration_EDRAS)  <- N_age
dim(c_migration_IDRAS)  <- N_age
dim(c_migration_EDV1)    <- N_age
dim(c_migration_IDV1)    <- N_age
dim(c_migration_RDV1)    <- N_age
dim(c_migration_EDRWV1) <- N_age
dim(c_migration_IDRWV1) <- N_age
dim(c_migration_EDRAV1) <- N_age
dim(c_migration_IDRAV1) <- N_age
dim(c_migration_EDV2)    <- N_age
dim(c_migration_IDV2)    <- N_age
dim(c_migration_RDV2)    <- N_age
dim(c_migration_EDRWV2) <- N_age
dim(c_migration_IDRWV2) <- N_age
dim(c_migration_EDRAV2) <- N_age
dim(c_migration_IDRAV2) <- N_age

dim(c_migration_EOS)     <- N_age
dim(c_migration_IOS)     <- N_age
dim(c_migration_ROS)     <- N_age
dim(c_migration_EORWS)  <- N_age
dim(c_migration_IORWS)  <- N_age
dim(c_migration_EORAS)  <- N_age
dim(c_migration_IORAS)  <- N_age
dim(c_migration_EORDS)  <- N_age
dim(c_migration_IORDS)  <- N_age
dim(c_migration_EOV1)    <- N_age
dim(c_migration_IOV1)    <- N_age
dim(c_migration_ROV1)    <- N_age
dim(c_migration_EORWV1) <- N_age
dim(c_migration_IORWV1) <- N_age
dim(c_migration_EORAV1) <- N_age
dim(c_migration_IORAV1) <- N_age
dim(c_migration_EORDV1) <- N_age
dim(c_migration_IORDV1) <- N_age
dim(c_migration_EOV2)    <- N_age
dim(c_migration_IOV2)    <- N_age
dim(c_migration_ROV2)    <- N_age
dim(c_migration_EORWV2) <- N_age
dim(c_migration_IORWV2) <- N_age
dim(c_migration_EORAV2) <- N_age
dim(c_migration_IORAV2) <- N_age
dim(c_migration_EORDV2) <- N_age
dim(c_migration_IORDV2) <- N_age

dim(c_migration_EXS)  <- N_age
dim(c_migration_IXS)  <- N_age
dim(c_migration_RXS)  <- N_age
dim(c_migration_EXV1) <- N_age
dim(c_migration_IXV1) <- N_age
dim(c_migration_RXV1) <- N_age
dim(c_migration_EXV2) <- N_age
dim(c_migration_IXV2) <- N_age
dim(c_migration_RXV2) <- N_age

dim(c_migration_EORXS)  <- N_age
dim(c_migration_IORXS)  <- N_age
dim(c_migration_EORXV1) <- N_age
dim(c_migration_IORXV1) <- N_age
dim(c_migration_EORXV2) <- N_age
dim(c_migration_IORXV2) <- N_age

dim(c_migration_EXROS)  <- N_age
dim(c_migration_IXROS)  <- N_age
dim(c_migration_EXROV1) <- N_age
dim(c_migration_IXROV1) <- N_age
dim(c_migration_EXROV2) <- N_age
dim(c_migration_IXROV2) <- N_age

dim(c_migration_ROXS)   <- N_age
dim(c_migration_ROXV1)  <- N_age
dim(c_migration_ROXV2)  <- N_age

# ----- 0.4.5.9. Deaths -----

dim(c_death_S)  <- N_age
dim(c_death_V1) <- N_age
dim(c_death_V2) <- N_age

dim(c_death_EWS)  <- N_age
dim(c_death_IWS)  <- N_age
dim(c_death_RWS)  <- N_age
dim(c_death_EWV1) <- N_age
dim(c_death_IWV1) <- N_age
dim(c_death_RWV1) <- N_age
dim(c_death_EWV2) <- N_age
dim(c_death_IWV2) <- N_age
dim(c_death_RWV2) <- N_age

dim(c_death_EAS)     <- N_age
dim(c_death_IAS)     <- N_age
dim(c_death_RAS)     <- N_age
dim(c_death_EARWS)  <- N_age
dim(c_death_IARWS)  <- N_age
dim(c_death_EAV1)    <- N_age
dim(c_death_IAV1)    <- N_age
dim(c_death_RAV1)    <- N_age
dim(c_death_EARWV1) <- N_age
dim(c_death_IARWV1) <- N_age
dim(c_death_EAV2)    <- N_age
dim(c_death_IAV2)    <- N_age
dim(c_death_RAV2)    <- N_age
dim(c_death_EARWV2) <- N_age
dim(c_death_IARWV2) <- N_age

dim(c_death_EDS)     <- N_age
dim(c_death_IDS)     <- N_age
dim(c_death_RDS)     <- N_age
dim(c_death_EDRWS)  <- N_age
dim(c_death_IDRWS)  <- N_age
dim(c_death_EDRAS)  <- N_age
dim(c_death_IDRAS)  <- N_age
dim(c_death_EDV1)    <- N_age
dim(c_death_IDV1)    <- N_age
dim(c_death_RDV1)    <- N_age
dim(c_death_EDRWV1) <- N_age
dim(c_death_IDRWV1) <- N_age
dim(c_death_EDRAV1) <- N_age
dim(c_death_IDRAV1) <- N_age
dim(c_death_EDV2)    <- N_age
dim(c_death_IDV2)    <- N_age
dim(c_death_RDV2)    <- N_age
dim(c_death_EDRWV2) <- N_age
dim(c_death_IDRWV2) <- N_age
dim(c_death_EDRAV2) <- N_age
dim(c_death_IDRAV2) <- N_age

dim(c_death_EOS)     <- N_age
dim(c_death_IOS)     <- N_age
dim(c_death_ROS)     <- N_age
dim(c_death_EORWS)  <- N_age
dim(c_death_IORWS)  <- N_age
dim(c_death_EORAS)  <- N_age
dim(c_death_IORAS)  <- N_age
dim(c_death_EORDS)  <- N_age
dim(c_death_IORDS)  <- N_age
dim(c_death_EOV1)    <- N_age
dim(c_death_IOV1)    <- N_age
dim(c_death_ROV1)    <- N_age
dim(c_death_EORWV1) <- N_age
dim(c_death_IORWV1) <- N_age
dim(c_death_EORAV1) <- N_age
dim(c_death_IORAV1) <- N_age
dim(c_death_EORDV1) <- N_age
dim(c_death_IORDV1) <- N_age
dim(c_death_EOV2)    <- N_age
dim(c_death_IOV2)    <- N_age
dim(c_death_ROV2)    <- N_age
dim(c_death_EORWV2) <- N_age
dim(c_death_IORWV2) <- N_age
dim(c_death_EORAV2) <- N_age
dim(c_death_IORAV2) <- N_age
dim(c_death_EORDV2) <- N_age
dim(c_death_IORDV2) <- N_age

dim(c_death_EXS)  <- N_age
dim(c_death_IXS)  <- N_age
dim(c_death_RXS)  <- N_age
dim(c_death_EXV1) <- N_age
dim(c_death_IXV1) <- N_age
dim(c_death_RXV1) <- N_age
dim(c_death_EXV2) <- N_age
dim(c_death_IXV2) <- N_age
dim(c_death_RXV2) <- N_age

dim(c_death_EORXS)  <- N_age
dim(c_death_IORXS)  <- N_age
dim(c_death_EORXV1) <- N_age
dim(c_death_IORXV1) <- N_age
dim(c_death_EORXV2) <- N_age
dim(c_death_IORXV2) <- N_age

dim(c_death_EXROS)  <- N_age
dim(c_death_IXROS)  <- N_age
dim(c_death_EXROV1) <- N_age
dim(c_death_IXROV1) <- N_age
dim(c_death_EXROV2) <- N_age
dim(c_death_IXROV2) <- N_age

dim(c_death_ROXS)   <- N_age
dim(c_death_ROXV1)  <- N_age
dim(c_death_ROXV2)  <- N_age

dim(c_death_inf_IWS) <- N_age
dim(c_death_inf_IAS) <- N_age
dim(c_death_inf_IDS) <- N_age
dim(c_death_inf_IOS) <- N_age
dim(c_death_inf_IXS)  <- N_age

dim(c_death_inf_IWV1) <- N_age
dim(c_death_inf_IAV1) <- N_age
dim(c_death_inf_IDV1) <- N_age
dim(c_death_inf_IOV1) <- N_age
dim(c_death_inf_IXV1)  <- N_age

dim(c_death_inf_IWV2) <- N_age
dim(c_death_inf_IAV2) <- N_age
dim(c_death_inf_IDV2) <- N_age
dim(c_death_inf_IOV2) <- N_age
dim(c_death_inf_IXV2)  <- N_age

dim(c_death_inf_IARWS) <- N_age
dim(c_death_inf_IDRWS) <- N_age
dim(c_death_inf_IDRAS) <- N_age
dim(c_death_inf_IORWS) <- N_age
dim(c_death_inf_IORAS) <- N_age
dim(c_death_inf_IORDS) <- N_age

dim(c_death_inf_IARWV1) <- N_age
dim(c_death_inf_IDRWV1) <- N_age
dim(c_death_inf_IDRAV1) <- N_age
dim(c_death_inf_IORWV1) <- N_age
dim(c_death_inf_IORAV1) <- N_age
dim(c_death_inf_IORDV1) <- N_age

dim(c_death_inf_IARWV2) <- N_age
dim(c_death_inf_IDRWV2) <- N_age
dim(c_death_inf_IDRAV2) <- N_age
dim(c_death_inf_IORWV2) <- N_age
dim(c_death_inf_IORAV2) <- N_age
dim(c_death_inf_IORDV2) <- N_age

dim(c_death_inf_IXROS)  <- N_age
dim(c_death_inf_IXROV1) <- N_age
dim(c_death_inf_IXROV2) <- N_age

dim(c_death_inf_IORXS)  <- N_age
dim(c_death_inf_IORXV1) <- N_age
dim(c_death_inf_IORXV2) <- N_age

# ----- 0.4.5.10. Aging -----

dim(c_aging_S)  <- N_age
dim(c_aging_V1) <- N_age
dim(c_aging_V2) <- N_age

dim(c_aging_EWS)  <- N_age
dim(c_aging_IWS)  <- N_age
dim(c_aging_RWS)  <- N_age
dim(c_aging_EWV1) <- N_age
dim(c_aging_IWV1) <- N_age
dim(c_aging_RWV1) <- N_age
dim(c_aging_EWV2) <- N_age
dim(c_aging_IWV2) <- N_age
dim(c_aging_RWV2) <- N_age

dim(c_aging_EAS)     <- N_age
dim(c_aging_IAS)     <- N_age
dim(c_aging_RAS)     <- N_age
dim(c_aging_EARWS)  <- N_age
dim(c_aging_IARWS)  <- N_age
dim(c_aging_EAV1)    <- N_age
dim(c_aging_IAV1)    <- N_age
dim(c_aging_RAV1)    <- N_age
dim(c_aging_EARWV1) <- N_age
dim(c_aging_IARWV1) <- N_age
dim(c_aging_EAV2)    <- N_age
dim(c_aging_IAV2)    <- N_age
dim(c_aging_RAV2)    <- N_age
dim(c_aging_EARWV2) <- N_age
dim(c_aging_IARWV2) <- N_age

dim(c_aging_EDS)     <- N_age
dim(c_aging_IDS)     <- N_age
dim(c_aging_RDS)     <- N_age
dim(c_aging_EDRWS)  <- N_age
dim(c_aging_IDRWS)  <- N_age
dim(c_aging_EDRAS)  <- N_age
dim(c_aging_IDRAS)  <- N_age
dim(c_aging_EDV1)    <- N_age
dim(c_aging_IDV1)    <- N_age
dim(c_aging_RDV1)    <- N_age
dim(c_aging_EDRWV1) <- N_age
dim(c_aging_IDRWV1) <- N_age
dim(c_aging_EDRAV1) <- N_age
dim(c_aging_IDRAV1) <- N_age
dim(c_aging_EDV2)    <- N_age
dim(c_aging_IDV2)    <- N_age
dim(c_aging_RDV2)    <- N_age
dim(c_aging_EDRWV2) <- N_age
dim(c_aging_IDRWV2) <- N_age
dim(c_aging_EDRAV2) <- N_age
dim(c_aging_IDRAV2) <- N_age

dim(c_aging_EOS)     <- N_age
dim(c_aging_IOS)     <- N_age
dim(c_aging_ROS)     <- N_age
dim(c_aging_EORWS)  <- N_age
dim(c_aging_IORWS)  <- N_age
dim(c_aging_EORAS)  <- N_age
dim(c_aging_IORAS)  <- N_age
dim(c_aging_EORDS)  <- N_age
dim(c_aging_IORDS)  <- N_age
dim(c_aging_EOV1)    <- N_age
dim(c_aging_IOV1)    <- N_age
dim(c_aging_ROV1)    <- N_age
dim(c_aging_EORWV1) <- N_age
dim(c_aging_IORWV1) <- N_age
dim(c_aging_EORAV1) <- N_age
dim(c_aging_IORAV1) <- N_age
dim(c_aging_EORDV1) <- N_age
dim(c_aging_IORDV1) <- N_age
dim(c_aging_EOV2)    <- N_age
dim(c_aging_IOV2)    <- N_age
dim(c_aging_ROV2)    <- N_age
dim(c_aging_EORWV2) <- N_age
dim(c_aging_IORWV2) <- N_age
dim(c_aging_EORAV2) <- N_age
dim(c_aging_IORAV2) <- N_age
dim(c_aging_EORDV2) <- N_age
dim(c_aging_IORDV2) <- N_age

dim(c_aging_EXS)  <- N_age
dim(c_aging_IXS)  <- N_age
dim(c_aging_RXS)  <- N_age
dim(c_aging_EXV1) <- N_age
dim(c_aging_IXV1) <- N_age
dim(c_aging_RXV1) <- N_age
dim(c_aging_EXV2) <- N_age
dim(c_aging_IXV2) <- N_age
dim(c_aging_RXV2) <- N_age

dim(c_aging_EORXS)  <- N_age
dim(c_aging_IORXS)  <- N_age
dim(c_aging_EORXV1) <- N_age
dim(c_aging_IORXV1) <- N_age
dim(c_aging_EORXV2) <- N_age
dim(c_aging_IORXV2) <- N_age

dim(c_aging_EXROS)  <- N_age
dim(c_aging_IXROS)  <- N_age
dim(c_aging_EXROV1) <- N_age
dim(c_aging_IXROV1) <- N_age
dim(c_aging_EXROV2) <- N_age
dim(c_aging_IXROV2) <- N_age

dim(c_aging_ROXS)   <- N_age
dim(c_aging_ROXV1)  <- N_age
dim(c_aging_ROXV2)  <- N_age

# ----- 0.4.6. Initial Exposures

dim(c_S_exp_X)  <- N_age
dim(c_S_exp_W) <- N_age

# ----- 0.4.6. Temporary Compartments -----

dim(Sa) <- N_age
dim(Sb) <- N_age
dim(Sc) <- N_age
dim(Sd) <- N_age
dim(Se) <- N_age
dim(Sf) <- N_age
dim(Sg) <- N_age

dim(V1a) <- N_age
dim(V1b) <- N_age
dim(V1c) <- N_age
dim(V1d) <- N_age
dim(V1e) <- N_age
dim(V1f) <- N_age
dim(V1g) <- N_age
dim(V1h) <- N_age


dim(V2a) <- N_age
dim(V2b) <- N_age
dim(V2c) <- N_age
dim(V2d) <- N_age
dim(V2e) <- N_age
dim(V2f) <- N_age
dim(V2g) <- N_age

dim(EWSa) <- N_age
dim(EWSb) <- N_age
dim(EWSc) <- N_age

dim(IWSa) <- N_age
dim(IWSb) <- N_age
dim(IWSc) <- N_age
dim(IWSd) <- N_age

dim(RWSa) <- N_age
dim(RWSb) <- N_age
dim(RWSc) <- N_age
dim(RWSd) <- N_age
dim(RWSe) <- N_age
dim(RWSf) <- N_age

dim(EWV1a) <- N_age
dim(EWV1b) <- N_age
dim(EWV1c) <- N_age
dim(EWV1d) <- N_age

dim(IWV1a) <- N_age
dim(IWV1b) <- N_age
dim(IWV1c) <- N_age
dim(IWV1d) <- N_age
dim(IWV1e) <- N_age

dim(RWV1a) <- N_age
dim(RWV1b) <- N_age
dim(RWV1c) <- N_age
dim(RWV1d) <- N_age
dim(RWV1e) <- N_age
dim(RWV1f) <- N_age
dim(RWV1g) <- N_age

dim(EWV2a) <- N_age
dim(EWV2b) <- N_age
dim(EWV2c) <- N_age

dim(IWV2a) <- N_age
dim(IWV2b) <- N_age
dim(IWV2c) <- N_age
dim(IWV2d) <- N_age

dim(RWV2a) <- N_age
dim(RWV2b) <- N_age
dim(RWV2c) <- N_age
dim(RWV2d) <- N_age
dim(RWV2e) <- N_age
dim(RWV2f) <- N_age

dim(EASa) <- N_age
dim(EASb) <- N_age
dim(EASc) <- N_age

dim(IASa) <- N_age
dim(IASb) <- N_age
dim(IASc) <- N_age
dim(IASd) <- N_age

dim(RASa) <- N_age
dim(RASb) <- N_age
dim(RASc) <- N_age
dim(RASd) <- N_age
dim(RASe) <- N_age

dim(EARWSa) <- N_age
dim(EARWSb) <- N_age
dim(EARWSc) <- N_age

dim(IARWSa) <- N_age
dim(IARWSb) <- N_age
dim(IARWSc) <- N_age
dim(IARWSd) <- N_age

dim(EAV1a) <- N_age
dim(EAV1b) <- N_age
dim(EAV1c) <- N_age
dim(EAV1d) <- N_age

dim(IAV1a) <- N_age
dim(IAV1b) <- N_age
dim(IAV1c) <- N_age
dim(IAV1d) <- N_age
dim(IAV1e) <- N_age

dim(RAV1a) <- N_age
dim(RAV1b) <- N_age
dim(RAV1c) <- N_age
dim(RAV1d) <- N_age
dim(RAV1e) <- N_age
dim(RAV1f) <- N_age

dim(EARWV1a) <- N_age
dim(EARWV1b) <- N_age
dim(EARWV1c) <- N_age
dim(EARWV1d) <- N_age

dim(IARWV1a) <- N_age
dim(IARWV1b) <- N_age
dim(IARWV1c) <- N_age
dim(IARWV1d) <- N_age
dim(IARWV1e) <- N_age

dim(EAV2a) <- N_age
dim(EAV2b) <- N_age
dim(EAV2c) <- N_age

dim(IAV2a) <- N_age
dim(IAV2b) <- N_age
dim(IAV2c) <- N_age
dim(IAV2d) <- N_age

dim(RAV2a) <- N_age
dim(RAV2b) <- N_age
dim(RAV2c) <- N_age
dim(RAV2d) <- N_age
dim(RAV2e) <- N_age

dim(EARWV2a) <- N_age
dim(EARWV2b) <- N_age
dim(EARWV2c) <- N_age

dim(IARWV2a) <- N_age
dim(IARWV2b) <- N_age
dim(IARWV2c) <- N_age
dim(IARWV2d) <- N_age

dim(EDSa) <- N_age
dim(EDSb) <- N_age
dim(EDSc) <- N_age

dim(IDSa) <- N_age
dim(IDSb) <- N_age
dim(IDSc) <- N_age
dim(IDSd) <- N_age

dim(RDSa) <- N_age
dim(RDSb) <- N_age
dim(RDSc) <- N_age
dim(RDSd) <- N_age

dim(EDRWSa) <- N_age
dim(EDRWSb) <- N_age
dim(EDRWSc) <- N_age

dim(IDRWSa) <- N_age
dim(IDRWSb) <- N_age
dim(IDRWSc) <- N_age
dim(IDRWSd) <- N_age

dim(EDRASa) <- N_age
dim(EDRASb) <- N_age
dim(EDRASc) <- N_age

dim(IDRASa) <- N_age
dim(IDRASb) <- N_age
dim(IDRASc) <- N_age
dim(IDRASd) <- N_age

dim(EDV1a) <- N_age
dim(EDV1b) <- N_age
dim(EDV1c) <- N_age
dim(EDV1d) <- N_age

dim(IDV1a) <- N_age
dim(IDV1b) <- N_age
dim(IDV1c) <- N_age
dim(IDV1d) <- N_age
dim(IDV1e) <- N_age

dim(RDV1a) <- N_age
dim(RDV1b) <- N_age
dim(RDV1c) <- N_age
dim(RDV1d) <- N_age
dim(RDV1e) <- N_age

dim(EDRWV1a) <- N_age
dim(EDRWV1b) <- N_age
dim(EDRWV1c) <- N_age
dim(EDRWV1d) <- N_age

dim(IDRWV1a) <- N_age
dim(IDRWV1b) <- N_age
dim(IDRWV1c) <- N_age
dim(IDRWV1d) <- N_age
dim(IDRWV1e) <- N_age

dim(EDRAV1a) <- N_age
dim(EDRAV1b) <- N_age
dim(EDRAV1c) <- N_age
dim(EDRAV1d) <- N_age

dim(IDRAV1a) <- N_age
dim(IDRAV1b) <- N_age
dim(IDRAV1c) <- N_age
dim(IDRAV1d) <- N_age
dim(IDRAV1e) <- N_age

dim(EDV2a) <- N_age
dim(EDV2b) <- N_age
dim(EDV2c) <- N_age

dim(IDV2a) <- N_age
dim(IDV2b) <- N_age
dim(IDV2c) <- N_age
dim(IDV2d) <- N_age

dim(RDV2a) <- N_age
dim(RDV2b) <- N_age
dim(RDV2c) <- N_age
dim(RDV2d) <- N_age

dim(EDRWV2a) <- N_age
dim(EDRWV2b) <- N_age
dim(EDRWV2c) <- N_age

dim(IDRWV2a) <- N_age
dim(IDRWV2b) <- N_age
dim(IDRWV2c) <- N_age
dim(IDRWV2d) <- N_age

dim(EDRAV2a) <- N_age
dim(EDRAV2b) <- N_age
dim(EDRAV2c) <- N_age

dim(IDRAV2a) <- N_age
dim(IDRAV2b) <- N_age
dim(IDRAV2c) <- N_age
dim(IDRAV2d) <- N_age

dim(EOSa) <- N_age
dim(EOSb) <- N_age
dim(EOSc) <- N_age

dim(IOSa) <- N_age
dim(IOSb) <- N_age
dim(IOSc) <- N_age
dim(IOSd) <- N_age

dim(ROSa) <- N_age
dim(ROSb) <- N_age
dim(ROSc) <- N_age
dim(ROSd) <- N_age

dim(EORWSa) <- N_age
dim(EORWSb) <- N_age
dim(EORWSc) <- N_age

dim(IORWSa) <- N_age
dim(IORWSb) <- N_age
dim(IORWSc) <- N_age
dim(IORWSd) <- N_age

dim(EORASa) <- N_age
dim(EORASb) <- N_age
dim(EORASc) <- N_age

dim(IORASa) <- N_age
dim(IORASb) <- N_age
dim(IORASc) <- N_age
dim(IORASd) <- N_age

dim(EORDSa) <- N_age
dim(EORDSb) <- N_age
dim(EORDSc) <- N_age

dim(IORDSa) <- N_age
dim(IORDSb) <- N_age
dim(IORDSc) <- N_age
dim(IORDSd) <- N_age

dim(EOV1a) <- N_age
dim(EOV1b) <- N_age
dim(EOV1c) <- N_age
dim(EOV1d) <- N_age

dim(IOV1a) <- N_age
dim(IOV1b) <- N_age
dim(IOV1c) <- N_age
dim(IOV1d) <- N_age
dim(IOV1e) <- N_age

dim(ROV1a) <- N_age
dim(ROV1b) <- N_age
dim(ROV1c) <- N_age
dim(ROV1d) <- N_age
dim(ROV1e) <- N_age

dim(EORWV1a) <- N_age
dim(EORWV1b) <- N_age
dim(EORWV1c) <- N_age
dim(EORWV1d) <- N_age

dim(IORWV1a) <- N_age
dim(IORWV1b) <- N_age
dim(IORWV1c) <- N_age
dim(IORWV1d) <- N_age
dim(IORWV1e) <- N_age

dim(EORAV1a) <- N_age
dim(EORAV1b) <- N_age
dim(EORAV1c) <- N_age
dim(EORAV1d) <- N_age

dim(IORAV1a) <- N_age
dim(IORAV1b) <- N_age
dim(IORAV1c) <- N_age
dim(IORAV1d) <- N_age
dim(IORAV1e) <- N_age

dim(EORDV1a) <- N_age
dim(EORDV1b) <- N_age
dim(EORDV1c) <- N_age
dim(EORDV1d) <- N_age

dim(IORDV1a) <- N_age
dim(IORDV1b) <- N_age
dim(IORDV1c) <- N_age
dim(IORDV1d) <- N_age
dim(IORDV1e) <- N_age

dim(EOV2a) <- N_age
dim(EOV2b) <- N_age
dim(EOV2c) <- N_age

dim(IOV2a) <- N_age
dim(IOV2b) <- N_age
dim(IOV2c) <- N_age
dim(IOV2d) <- N_age

dim(ROV2a) <- N_age
dim(ROV2b) <- N_age
dim(ROV2c) <- N_age
dim(ROV2d) <- N_age

dim(EORWV2a) <- N_age
dim(EORWV2b) <- N_age
dim(EORWV2c) <- N_age

dim(IORWV2a) <- N_age
dim(IORWV2b) <- N_age
dim(IORWV2c) <- N_age
dim(IORWV2d) <- N_age

dim(EORAV2a) <- N_age
dim(EORAV2b) <- N_age
dim(EORAV2c) <- N_age

dim(IORAV2a) <- N_age
dim(IORAV2b) <- N_age
dim(IORAV2c) <- N_age
dim(IORAV2d) <- N_age

dim(EORDV2a) <- N_age
dim(EORDV2b) <- N_age
dim(EORDV2c) <- N_age

dim(IORDV2a) <- N_age
dim(IORDV2b) <- N_age
dim(IORDV2c) <- N_age
dim(IORDV2d) <- N_age

dim(EXSa) <- N_age
dim(EXSb) <- N_age
dim(EXSc) <- N_age

dim(IXSa) <- N_age
dim(IXSb) <- N_age
dim(IXSc) <- N_age
dim(IXSd) <- N_age

dim(RXSa) <- N_age
dim(RXSb) <- N_age
dim(RXSc) <- N_age
dim(RXSd) <- N_age

dim(EXV1a) <- N_age
dim(EXV1b) <- N_age
dim(EXV1c) <- N_age
dim(EXV1d) <- N_age

dim(IXV1a) <- N_age
dim(IXV1b) <- N_age
dim(IXV1c) <- N_age
dim(IXV1d) <- N_age
dim(IXV1e) <- N_age

dim(RXV1a) <- N_age
dim(RXV1b) <- N_age
dim(RXV1c) <- N_age
dim(RXV1d) <- N_age
dim(RXV1e) <- N_age

dim(EXV2a) <- N_age
dim(EXV2b) <- N_age
dim(EXV2c) <- N_age

dim(IXV2a) <- N_age
dim(IXV2b) <- N_age
dim(IXV2c) <- N_age
dim(IXV2d) <- N_age

dim(RXV2a) <- N_age
dim(RXV2b) <- N_age
dim(RXV2c) <- N_age
dim(RXV2d) <- N_age

dim(EORXSa) <- N_age
dim(EORXSb) <- N_age
dim(EORXSc) <- N_age

dim(IORXSa) <- N_age
dim(IORXSb) <- N_age
dim(IORXSc) <- N_age
dim(IORXSd) <- N_age

dim(ROXSa) <- N_age
dim(ROXSb) <- N_age
dim(ROXSc) <- N_age

dim(EORXV1a) <- N_age
dim(EORXV1b) <- N_age
dim(EORXV1c) <- N_age
dim(EORXV1d) <- N_age

dim(IORXV1a) <- N_age
dim(IORXV1b) <- N_age
dim(IORXV1c) <- N_age
dim(IORXV1d) <- N_age
dim(IORXV1e) <- N_age

dim(ROXV1a) <- N_age
dim(ROXV1b) <- N_age
dim(ROXV1c) <- N_age
dim(ROXV1d) <- N_age

dim(EORXV2a) <- N_age
dim(EORXV2b) <- N_age
dim(EORXV2c) <- N_age

dim(IORXV2a) <- N_age
dim(IORXV2b) <- N_age
dim(IORXV2c) <- N_age
dim(IORXV2d) <- N_age

dim(ROXV2a) <- N_age
dim(ROXV2b) <- N_age
dim(ROXV2c) <- N_age

dim(EXROSa) <- N_age
dim(EXROSb) <- N_age
dim(EXROSc) <- N_age

dim(IXROSa) <- N_age
dim(IXROSb) <- N_age
dim(IXROSc) <- N_age
dim(IXROSd) <- N_age

dim(EXROV1a) <- N_age
dim(EXROV1b) <- N_age
dim(EXROV1c) <- N_age
dim(EXROV1d) <- N_age

dim(IXROV1a) <- N_age
dim(IXROV1b) <- N_age
dim(IXROV1c) <- N_age
dim(IXROV1d) <- N_age
dim(IXROV1e) <- N_age

dim(EXROV2a) <- N_age
dim(EXROV2b) <- N_age
dim(EXROV2c) <- N_age

dim(IXROV2a) <- N_age
dim(IXROV2b) <- N_age
dim(IXROV2c) <- N_age
dim(IXROV2d) <- N_age

# ----- 0.4.6. In/Out -----

dim(c_S_in)         <- N_age
dim(c_S_out)        <- N_age
dim(c_V1_in)        <- N_age
dim(c_V1_out)       <- N_age
dim(c_V2_in)        <- N_age
dim(c_V2_out)       <- N_age
dim(c_EWS_in)      <- N_age
dim(c_EWS_out)     <- N_age
dim(c_IWS_in)      <- N_age
dim(c_IWS_out)     <- N_age
dim(c_RWS_in)      <- N_age
dim(c_RWS_out)     <- N_age
dim(c_EWV1_in)     <- N_age
dim(c_EWV1_out)    <- N_age
dim(c_IWV1_in)     <- N_age
dim(c_IWV1_out)    <- N_age
dim(c_RWV1_in)     <- N_age
dim(c_RWV1_out)    <- N_age
dim(c_EWV2_in)     <- N_age
dim(c_EWV2_out)    <- N_age
dim(c_IWV2_in)     <- N_age
dim(c_IWV2_out)    <- N_age
dim(c_RWV2_in)     <- N_age
dim(c_RWV2_out)    <- N_age
dim(c_EAS_in)      <- N_age
dim(c_EAS_out)     <- N_age
dim(c_IAS_in)      <- N_age
dim(c_IAS_out)     <- N_age
dim(c_RAS_in)      <- N_age
dim(c_RAS_out)     <- N_age
dim(c_EARWS_in)   <- N_age
dim(c_EARWS_out)  <- N_age
dim(c_IARWS_in)   <- N_age
dim(c_IARWS_out)  <- N_age
dim(c_EAV1_in)     <- N_age
dim(c_EAV1_out)    <- N_age
dim(c_IAV1_in)     <- N_age
dim(c_IAV1_out)    <- N_age
dim(c_RAV1_in)     <- N_age
dim(c_RAV1_out)    <- N_age
dim(c_EARWV1_in)  <- N_age
dim(c_EARWV1_out) <- N_age
dim(c_IARWV1_in)  <- N_age
dim(c_IARWV1_out) <- N_age
dim(c_EAV2_in)     <- N_age
dim(c_EAV2_out)    <- N_age
dim(c_IAV2_in)     <- N_age
dim(c_IAV2_out)    <- N_age
dim(c_RAV2_in)     <- N_age
dim(c_RAV2_out)    <- N_age
dim(c_EARWV2_in)  <- N_age
dim(c_EARWV2_out) <- N_age
dim(c_IARWV2_in)  <- N_age
dim(c_IARWV2_out) <- N_age
dim(c_EDS_in)      <- N_age
dim(c_EDS_out)     <- N_age
dim(c_IDS_in)      <- N_age
dim(c_IDS_out)     <- N_age
dim(c_RDS_in)      <- N_age
dim(c_RDS_out)     <- N_age
dim(c_EDRWS_in)   <- N_age
dim(c_EDRWS_out)  <- N_age
dim(c_IDRWS_in)   <- N_age
dim(c_IDRWS_out)  <- N_age
dim(c_EDRAS_in)   <- N_age
dim(c_EDRAS_out)  <- N_age
dim(c_IDRAS_in)   <- N_age
dim(c_IDRAS_out)  <- N_age
dim(c_EDV1_in)     <- N_age
dim(c_EDV1_out)    <- N_age
dim(c_IDV1_in)     <- N_age
dim(c_IDV1_out)    <- N_age
dim(c_RDV1_in)     <- N_age
dim(c_RDV1_out)    <- N_age
dim(c_EDRWV1_in)  <- N_age
dim(c_EDRWV1_out) <- N_age
dim(c_IDRWV1_in)  <- N_age
dim(c_IDRWV1_out) <- N_age
dim(c_EDRAV1_in)  <- N_age
dim(c_EDRAV1_out) <- N_age
dim(c_IDRAV1_in)  <- N_age
dim(c_IDRAV1_out) <- N_age
dim(c_EDV2_in)     <- N_age
dim(c_EDV2_out)    <- N_age
dim(c_IDV2_in)     <- N_age
dim(c_IDV2_out)    <- N_age
dim(c_RDV2_in)     <- N_age
dim(c_RDV2_out)    <- N_age
dim(c_EDRWV2_in)  <- N_age
dim(c_EDRWV2_out) <- N_age
dim(c_IDRWV2_in)  <- N_age
dim(c_IDRWV2_out) <- N_age
dim(c_EDRAV2_in)  <- N_age
dim(c_EDRAV2_out) <- N_age
dim(c_IDRAV2_in)  <- N_age
dim(c_IDRAV2_out) <- N_age
dim(c_EOS_in)      <- N_age
dim(c_EOS_out)     <- N_age
dim(c_IOS_in)      <- N_age
dim(c_IOS_out)     <- N_age
dim(c_ROS_in)      <- N_age
dim(c_ROS_out)     <- N_age
dim(c_EORWS_in)   <- N_age
dim(c_EORWS_out)  <- N_age
dim(c_IORWS_in)   <- N_age
dim(c_IORWS_out)  <- N_age
dim(c_EORAS_in)   <- N_age
dim(c_EORAS_out)  <- N_age
dim(c_IORAS_in)   <- N_age
dim(c_IORAS_out)  <- N_age
dim(c_EORDS_in)   <- N_age
dim(c_EORDS_out)  <- N_age
dim(c_IORDS_in)   <- N_age
dim(c_IORDS_out)  <- N_age
dim(c_EOV1_in)     <- N_age
dim(c_EOV1_out)    <- N_age
dim(c_IOV1_in)     <- N_age
dim(c_IOV1_out)    <- N_age
dim(c_ROV1_in)     <- N_age
dim(c_ROV1_out)    <- N_age
dim(c_EORWV1_in)  <- N_age
dim(c_EORWV1_out) <- N_age
dim(c_IORWV1_in)  <- N_age
dim(c_IORWV1_out) <- N_age
dim(c_EORAV1_in)  <- N_age
dim(c_EORAV1_out) <- N_age
dim(c_IORAV1_in)  <- N_age
dim(c_IORAV1_out) <- N_age
dim(c_EORDV1_in)  <- N_age
dim(c_EORDV1_out) <- N_age
dim(c_IORDV1_in)  <- N_age
dim(c_IORDV1_out) <- N_age
dim(c_EOV2_in)     <- N_age
dim(c_EOV2_out)    <- N_age
dim(c_IOV2_in)     <- N_age
dim(c_IOV2_out)    <- N_age
dim(c_ROV2_in)     <- N_age
dim(c_ROV2_out)    <- N_age
dim(c_EORWV2_in)  <- N_age
dim(c_EORWV2_out) <- N_age
dim(c_IORWV2_in)  <- N_age
dim(c_IORWV2_out) <- N_age
dim(c_EORAV2_in)  <- N_age
dim(c_EORAV2_out) <- N_age
dim(c_IORAV2_in)  <- N_age
dim(c_IORAV2_out) <- N_age
dim(c_EORDV2_in)  <- N_age
dim(c_EORDV2_out) <- N_age
dim(c_IORDV2_in)  <- N_age
dim(c_IORDV2_out) <- N_age
dim(c_EXS_in)       <- N_age
dim(c_EXS_out)      <- N_age
dim(c_IXS_in)       <- N_age
dim(c_IXS_out)      <- N_age
dim(c_RXS_in)       <- N_age
dim(c_RXS_out)      <- N_age
dim(c_EXV1_in)      <- N_age
dim(c_EXV1_out)     <- N_age
dim(c_IXV1_in)      <- N_age
dim(c_IXV1_out)     <- N_age
dim(c_RXV1_in)      <- N_age
dim(c_RXV1_out)     <- N_age
dim(c_EXV2_in)      <- N_age
dim(c_EXV2_out)     <- N_age
dim(c_IXV2_in)      <- N_age
dim(c_IXV2_out)     <- N_age
dim(c_RXV2_in)      <- N_age
dim(c_RXV2_out)     <- N_age
dim(c_EORXS_in)    <- N_age
dim(c_EORXS_out)   <- N_age
dim(c_IORXS_in)    <- N_age
dim(c_IORXS_out)   <- N_age
dim(c_EORXV1_in)   <- N_age
dim(c_EORXV1_out)  <- N_age
dim(c_IORXV1_in)   <- N_age
dim(c_IORXV1_out)  <- N_age
dim(c_EORXV2_in)   <- N_age
dim(c_EORXV2_out)  <- N_age
dim(c_IORXV2_in)   <- N_age
dim(c_IORXV2_out)  <- N_age
dim(c_EXROS_in)    <- N_age
dim(c_EXROS_out)   <- N_age
dim(c_IXROS_in)    <- N_age
dim(c_IXROS_out)   <- N_age
dim(c_EXROV1_in)   <- N_age
dim(c_EXROV1_out)  <- N_age
dim(c_IXROV1_in)   <- N_age
dim(c_IXROV1_out)  <- N_age
dim(c_EXROV2_in)   <- N_age
dim(c_EXROV2_out)  <- N_age
dim(c_IXROV2_in)   <- N_age
dim(c_IXROV2_out)  <- N_age
dim(c_ROXS_in)     <- N_age
dim(c_ROXS_out)    <- N_age
dim(c_ROXV1_in)    <- N_age
dim(c_ROXV1_out)   <- N_age
dim(c_ROXV2_in)    <- N_age
dim(c_ROXV2_out)   <- N_age

dim(I2_tot) <- N_age
dim(IX_tot) <- N_age

# ------------------------------------------------------------------------------
# ----- 1. Transitions ---------------------------------------------------------
# ------------------------------------------------------------------------------
# ----- 1.1. Infection Contacts ------------------------------------------------
# ----- 1.1.1 Lockdowns -----

I2_tot[] <-
  IWS[i] + IWV1[i] + IWV2[i] +
  IAS[i] + IAV1[i] + IAV2[i] +
  IDS[i] + IDV1[i] + IDV2[i] +
  IOS[i] + IOV1[i] + IOV2[i] +
  IORXS[i] +  IORXV1[i] +  IORXV2[i] +
  IARWS[i] +  IARWV1[i] +  IARWV2[i] +
  IDRWS[i] +  IDRWV1[i] +  IDRWV2[i] +
  IDRAS[i] +  IDRAV1[i] +  IDRAV2[i] +
  IORWS[i] +  IORWV1[i] +  IORWV2[i] +
  IORAS[i] +  IORAV1[i] +  IORAV2[i] +
  IORDS[i] +  IORDV1[i] +  IORDV2[i]

IX_tot[] <-
  IXS[i] + IXV1[i] + IXV2[i] +
  IXROS[i] + IXROV1[i] + IXROV2[i]

contacts_2[,] <- contact_matrices[i,j,step] * I2_tot[i]
contacts_X[,] <- contact_matrices[i,j,step] * IX_tot[i]

# ----- 1.2. Mean Rates --------------------------------------------------------
# ----- 1.2.1. Naive to Exposed -----
u_S_exposed_W[] <- u_trans_W[i] * prevalence_W[step] * sum(contacts_2[,i])
u_S_exposed_A[] <- u_trans_A[i] * prevalence_A[step] * sum(contacts_2[,i])
u_S_exposed_D[] <- u_trans_D[i] * prevalence_D[step] * sum(contacts_2[,i])
u_S_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i])
u_S_exposed_X[]  <- u_trans_X[i] * sum(contacts_X[,i])

u_V1_exposed_W[] <- u_trans_W[i] * prevalence_W[step] * sum(contacts_2[,i]) * (1-immunity_inf_W[i,1])
u_V1_exposed_A[] <- u_trans_A[i] * prevalence_A[step] * sum(contacts_2[,i]) * (1-immunity_inf_A[i,1])
u_V1_exposed_D[] <- u_trans_D[i] * prevalence_D[step] * sum(contacts_2[,i]) * (1-immunity_inf_D[i,1])
u_V1_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,1])
u_V1_exposed_X[]  <- u_trans_X[i] * sum(contacts_X[,i]) * (1-immunity_inf_X[i,1])

u_V2_exposed_W[] <- u_trans_W[i] * prevalence_W[step] * sum(contacts_2[,i]) * (1-immunity_inf_W[i,2])
u_V2_exposed_A[] <- u_trans_A[i] * prevalence_A[step] * sum(contacts_2[,i]) * (1-immunity_inf_A[i,2])
u_V2_exposed_D[] <- u_trans_D[i] * prevalence_D[step] * sum(contacts_2[,i]) * (1-immunity_inf_D[i,2])
u_V2_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,2])
u_V2_exposed_X[]  <- u_trans_X[i] * sum(contacts_X[,i]) * (1-immunity_inf_X[i,2])

# ----- 1.2.2. Recovered to Exposed -----
u_RWS_exposed_A[] <- u_trans_A[i] * prevalence_A[step] * sum(contacts_2[,i]) * (1-immunity_inf_A[i,4])
u_RWS_exposed_D[] <- u_trans_D[i] * prevalence_D[step] * sum(contacts_2[,i]) * (1-immunity_inf_D[i,4])
u_RAS_exposed_D[] <- u_trans_D[i] * prevalence_D[step] * sum(contacts_2[,i]) * (1-immunity_inf_D[i,5])
u_RWS_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,4])
u_RAS_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,5])
u_RDS_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,6])

u_RWV1_exposed_A[] <- u_trans_A[i] * prevalence_A[step] * sum(contacts_2[,i]) * (1-immunity_inf_A[i,8])
u_RWV1_exposed_D[] <- u_trans_D[i] * prevalence_D[step] * sum(contacts_2[,i]) * (1-immunity_inf_D[i,8])
u_RAV1_exposed_D[] <- u_trans_D[i] * prevalence_D[step] * sum(contacts_2[,i]) * (1-immunity_inf_D[i,9])
u_RWV1_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,8])
u_RAV1_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,9])
u_RDV1_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,10])

u_RWV2_exposed_A[] <- u_trans_A[i] * prevalence_A[step] * sum(contacts_2[,i]) * (1-immunity_inf_A[i,12])
u_RWV2_exposed_D[] <- u_trans_D[i] * prevalence_D[step] * sum(contacts_2[,i]) * (1-immunity_inf_D[i,12])
u_RAV2_exposed_D[] <- u_trans_D[i] * prevalence_D[step] * sum(contacts_2[,i]) * (1-immunity_inf_D[i,13])
u_RWV2_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,12])
u_RAV2_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,13])
u_RDV2_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_O[i,14])

u_RXS_exposed_O[]  <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_X[i,7])
u_RXV1_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_X[i,11])
u_RXV2_exposed_O[] <- u_trans_O[i] * prevalence_O[step] * sum(contacts_2[,i]) * (1-immunity_inf_X[i,15])

u_ROS_exposed_X[] <- u_trans_X[i] * sum(contacts_X[,i]) * (1-immunity_inf_X[i,7])
u_ROV1_exposed_X[] <- u_trans_X[i] * sum(contacts_X[,i]) * (1-immunity_inf_X[i,11])
u_ROV2_exposed_X[] <- u_trans_X[i] * sum(contacts_X[,i]) * (1-immunity_inf_X[i,15])

# ----- 1.2.3. Infection Deaths -----

u_death_W_S[] <- u_death_W[i]
u_death_A_S[] <- u_death_A[i]
u_death_D_S[] <- u_death_D[i]
u_death_O_S[] <- u_death_O[i]
u_death_X_S[] <- u_death_X[i]

u_death_W_V1[] <- u_death_W[i] * (1-immunity_vir_W[i,1])
u_death_A_V1[] <- u_death_A[i] * (1-immunity_vir_A[i,1])
u_death_D_V1[] <- u_death_D[i] * (1-immunity_vir_D[i,1])
u_death_O_V1[] <- u_death_O[i] * (1-immunity_vir_O[i,1])
u_death_X_V1[] <- u_death_X[i] * (1-immunity_vir_X[i,1])

u_death_W_V2[] <- u_death_W[i] * (1-immunity_vir_W[i,2])
u_death_A_V2[] <- u_death_A[i] * (1-immunity_vir_A[i,2])
u_death_D_V2[] <- u_death_D[i] * (1-immunity_vir_D[i,2])
u_death_O_V2[] <- u_death_O[i] * (1-immunity_vir_O[i,2])
u_death_X_V2[] <- u_death_X[i] * (1-immunity_vir_X[i,2])

u_death_A_RWS[] <- u_death_A[i] * (1-immunity_vir_A[i,4])
u_death_D_RWS[] <- u_death_D[i] * (1-immunity_vir_D[i,4])
u_death_D_RAS[] <- u_death_D[i] * (1-immunity_vir_D[i,5])
u_death_O_RWS[] <- u_death_O[i] * (1-immunity_vir_O[i,4])
u_death_O_RAS[] <- u_death_O[i] * (1-immunity_vir_O[i,5])
u_death_O_RDS[] <- u_death_O[i] * (1-immunity_vir_O[i,6])

u_death_A_RWV1[] <- u_death_A[i] * (1-immunity_vir_A[i,8])
u_death_D_RWV1[] <- u_death_D[i] * (1-immunity_vir_D[i,8])
u_death_D_RAV1[] <- u_death_D[i] * (1-immunity_vir_D[i,9])
u_death_O_RWV1[] <- u_death_O[i] * (1-immunity_vir_O[i,8])
u_death_O_RAV1[] <- u_death_O[i] * (1-immunity_vir_O[i,9])
u_death_O_RDV1[] <- u_death_O[i] * (1-immunity_vir_O[i,10])

u_death_A_RWV2[] <- u_death_A[i] * (1-immunity_vir_A[i,12])
u_death_D_RWV2[] <- u_death_D[i] * (1-immunity_vir_D[i,12])
u_death_D_RAV2[] <- u_death_D[i] * (1-immunity_vir_D[i,13])
u_death_O_RWV2[] <- u_death_O[i] * (1-immunity_vir_O[i,12])
u_death_O_RAV2[] <- u_death_O[i] * (1-immunity_vir_O[i,13])
u_death_O_RDV2[] <- u_death_O[i] * (1-immunity_vir_O[i,14])

u_death_X_ROS[] <- u_death_X[i] * (1-immunity_vir_X[i,7])
u_death_X_ROV1[] <- u_death_X[i] * (1-immunity_vir_X[i,11])
u_death_X_ROV2[] <- u_death_X[i] * (1-immunity_vir_X[i,15])

u_death_O_RXS[] <- u_death_O[i] * (1-immunity_vir_X[i,7])
u_death_O_RXV1[] <- u_death_O[i] * (1-immunity_vir_X[i,11])
u_death_O_RXV2[] <- u_death_O[i] * (1-immunity_vir_X[i,15])

# ----- 1.3. Transition Probabilities ------------------------------------------
# ----- 1.3.1. Births/Migrations -----

p_birth[] <- 1 - exp(-u_birth[i])

p_migration_S[] <- 1 - exp(-(S[i]/N[i]))
p_migration_V1[] <- 1 - exp(-(V1[i]/N[i]))
p_migration_V2[] <- 1 - exp(-(V2[i]/N[i]))

p_migration_EWS[] <- 1 - exp(-(EWS[i]/N[i]))
p_migration_IWS[] <- 1 - exp(-(IWS[i]/N[i]))
p_migration_RWS[] <- 1 - exp(-(RWS[i]/N[i]))
p_migration_EWV1[] <- 1 - exp(-(EWV1[i]/N[i]))
p_migration_IWV1[] <- 1 - exp(-(IWV1[i]/N[i]))
p_migration_RWV1[] <- 1 - exp(-(RWV1[i]/N[i]))
p_migration_EWV2[] <- 1 - exp(-(EWV2[i]/N[i]))
p_migration_IWV2[] <- 1 - exp(-(IWV2[i]/N[i]))
p_migration_RWV2[] <- 1 - exp(-(RWV2[i]/N[i]))

p_migration_EAS[] <- 1 - exp(-(EAS[i]/N[i]))
p_migration_IAS[] <- 1 - exp(-(IAS[i]/N[i]))
p_migration_RAS[] <- 1 - exp(-(RAS[i]/N[i]))
p_migration_EARWS[] <- 1 - exp(-(EARWS[i]/N[i]))
p_migration_IARWS[] <- 1 - exp(-(IARWS[i]/N[i]))
p_migration_EAV1[] <- 1 - exp(-(EAV1[i]/N[i]))
p_migration_IAV1[] <- 1 - exp(-(IAV1[i]/N[i]))
p_migration_RAV1[] <- 1 - exp(-(RAV1[i]/N[i]))
p_migration_EARWV1[] <- 1 - exp(-(EARWV1[i]/N[i]))
p_migration_IARWV1[] <- 1 - exp(-(IARWV1[i]/N[i]))
p_migration_EAV2[] <- 1 - exp(-(EAV2[i]/N[i]))
p_migration_IAV2[] <- 1 - exp(-(IAV2[i]/N[i]))
p_migration_RAV2[] <- 1 - exp(-(RAV2[i]/N[i]))
p_migration_EARWV2[] <- 1 - exp(-(EARWV2[i]/N[i]))
p_migration_IARWV2[] <- 1 - exp(-(IARWV2[i]/N[i]))

p_migration_EDS[] <- 1 - exp(-(EDS[i]/N[i]))
p_migration_IDS[] <- 1 - exp(-(IDS[i]/N[i]))
p_migration_RDS[] <- 1 - exp(-(RDS[i]/N[i]))
p_migration_EDRWS[] <- 1 - exp(-(EDRWS[i]/N[i]))
p_migration_IDRWS[] <- 1 - exp(-(IDRWS[i]/N[i]))
p_migration_EDRAS[] <- 1 - exp(-(EDRAS[i]/N[i]))
p_migration_IDRAS[] <- 1 - exp(-(IDRAS[i]/N[i]))
p_migration_EDV1[] <- 1 - exp(-(EDV1[i]/N[i]))
p_migration_IDV1[] <- 1 - exp(-(IDV1[i]/N[i]))
p_migration_RDV1[] <- 1 - exp(-(RDV1[i]/N[i]))
p_migration_EDRWV1[] <- 1 - exp(-(EDRWV1[i]/N[i]))
p_migration_IDRWV1[] <- 1 - exp(-(IDRWV1[i]/N[i]))
p_migration_EDRAV1[] <- 1 - exp(-(EDRAV1[i]/N[i]))
p_migration_IDRAV1[] <- 1 - exp(-(IDRAV1[i]/N[i]))
p_migration_EDV2[] <- 1 - exp(-(EDV2[i]/N[i]))
p_migration_IDV2[] <- 1 - exp(-(IDV2[i]/N[i]))
p_migration_RDV2[] <- 1 - exp(-(RDV2[i]/N[i]))
p_migration_EDRWV2[] <- 1 - exp(-(EDRWV2[i]/N[i]))
p_migration_IDRWV2[] <- 1 - exp(-(IDRWV2[i]/N[i]))
p_migration_EDRAV2[] <- 1 - exp(-(EDRAV2[i]/N[i]))
p_migration_IDRAV2[] <- 1 - exp(-(IDRAV2[i]/N[i]))

p_migration_EOS[] <- 1 - exp(-(EOS[i]/N[i]))
p_migration_IOS[] <- 1 - exp(-(IOS[i]/N[i]))
p_migration_ROS[] <- 1 - exp(-(ROS[i]/N[i]))
p_migration_EORWS[] <- 1 - exp(-(EORWS[i]/N[i]))
p_migration_IORWS[] <- 1 - exp(-(IORWS[i]/N[i]))
p_migration_EORAS[] <- 1 - exp(-(EORAS[i]/N[i]))
p_migration_IORAS[] <- 1 - exp(-(IORAS[i]/N[i]))
p_migration_EORDS[] <- 1 - exp(-(EORDS[i]/N[i]))
p_migration_IORDS[] <- 1 - exp(-(IORDS[i]/N[i]))
p_migration_EOV1[] <- 1 - exp(-(EOV1[i]/N[i]))
p_migration_IOV1[] <- 1 - exp(-(IOV1[i]/N[i]))
p_migration_ROV1[] <- 1 - exp(-(ROV1[i]/N[i]))
p_migration_EORWV1[] <- 1 - exp(-(EORWV1[i]/N[i]))
p_migration_IORWV1[] <- 1 - exp(-(IORWV1[i]/N[i]))
p_migration_EORAV1[] <- 1 - exp(-(EORAV1[i]/N[i]))
p_migration_IORAV1[] <- 1 - exp(-(IORAV1[i]/N[i]))
p_migration_EORDV1[] <- 1 - exp(-(EORDV1[i]/N[i]))
p_migration_IORDV1[] <- 1 - exp(-(IORDV1[i]/N[i]))
p_migration_EOV2[] <- 1 - exp(-(EOV2[i]/N[i]))
p_migration_IOV2[] <- 1 - exp(-(IOV2[i]/N[i]))
p_migration_ROV2[] <- 1 - exp(-(ROV2[i]/N[i]))
p_migration_EORWV2[] <- 1 - exp(-(EORWV2[i]/N[i]))
p_migration_IORWV2[] <- 1 - exp(-(IORWV2[i]/N[i]))
p_migration_EORAV2[] <- 1 - exp(-(EORAV2[i]/N[i]))
p_migration_IORAV2[] <- 1 - exp(-(IORAV2[i]/N[i]))
p_migration_EORDV2[] <- 1 - exp(-(EORDV2[i]/N[i]))
p_migration_IORDV2[] <- 1 - exp(-(IORDV2[i]/N[i]))

p_migration_EXS[] <- 1 - exp(-(EXS[i]/N[i]))
p_migration_IXS[] <- 1 - exp(-(IXS[i]/N[i]))
p_migration_RXS[] <- 1 - exp(-(RXS[i]/N[i]))
p_migration_EXV1[] <- 1 - exp(-(EXV1[i]/N[i]))
p_migration_IXV1[] <- 1 - exp(-(IXV1[i]/N[i]))
p_migration_RXV1[] <- 1 - exp(-(RXV1[i]/N[i]))
p_migration_EXV2[] <- 1 - exp(-(EXV2[i]/N[i]))
p_migration_IXV2[] <- 1 - exp(-(IXV2[i]/N[i]))
p_migration_RXV2[] <- 1 - exp(-(RXV2[i]/N[i]))

p_migration_EORXS[] <- 1 - exp(-(EORXS[i]/N[i]))
p_migration_IORXS[] <- 1 - exp(-(IORXS[i]/N[i]))
p_migration_EORXV1[] <- 1 - exp(-(EORXV1[i]/N[i]))
p_migration_IORXV1[] <- 1 - exp(-(IORXV1[i]/N[i]))
p_migration_EORXV2[] <- 1 - exp(-(EORXV2[i]/N[i]))
p_migration_IORXV2[] <- 1 - exp(-(IORXV2[i]/N[i]))

p_migration_EXROS[] <- 1 - exp(-(EXROS[i]/N[i]))
p_migration_IXROS[] <- 1 - exp(-(IXROS[i]/N[i]))
p_migration_EXROV1[] <- 1 - exp(-(EXROV1[i]/N[i]))
p_migration_IXROV1[] <- 1 - exp(-(IXROV1[i]/N[i]))
p_migration_EXROV2[] <- 1 - exp(-(EXROV2[i]/N[i]))
p_migration_IXROV2[] <- 1 - exp(-(IXROV2[i]/N[i]))

p_migration_ROXS[] <- 1 - exp(-(ROXS[i]/N[i]))
p_migration_ROXV1[] <- 1 - exp(-(ROXV1[i]/N[i]))
p_migration_ROXV2[] <- 1 - exp(-(ROXV2[i]/N[i]))

# ----- 1.3.2. Deaths -----

p_death[] <- 1 - exp(-u_death[i])

p_death_IWS[] <- 1 - exp(-u_death_W_S[i])
p_death_IAS[] <- 1 - exp(-u_death_A_S[i])
p_death_IDS[] <- 1 - exp(-u_death_D_S[i])
p_death_IOS[] <- 1 - exp(-u_death_O_S[i])
p_death_IXS[]  <- 1 - exp(-u_death_X_S[i])

p_death_IWV1[] <- 1 - exp(-u_death_W_V1[i])
p_death_IAV1[] <- 1 - exp(-u_death_A_V1[i])
p_death_IDV1[] <- 1 - exp(-u_death_D_V1[i])
p_death_IOV1[] <- 1 - exp(-u_death_O_V1[i])
p_death_IXV1[]  <- 1 - exp(-u_death_X_V1[i])

p_death_IWV2[] <- 1 - exp(-u_death_W_V2[i])
p_death_IAV2[] <- 1 - exp(-u_death_A_V2[i])
p_death_IDV2[] <- 1 - exp(-u_death_D_V2[i])
p_death_IOV2[] <- 1 - exp(-u_death_O_V2[i])
p_death_IXV2[]  <- 1 - exp(-u_death_X_V2[i])

p_death_IARWS[] <- 1 - exp(-u_death_A_RWS[i])
p_death_IDRWS[] <- 1 - exp(-u_death_D_RWS[i])
p_death_IDRAS[] <- 1 - exp(-u_death_D_RAS[i])
p_death_IORWS[] <- 1 - exp(-u_death_O_RWS[i])
p_death_IORAS[] <- 1 - exp(-u_death_O_RAS[i])
p_death_IORDS[] <- 1 - exp(-u_death_O_RDS[i])

p_death_IARWV1[] <- 1 - exp(-u_death_A_RWV1[i])
p_death_IDRWV1[] <- 1 - exp(-u_death_D_RWV1[i])
p_death_IDRAV1[] <- 1 - exp(-u_death_D_RAV1[i])
p_death_IORWV1[] <- 1 - exp(-u_death_O_RWV1[i])
p_death_IORAV1[] <- 1 - exp(-u_death_O_RAV1[i])
p_death_IORDV1[] <- 1 - exp(-u_death_O_RDV1[i])

p_death_IARWV2[] <- 1 - exp(-u_death_A_RWV2[i])
p_death_IDRWV2[] <- 1 - exp(-u_death_D_RWV2[i])
p_death_IDRAV2[] <- 1 - exp(-u_death_D_RAV2[i])
p_death_IORWV2[] <- 1 - exp(-u_death_O_RWV2[i])
p_death_IORAV2[] <- 1 - exp(-u_death_O_RAV2[i])
p_death_IORDV2[] <- 1 - exp(-u_death_O_RDV2[i])

p_death_IXROS[]  <- 1 - exp(-u_death_X_ROS[i])
p_death_IXROV1[] <- 1 - exp(-u_death_X_ROV1[i])
p_death_IXROV2[] <- 1 - exp(-u_death_X_ROV2[i])

p_death_IORXS[]  <- 1 - exp(-u_death_O_RXS[i])
p_death_IORXV1[] <- 1 - exp(-u_death_O_RXV1[i])
p_death_IORXV2[] <- 1 - exp(-u_death_O_RXV2[i])

# ----- 1.3.3. Naive to Exposed -----

p_S_to_EWS[] <- 1 - exp(-u_S_exposed_W[i])
p_S_to_EAS[] <- 1 - exp(-u_S_exposed_A[i])
p_S_to_EDS[] <- 1 - exp(-u_S_exposed_D[i])
p_S_to_EOS[] <- 1 - exp(-u_S_exposed_O[i])
p_S_to_EXS[] <- 1 - exp(-u_S_exposed_X[i])

p_V1_to_EWV1[] <- 1 - exp(-u_V1_exposed_W[i])
p_V1_to_EAV1[] <- 1 - exp(-u_V1_exposed_A[i])
p_V1_to_EDV1[] <- 1 - exp(-u_V1_exposed_D[i])
p_V1_to_EOV1[] <- 1 - exp(-u_V1_exposed_O[i])
p_V1_to_EXV1[] <- 1 - exp(-u_V1_exposed_X[i])

p_V2_to_EWV2[] <- 1 - exp(-u_V2_exposed_W[i])
p_V2_to_EAV2[] <- 1 - exp(-u_V2_exposed_A[i])
p_V2_to_EDV2[] <- 1 - exp(-u_V2_exposed_D[i])
p_V2_to_EOV2[] <- 1 - exp(-u_V2_exposed_O[i])
p_V2_to_EXV2[] <- 1 - exp(-u_V2_exposed_X[i])

# ----- 1.3.3. Recovered to Exposed -----

p_RWS_to_EARWS[] <- 1 - exp(-u_RWS_exposed_A[i])
p_RWS_to_EDRWS[] <- 1 - exp(-u_RWS_exposed_D[i])
p_RAS_to_EDRAS[] <- 1 - exp(-u_RAS_exposed_D[i])
p_RWS_to_EORWS[] <- 1 - exp(-u_RWS_exposed_O[i])
p_RAS_to_EORAS[] <- 1 - exp(-u_RAS_exposed_O[i])
p_RDS_to_EORDS[] <- 1 - exp(-u_RDS_exposed_O[i])

p_RWV1_to_EARWV1[] <- 1 - exp(-u_RWV1_exposed_A[i])
p_RWV1_to_EDRWV1[] <- 1 - exp(-u_RWV1_exposed_D[i])
p_RAV1_to_EDRAV1[] <- 1 - exp(-u_RAV1_exposed_D[i])
p_RWV1_to_EORWV1[] <- 1 - exp(-u_RWV1_exposed_O[i])
p_RAV1_to_EORAV1[] <- 1 - exp(-u_RAV1_exposed_O[i])
p_RDV1_to_EORDV1[] <- 1 - exp(-u_RDV1_exposed_O[i])

p_RWV2_to_EARWV2[] <- 1 - exp(-u_RWV2_exposed_A[i])
p_RWV2_to_EDRWV2[] <- 1 - exp(-u_RWV2_exposed_D[i])
p_RAV2_to_EDRAV2[] <- 1 - exp(-u_RAV2_exposed_D[i])
p_RWV2_to_EORWV2[] <- 1 - exp(-u_RWV2_exposed_O[i])
p_RAV2_to_EORAV2[] <- 1 - exp(-u_RAV2_exposed_O[i])
p_RDV2_to_EORDV2[] <- 1 - exp(-u_RDV2_exposed_O[i])

p_RXS_to_EORXS[] <- 1 - exp(-u_RXS_exposed_O[i])
p_RXV1_to_EORXV1[] <- 1 - exp(-u_RXV1_exposed_O[i])
p_RXV2_to_EORXV2[] <- 1 - exp(-u_RXV2_exposed_O[i])

p_ROS_to_EXROS[] <- 1 - exp(-u_ROS_exposed_X[i])
p_ROV1_to_EXROV1[] <- 1 - exp(-u_ROV1_exposed_X[i])
p_ROV2_to_EXROV2[] <- 1 - exp(-u_ROV2_exposed_X[i])

# ----- 1.3.4. Exposed to Infected -----

p_incub_W[] <- 1 - exp(-u_incub_W[i])
p_incub_A[] <- 1 - exp(-u_incub_A[i])
p_incub_D[] <- 1 - exp(-u_incub_D[i])
p_incub_O[] <- 1 - exp(-u_incub_O[i])
p_incub_X[] <- 1 - exp(-u_incub_X[i])

# ----- 1.3.5. Infected to Recovered -----

p_recov_W[] <- 1 - exp(-u_recov_W[i])
p_recov_A[] <- 1 - exp(-u_recov_A[i])
p_recov_D[] <- 1 - exp(-u_recov_D[i])
p_recov_O[] <- 1 - exp(-u_recov_O[i])
p_recov_X[] <- 1 - exp(-u_recov_X[i])

# ----- 1.3.6. Recovered to Naive -----

p_waning_W[] <- 1 - exp(-u_waning_W[i])
p_waning_A[] <- 1 - exp(-u_waning_A[i])
p_waning_D[] <- 1 - exp(-u_waning_D[i])
p_waning_O[] <- 1 - exp(-u_waning_O[i])
p_waning_X[] <- 1 - exp(-u_waning_X[i])

# ----- 1.3.7. Vaccination -----

p_S_to_V1[]  <- 1 - exp(-(u_V1[step,i]+u_V2[step,i]))
p_V1_to_V2[] <- 1 - exp(-(u_V2[step,i]))

# ----- 1.3.10. Vaccination Waning -----

p_V1_to_S[]  <- 1 - exp(-u_waning_V1[i])
p_V2_to_S[] <- 1 - exp(-u_waning_V2[i])

# ----- 1.4. Transitions -------------------------------------------------------
# ----- 1.4.1. S Transitions -----

# S -> EXS
c_S_to_EXS[] <- rbinom(S[i],p_S_to_EXS[i])
Sa[] <- S[i] - c_S_to_EXS[i]

# S -> EWS
c_S_to_EWS[] <- rbinom(Sa[i],p_S_to_EWS[i])
Sb[] <- Sa[i] - c_S_to_EWS[i]

# S -> EAS
c_S_to_EAS[] <- rbinom(Sb[i],p_S_to_EAS[i])
Sc[] <- Sb[i] - c_S_to_EAS[i]

# S -> EDS
c_S_to_EDS[] <- rbinom(Sc[i],p_S_to_EDS[i])
Sd[] <- Sc[i] - c_S_to_EDS[i]

# S -> EOS
c_S_to_EOS[] <- rbinom(Sd[i],p_S_to_EOS[i])
Se[] <- Sd[i] - c_S_to_EOS[i]

# S -> V1
c_S_to_V1[] <-  rbinom(Se[i],p_S_to_V1[i])
Sf[] <- Se[i] - c_S_to_V1[i]

# Deaths
c_death_S[] <- rbinom(Sf[i],p_death[i])
Sg[] <- Sf[i] - c_death_S[i]

# Aging
c_aging_S[] <- floor(Sg[i] * u_aging)

# ----- 1.4.2. V1 Transitions -----

# V1 -> EXV1
c_V1_to_EXV1[] <- rbinom(V1[i],p_V1_to_EXV1[i])
V1a[] <- V1[i] - c_V1_to_EXV1[i]
# V1 -> EWV1
c_V1_to_EWV1[] <- rbinom(V1a[i],p_V1_to_EWV1[i])
V1b[] <- V1a[i] - c_V1_to_EWV1[i]

# V1 -> EAV1
c_V1_to_EAV1[] <- rbinom(V1b[i],p_V1_to_EAV1[i])
V1c[] <- V1b[i] - c_V1_to_EAV1[i]

# V1 -> EDV1
c_V1_to_EDV1[] <- rbinom(V1c[i],p_V1_to_EDV1[i])
V1d[] <- V1c[i] - c_V1_to_EDV1[i]

# V1 -> EOV1
c_V1_to_EOV1[] <- rbinom(V1d[i],p_V1_to_EOV1[i])
V1e[] <- V1d[i] - c_V1_to_EOV1[i]

# V1 -> V2
c_V1_to_V2[] <- rbinom(V1e[i],p_V1_to_V2[i])
V1f[] <- V1e[i] - c_V1_to_V2[i]

# V1 -> S
c_V1_to_S[] <- rbinom(V1f[i],p_V1_to_S[i])
V1g[] <- V1f[i] - c_V1_to_S[i]

# Death
c_death_V1[] <- rbinom(V1g[i],p_death[i])
V1h[] <- V1g[i] - c_death_V1[i]

# Aging
c_aging_V1[] <- floor(V1h[i] * u_aging)

# ----- 1.4.2. V2 Transitions -----

# V2 -> EXV2
c_V2_to_EXV2[] <- rbinom(V2[i],p_V2_to_EXV2[i])
V2a[] <- V2[i] - c_V2_to_EXV2[i]
# V2 -> EWV2
c_V2_to_EWV2[] <- rbinom(V2a[i],p_V2_to_EWV2[i])
V2b[] <- V2a[i] - c_V2_to_EWV2[i]

# V2 -> EAV2
c_V2_to_EAV2[] <- rbinom(V2b[i],p_V2_to_EAV2[i])
V2c[] <- V2b[i] - c_V2_to_EAV2[i]

# V2 -> EDV2
c_V2_to_EDV2[] <- rbinom(V2c[i],p_V2_to_EDV2[i])
V2d[] <- V2c[i] - c_V2_to_EDV2[i]

# V2 -> EOV2
c_V2_to_EOV2[] <- rbinom(V2d[i],p_V2_to_EOV2[i])
V2e[] <- V2d[i] - c_V2_to_EOV2[i]

# V2 -> S
c_V2_to_S[] <- rbinom(V2e[i],p_V2_to_S[i])
V2f[] <- V2e[i] - c_V2_to_S[i]

# Death
c_death_V2[] <- rbinom(V2f[i],p_death[i])
V2g[] <- V2f[i] - c_death_V2[i]

# Aging
c_aging_V2[] <- floor(V2g[i] * u_aging)

# ----- 1.4.3. EWS Transitions -----

# EWS -> IWS
c_EWS_to_IWS[] <- rbinom(EWS[i],p_incub_W[i])
EWSa[] <- EWS[i] - c_EWS_to_IWS[i]
  
# EWS -> EWV1
c_EWS_to_EWV1[] <- rbinom(EWSa[i],p_S_to_V1[i])
EWSb[] <- EWSa[i] - c_EWS_to_EWV1[i]

# Death
c_death_EWS[] <- rbinom(EWSb[i],p_death[i])
EWSc[] <- EWSb[i] - c_death_EWS[i]

# Aging
c_aging_EWS[] <- floor(EWSc[i] * u_aging)

# ----- 1.4.3. EWV1 Transitions -----

# EWV1 -> IWV1
c_EWV1_to_IWV1[] <- rbinom(EWV1[i],p_incub_W[i])
EWV1a[] <- EWV1[i] - c_EWV1_to_IWV1[i]

# EWV1 -> EWV2
c_EWV1_to_EWV2[] <- rbinom(EWV1a[i],p_V1_to_V2[i])
EWV1b[] <- EWV1a[i] - c_EWV1_to_EWV2[i]

# EWV1 -> EWS
c_EWV1_to_EWS[] <- rbinom(EWV1b[i],p_V1_to_S[i])
EWV1c[] <- EWV1b[i] - c_EWV1_to_EWS[i]

# Death
c_death_EWV1[] <- rbinom(EWV1c[i],p_death[i])
EWV1d[] <- EWV1c[i] - c_death_EWV1[i]

# Aging
c_aging_EWV1[] <- floor(EWV1d[i] * u_aging)

# ----- 1.4.3. EWV2 Transitions -----

# EWV2 -> IWV2
c_EWV2_to_IWV2[] <- rbinom(EWV2[i],p_incub_W[i])
EWV2a[] <- EWV2[i] - c_EWV2_to_IWV2[i]

# EWV2 -> EWS
c_EWV2_to_EWS[] <- rbinom(EWV2a[i],p_V2_to_S[i])
EWV2b[] <- EWV2a[i] - c_EWV2_to_EWS[i]

# Death
c_death_EWV2[] <- rbinom(EWV2b[i],p_death[i])
EWV2c[] <- EWV2b[i] - c_death_EWV2[i]

# Aging
c_aging_EWV2[] <- floor(EWV2c[i] * u_aging)

# ----- 1.4.3. EAS Transitions -----

# EAS -> IAS
c_EAS_to_IAS[] <- rbinom(EAS[i],p_incub_A[i])
EASa[] <- EAS[i] - c_EAS_to_IAS[i]

# EAS -> EAV1
c_EAS_to_EAV1[] <- rbinom(EASa[i],p_S_to_V1[i])
EASb[] <- EASa[i] - c_EAS_to_EAV1[i]

# Death
c_death_EAS[] <- rbinom(EASb[i],p_death[i])
EASc[] <- EASb[i] - c_death_EAS[i]

# Aging
c_aging_EAS[] <- floor(EASc[i] * u_aging)

# ----- 1.4.3. EAV1 Transitions -----

# EAV1 -> IAV1
c_EAV1_to_IAV1[] <- rbinom(EAV1[i],p_incub_A[i])
EAV1a[] <- EAV1[i] - c_EAV1_to_IAV1[i]

# EAV1 -> EAV2
c_EAV1_to_EAV2[] <- rbinom(EAV1a[i],p_V1_to_V2[i])
EAV1b[] <- EAV1a[i] - c_EAV1_to_EAV2[i]

# EAV1 -> EAS
c_EAV1_to_EAS[] <- rbinom(EAV1b[i],p_V1_to_S[i])
EAV1c[] <- EAV1b[i] - c_EAV1_to_EAS[i]

# Death
c_death_EAV1[] <- rbinom(EAV1c[i],p_death[i])
EAV1d[] <- EAV1c[i] - c_death_EAV1[i]

# Aging
c_aging_EAV1[] <- floor(EAV1d[i] * u_aging)

# ----- 1.4.3. EAV2 Transitions -----

# EAV2 -> IAV2
c_EAV2_to_IAV2[] <- rbinom(EAV2[i],p_incub_A[i])
EAV2a[] <- EAV2[i] - c_EAV2_to_IAV2[i]

# EAV2 -> EAS
c_EAV2_to_EAS[] <- rbinom(EAV2a[i],p_V2_to_S[i])
EAV2b[] <- EAV2a[i] - c_EAV2_to_EAS[i]

# Death
c_death_EAV2[] <- rbinom(EAV2b[i],p_death[i])
EAV2c[] <- EAV2b[i] - c_death_EAV2[i]

# Aging
c_aging_EAV2[] <- floor(EAV2c[i] * u_aging)

# ----- 1.4.3. EDS Transitions -----

# EDS -> IDS
c_EDS_to_IDS[] <- rbinom(EDS[i],p_incub_D[i])
EDSa[] <- EDS[i] - c_EDS_to_IDS[i]

# EDS -> EDV1
c_EDS_to_EDV1[] <- rbinom(EDSa[i],p_S_to_V1[i])
EDSb[] <- EDSa[i] - c_EDS_to_EDV1[i]

# Death
c_death_EDS[] <- rbinom(EDSb[i],p_death[i])
EDSc[] <- EDSb[i] - c_death_EDS[i]

# Aging
c_aging_EDS[] <- floor(EDSc[i] * u_aging)

# ----- 1.4.3. EDV1 Transitions -----

# EDV1 -> IDV1
c_EDV1_to_IDV1[] <- rbinom(EDV1[i],p_incub_D[i])
EDV1a[] <- EDV1[i] - c_EDV1_to_IDV1[i]

# EDV1 -> EDV2
c_EDV1_to_EDV2[] <- rbinom(EDV1a[i],p_V1_to_V2[i])
EDV1b[] <- EDV1a[i] - c_EDV1_to_EDV2[i]

# EDV1 -> EDS
c_EDV1_to_EDS[] <- rbinom(EDV1b[i],p_V1_to_S[i])
EDV1c[] <- EDV1b[i] - c_EDV1_to_EDS[i]

# Death
c_death_EDV1[] <- rbinom(EDV1c[i],p_death[i])
EDV1d[] <- EDV1c[i] - c_death_EDV1[i]

# Aging
c_aging_EDV1[] <- floor(EDV1d[i] * u_aging)

# ----- 1.4.3. EDV2 Transitions -----

# EDV2 -> IDV2
c_EDV2_to_IDV2[] <- rbinom(EDV2[i],p_incub_D[i])
EDV2a[] <- EDV2[i] - c_EDV2_to_IDV2[i]

# EDV2 -> EDS
c_EDV2_to_EDS[] <- rbinom(EDV2a[i],p_V2_to_S[i])
EDV2b[] <- EDV2a[i] - c_EDV2_to_EDS[i]

# Death
c_death_EDV2[] <- rbinom(EDV2b[i],p_death[i])
EDV2c[] <- EDV2b[i] - c_death_EDV2[i]

# Aging
c_aging_EDV2[] <- floor(EDV2c[i] * u_aging)

# ----- 1.4.3. EOS Transitions -----

# EOS -> IOS
c_EOS_to_IOS[] <- rbinom(EOS[i],p_incub_O[i])
EOSa[] <- EOS[i] - c_EOS_to_IOS[i]

# EOS -> EOV1
c_EOS_to_EOV1[] <- rbinom(EOSa[i],p_S_to_V1[i])
EOSb[] <- EOSa[i] - c_EOS_to_EOV1[i]

# Death
c_death_EOS[] <- rbinom(EOSb[i],p_death[i])
EOSc[] <- EOSb[i] - c_death_EOS[i]

# Aging
c_aging_EOS[] <- floor(EOSc[i] * u_aging)

# ----- 1.4.3. EOV1 Transitions -----

# EOV1 -> IOV1
c_EOV1_to_IOV1[] <- rbinom(EOV1[i],p_incub_O[i])
EOV1a[] <- EOV1[i] - c_EOV1_to_IOV1[i]

# EOV1 -> EOV2
c_EOV1_to_EOV2[] <- rbinom(EOV1a[i],p_V1_to_V2[i])
EOV1b[] <- EOV1a[i] - c_EOV1_to_EOV2[i]

# EOV1 -> EOS
c_EOV1_to_EOS[] <- rbinom(EOV1b[i],p_V1_to_S[i])
EOV1c[] <- EOV1b[i] - c_EOV1_to_EOS[i]

# Death
c_death_EOV1[] <- rbinom(EOV1c[i],p_death[i])
EOV1d[] <- EOV1c[i] - c_death_EOV1[i]

# Aging
c_aging_EOV1[] <- floor(EOV1d[i] * u_aging)

# ----- 1.4.3. EOV2 Transitions -----

# EOV2 -> IOV2
c_EOV2_to_IOV2[] <- rbinom(EOV2[i],p_incub_O[i])
EOV2a[] <- EOV2[i] - c_EOV2_to_IOV2[i]

# EOV2 -> EOS
c_EOV2_to_EOS[] <- rbinom(EOV2a[i],p_V2_to_S[i])
EOV2b[] <- EOV2a[i] - c_EOV2_to_EOS[i]

# Death
c_death_EOV2[] <- rbinom(EOV2b[i],p_death[i])
EOV2c[] <- EOV2b[i] - c_death_EOV2[i]

# Aging
c_aging_EOV2[] <- floor(EOV2c[i] * u_aging)

# ----- 1.4.11. IWS Transitions -----

# IWS -> RWS
c_IWS_to_RWS[] <- rbinom(IWS[i],p_recov_W[i])
IWSa[] <- IWS[i] - c_IWS_to_RWS[i]

# Infection Deaths
c_death_inf_IWS[] <- rbinom(IWSa[i],p_death_IWS[i])
IWSb[] <- IWSa[i] - c_death_inf_IWS[i]

# IWS -> IWV1
c_IWS_to_IWV1[] <- rbinom(IWSb[i],p_S_to_V1[i])
IWSc[] <- IWSb[i] - c_IWS_to_IWV1[i]

# Death
c_death_IWS[] <- rbinom(IWSc[i],p_death[i])
IWSd[] <- IWSc[i] - c_death_IWS[i]

# Aging
c_aging_IWS[] <- floor(IWSd[i] * u_aging)

# ----- 1.4.11. IWV1 Transitions -----

# IWV1 -> RWV1
c_IWV1_to_RWV1[] <- rbinom(IWV1[i],p_recov_W[i])
IWV1a[] <- IWV1[i] - c_IWV1_to_RWV1[i]

# Infection Deaths
c_death_inf_IWV1[] <- rbinom(IWV1a[i],p_death_IWV1[i])
IWV1b[] <- IWV1a[i] - c_death_inf_IWV1[i]

# IWV1 -> IWV2
c_IWV1_to_IWV2[] <- rbinom(IWV1b[i],p_V1_to_V2[i])
IWV1c[] <- IWV1b[i] - c_IWV1_to_IWV2[i]

# IWV1 -> IWS
c_IWV1_to_IWS[] <- rbinom(IWV1c[i],p_V1_to_S[i])
IWV1d[] <- IWV1c[i] - c_IWV1_to_IWS[i]

# Death
c_death_IWV1[] <- rbinom(IWV1d[i],p_death[i])
IWV1e[] <- IWV1d[i] - c_death_IWV1[i]

# Aging
c_aging_IWV1[] <- floor(IWV1e[i] * u_aging)

# ----- 1.4.11. IWV2 Transitions -----

# IWV2 -> RWV2
c_IWV2_to_RWV2[] <- rbinom(IWV2[i],p_recov_W[i])
IWV2a[] <- IWV2[i] - c_IWV2_to_RWV2[i]

# Infection Deaths
c_death_inf_IWV2[] <- rbinom(IWV2a[i],p_death_IWV2[i])
IWV2b[] <- IWV2a[i] - c_death_inf_IWV2[i]

# IWV2 -> IWS
c_IWV2_to_IWS[] <- rbinom(IWV2b[i],p_V2_to_S[i])
IWV2c[] <- IWV2b[i] - c_IWV2_to_IWS[i]

# Death
c_death_IWV2[] <- rbinom(IWV2c[i],p_death[i])
IWV2d[] <- IWV2c[i] - c_death_IWV2[i]

# Aging
c_aging_IWV2[] <- floor(IWV2d[i] * u_aging)

# ----- 1.4.11. IAS Transitions -----

# IAS -> RAS
c_IAS_to_RAS[] <- rbinom(IAS[i],p_recov_A[i])
IASa[] <- IAS[i] - c_IAS_to_RAS[i]

# Infection Deaths
c_death_inf_IAS[] <- rbinom(IASa[i],p_death_IAS[i])
IASb[] <- IASa[i] - c_death_inf_IAS[i]

# IAS -> IAV1
c_IAS_to_IAV1[] <- rbinom(IASb[i],p_S_to_V1[i])
IASc[] <- IASb[i] - c_IAS_to_IAV1[i]

# Death
c_death_IAS[] <- rbinom(IASc[i],p_death[i])
IASd[] <- IASc[i] - c_death_IAS[i]

# Aging
c_aging_IAS[] <- floor(IASd[i] * u_aging)

# ----- 1.4.11. IAV1 Transitions -----

# IAV1 -> RAV1
c_IAV1_to_RAV1[] <- rbinom(IAV1[i],p_recov_A[i])
IAV1a[] <- IAV1[i] - c_IAV1_to_RAV1[i]

# Infection Deaths
c_death_inf_IAV1[] <- rbinom(IAV1a[i],p_death_IAV1[i])
IAV1b[] <- IAV1a[i] - c_death_inf_IAV1[i]

# IAV1 -> IAV2
c_IAV1_to_IAV2[] <- rbinom(IAV1b[i],p_V1_to_V2[i])
IAV1c[] <- IAV1b[i] - c_IAV1_to_IAV2[i]

# IAV1 -> IAS
c_IAV1_to_IAS[] <- rbinom(IAV1c[i],p_V1_to_S[i])
IAV1d[] <- IAV1c[i] - c_IAV1_to_IAS[i]

# Death
c_death_IAV1[] <- rbinom(IAV1d[i],p_death[i])
IAV1e[] <- IAV1d[i] - c_death_IAV1[i]

# Aging
c_aging_IAV1[] <- floor(IAV1e[i] * u_aging)

# ----- 1.4.11. IAV2 Transitions -----

# IAV2 -> RAV2
c_IAV2_to_RAV2[] <- rbinom(IAV2[i],p_recov_A[i])
IAV2a[] <- IAV2[i] - c_IAV2_to_RAV2[i]

# Infection Deaths
c_death_inf_IAV2[] <- rbinom(IAV2a[i],p_death_IAV2[i])
IAV2b[] <- IAV2a[i] - c_death_inf_IAV2[i]

# IAV2 -> IAS
c_IAV2_to_IAS[] <- rbinom(IAV2b[i],p_V2_to_S[i])
IAV2c[] <- IAV2b[i] - c_IAV2_to_IAS[i]

# Death
c_death_IAV2[] <- rbinom(IAV2c[i],p_death[i])
IAV2d[] <- IAV2c[i] - c_death_IAV2[i]

# Aging
c_aging_IAV2[] <- floor(IAV2d[i] * u_aging)

# ----- 1.4.11. IDS Transitions -----

# IDS -> RDS
c_IDS_to_RDS[] <- rbinom(IDS[i],p_recov_D[i])
IDSa[] <- IDS[i] - c_IDS_to_RDS[i]

# Infection Deaths
c_death_inf_IDS[] <- rbinom(IDSa[i],p_death_IDS[i])
IDSb[] <- IDSa[i] - c_death_inf_IDS[i]

# IDS -> IDV1
c_IDS_to_IDV1[] <- rbinom(IDSb[i],p_S_to_V1[i])
IDSc[] <- IDSb[i] - c_IDS_to_IDV1[i]

# Death
c_death_IDS[] <- rbinom(IDSc[i],p_death[i])
IDSd[] <- IDSc[i] - c_death_IDS[i]

# Aging
c_aging_IDS[] <- floor(IDSd[i] * u_aging)

# ----- 1.4.11. IDV1 Transitions -----

# IDV1 -> RDV1
c_IDV1_to_RDV1[] <- rbinom(IDV1[i],p_recov_D[i])
IDV1a[] <- IDV1[i] - c_IDV1_to_RDV1[i]

# Infection Deaths
c_death_inf_IDV1[] <- rbinom(IDV1a[i],p_death_IDV1[i])
IDV1b[] <- IDV1a[i] - c_death_inf_IDV1[i]

# IDV1 -> IDV2
c_IDV1_to_IDV2[] <- rbinom(IDV1b[i],p_V1_to_V2[i])
IDV1c[] <- IDV1b[i] - c_IDV1_to_IDV2[i]

# IDV1 -> IDS
c_IDV1_to_IDS[] <- rbinom(IDV1c[i],p_V1_to_S[i])
IDV1d[] <- IDV1c[i] - c_IDV1_to_IDS[i]

# Death
c_death_IDV1[] <- rbinom(IDV1d[i],p_death[i])
IDV1e[] <- IDV1d[i] - c_death_IDV1[i]

# Aging
c_aging_IDV1[] <- floor(IDV1e[i] * u_aging)

# ----- 1.4.11. IDV2 Transitions -----

# IDV2 -> RDV2
c_IDV2_to_RDV2[] <- rbinom(IDV2[i],p_recov_D[i])
IDV2a[] <- IDV2[i] - c_IDV2_to_RDV2[i]

# Infection Deaths
c_death_inf_IDV2[] <- rbinom(IDV2a[i],p_death_IDV2[i])
IDV2b[] <- IDV2a[i] - c_death_inf_IDV2[i]

# IDV2 -> IDS
c_IDV2_to_IDS[] <- rbinom(IDV2b[i],p_V2_to_S[i])
IDV2c[] <- IDV2b[i] - c_IDV2_to_IDS[i]

# Death
c_death_IDV2[] <- rbinom(IDV2c[i],p_death[i])
IDV2d[] <- IDV2c[i] - c_death_IDV2[i]

# Aging
c_aging_IDV2[] <- floor(IDV2d[i] * u_aging)

# ----- 1.4.11. IOS Transitions -----

# IOS -> ROS
c_IOS_to_ROS[] <- rbinom(IOS[i],p_recov_O[i])
IOSa[] <- IOS[i] - c_IOS_to_ROS[i]

# Infection Deaths
c_death_inf_IOS[] <- rbinom(IOSa[i],p_death_IOS[i])
IOSb[] <- IOSa[i] - c_death_inf_IOS[i]

# IOS -> IOV1
c_IOS_to_IOV1[] <- rbinom(IOSb[i],p_S_to_V1[i])
IOSc[] <- IOSb[i] - c_IOS_to_IOV1[i]

# Death
c_death_IOS[] <- rbinom(IOSc[i],p_death[i])
IOSd[] <- IOSc[i] - c_death_IOS[i]

# Aging
c_aging_IOS[] <- floor(IOSd[i] * u_aging)

# ----- 1.4.11. IOV1 Transitions -----

# IOV1 -> ROV1
c_IOV1_to_ROV1[] <- rbinom(IOV1[i],p_recov_O[i])
IOV1a[] <- IOV1[i] - c_IOV1_to_ROV1[i]

# Infection Deaths
c_death_inf_IOV1[] <- rbinom(IOV1a[i],p_death_IOV1[i])
IOV1b[] <- IOV1a[i] - c_death_inf_IOV1[i]

# IOV1 -> IOV2
c_IOV1_to_IOV2[] <- rbinom(IOV1b[i],p_V1_to_V2[i])
IOV1c[] <- IOV1b[i] - c_IOV1_to_IOV2[i]

# IOV1 -> IOS
c_IOV1_to_IOS[] <- rbinom(IOV1c[i],p_V1_to_S[i])
IOV1d[] <- IOV1c[i] - c_IOV1_to_IOS[i]

# Death
c_death_IOV1[] <- rbinom(IOV1d[i],p_death[i])
IOV1e[] <- IOV1d[i] - c_death_IOV1[i]

# Aging
c_aging_IOV1[] <- floor(IOV1e[i] * u_aging)

# ----- 1.4.11. IOV2 Transitions -----

# IOV2 -> ROV2
c_IOV2_to_ROV2[] <- rbinom(IOV2[i],p_recov_O[i])
IOV2a[] <- IOV2[i] - c_IOV2_to_ROV2[i]

# Infection Deaths
c_death_inf_IOV2[] <- rbinom(IOV2a[i],p_death_IOV2[i])
IOV2b[] <- IOV2a[i] - c_death_inf_IOV2[i]

# IOV2 -> IOS
c_IOV2_to_IOS[] <- rbinom(IOV2b[i],p_V2_to_S[i])
IOV2c[] <- IOV2b[i] - c_IOV2_to_IOS[i]

# Death
c_death_IOV2[] <- rbinom(IOV2c[i],p_death[i])
IOV2d[] <- IOV2c[i] - c_death_IOV2[i]

# Aging
c_aging_IOV2[] <- floor(IOV2d[i] * u_aging)

# ----- 1.4.19. RWS Transitions -----

# RWS -> EARWS
c_RWS_to_EARWS[] <- rbinom(RWS[i],p_RWS_to_EARWS[i])
RWSa[] <- RWS[i] - c_RWS_to_EARWS[i]

# RWS -> EDRWS
c_RWS_to_EDRWS[] <- rbinom(RWSa[i],p_RWS_to_EDRWS[i])
RWSb[] <- RWSa[i] - c_RWS_to_EDRWS[i]

# RWS -> EORWS
c_RWS_to_EORWS[] <- rbinom(RWSb[i],p_RWS_to_EORWS[i])
RWSc[] <- RWSb[i] - c_RWS_to_EORWS[i]

# RWS -> RWV1
c_RWS_to_RWV1[] <- rbinom(RWSc[i],p_S_to_V1[i])
RWSd[] <- RWSc[i] - c_RWS_to_RWV1[i]

# RWS -> S
c_RWS_to_S[] <- rbinom(RWSd[i],p_waning_W[i])
RWSe[] <- RWSd[i] - c_RWS_to_S[i]

# Death
c_death_RWS[] <- rbinom(RWSe[i],p_death[i])
RWSf[] <- RWSe[i] - c_death_RWS[i]

# Aging
c_aging_RWS[] <- floor(RWSf[i] * u_aging)

# ----- 1.4.19. RWV1 Transitions -----

# RWV1 -> EARWV1
c_RWV1_to_EARWV1[] <- rbinom(RWV1[i],p_RWV1_to_EARWV1[i])
RWV1a[] <- RWV1[i] - c_RWV1_to_EARWV1[i]

# RWV1 -> EDRWV1
c_RWV1_to_EDRWV1[] <- rbinom(RWV1a[i],p_RWV1_to_EDRWV1[i])
RWV1b[] <- RWV1a[i] - c_RWV1_to_EDRWV1[i]

# RWV1 -> EORWV1
c_RWV1_to_EORWV1[] <- rbinom(RWV1b[i],p_RWV1_to_EORWV1[i])
RWV1c[] <- RWV1b[i] - c_RWV1_to_EORWV1[i]

# RWV1 -> RWV2
c_RWV1_to_RWV2[] <- rbinom(RWV1c[i],p_V1_to_V2[i])
RWV1d[] <- RWV1c[i] - c_RWV1_to_RWV2[i]

# RWV1 -> V1
c_RWV1_to_V1[] <- rbinom(RWV1d[i],p_waning_W[i])
RWV1e[] <- RWV1d[i] - c_RWV1_to_V1[i]

# RWV1 -> RWS
c_RWV1_to_RWS[] <- rbinom(RWV1e[i],p_V1_to_S[i])
RWV1f[] <- RWV1e[i] - c_RWV1_to_RWS[i]

# Death
c_death_RWV1[] <- rbinom(RWV1f[i],p_death[i])
RWV1g[] <- RWV1f[i] - c_death_RWV1[i]

# Aging
c_aging_RWV1[] <- floor(RWV1g[i] * u_aging)

# ----- 1.4.19. RWV2 Transitions -----

# RWV2 -> EARWV2
c_RWV2_to_EARWV2[] <- rbinom(RWV2[i],p_RWV2_to_EARWV2[i])
RWV2a[] <- RWV2[i] - c_RWV2_to_EARWV2[i]

# RWV2 -> EDRWV2
c_RWV2_to_EDRWV2[] <- rbinom(RWV2a[i],p_RWV2_to_EDRWV2[i])
RWV2b[] <- RWV2a[i] - c_RWV2_to_EDRWV2[i]

# RWV2 -> EORWV2
c_RWV2_to_EORWV2[] <- rbinom(RWV2b[i],p_RWV2_to_EORWV2[i])
RWV2c[] <- RWV2b[i] - c_RWV2_to_EORWV2[i]

# RWV2 -> V2
c_RWV2_to_V2[] <- rbinom(RWV2c[i],p_waning_W[i])
RWV2d[] <- RWV2c[i] - c_RWV2_to_V2[i]

# RWV2 -> RWS
c_RWV2_to_RWS[] <- rbinom(RWV2d[i],p_V2_to_S[i])
RWV2e[] <- RWV2d[i] - c_RWV2_to_RWS[i]

# Death
c_death_RWV2[] <- rbinom(RWV2e[i],p_death[i])
RWV2f[] <- RWV2e[i] - c_death_RWV2[i]

# Aging
c_aging_RWV2[] <- floor(RWV2f[i] * u_aging)

# ----- 1.4.19. RAS Transitions -----

# RAS -> EDRAS
c_RAS_to_EDRAS[] <- rbinom(RAS[i],p_RAS_to_EDRAS[i])
RASa[] <- RAS[i] - c_RAS_to_EDRAS[i]

# RAS -> EORAS
c_RAS_to_EORAS[] <- rbinom(RASa[i],p_RAS_to_EORAS[i])
RASb[] <- RASa[i] - c_RAS_to_EORAS[i]

# RAS -> RAV1
c_RAS_to_RAV1[] <- rbinom(RASb[i],p_S_to_V1[i])
RASc[] <- RASb[i] - c_RAS_to_RAV1[i]

# RAS -> S
c_RAS_to_S[] <- rbinom(RASc[i],p_waning_A[i])
RASd[] <- RASc[i] - c_RAS_to_S[i]

# Death
c_death_RAS[] <- rbinom(RASd[i],p_death[i])
RASe[] <- RASd[i] - c_death_RAS[i]

# Aging
c_aging_RAS[] <- floor(RASe[i] * u_aging)

# ----- 1.4.19. RAV1 Transitions -----

# RAV1 -> EDRAV1
c_RAV1_to_EDRAV1[] <- rbinom(RAV1[i],p_RAV1_to_EDRAV1[i])
RAV1a[] <- RAV1[i] - c_RAV1_to_EDRAV1[i]

# RAV1 -> EORAV1
c_RAV1_to_EORAV1[] <- rbinom(RAV1a[i],p_RAV1_to_EORAV1[i])
RAV1b[] <- RAV1a[i] - c_RAV1_to_EORAV1[i]

# RAV1 -> RAV2
c_RAV1_to_RAV2[] <- rbinom(RAV1b[i],p_V1_to_V2[i])
RAV1c[] <- RAV1b[i] - c_RAV1_to_RAV2[i]

# RAV1 -> V1
c_RAV1_to_V1[] <- rbinom(RAV1c[i],p_waning_A[i])
RAV1d[] <- RAV1c[i] - c_RAV1_to_V1[i]

# RAV1 -> RAS
c_RAV1_to_RAS[] <- rbinom(RAV1d[i],p_V2_to_S[i])
RAV1e[] <- RAV1d[i] - c_RAV1_to_RAS[i]

# Death
c_death_RAV1[] <- rbinom(RAV1e[i],p_death[i])
RAV1f[] <- RAV1e[i] - c_death_RAV1[i]

# Aging
c_aging_RAV1[] <- floor(RAV1f[i] * u_aging)

# ----- 1.4.19. RAV2 Transitions -----

# RAV2 -> EDRAV2
c_RAV2_to_EDRAV2[] <- rbinom(RAV2[i],p_RAV2_to_EDRAV2[i])
RAV2a[] <- RAV2[i] - c_RAV2_to_EDRAV2[i]

# RAV2 -> EORAV2
c_RAV2_to_EORAV2[] <- rbinom(RAV2a[i],p_RAV2_to_EORAV2[i])
RAV2b[] <- RAV2a[i] - c_RAV2_to_EORAV2[i]

# RAV2 -> V2
c_RAV2_to_V2[] <- rbinom(RAV2b[i],p_waning_A[i])
RAV2c[] <- RAV2b[i] - c_RAV2_to_V2[i]

# RAV2 -> RAS
c_RAV2_to_RAS[] <- rbinom(RAV2c[i],p_V2_to_S[i])
RAV2d[] <- RAV2c[i] - c_RAV2_to_RAS[i]

# Death
c_death_RAV2[] <- rbinom(RAV2d[i],p_death[i])
RAV2e[] <- RAV2d[i] - c_death_RAV2[i]

# Aging
c_aging_RAV2[] <- floor(RAV2e[i] * u_aging)

# ----- 1.4.19. RDS Transitions -----

# RDS -> EORDS
c_RDS_to_EORDS[] <- rbinom(RDS[i],p_RDS_to_EORDS[i])
RDSa[] <- RDS[i] - c_RDS_to_EORDS[i]

# RDS -> RDV1
c_RDS_to_RDV1[] <- rbinom(RDSa[i],p_S_to_V1[i])
RDSb[] <- RDSa[i] - c_RDS_to_RDV1[i]

# RDS -> S
c_RDS_to_S[] <- rbinom(RDSb[i],p_waning_D[i])
RDSc[] <- RDSb[i] - c_RDS_to_S[i]

# Death
c_death_RDS[] <- rbinom(RDSc[i],p_death[i])
RDSd[] <- RDSc[i] - c_death_RDS[i]

# Aging
c_aging_RDS[] <- floor(RDSd[i] * u_aging)

# ----- 1.4.19. RDV1 Transitions -----

# RDV1 -> EORDV1
c_RDV1_to_EORDV1[] <- rbinom(RDV1[i],p_RDV1_to_EORDV1[i])
RDV1a[] <- RDV1[i] - c_RDV1_to_EORDV1[i]

# RDV1 -> RDV2
c_RDV1_to_RDV2[] <- rbinom(RDV1a[i],p_V1_to_V2[i])
RDV1b[] <- RDV1a[i] - c_RDV1_to_RDV2[i]

# RDV1 -> V1
c_RDV1_to_V1[] <- rbinom(RDV1b[i],p_waning_D[i])
RDV1c[] <- RDV1b[i] - c_RDV1_to_V1[i]

# RDV1 -> RDS
c_RDV1_to_RDS[] <- rbinom(RDV1c[i],p_V1_to_S[i])
RDV1d[] <- RDV1c[i] - c_RDV1_to_RDS[i]

# Death
c_death_RDV1[] <- rbinom(RDV1d[i],p_death[i])
RDV1e[] <- RDV1d[i] - c_death_RDV1[i]

# Aging
c_aging_RDV1[] <- floor(RDV1e[i] * u_aging)

# ----- 1.4.19. RDV2 Transitions -----

# RDV2 -> EORDV2
c_RDV2_to_EORDV2[] <- rbinom(RDV2[i],p_RDV2_to_EORDV2[i])
RDV2a[] <- RDV2[i] - c_RDV2_to_EORDV2[i]

# RDV2 -> V2
c_RDV2_to_V2[] <- rbinom(RDV2a[i],p_waning_D[i])
RDV2b[] <- RDV2a[i] - c_RDV2_to_V2[i]

# RDV2 -> RDS
c_RDV2_to_RDS[] <- rbinom(RDV2b[i],p_V2_to_S[i])
RDV2c[] <- RDV2b[i] - c_RDV2_to_RDS[i]

# Death
c_death_RDV2[] <- rbinom(RDV2c[i],p_death[i])
RDV2d[] <- RDV2c[i] - c_death_RDV2[i]

# Aging
c_aging_RDV2[] <- floor(RDV2d[i] * u_aging)

# ----- 1.4.19. ROS Transitions -----

# ROS -> EXROS
c_ROS_to_EXROS[] <- rbinom(ROS[i],p_ROS_to_EXROS[i])
ROSa[] <- ROS[i] - c_ROS_to_EXROS[i]

# ROS -> ROV1
c_ROS_to_ROV1[] <- rbinom(ROSa[i],p_S_to_V1[i])
ROSb[] <- ROSa[i] - c_ROS_to_ROV1[i]

# ROS -> S
c_ROS_to_S[] <- rbinom(ROSb[i],p_waning_O[i])
ROSc[] <- ROSb[i] - c_ROS_to_S[i]

# Death
c_death_ROS[] <- rbinom(ROSc[i],p_death[i])
ROSd[] <- ROSc[i] - c_death_ROS[i]

# Aging
c_aging_ROS[] <- floor(ROSd[i] * u_aging)

# ----- 1.4.19. ROV1 Transitions -----

# ROV1 -> EXROV1
c_ROV1_to_EXROV1[] <- rbinom(ROV1[i],p_ROV1_to_EXROV1[i])
ROV1a[] <- ROV1[i] - c_ROV1_to_EXROV1[i]

# ROV1 -> ROV2
c_ROV1_to_ROV2[] <- rbinom(ROV1a[i],p_V1_to_V2[i])
ROV1b[] <- ROV1a[i] - c_ROV1_to_ROV2[i]

# ROV1 -> V1
c_ROV1_to_V1[] <- rbinom(ROV1b[i],p_waning_O[i])
ROV1c[] <- ROV1b[i] - c_ROV1_to_V1[i]

# ROV1 -> ROS
c_ROV1_to_ROS[] <- rbinom(ROV1c[i],p_V2_to_S[i])
ROV1d[] <- ROV1c[i] - c_ROV1_to_ROS[i]

# Death
c_death_ROV1[] <- rbinom(ROV1d[i],p_death[i])
ROV1e[] <- ROV1d[i] - c_death_ROV1[i]

# Aging
c_aging_ROV1[] <- floor(ROV1e[i] * u_aging)

# ----- 1.4.19. ROV2 Transitions -----

# ROV2 -> EXROV2
c_ROV2_to_EXROV2[] <- rbinom(ROV2[i],p_ROV2_to_EXROV2[i])
ROV2a[] <- ROV2[i] - c_ROV2_to_EXROV2[i]

# ROV2 -> V2
c_ROV2_to_V2[] <- rbinom(ROV2a[i],p_waning_O[i])
ROV2b[] <- ROV2a[i] - c_ROV2_to_V2[i]

# ROV2 -> ROS
c_ROV2_to_ROS[] <- rbinom(ROV2b[i],p_V2_to_S[i])
ROV2c[] <- ROV2b[i] - c_ROV2_to_ROS[i]

# Death
c_death_ROV2[] <- rbinom(ROV2c[i],p_death[i])
ROV2d[] <- ROV2c[i] - c_death_ROV2[i]

# Aging
c_aging_ROV2[] <- floor(ROV2d[i] * u_aging)

# ----- 1.4.19. EARWS Transitions -----

# EARWS -> IARWS
c_EARWS_to_IARWS[] <- rbinom(EARWS[i],p_incub_A[i])
EARWSa[] <- EARWS[i] - c_EARWS_to_IARWS[i]

# EARWS -> EARWV1
c_EARWS_to_EARWV1[] <- rbinom(EARWSa[i],p_S_to_V1[i])
EARWSb[] <- EARWSa[i] - c_EARWS_to_EARWV1[i]

# Death
c_death_EARWS[] <- rbinom(EARWSb[i],p_death[i])
EARWSc[] <- EARWSb[i] - c_death_EARWS[i]

# Aging
c_aging_EARWS[] <- floor(EARWSc[i] * u_aging)

# ----- 1.4.19. EARWV1 Transitions -----

# EARWV1 -> IARWV1
c_EARWV1_to_IARWV1[] <- rbinom(EARWV1[i],p_incub_A[i])
EARWV1a[] <- EARWV1[i] - c_EARWV1_to_IARWV1[i]

# EARWV1 -> EARWV2
c_EARWV1_to_EARWV2[] <- rbinom(EARWV1a[i],p_V1_to_V2[i])
EARWV1b[] <- EARWV1a[i] - c_EARWV1_to_EARWV2[i]

# EARWV1 -> EARWS
c_EARWV1_to_EARWS[] <- rbinom(EARWV1b[i],p_V1_to_S[i])
EARWV1c[] <- EARWV1b[i] - c_EARWV1_to_EARWS[i]

# Death
c_death_EARWV1[] <- rbinom(EARWV1c[i],p_death[i])
EARWV1d[] <- EARWV1c[i] - c_death_EARWV1[i]

# Aging
c_aging_EARWV1[] <- floor(EARWV1d[i] * u_aging)

# ----- 1.4.19. EARWV2 Transitions -----

# EARWV2 -> IARWV2
c_EARWV2_to_IARWV2[] <- rbinom(EARWV2[i],p_incub_A[i])
EARWV2a[] <- EARWV2[i] - c_EARWV2_to_IARWV2[i]

# EARWV2 -> EARWS
c_EARWV2_to_EARWS[] <- rbinom(EARWV2a[i],p_V2_to_S[i])
EARWV2b[] <- EARWV2a[i] - c_EARWV2_to_EARWS[i]

# Death
c_death_EARWV2[] <- rbinom(EARWV2b[i],p_death[i])
EARWV2c[] <- EARWV2b[i] - c_death_EARWV2[i]

# Aging
c_aging_EARWV2[] <- floor(EARWV2c[i] * u_aging)

# ----- 1.4.19. EDRWS Transitions -----

# EDRWS -> IDRWS
c_EDRWS_to_IDRWS[] <- rbinom(EDRWS[i],p_incub_D[i])
EDRWSa[] <- EDRWS[i] - c_EDRWS_to_IDRWS[i]

# EDRWS -> EDRWV1
c_EDRWS_to_EDRWV1[] <- rbinom(EDRWSa[i],p_S_to_V1[i])
EDRWSb[] <- EDRWSa[i] - c_EDRWS_to_EDRWV1[i]

# Death
c_death_EDRWS[] <- rbinom(EDRWSb[i],p_death[i])
EDRWSc[] <- EDRWSb[i] - c_death_EDRWS[i]

# Aging
c_aging_EDRWS[] <- floor(EDRWSc[i] * u_aging)

# ----- 1.4.19. EDRWV1 Transitions -----

# EDRWV1 -> IDRWV1
c_EDRWV1_to_IDRWV1[] <- rbinom(EDRWV1[i],p_incub_D[i])
EDRWV1a[] <- EDRWV1[i] - c_EDRWV1_to_IDRWV1[i]

# EDRWV1 -> EDRWV2
c_EDRWV1_to_EDRWV2[] <- rbinom(EDRWV1a[i],p_V1_to_V2[i])
EDRWV1b[] <- EDRWV1a[i] - c_EDRWV1_to_EDRWV2[i]

# EDRWV1 -> EDRWS
c_EDRWV1_to_EDRWS[] <- rbinom(EDRWV1b[i],p_V1_to_S[i])
EDRWV1c[] <- EDRWV1b[i] - c_EDRWV1_to_EDRWS[i]

# Death
c_death_EDRWV1[] <- rbinom(EDRWV1c[i],p_death[i])
EDRWV1d[] <- EDRWV1c[i] - c_death_EDRWV1[i]

# Aging
c_aging_EDRWV1[] <- floor(EDRWV1d[i] * u_aging)

# ----- 1.4.19. EDRWV2 Transitions -----

# EDRWV2 -> IDRWV2
c_EDRWV2_to_IDRWV2[] <- rbinom(EDRWV2[i],p_incub_D[i])
EDRWV2a[] <- EDRWV2[i] - c_EDRWV2_to_IDRWV2[i]

# EDRWV2 -> EDRWS
c_EDRWV2_to_EDRWS[] <- rbinom(EDRWV2a[i],p_V2_to_S[i])
EDRWV2b[] <- EDRWV2a[i] - c_EDRWV2_to_EDRWS[i]

# Death
c_death_EDRWV2[] <- rbinom(EDRWV2b[i],p_death[i])
EDRWV2c[] <- EDRWV2b[i] - c_death_EDRWV2[i]

# Aging
c_aging_EDRWV2[] <- floor(EDRWV2c[i] * u_aging)

# ----- 1.4.19. EORWS Transitions -----

# EORWS -> IORWS
c_EORWS_to_IORWS[] <- rbinom(EORWS[i],p_incub_O[i])
EORWSa[] <- EORWS[i] - c_EORWS_to_IORWS[i]

# EORWS -> EORWV1
c_EORWS_to_EORWV1[] <- rbinom(EORWSa[i],p_S_to_V1[i])
EORWSb[] <- EORWSa[i] - c_EORWS_to_EORWV1[i]

# Death
c_death_EORWS[] <- rbinom(EORWSb[i],p_death[i])
EORWSc[] <- EORWSb[i] - c_death_EORWS[i]

# Aging
c_aging_EORWS[] <- floor(EORWSc[i] * u_aging)

# ----- 1.4.19. EORWV1 Transitions -----

# EORWV1 -> IORWV1
c_EORWV1_to_IORWV1[] <- rbinom(EORWV1[i],p_incub_O[i])
EORWV1a[] <- EORWV1[i] - c_EORWV1_to_IORWV1[i]

# EORWV1 -> EORWV2
c_EORWV1_to_EORWV2[] <- rbinom(EORWV1a[i],p_V1_to_V2[i])
EORWV1b[] <- EORWV1a[i] - c_EORWV1_to_EORWV2[i]

# EORWV1 -> EORWS
c_EORWV1_to_EORWS[] <- rbinom(EORWV1b[i],p_V1_to_S[i])
EORWV1c[] <- EORWV1b[i] - c_EORWV1_to_EORWS[i]

# Death
c_death_EORWV1[] <- rbinom(EORWV1c[i],p_death[i])
EORWV1d[] <- EORWV1c[i] - c_death_EORWV1[i]

# Aging
c_aging_EORWV1[] <- floor(EORWV1d[i] * u_aging)

# ----- 1.4.19. EORWV2 Transitions -----

# EORWV2 -> IORWV2
c_EORWV2_to_IORWV2[] <- rbinom(EORWV2[i],p_incub_O[i])
EORWV2a[] <- EORWV2[i] - c_EORWV2_to_IORWV2[i]

# EORWV2 -> EORWS
c_EORWV2_to_EORWS[] <- rbinom(EORWV2a[i],p_V2_to_S[i])
EORWV2b[] <- EORWV2a[i] - c_EORWV2_to_EORWS[i]

# Death
c_death_EORWV2[] <- rbinom(EORWV2b[i],p_death[i])
EORWV2c[] <- EORWV2b[i] - c_death_EORWV2[i]

# Aging
c_aging_EORWV2[] <- floor(EORWV2c[i] * u_aging)

# ----- 1.4.19. EDRAS Transitions -----

# EDRAS -> IDRAS
c_EDRAS_to_IDRAS[] <- rbinom(EDRAS[i],p_incub_D[i])
EDRASa[] <- EDRAS[i] - c_EDRAS_to_IDRAS[i]

# EDRAS -> EDRAV1
c_EDRAS_to_EDRAV1[] <- rbinom(EDRASa[i],p_S_to_V1[i])
EDRASb[] <- EDRASa[i] - c_EDRAS_to_EDRAV1[i]

# Death
c_death_EDRAS[] <- rbinom(EDRASb[i],p_death[i])
EDRASc[] <- EDRASb[i] - c_death_EDRAS[i]

# Aging
c_aging_EDRAS[] <- floor(EDRASc[i] * u_aging)

# ----- 1.4.19. EDRAV1 Transitions -----

# EDRAV1 -> IDRAV1
c_EDRAV1_to_IDRAV1[] <- rbinom(EDRAV1[i],p_incub_D[i])
EDRAV1a[] <- EDRAV1[i] - c_EDRAV1_to_IDRAV1[i]

# EDRAV1 -> EDRAV2
c_EDRAV1_to_EDRAV2[] <- rbinom(EDRAV1a[i],p_V1_to_V2[i])
EDRAV1b[] <- EDRAV1a[i] - c_EDRAV1_to_EDRAV2[i]

# EDRAV1 -> EDRAS
c_EDRAV1_to_EDRAS[] <- rbinom(EDRAV1b[i],p_V1_to_S[i])
EDRAV1c[] <- EDRAV1b[i] - c_EDRAV1_to_EDRAS[i]

# Death
c_death_EDRAV1[] <- rbinom(EDRAV1c[i],p_death[i])
EDRAV1d[] <- EDRAV1c[i] - c_death_EDRAV1[i]

# Aging
c_aging_EDRAV1[] <- floor(EDRAV1d[i] * u_aging)

# ----- 1.4.19. EDRAV2 Transitions -----

# EDRAV2 -> IDRAV2
c_EDRAV2_to_IDRAV2[] <- rbinom(EDRAV2[i],p_incub_D[i])
EDRAV2a[] <- EDRAV2[i] - c_EDRAV2_to_IDRAV2[i]

# EDRAV2 -> EDRAS
c_EDRAV2_to_EDRAS[] <- rbinom(EDRAV2a[i],p_V2_to_S[i])
EDRAV2b[] <- EDRAV2a[i] - c_EDRAV2_to_EDRAS[i]

# Death
c_death_EDRAV2[] <- rbinom(EDRAV2b[i],p_death[i])
EDRAV2c[] <- EDRAV2b[i] - c_death_EDRAV2[i]

# Aging
c_aging_EDRAV2[] <- floor(EDRAV2c[i] * u_aging)

# ----- 1.4.19. EORAS Transitions -----

# EORAS -> IORAS
c_EORAS_to_IORAS[] <- rbinom(EORAS[i],p_incub_O[i])
EORASa[] <- EORAS[i] - c_EORAS_to_IORAS[i]

# EORAS -> EORAV1
c_EORAS_to_EORAV1[] <- rbinom(EORASa[i],p_S_to_V1[i])
EORASb[] <- EORASa[i] - c_EORAS_to_EORAV1[i]

# Death
c_death_EORAS[] <- rbinom(EORASb[i],p_death[i])
EORASc[] <- EORASb[i] - c_death_EORAS[i]

# Aging
c_aging_EORAS[] <- floor(EORASc[i] * u_aging)

# ----- 1.4.19. EORAV1 Transitions -----

# EORAV1 -> IORAV1
c_EORAV1_to_IORAV1[] <- rbinom(EORAV1[i],p_incub_O[i])
EORAV1a[] <- EORAV1[i] - c_EORAV1_to_IORAV1[i]

# EORAV1 -> EORAV2
c_EORAV1_to_EORAV2[] <- rbinom(EORAV1a[i],p_V1_to_V2[i])
EORAV1b[] <- EORAV1a[i] - c_EORAV1_to_EORAV2[i]

# EORAV1 -> EORAS
c_EORAV1_to_EORAS[] <- rbinom(EORAV1b[i],p_V1_to_S[i])
EORAV1c[] <- EORAV1b[i] - c_EORAV1_to_EORAS[i]

# Death
c_death_EORAV1[] <- rbinom(EORAV1c[i],p_death[i])
EORAV1d[] <- EORAV1c[i] - c_death_EORAV1[i]

# Aging
c_aging_EORAV1[] <- floor(EORAV1d[i] * u_aging)

# ----- 1.4.19. EORAV2 Transitions -----

# EORAV2 -> IORAV2
c_EORAV2_to_IORAV2[] <- rbinom(EORAV2[i],p_incub_O[i])
EORAV2a[] <- EORAV2[i] - c_EORAV2_to_IORAV2[i]

# EORAV2 -> EORAS
c_EORAV2_to_EORAS[] <- rbinom(EORAV2a[i],p_V2_to_S[i])
EORAV2b[] <- EORAV2a[i] - c_EORAV2_to_EORAS[i]

# Death
c_death_EORAV2[] <- rbinom(EORAV2b[i],p_death[i])
EORAV2c[] <- EORAV2b[i] - c_death_EORAV2[i]

# Aging
c_aging_EORAV2[] <- floor(EORAV2c[i] * u_aging)

# ----- 1.4.19. EORDS Transitions -----

# EORDS -> IORDS
c_EORDS_to_IORDS[] <- rbinom(EORDS[i],p_incub_O[i])
EORDSa[] <- EORDS[i] - c_EORDS_to_IORDS[i]

# EORDS -> EORDV1
c_EORDS_to_EORDV1[] <- rbinom(EORDSa[i],p_S_to_V1[i])
EORDSb[] <- EORDSa[i] - c_EORDS_to_EORDV1[i]

# Death
c_death_EORDS[] <- rbinom(EORDSb[i],p_death[i])
EORDSc[] <- EORDSb[i] - c_death_EORDS[i]

# Aging
c_aging_EORDS[] <- floor(EORDSc[i] * u_aging)

# ----- 1.4.19. EORDV1 Transitions -----

# EORDV1 -> IORDV1
c_EORDV1_to_IORDV1[] <- rbinom(EORDV1[i],p_incub_O[i])
EORDV1a[] <- EORDV1[i] - c_EORDV1_to_IORDV1[i]

# EORDV1 -> EORDV2
c_EORDV1_to_EORDV2[] <- rbinom(EORDV1a[i],p_V1_to_V2[i])
EORDV1b[] <- EORDV1a[i] - c_EORDV1_to_EORDV2[i]

# EORDV1 -> EORDS
c_EORDV1_to_EORDS[] <- rbinom(EORDV1b[i],p_V1_to_S[i])
EORDV1c[] <- EORDV1b[i] - c_EORDV1_to_EORDS[i]

# Death
c_death_EORDV1[] <- rbinom(EORDV1c[i],p_death[i])
EORDV1d[] <- EORDV1c[i] - c_death_EORDV1[i]

# Aging
c_aging_EORDV1[] <- floor(EORDV1d[i] * u_aging)

# ----- 1.4.19. EORDV2 Transitions -----

# EORDV2 -> IORDV2
c_EORDV2_to_IORDV2[] <- rbinom(EORDV2[i],p_incub_O[i])
EORDV2a[] <- EORDV2[i] - c_EORDV2_to_IORDV2[i]

# EORDV2 -> EORDS
c_EORDV2_to_EORDS[] <- rbinom(EORDV2a[i],p_V2_to_S[i])
EORDV2b[] <- EORDV2a[i] - c_EORDV2_to_EORDS[i]

# Death
c_death_EORDV2[] <- rbinom(EORDV2b[i],p_death[i])
EORDV2c[] <- EORDV2b[i] - c_death_EORDV2[i]

# Aging
c_aging_EORDV2[] <- floor(EORDV2c[i] * u_aging)

# ----- 1.4.11. IARWS Transitions -----

# IARWS -> RAS
c_IARWS_to_RAS[] <- rbinom(IARWS[i],p_recov_A[i])
IARWSa[] <- IARWS[i] - c_IARWS_to_RAS[i]

# Infection Deaths
c_death_inf_IARWS[] <- rbinom(IARWSa[i],p_death_IARWS[i])
IARWSb[] <- IARWSa[i] - c_death_inf_IARWS[i]

# IARWS -> IARWV1
c_IARWS_to_IARWV1[] <- rbinom(IARWSb[i],p_S_to_V1[i])
IARWSc[] <- IARWSb[i] - c_IARWS_to_IARWV1[i]

# Death
c_death_IARWS[] <- rbinom(IARWSc[i],p_death[i])
IARWSd[] <- IARWSc[i] - c_death_IARWS[i]

# Aging
c_aging_IARWS[] <- floor(IARWSd[i] * u_aging)

# ----- 1.4.11. IARWV1 Transitions -----

# IARWV1 -> RAV1
c_IARWV1_to_RAV1[] <- rbinom(IARWV1[i],p_recov_A[i])
IARWV1a[] <- IARWV1[i] - c_IARWV1_to_RAV1[i]

# Infection Deaths
c_death_inf_IARWV1[] <- rbinom(IARWV1a[i],p_death_IARWV1[i])
IARWV1b[] <- IARWV1a[i] - c_death_inf_IARWV1[i]

# IARWV1 -> IARWV2
c_IARWV1_to_IARWV2[] <- rbinom(IARWV1b[i],p_V1_to_V2[i])
IARWV1c[] <- IARWV1b[i] - c_IARWV1_to_IARWV2[i]

# IARWV1 -> IARWS
c_IARWV1_to_IARWS[] <- rbinom(IARWV1c[i],p_V1_to_S[i])
IARWV1d[] <- IARWV1c[i] - c_IARWV1_to_IARWS[i]

# Death
c_death_IARWV1[] <- rbinom(IARWV1d[i],p_death[i])
IARWV1e[] <- IARWV1d[i] - c_death_IARWV1[i]

# Aging
c_aging_IARWV1[] <- floor(IARWV1e[i] * u_aging)

# ----- 1.4.11. IARWV2 Transitions -----

# IARWV2 -> RAV2
c_IARWV2_to_RAV2[] <- rbinom(IARWV2[i],p_recov_A[i])
IARWV2a[] <- IARWV2[i] - c_IARWV2_to_RAV2[i]

# Infection Deaths
c_death_inf_IARWV2[] <- rbinom(IARWV2a[i],p_death_IARWV2[i])
IARWV2b[] <- IARWV2a[i] - c_death_inf_IARWV2[i]

# IARWV2 -> IARWS
c_IARWV2_to_IARWS[] <- rbinom(IARWV2b[i],p_V2_to_S[i])
IARWV2c[] <- IARWV2b[i] - c_IARWV2_to_IARWS[i]

# Death
c_death_IARWV2[] <- rbinom(IARWV2c[i],p_death[i])
IARWV2d[] <- IARWV2c[i] - c_death_IARWV2[i]

# Aging
c_aging_IARWV2[] <- floor(IARWV2d[i] * u_aging)

# ----- 1.4.11. IDRWS Transitions -----

# IDRWS -> RDS
c_IDRWS_to_RDS[] <- rbinom(IDRWS[i],p_recov_D[i])
IDRWSa[] <- IDRWS[i] - c_IDRWS_to_RDS[i]

# Infection Deaths
c_death_inf_IDRWS[] <- rbinom(IDRWSa[i],p_death_IDRWS[i])
IDRWSb[] <- IDRWSa[i] - c_death_inf_IDRWS[i]

# IDRWS -> IDRWV1
c_IDRWS_to_IDRWV1[] <- rbinom(IDRWSb[i],p_S_to_V1[i])
IDRWSc[] <- IDRWSb[i] - c_IDRWS_to_IDRWV1[i]

# Death
c_death_IDRWS[] <- rbinom(IDRWSc[i],p_death[i])
IDRWSd[] <- IDRWSc[i] - c_death_IDRWS[i]

# Aging
c_aging_IDRWS[] <- floor(IDRWSd[i] * u_aging)

# ----- 1.4.11. IDRWV1 Transitions -----

# IDRWV1 -> RDV1
c_IDRWV1_to_RDV1[] <- rbinom(IDRWV1[i],p_recov_D[i])
IDRWV1a[] <- IDRWV1[i] - c_IDRWV1_to_RDV1[i]

# Infection Deaths
c_death_inf_IDRWV1[] <- rbinom(IDRWV1a[i],p_death_IDRWV1[i])
IDRWV1b[] <- IDRWV1a[i] - c_death_inf_IDRWV1[i]

# IDRWV1 -> IDRWV2
c_IDRWV1_to_IDRWV2[] <- rbinom(IDRWV1b[i],p_V1_to_V2[i])
IDRWV1c[] <- IDRWV1b[i] - c_IDRWV1_to_IDRWV2[i]

# IDRWV1 -> IDRWS
c_IDRWV1_to_IDRWS[] <- rbinom(IDRWV1c[i],p_V1_to_S[i])
IDRWV1d[] <- IDRWV1c[i] - c_IDRWV1_to_IDRWS[i]

# Death
c_death_IDRWV1[] <- rbinom(IDRWV1d[i],p_death[i])
IDRWV1e[] <- IDRWV1d[i] - c_death_IDRWV1[i]

# Aging
c_aging_IDRWV1[] <- floor(IDRWV1e[i] * u_aging)

# ----- 1.4.11. IDRWV2 Transitions -----

# IDRWV2 -> RDV2
c_IDRWV2_to_RDV2[] <- rbinom(IDRWV2[i],p_recov_D[i])
IDRWV2a[] <- IDRWV2[i] - c_IDRWV2_to_RDV2[i]

# Infection Deaths
c_death_inf_IDRWV2[] <- rbinom(IDRWV2a[i],p_death_IDRWV2[i])
IDRWV2b[] <- IDRWV2a[i] - c_death_inf_IDRWV2[i]

# IDRWV2 -> IDRWS
c_IDRWV2_to_IDRWS[] <- rbinom(IDRWV2b[i],p_V2_to_S[i])
IDRWV2c[] <- IDRWV2b[i] - c_IDRWV2_to_IDRWS[i]

# Death
c_death_IDRWV2[] <- rbinom(IDRWV2c[i],p_death[i])
IDRWV2d[] <- IDRWV2c[i] - c_death_IDRWV2[i]

# Aging
c_aging_IDRWV2[] <- floor(IDRWV2d[i] * u_aging)

# ----- 1.4.11. IORWS Transitions -----

# IORWS -> ROS
c_IORWS_to_ROS[] <- rbinom(IORWS[i],p_recov_O[i])
IORWSa[] <- IORWS[i] - c_IORWS_to_ROS[i]

# Infection Deaths
c_death_inf_IORWS[] <- rbinom(IORWSa[i],p_death_IORWS[i])
IORWSb[] <- IORWSa[i] - c_death_inf_IORWS[i]

# IORWS -> IORWV1
c_IORWS_to_IORWV1[] <- rbinom(IORWSb[i],p_S_to_V1[i])
IORWSc[] <- IORWSb[i] - c_IORWS_to_IORWV1[i]

# Death
c_death_IORWS[] <- rbinom(IORWSc[i],p_death[i])
IORWSd[] <- IORWSc[i] - c_death_IORWS[i]

# Aging
c_aging_IORWS[] <- floor(IORWSd[i] * u_aging)

# ----- 1.4.11. IORWV1 Transitions -----

# IORWV1 -> ROV1
c_IORWV1_to_ROV1[] <- rbinom(IORWV1[i],p_recov_O[i])
IORWV1a[] <- IORWV1[i] - c_IORWV1_to_ROV1[i]

# Infection Deaths
c_death_inf_IORWV1[] <- rbinom(IORWV1a[i],p_death_IORWV1[i])
IORWV1b[] <- IORWV1a[i] - c_death_inf_IORWV1[i]

# IORWV1 -> IORWV2
c_IORWV1_to_IORWV2[] <- rbinom(IORWV1b[i],p_V1_to_V2[i])
IORWV1c[] <- IORWV1b[i] - c_IORWV1_to_IORWV2[i]

# IORWV1 -> IORWS
c_IORWV1_to_IORWS[] <- rbinom(IORWV1c[i],p_V1_to_S[i])
IORWV1d[] <- IORWV1c[i] - c_IORWV1_to_IORWS[i]

# Death
c_death_IORWV1[] <- rbinom(IORWV1d[i],p_death[i])
IORWV1e[] <- IORWV1d[i] - c_death_IORWV1[i]

# Aging
c_aging_IORWV1[] <- floor(IORWV1e[i] * u_aging)

# ----- 1.4.11. IORWV2 Transitions -----

# IORWV2 -> ROV2
c_IORWV2_to_ROV2[] <- rbinom(IORWV2[i],p_recov_O[i])
IORWV2a[] <- IORWV2[i] - c_IORWV2_to_ROV2[i]

# Infection Deaths
c_death_inf_IORWV2[] <- rbinom(IORWV2a[i],p_death_IORWV2[i])
IORWV2b[] <- IORWV2a[i] - c_death_inf_IORWV2[i]

# IORWV2 -> IORWS
c_IORWV2_to_IORWS[] <- rbinom(IORWV2b[i],p_V2_to_S[i])
IORWV2c[] <- IORWV2b[i] - c_IORWV2_to_IORWS[i]

# Death
c_death_IORWV2[] <- rbinom(IORWV2c[i],p_death[i])
IORWV2d[] <- IORWV2c[i] - c_death_IORWV2[i]

# Aging
c_aging_IORWV2[] <- floor(IORWV2d[i] * u_aging)

# ----- 1.4.11. IDRAS Transitions -----

# IDRAS -> RDS
c_IDRAS_to_RDS[] <- rbinom(IDRAS[i],p_recov_D[i])
IDRASa[] <- IDRAS[i] - c_IDRAS_to_RDS[i]

# Infection Deaths
c_death_inf_IDRAS[] <- rbinom(IDRASa[i],p_death_IDRAS[i])
IDRASb[] <- IDRASa[i] - c_death_inf_IDRAS[i]

# IDRAS -> IDRAV1
c_IDRAS_to_IDRAV1[] <- rbinom(IDRASb[i],p_S_to_V1[i])
IDRASc[] <- IDRASb[i] - c_IDRAS_to_IDRAV1[i]

# Death
c_death_IDRAS[] <- rbinom(IDRASc[i],p_death[i])
IDRASd[] <- IDRASc[i] - c_death_IDRAS[i]

# Aging
c_aging_IDRAS[] <- floor(IDRASd[i] * u_aging)

# ----- 1.4.11. IDRAV1 Transitions -----

# IDRAV1 -> RDV1
c_IDRAV1_to_RDV1[] <- rbinom(IDRAV1[i],p_recov_D[i])
IDRAV1a[] <- IDRAV1[i] - c_IDRAV1_to_RDV1[i]

# Infection Deaths
c_death_inf_IDRAV1[] <- rbinom(IDRAV1a[i],p_death_IDRAV1[i])
IDRAV1b[] <- IDRAV1a[i] - c_death_inf_IDRAV1[i]

# IDRAV1 -> IDRAV2
c_IDRAV1_to_IDRAV2[] <- rbinom(IDRAV1b[i],p_V1_to_V2[i])
IDRAV1c[] <- IDRAV1b[i] - c_IDRAV1_to_IDRAV2[i]

# IDRAV1 -> IDRAS
c_IDRAV1_to_IDRAS[] <- rbinom(IDRAV1c[i],p_V1_to_S[i])
IDRAV1d[] <- IDRAV1c[i] - c_IDRAV1_to_IDRAS[i]

# Death
c_death_IDRAV1[] <- rbinom(IDRAV1d[i],p_death[i])
IDRAV1e[] <- IDRAV1d[i] - c_death_IDRAV1[i]

# Aging
c_aging_IDRAV1[] <- floor(IDRAV1e[i] * u_aging)

# ----- 1.4.11. IDRAV2 Transitions -----

# IDRAV2 -> RDV2
c_IDRAV2_to_RDV2[] <- rbinom(IDRAV2[i],p_recov_D[i])
IDRAV2a[] <- IDRAV2[i] - c_IDRAV2_to_RDV2[i]

# Infection Deaths
c_death_inf_IDRAV2[] <- rbinom(IDRAV2a[i],p_death_IDRAV2[i])
IDRAV2b[] <- IDRAV2a[i] - c_death_inf_IDRAV2[i]

# IDRAV2 -> IDRAS
c_IDRAV2_to_IDRAS[] <- rbinom(IDRAV2b[i],p_V2_to_S[i])
IDRAV2c[] <- IDRAV2b[i] - c_IDRAV2_to_IDRAS[i]

# Death
c_death_IDRAV2[] <- rbinom(IDRAV2c[i],p_death[i])
IDRAV2d[] <- IDRAV2c[i] - c_death_IDRAV2[i]

# Aging
c_aging_IDRAV2[] <- floor(IDRAV2d[i] * u_aging)

# ----- 1.4.11. IORAS Transitions -----

# IORAS -> ROS
c_IORAS_to_ROS[] <- rbinom(IORAS[i],p_recov_O[i])
IORASa[] <- IORAS[i] - c_IORAS_to_ROS[i]

# Infection Deaths
c_death_inf_IORAS[] <- rbinom(IORASa[i],p_death_IORAS[i])
IORASb[] <- IORASa[i] - c_death_inf_IORAS[i]

# IORAS -> IORAV1
c_IORAS_to_IORAV1[] <- rbinom(IORASb[i],p_S_to_V1[i])
IORASc[] <- IORASb[i] - c_IORAS_to_IORAV1[i]

# Death
c_death_IORAS[] <- rbinom(IORASc[i],p_death[i])
IORASd[] <- IORASc[i] - c_death_IORAS[i]

# Aging
c_aging_IORAS[] <- floor(IORASd[i] * u_aging)

# ----- 1.4.11. IORAV1 Transitions -----

# IORAV1 -> ROV1
c_IORAV1_to_ROV1[] <- rbinom(IORAV1[i],p_recov_O[i])
IORAV1a[] <- IORAV1[i] - c_IORAV1_to_ROV1[i]

# Infection Deaths
c_death_inf_IORAV1[] <- rbinom(IORAV1a[i],p_death_IORAV1[i])
IORAV1b[] <- IORAV1a[i] - c_death_inf_IORAV1[i]

# IORAV1 -> IORAV2
c_IORAV1_to_IORAV2[] <- rbinom(IORAV1b[i],p_V1_to_V2[i])
IORAV1c[] <- IORAV1b[i] - c_IORAV1_to_IORAV2[i]

# IORAV1 -> IORAS
c_IORAV1_to_IORAS[] <- rbinom(IORAV1c[i],p_V1_to_S[i])
IORAV1d[] <- IORAV1c[i] - c_IORAV1_to_IORAS[i]

# Death
c_death_IORAV1[] <- rbinom(IORAV1d[i],p_death[i])
IORAV1e[] <- IORAV1d[i] - c_death_IORAV1[i]

# Aging
c_aging_IORAV1[] <- floor(IORAV1e[i] * u_aging)

# ----- 1.4.11. IORAV2 Transitions -----

# IORAV2 -> ROV2
c_IORAV2_to_ROV2[] <- rbinom(IORAV2[i],p_recov_O[i])
IORAV2a[] <- IORAV2[i] - c_IORAV2_to_ROV2[i]

# Infection Deaths
c_death_inf_IORAV2[] <- rbinom(IORAV2a[i],p_death_IORAV2[i])
IORAV2b[] <- IORAV2a[i] - c_death_inf_IORAV2[i]

# IORAV2 -> IORAS
c_IORAV2_to_IORAS[] <- rbinom(IORAV2b[i],p_V2_to_S[i])
IORAV2c[] <- IORAV2b[i] - c_IORAV2_to_IORAS[i]

# Death
c_death_IORAV2[] <- rbinom(IORAV2c[i],p_death[i])
IORAV2d[] <- IORAV2c[i] - c_death_IORAV2[i]

# Aging
c_aging_IORAV2[] <- floor(IORAV2d[i] * u_aging)

# ----- 1.4.11. IORDS Transitions -----

# IORDS -> ROS
c_IORDS_to_ROS[] <- rbinom(IORDS[i],p_recov_O[i])
IORDSa[] <- IORDS[i] - c_IORDS_to_ROS[i]

# Infection Deaths
c_death_inf_IORDS[] <- rbinom(IORDSa[i],p_death_IORDS[i])
IORDSb[] <- IORDSa[i] - c_death_inf_IORDS[i]

# IORDS -> IORDV1
c_IORDS_to_IORDV1[] <- rbinom(IORDSb[i],p_S_to_V1[i])
IORDSc[] <- IORDSb[i] - c_IORDS_to_IORDV1[i]

# Death
c_death_IORDS[] <- rbinom(IORDSc[i],p_death[i])
IORDSd[] <- IORDSc[i] - c_death_IORDS[i]

# Aging
c_aging_IORDS[] <- floor(IORDSd[i] * u_aging)

# ----- 1.4.11. IORDV1 Transitions -----

# IORDV1 -> ROV1
c_IORDV1_to_ROV1[] <- rbinom(IORDV1[i],p_recov_O[i])
IORDV1a[] <- IORDV1[i] - c_IORDV1_to_ROV1[i]

# Infection Deaths
c_death_inf_IORDV1[] <- rbinom(IORDV1a[i],p_death_IORDV1[i])
IORDV1b[] <- IORDV1a[i] - c_death_inf_IORDV1[i]

# IORDV1 -> IORDV2
c_IORDV1_to_IORDV2[] <- rbinom(IORDV1b[i],p_V1_to_V2[i])
IORDV1c[] <- IORDV1b[i] - c_IORDV1_to_IORDV2[i]

# IORDV1 -> IORDS
c_IORDV1_to_IORDS[] <- rbinom(IORDV1c[i],p_V1_to_S[i])
IORDV1d[] <- IORDV1c[i] - c_IORDV1_to_IORDS[i]

# Death
c_death_IORDV1[] <- rbinom(IORDV1d[i],p_death[i])
IORDV1e[] <- IORDV1d[i] - c_death_IORDV1[i]

# Aging
c_aging_IORDV1[] <- floor(IORDV1e[i] * u_aging)

# ----- 1.4.11. IORDV2 Transitions -----

# IORDV2 -> ROV2
c_IORDV2_to_ROV2[] <- rbinom(IORDV2[i],p_recov_O[i])
IORDV2a[] <- IORDV2[i] - c_IORDV2_to_ROV2[i]

# Infection Deaths
c_death_inf_IORDV2[] <- rbinom(IORDV2a[i],p_death_IORDV2[i])
IORDV2b[] <- IORDV2a[i] - c_death_inf_IORDV2[i]

# IORDV2 -> IORDS
c_IORDV2_to_IORDS[] <- rbinom(IORDV2b[i],p_V2_to_S[i])
IORDV2c[] <- IORDV2b[i] - c_IORDV2_to_IORDS[i]

# Death
c_death_IORDV2[] <- rbinom(IORDV2c[i],p_death[i])
IORDV2d[] <- IORDV2c[i] - c_death_IORDV2[i]

# Aging
c_aging_IORDV2[] <- floor(IORDV2d[i] * u_aging)

# ----- 1.4.3. EXS Transitions -----

# EXS -> IXS
c_EXS_to_IXS[] <- rbinom(EXS[i],p_incub_X[i])
EXSa[] <- EXS[i] - c_EXS_to_IXS[i]

# EXS -> EXV1
c_EXS_to_EXV1[] <- rbinom(EXSa[i],p_S_to_V1[i])
EXSb[] <- EXSa[i] - c_EXS_to_EXV1[i]

# Death
c_death_EXS[] <- rbinom(EXSb[i],p_death[i])
EXSc[] <- EXSb[i] - c_death_EXS[i]

# Aging
c_aging_EXS[] <- floor(EXSc[i] * u_aging)

# ----- 1.4.3. EXV1 Transitions -----

# EXV1 -> IXV1
c_EXV1_to_IXV1[] <- rbinom(EXV1[i],p_incub_X[i])
EXV1a[] <- EXV1[i] - c_EXV1_to_IXV1[i]

# EXV1 -> EXV2
c_EXV1_to_EXV2[] <- rbinom(EXV1a[i],p_V1_to_V2[i])
EXV1b[] <- EXV1a[i] - c_EXV1_to_EXV2[i]

# EXV1 -> EXS
c_EXV1_to_EXS[] <- rbinom(EXV1b[i],p_V1_to_S[i])
EXV1c[] <- EXV1b[i] - c_EXV1_to_EXS[i]

# Death
c_death_EXV1[] <- rbinom(EXV1c[i],p_death[i])
EXV1d[] <- EXV1c[i] - c_death_EXV1[i]

# Aging
c_aging_EXV1[] <- floor(EXV1d[i] * u_aging)

# ----- 1.4.3. EXV2 Transitions -----

# EXV2 -> IXV2
c_EXV2_to_IXV2[] <- rbinom(EXV2[i],p_incub_X[i])
EXV2a[] <- EXV2[i] - c_EXV2_to_IXV2[i]

# EXV2 -> EXS
c_EXV2_to_EXS[] <- rbinom(EXV2a[i],p_V2_to_S[i])
EXV2b[] <- EXV2a[i] - c_EXV2_to_EXS[i]

# Death
c_death_EXV2[] <- rbinom(EXV2b[i],p_death[i])
EXV2c[] <- EXV2b[i] - c_death_EXV2[i]

# Aging
c_aging_EXV2[] <- floor(EXV2c[i] * u_aging)

# ----- 1.4.11. IXS Transitions -----

# IXS -> RXS
c_IXS_to_RXS[] <- rbinom(IXS[i],p_recov_X[i])
IXSa[] <- IXS[i] - c_IXS_to_RXS[i]

# Infection Deaths
c_death_inf_IXS[] <- rbinom(IXSa[i],p_death_IXS[i])
IXSb[] <- IXSa[i] - c_death_inf_IXS[i]

# IXS -> IXV1
c_IXS_to_IXV1[] <- rbinom(IXSb[i],p_S_to_V1[i])
IXSc[] <- IXSb[i] - c_IXS_to_IXV1[i]

# Death
c_death_IXS[] <- rbinom(IXSc[i],p_death[i])
IXSd[] <- IXSc[i] - c_death_IXS[i]

# Aging
c_aging_IXS[] <- floor(IXSd[i] * u_aging)

# ----- 1.4.11. IXV1 Transitions -----

# IXV1 -> RXV1
c_IXV1_to_RXV1[] <- rbinom(IXV1[i],p_recov_X[i])
IXV1a[] <- IXV1[i] - c_IXV1_to_RXV1[i]

# Infection Deaths
c_death_inf_IXV1[] <- rbinom(IXV1a[i],p_death_IXV1[i])
IXV1b[] <- IXV1a[i] - c_death_inf_IXV1[i]

# IXV1 -> IXV2
c_IXV1_to_IXV2[] <- rbinom(IXV1b[i],p_V1_to_V2[i])
IXV1c[] <- IXV1b[i] - c_IXV1_to_IXV2[i]

# IXV1 -> IXS
c_IXV1_to_IXS[] <- rbinom(IXV1c[i],p_V1_to_S[i])
IXV1d[] <- IXV1c[i] - c_IXV1_to_IXS[i]

# Death
c_death_IXV1[] <- rbinom(IXV1d[i],p_death[i])
IXV1e[] <- IXV1d[i] - c_death_IXV1[i]

# Aging
c_aging_IXV1[] <- floor(IXV1e[i] * u_aging)

# ----- 1.4.11. IXV2 Transitions -----

# IXV2 -> RXV2
c_IXV2_to_RXV2[] <- rbinom(IXV2[i],p_recov_X[i])
IXV2a[] <- IXV2[i] - c_IXV2_to_RXV2[i]

# Infection Deaths
c_death_inf_IXV2[] <- rbinom(IXV2a[i],p_death_IXV2[i])
IXV2b[] <- IXV2a[i] - c_death_inf_IXV2[i]

# IXV2 -> IXS
c_IXV2_to_IXS[] <- rbinom(IXV2b[i],p_V2_to_S[i])
IXV2c[] <- IXV2b[i] - c_IXV2_to_IXS[i]

# Death
c_death_IXV2[] <- rbinom(IXV2c[i],p_death[i])
IXV2d[] <- IXV2c[i] - c_death_IXV2[i]

# Aging
c_aging_IXV2[] <- floor(IXV2d[i] * u_aging)

# ----- 1.4.19. RXS Transitions -----

# RXS -> EORXS
c_RXS_to_EORXS[] <- rbinom(RXS[i],p_RXS_to_EORXS[i])
RXSa[] <- RXS[i] - c_RXS_to_EORXS[i]

# RXS -> RXV1
c_RXS_to_RXV1[] <- rbinom(RXSa[i],p_S_to_V1[i])
RXSb[] <- RXSa[i] - c_RXS_to_RXV1[i]

# RXS -> S
c_RXS_to_S[] <- rbinom(RXSb[i],p_waning_X[i])
RXSc[] <- RXSb[i] - c_RXS_to_S[i]

# Death
c_death_RXS[] <- rbinom(RXSc[i],p_death[i])
RXSd[] <- RXSc[i] - c_death_RXS[i]

# Aging
c_aging_RXS[] <- floor(RXSd[i] * u_aging)

# ----- 1.4.19. RXV1 Transitions -----

# RXV1 -> EORXV1
c_RXV1_to_EORXV1[] <- rbinom(RXV1[i],p_RXV1_to_EORXV1[i])
RXV1a[] <- RXV1[i] - c_RXV1_to_EORXV1[i]

# RXV1 -> RXV1
c_RXV1_to_RXV2[] <- rbinom(RXV1a[i],p_V1_to_V2[i])
RXV1b[] <- RXV1a[i] - c_RXV1_to_RXV2[i]

# RXV1 -> V1
c_RXV1_to_V1[] <- rbinom(RXV1b[i],p_waning_X[i])
RXV1c[] <- RXV1b[i] - c_RXV1_to_V1[i]

# RXV1 -> RXS
c_RXV1_to_RXS[] <- rbinom(RXV1c[i],p_V1_to_S[i])
RXV1d[] <- RXV1c[i] - c_RXV1_to_RXS[i]

# Death
c_death_RXV1[] <- rbinom(RXV1d[i],p_death[i])
RXV1e[] <- RXV1d[i] - c_death_RXV1[i]

# Aging
c_aging_RXV1[] <- floor(RXV1e[i] * u_aging)

# ----- 1.4.19. RXV2 Transitions -----

# RXV2 -> EORXV2
c_RXV2_to_EORXV2[] <- rbinom(RXV2[i],p_RXV2_to_EORXV2[i])
RXV2a[] <- RXV2[i] - c_RXV2_to_EORXV2[i]

# RXV2 -> V2
c_RXV2_to_V2[] <- rbinom(RXV2a[i],p_waning_X[i])
RXV2b[] <- RXV2a[i] - c_RXV2_to_V2[i]

# RXV2 -> RXS
c_RXV2_to_RXS[] <- rbinom(RXV2b[i],p_V2_to_S[i])
RXV2c[] <- RXV2b[i] - c_RXV2_to_RXS[i]

# Death
c_death_RXV2[] <- rbinom(RXV2c[i],p_death[i])
RXV2d[] <- RXV2c[i] - c_death_RXV2[i]

# Aging
c_aging_RXV2[] <- floor(RXV2d[i] * u_aging)

# ----- 1.4.3. EORXS Transitions -----

# EORXS -> IORXS
c_EORXS_to_IORXS[] <- rbinom(EORXS[i],p_incub_O[i])
EORXSa[] <- EORXS[i] - c_EORXS_to_IORXS[i]

# EORXS -> EORXV1
c_EORXS_to_EORXV1[] <- rbinom(EORXSa[i],p_S_to_V1[i])
EORXSb[] <- EORXSa[i] - c_EORXS_to_EORXV1[i]

# Death
c_death_EORXS[] <- rbinom(EORXSb[i],p_death[i])
EORXSc[] <- EORXSb[i] - c_death_EORXS[i]

# Aging
c_aging_EORXS[] <- floor(EORXSc[i] * u_aging)

# ----- 1.4.3. EORXV1 Transitions -----

# EORXV1 -> IORXV1
c_EORXV1_to_IORXV1[] <- rbinom(EORXV1[i],p_incub_O[i])
EORXV1a[] <- EORXV1[i] - c_EORXV1_to_IORXV1[i]

# EORXV1 -> EORXV2
c_EORXV1_to_EORXV2[] <- rbinom(EORXV1a[i],p_V1_to_V2[i])
EORXV1b[] <- EORXV1a[i] - c_EORXV1_to_EORXV2[i]

# EORXV1 -> EORXS
c_EORXV1_to_EORXS[] <- rbinom(EORXV1b[i],p_V1_to_S[i])
EORXV1c[] <- EORXV1b[i] - c_EORXV1_to_EORXS[i]

# Death
c_death_EORXV1[] <- rbinom(EORXV1c[i],p_death[i])
EORXV1d[] <- EORXV1c[i] - c_death_EORXV1[i]

# Aging
c_aging_EORXV1[] <- floor(EORXV1d[i] * u_aging)

# ----- 1.4.3. EORXV2 Transitions -----

# EORXV2 -> IORXV2
c_EORXV2_to_IORXV2[] <- rbinom(EORXV2[i],p_incub_O[i])
EORXV2a[] <- EORXV2[i] - c_EORXV2_to_IORXV2[i]

# EORXV2 -> EORXS
c_EORXV2_to_EORXS[] <- rbinom(EORXV2a[i],p_V2_to_S[i])
EORXV2b[] <- EORXV2a[i] - c_EORXV2_to_EORXS[i]

# Death
c_death_EORXV2[] <- rbinom(EORXV2b[i],p_death[i])
EORXV2c[] <- EORXV2b[i] - c_death_EORXV2[i]

# Aging
c_aging_EORXV2[] <- floor(EORXV2c[i] * u_aging)

# ----- 1.4.11. IORXS Transitions -----

# IORXS -> ROXS
c_IORXS_to_ROXS[] <- rbinom(IORXS[i],p_recov_O[i])
IORXSa[] <- IORXS[i] - c_IORXS_to_ROXS[i]

# Infection Deaths
c_death_inf_IORXS[] <- rbinom(IORXSa[i],p_death_IORXS[i])
IORXSb[] <- IORXSa[i] - c_death_inf_IORXS[i]

# IORXS -> IORXV1
c_IORXS_to_IORXV1[] <- rbinom(IORXSb[i],p_S_to_V1[i])
IORXSc[] <- IORXSb[i] - c_IORXS_to_IORXV1[i]

# Death
c_death_IORXS[] <- rbinom(IORXSc[i],p_death[i])
IORXSd[] <- IORXSc[i] - c_death_IORXS[i]

# Aging
c_aging_IORXS[] <- floor(IORXSd[i] * u_aging)

# ----- 1.4.11. IORXV1 Transitions -----

# IORXV1 -> ROXV1
c_IORXV1_to_ROXV1[] <- rbinom(IORXV1[i],p_recov_O[i])
IORXV1a[] <- IORXV1[i] - c_IORXV1_to_ROXV1[i]

# Infection Deaths
c_death_inf_IORXV1[] <- rbinom(IORXV1a[i],p_death_IORXV1[i])
IORXV1b[] <- IORXV1a[i] - c_death_inf_IORXV1[i]

# IORXV1 -> IORXV2
c_IORXV1_to_IORXV2[] <- rbinom(IORXV1b[i],p_V1_to_V2[i])
IORXV1c[] <- IORXV1b[i] - c_IORXV1_to_IORXV2[i]

# IORXV1 -> IORXS
c_IORXV1_to_IORXS[] <- rbinom(IORXV1c[i],p_V1_to_S[i])
IORXV1d[] <- IORXV1c[i] - c_IORXV1_to_IORXS[i]

# Death
c_death_IORXV1[] <- rbinom(IORXV1d[i],p_death[i])
IORXV1e[] <- IORXV1d[i] - c_death_IORXV1[i]

# Aging
c_aging_IORXV1[] <- floor(IORXV1e[i] * u_aging)

# ----- 1.4.11. IORXV2 Transitions -----

# IORXV2 -> ROXV2
c_IORXV2_to_ROXV2[] <- rbinom(IORXV2[i],p_recov_O[i])
IORXV2a[] <- IORXV2[i] - c_IORXV2_to_ROXV2[i]

# Infection Deaths
c_death_inf_IORXV2[] <- rbinom(IORXV2a[i],p_death_IORXV2[i])
IORXV2b[] <- IORXV2a[i] - c_death_inf_IORXV2[i]

# IORXV2 -> IORXS
c_IORXV2_to_IORXS[] <- rbinom(IORXV2b[i],p_V2_to_S[i])
IORXV2c[] <- IORXV2b[i] - c_IORXV2_to_IORXS[i]

# Death
c_death_IORXV2[] <- rbinom(IORXV2c[i],p_death[i])
IORXV2d[] <- IORXV2c[i] - c_death_IORXV2[i]

# Aging
c_aging_IORXV2[] <- floor(IORXV2d[i] * u_aging)

# ----- 1.4.3. EXROS Transitions -----

# EXROS -> IXROS
c_EXROS_to_IXROS[] <- rbinom(EXROS[i],p_incub_X[i])
EXROSa[] <- EXROS[i] - c_EXROS_to_IXROS[i]

# EXROS -> EXROV1
c_EXROS_to_EXROV1[] <- rbinom(EXROSa[i],p_S_to_V1[i])
EXROSb[] <- EXROSa[i] - c_EXROS_to_EXROV1[i]

# Death
c_death_EXROS[] <- rbinom(EXROSb[i],p_death[i])
EXROSc[] <- EXROSb[i] - c_death_EXROS[i]

# Aging
c_aging_EXROS[] <- floor(EXROSc[i] * u_aging)

# ----- 1.4.3. EXROV1 Transitions -----

# EXROV1 -> IXROV1
c_EXROV1_to_IXROV1[] <- rbinom(EXROV1[i],p_incub_X[i])
EXROV1a[] <- EXROV1[i] - c_EXROV1_to_IXROV1[i]

# EXROV1 -> EXROV2
c_EXROV1_to_EXROV2[] <- rbinom(EXROV1a[i],p_V1_to_V2[i])
EXROV1b[] <- EXROV1a[i] - c_EXROV1_to_EXROV2[i]

# EXROV1 -> EXROS
c_EXROV1_to_EXROS[] <- rbinom(EXROV1b[i],p_V1_to_S[i])
EXROV1c[] <- EXROV1b[i] - c_EXROV1_to_EXROS[i]

# Death
c_death_EXROV1[] <- rbinom(EXROV1c[i],p_death[i])
EXROV1d[] <- EXROV1c[i] - c_death_EXROV1[i]

# Aging
c_aging_EXROV1[] <- floor(EXROV1d[i] * u_aging)

# ----- 1.4.3. EXROV2 Transitions -----

# EXROV2 -> IXROV2
c_EXROV2_to_IXROV2[] <- rbinom(EXROV2[i],p_incub_X[i])
EXROV2a[] <- EXROV2[i] - c_EXROV2_to_IXROV2[i]

# EXROV2 -> EXROS
c_EXROV2_to_EXROS[] <- rbinom(EXROV2a[i],p_V2_to_S[i])
EXROV2b[] <- EXROV2a[i] - c_EXROV2_to_EXROS[i]

# Death
c_death_EXROV2[] <- rbinom(EXROV2b[i],p_death[i])
EXROV2c[] <- EXROV2b[i] - c_death_EXROV2[i]

# Aging
c_aging_EXROV2[] <- floor(EXROV2c[i] * u_aging)

# ----- 1.4.11. IXROS Transitions -----

# IXROS -> ROXS
c_IXROS_to_ROXS[] <- rbinom(IXROS[i],p_recov_X[i])
IXROSa[] <- IXROS[i] - c_IXROS_to_ROXS[i]

# Infection Deaths
c_death_inf_IXROS[] <- rbinom(IXROSa[i],p_death_IXROS[i])
IXROSb[] <- IXROSa[i] - c_death_inf_IXROS[i]

# IXROS -> IXROV1
c_IXROS_to_IXROV1[] <- rbinom(IXROSb[i],p_S_to_V1[i])
IXROSc[] <- IXROSb[i] - c_IXROS_to_IXROV1[i]

# Death
c_death_IXROS[] <- rbinom(IXROSc[i],p_death[i])
IXROSd[] <- IXROSc[i] - c_death_IXROS[i]

# Aging
c_aging_IXROS[] <- floor(IXROSd[i] * u_aging)

# ----- 1.4.11. IXROV1 Transitions -----

# IXROV1 -> ROXV1
c_IXROV1_to_ROXV1[] <- rbinom(IXROV1[i],p_recov_X[i])
IXROV1a[] <- IXROV1[i] - c_IXROV1_to_ROXV1[i]

# Infection Deaths
c_death_inf_IXROV1[] <- rbinom(IXROV1a[i],p_death_IXROV1[i])
IXROV1b[] <- IXROV1a[i] - c_death_inf_IXROV1[i]

# IXROV1 -> IXROV2
c_IXROV1_to_IXROV2[] <- rbinom(IXROV1b[i],p_V1_to_V2[i])
IXROV1c[] <- IXROV1b[i] - c_IXROV1_to_IXROV2[i]

# IXROV1 -> IXROS
c_IXROV1_to_IXROS[] <- rbinom(IXROV1c[i],p_V1_to_S[i])
IXROV1d[] <- IXROV1c[i] - c_IXROV1_to_IXROS[i]

# Death
c_death_IXROV1[] <- rbinom(IXROV1d[i],p_death[i])
IXROV1e[] <- IXROV1d[i] - c_death_IXROV1[i]

# Aging
c_aging_IXROV1[] <- floor(IXROV1e[i] * u_aging)

# ----- 1.4.11. IXROV2 Transitions -----

# IXROV2 -> ROXV2
c_IXROV2_to_ROXV2[] <- rbinom(IXROV2[i],p_recov_X[i])
IXROV2a[] <- IXROV2[i] - c_IXROV2_to_ROXV2[i]

# Infection Deaths
c_death_inf_IXROV2[] <- rbinom(IXROV2a[i],p_death_IXROV2[i])
IXROV2b[] <- IXROV2a[i] - c_death_inf_IXROV2[i]

# IXROV2 -> IXROS
c_IXROV2_to_IXROS[] <- rbinom(IXROV2b[i],p_V2_to_S[i])
IXROV2c[] <- IXROV2b[i] - c_IXROV2_to_IXROS[i]

# Death
c_death_IXROV2[] <- rbinom(IXROV2c[i],p_death[i])
IXROV2d[] <- IXROV2c[i] - c_death_IXROV2[i]

# Aging
c_aging_IXROV2[] <- floor(IXROV2d[i] * u_aging)

# ----- 1.4.19. ROXS Transitions -----

# ROXS -> ROXV1
c_ROXS_to_ROXV1[] <- rbinom(ROXS[i],p_S_to_V1[i])
ROXSa[] <- ROXS[i] - c_ROXS_to_ROXV1[i]

# ROXS -> S
c_ROXS_to_S[] <- rbinom(ROXSa[i],p_waning_X[i])
ROXSb[] <- ROXSa[i] - c_ROXS_to_S[i]

# Death
c_death_ROXS[] <- rbinom(ROXSb[i],p_death[i])
ROXSc[] <- ROXSb[i] - c_death_ROXS[i]

# Aging
c_aging_ROXS[] <- floor(ROXSc[i] * u_aging)

# ----- 1.4.19. ROXV1 Transitions -----

# ROXV1 -> ROXV2
c_ROXV1_to_ROXV2[] <- rbinom(ROXV1[i],p_V1_to_V2[i])
ROXV1a[] <- ROXV1[i] - c_ROXV1_to_ROXV2[i]

# ROXV1 -> V1
c_ROXV1_to_V1[] <- rbinom(ROXV1a[i],p_waning_X[i])
ROXV1b[] <- ROXV1a[i] - c_ROXV1_to_V1[i]

# ROXV1 -> ROXS
c_ROXV1_to_ROXS[] <- rbinom(ROXV1b[i],p_V1_to_S[i])
ROXV1c[] <- ROXV1b[i] - c_ROXV1_to_ROXS[i]

# Death
c_death_ROXV1[] <- rbinom(ROXV1c[i],p_death[i])
ROXV1d[] <- ROXV1c[i] - c_death_ROXV1[i]

# Aging
c_aging_ROXV1[] <- floor(ROXV1d[i] * u_aging)

# ----- 1.4.19. ROXV2 Transitions -----

# ROXV2 -> V2
c_ROXV2_to_V2[] <- rbinom(ROXV2[i],p_waning_X[i])
ROXV2a[] <- ROXV2[i] - c_ROXV2_to_V2[i]

# ROXV2 -> ROXS
c_ROXV2_to_ROXS[] <- rbinom(ROXV2a[i],p_V2_to_S[i])
ROXV2b[] <- ROXV2a[i] - c_ROXV2_to_ROXS[i]

# Death
c_death_ROXV2[] <- rbinom(ROXV2b[i],p_death[i])
ROXV2c[] <- ROXV2b[i] - c_death_ROXV2[i]

# Aging
c_aging_ROXV2[] <- floor(ROXV2c[i] * u_aging)

# ----- 1.5. Births and Migrations ---------------------------------------------

c_birth[] <- rbinom(N[i],p_birth[i])

c_migration_S[]  <- rbinom(n_migration[i],p_migration_S[i])
c_migration_V1[] <- rbinom(n_migration[i],p_migration_V1[i])
c_migration_V2[] <- rbinom(n_migration[i],p_migration_V2[i])

c_migration_EWS[]  <- rbinom(n_migration[i],p_migration_EWS[i])
c_migration_IWS[]  <- rbinom(n_migration[i],p_migration_IWS[i])
c_migration_RWS[]  <- rbinom(n_migration[i],p_migration_RWS[i])
c_migration_EWV1[] <- rbinom(n_migration[i],p_migration_EWV1[i])
c_migration_IWV1[] <- rbinom(n_migration[i],p_migration_IWV1[i])
c_migration_RWV1[] <- rbinom(n_migration[i],p_migration_RWV1[i])
c_migration_EWV2[] <- rbinom(n_migration[i],p_migration_EWV2[i])
c_migration_IWV2[] <- rbinom(n_migration[i],p_migration_IWV2[i])
c_migration_RWV2[] <- rbinom(n_migration[i],p_migration_RWV2[i])

c_migration_EAS[]     <- rbinom(n_migration[i],p_migration_EAS[i])
c_migration_IAS[]     <- rbinom(n_migration[i],p_migration_IAS[i])
c_migration_RAS[]     <- rbinom(n_migration[i],p_migration_RAS[i])
c_migration_EARWS[]  <- rbinom(n_migration[i],p_migration_EARWS[i])
c_migration_IARWS[]  <- rbinom(n_migration[i],p_migration_IARWS[i])
c_migration_EAV1[]    <- rbinom(n_migration[i],p_migration_EAV1[i])
c_migration_IAV1[]    <- rbinom(n_migration[i],p_migration_IAV1[i])
c_migration_RAV1[]    <- rbinom(n_migration[i],p_migration_RAV1[i])
c_migration_EARWV1[] <- rbinom(n_migration[i],p_migration_EARWV1[i])
c_migration_IARWV1[] <- rbinom(n_migration[i],p_migration_IARWV1[i])
c_migration_EAV2[]    <- rbinom(n_migration[i],p_migration_EAV2[i])
c_migration_IAV2[]    <- rbinom(n_migration[i],p_migration_IAV2[i])
c_migration_RAV2[]    <- rbinom(n_migration[i],p_migration_RAV2[i])
c_migration_EARWV2[] <- rbinom(n_migration[i],p_migration_EARWV2[i])
c_migration_IARWV2[] <- rbinom(n_migration[i],p_migration_IARWV2[i])

c_migration_EDS[]     <- rbinom(n_migration[i],p_migration_EDS[i])
c_migration_IDS[]     <- rbinom(n_migration[i],p_migration_IDS[i])
c_migration_RDS[]     <- rbinom(n_migration[i],p_migration_RDS[i])
c_migration_EDRWS[]  <- rbinom(n_migration[i],p_migration_EDRWS[i])
c_migration_IDRWS[]  <- rbinom(n_migration[i],p_migration_IDRWS[i])
c_migration_EDRAS[]  <- rbinom(n_migration[i],p_migration_EDRAS[i])
c_migration_IDRAS[]  <- rbinom(n_migration[i],p_migration_IDRAS[i])
c_migration_EDV1[]    <- rbinom(n_migration[i],p_migration_EDV1[i])
c_migration_IDV1[]    <- rbinom(n_migration[i],p_migration_IDV1[i])
c_migration_RDV1[]    <- rbinom(n_migration[i],p_migration_RDV1[i])
c_migration_EDRWV1[] <- rbinom(n_migration[i],p_migration_EDRWV1[i])
c_migration_IDRWV1[] <- rbinom(n_migration[i],p_migration_IDRWV1[i])
c_migration_EDRAV1[] <- rbinom(n_migration[i],p_migration_EDRAV1[i])
c_migration_IDRAV1[] <- rbinom(n_migration[i],p_migration_IDRAV1[i])
c_migration_EDV2[]    <- rbinom(n_migration[i],p_migration_EDV2[i])
c_migration_IDV2[]    <- rbinom(n_migration[i],p_migration_IDV2[i])
c_migration_RDV2[]    <- rbinom(n_migration[i],p_migration_RDV2[i])
c_migration_EDRWV2[] <- rbinom(n_migration[i],p_migration_EDRWV2[i])
c_migration_IDRWV2[] <- rbinom(n_migration[i],p_migration_IDRWV2[i])
c_migration_EDRAV2[] <- rbinom(n_migration[i],p_migration_EDRAV2[i])
c_migration_IDRAV2[] <- rbinom(n_migration[i],p_migration_IDRAV2[i])

c_migration_EOS[]     <- rbinom(n_migration[i],p_migration_EOS[i])
c_migration_IOS[]     <- rbinom(n_migration[i],p_migration_IOS[i])
c_migration_ROS[]     <- rbinom(n_migration[i],p_migration_ROS[i])
c_migration_EORWS[]  <- rbinom(n_migration[i],p_migration_EORWS[i])
c_migration_IORWS[]  <- rbinom(n_migration[i],p_migration_IORWS[i])
c_migration_EORAS[]  <- rbinom(n_migration[i],p_migration_EORAS[i])
c_migration_IORAS[]  <- rbinom(n_migration[i],p_migration_IORAS[i])
c_migration_EORDS[]  <- rbinom(n_migration[i],p_migration_EORDS[i])
c_migration_IORDS[]  <- rbinom(n_migration[i],p_migration_IORDS[i])
c_migration_EOV1[]    <- rbinom(n_migration[i],p_migration_EOV1[i])
c_migration_IOV1[]    <- rbinom(n_migration[i],p_migration_IOV1[i])
c_migration_ROV1[]    <- rbinom(n_migration[i],p_migration_ROV1[i])
c_migration_EORWV1[] <- rbinom(n_migration[i],p_migration_EORWV1[i])
c_migration_IORWV1[] <- rbinom(n_migration[i],p_migration_IORWV1[i])
c_migration_EORAV1[] <- rbinom(n_migration[i],p_migration_EORAV1[i])
c_migration_IORAV1[] <- rbinom(n_migration[i],p_migration_IORAV1[i])
c_migration_EORDV1[] <- rbinom(n_migration[i],p_migration_EORDV1[i])
c_migration_IORDV1[] <- rbinom(n_migration[i],p_migration_IORDV1[i])
c_migration_EOV2[]    <- rbinom(n_migration[i],p_migration_EOV2[i])
c_migration_IOV2[]    <- rbinom(n_migration[i],p_migration_IOV2[i])
c_migration_ROV2[]    <- rbinom(n_migration[i],p_migration_ROV2[i])
c_migration_EORWV2[] <- rbinom(n_migration[i],p_migration_EORWV2[i])
c_migration_IORWV2[] <- rbinom(n_migration[i],p_migration_IORWV2[i])
c_migration_EORAV2[] <- rbinom(n_migration[i],p_migration_EORAV2[i])
c_migration_IORAV2[] <- rbinom(n_migration[i],p_migration_IORAV2[i])
c_migration_EORDV2[] <- rbinom(n_migration[i],p_migration_EORDV2[i])
c_migration_IORDV2[] <- rbinom(n_migration[i],p_migration_IORDV2[i])

c_migration_EXS[]  <- rbinom(n_migration[i],p_migration_EXS[i])
c_migration_IXS[]  <- rbinom(n_migration[i],p_migration_IXS[i])
c_migration_RXS[]  <- rbinom(n_migration[i],p_migration_RXS[i])
c_migration_EXV1[] <- rbinom(n_migration[i],p_migration_EXV1[i])
c_migration_IXV1[] <- rbinom(n_migration[i],p_migration_IXV1[i])
c_migration_RXV1[] <- rbinom(n_migration[i],p_migration_RXV1[i])
c_migration_EXV2[] <- rbinom(n_migration[i],p_migration_EXV2[i])
c_migration_IXV2[] <- rbinom(n_migration[i],p_migration_IXV2[i])
c_migration_RXV2[] <- rbinom(n_migration[i],p_migration_RXV2[i])

c_migration_EORXS[]  <- rbinom(n_migration[i],p_migration_EORXS[i])
c_migration_IORXS[]  <- rbinom(n_migration[i],p_migration_IORXS[i])
c_migration_EORXV1[] <- rbinom(n_migration[i],p_migration_EORXV1[i])
c_migration_IORXV1[] <- rbinom(n_migration[i],p_migration_IORXV1[i])
c_migration_EORXV2[] <- rbinom(n_migration[i],p_migration_EORXV2[i])
c_migration_IORXV2[] <- rbinom(n_migration[i],p_migration_IORXV2[i])

c_migration_EXROS[]  <- rbinom(n_migration[i],p_migration_EXROS[i])
c_migration_IXROS[]  <- rbinom(n_migration[i],p_migration_IXROS[i])
c_migration_EXROV1[] <- rbinom(n_migration[i],p_migration_EXROV1[i])
c_migration_IXROV1[] <- rbinom(n_migration[i],p_migration_IXROV1[i])
c_migration_EXROV2[] <- rbinom(n_migration[i],p_migration_EXROV2[i])
c_migration_IXROV2[] <- rbinom(n_migration[i],p_migration_IXROV2[i])

c_migration_ROXS[]   <- rbinom(n_migration[i],p_migration_ROXS[i])
c_migration_ROXV1[]  <- rbinom(n_migration[i],p_migration_ROXV1[i])
c_migration_ROXV2[]  <- rbinom(n_migration[i],p_migration_ROXV2[i])

# ----- 1.6 Exposures ----------------------------------------------------------

c_S_exp_W[] <- amount_intro_W[i] * (step == time_intro_W)
c_S_exp_X[] <- amount_intro_X[i] * (step == time_intro_X)

# ------------------------------------------------------------------------------
# ----- 2. Compartment Updates -------------------------------------------------
# ------------------------------------------------------------------------------

# ----- 2.1. Process Trackers --------------------------------------------------

update(deaths_IX[]) <- deaths_IX[i] +
  c_death_inf_IXS[i] + c_death_inf_IXV1[i] + c_death_inf_IXV2[i] +
  c_death_inf_IXROS[i] + c_death_inf_IXROV1[i] + c_death_inf_IXROV2[i]

update(days_IXS[]) <- days_IXS[i] + IXS[i]
update(days_IXV1[]) <- days_IXV1[i] + IXV1[i]
update(days_IXV2[]) <- days_IXV2[i] + IXV2[i]
update(days_IXROS[]) <- days_IXROS[i] + IXROS[i]
update(days_IXROV1[]) <- days_IXROV1[i] + IXROV1[i]
update(days_IXROV2[]) <- days_IXROV2[i] + IXROV2[i]

# ----- 2.2. Composite Compartments --------------------------------------------

update(N[]) <- 
  S[i] + V1[i] + V2[i] +
  EWS[i] + IWS[i] + RWS[i] +
  EWV1[i] + IWV1[i] + RWV1[i] +
  EWV2[i] + IWV2[i] + RWV2[i] +
  EAS[i] + IAS[i] + RAS[i] +
  EARWS[i] + IARWS[i] +
  EAV1[i] + IAV1[i] + RAV1[i] +
  EARWV1[i] + IARWV1[i] +
  EAV2[i] + IAV2[i] + RAV2[i] +
  EARWV2[i] + IARWV2[i] +
  EDS[i] + IDS[i] + RDS[i] +
  EDRWS[i] + IDRWS[i] +
  EDRAS[i] + IDRAS[i] +
  EDV1[i] + IDV1[i] + RDV1[i] +
  EDRWV1[i] + IDRWV1[i] +
  EDRAV1[i] + IDRAV1[i] +
  EDV2[i] + IDV2[i] + RDV2[i] +
  EDRWV2[i] + IDRWV2[i] +
  EDRAV2[i] + IDRAV2[i] +
  EOS[i] + IOS[i] + ROS[i] +
  EORWS[i] + IORWS[i] +
  EORAS[i] + IORAS[i] +
  EORDS[i] + IORDS[i] +
  EOV1[i] + IOV1[i] + ROV1[i] +
  EORWV1[i] + IORWV1[i] +
  EORAV1[i] + IORAV1[i] +
  EORDV1[i] + IORDV1[i] +
  EOV2[i] + IOV2[i] + ROV2[i] +
  EORWV2[i] + IORWV2[i] +
  EORAV2[i] + IORAV2[i] +
  EORDV2[i] + IORDV2[i] +
  EXS[i] + IXS[i] + RXS[i] +
  EXV1[i] + IXV1[i] + RXV1[i] +
  EXV2[i] + IXV2[i] + RXV2[i] +
  EORXS[i] + IORXS[i] +
  EORXV1[i] + IORXV1[i] +
  EORXV2[i] + IORXV2[i] +
  EXROS[i] + IXROS[i] +
  EXROV1[i] + IXROV1[i] +
  EXROV2[i] + IXROV2[i] +
  ROXS[i] + ROXV1[i] + ROXV2[i]
  

# ----- 2.3. S Compartment ----------------------------------------------------

c_S_in[] <- c_RAS_to_S[i] + c_RDS_to_S[i] + c_ROS_to_S[i] + c_RWS_to_S[i] + c_ROXS_to_S[i] + c_RXS_to_S[i] + c_V1_to_S[i] + c_V2_to_S[i] + c_migration_S[i]

c_S_out[] <- c_S_to_EAS[i] + c_S_to_EDS[i] + c_S_to_EOS[i] + c_S_to_EWS[i] + c_S_to_EXS[i] + c_S_to_V1[i] + c_death_S[i] + c_S_exp_X[i] + c_S_exp_W[i]

update(S[1]) <- S[1] + c_S_in[1] - c_S_out[1] + sum(c_birth) - c_aging_S[1]
update(S[2]) <- S[2] + c_S_in[2] - c_S_out[2] + c_aging_S[1] - c_aging_S[2]
update(S[3]) <- S[3] + c_S_in[3] - c_S_out[3] + c_aging_S[2] - c_aging_S[3]
update(S[4]) <- S[4] + c_S_in[4] - c_S_out[4] + c_aging_S[3] - c_aging_S[4]
update(S[5]) <- S[5] + c_S_in[5] - c_S_out[5] + c_aging_S[4] - c_aging_S[5]
update(S[6]) <- S[6] + c_S_in[6] - c_S_out[6] + c_aging_S[5] - c_aging_S[6]
update(S[7]) <- S[7] + c_S_in[7] - c_S_out[7] + c_aging_S[6] - c_aging_S[7]
update(S[8]) <- S[8] + c_S_in[8] - c_S_out[8] + c_aging_S[7] - c_aging_S[8]
update(S[9]) <- S[9] + c_S_in[9] - c_S_out[9] + c_aging_S[8] - c_aging_S[9]
update(S[10]) <- S[10] + c_S_in[10] - c_S_out[10] + c_aging_S[9] - c_aging_S[10]
update(S[11]) <- S[11] + c_S_in[11] - c_S_out[11] + c_aging_S[10] - c_aging_S[11]
update(S[12]) <- S[12] + c_S_in[12] - c_S_out[12] + c_aging_S[11] - c_aging_S[12]
update(S[13]) <- S[13] + c_S_in[13] - c_S_out[13] + c_aging_S[12] - c_aging_S[13]
update(S[14]) <- S[14] + c_S_in[14] - c_S_out[14] + c_aging_S[13] - c_aging_S[14]
update(S[15]) <- S[15] + c_S_in[15] - c_S_out[15] + c_aging_S[14] - c_aging_S[15]
update(S[16]) <- S[16] + c_S_in[16] - c_S_out[16] + c_aging_S[15]

# ----- 2.4. V1 Compartment -----------------------------------------------------

c_V1_in[] <- c_RAV1_to_V1[i] + c_RDV1_to_V1[i] + c_ROV1_to_V1[i] + c_RWV1_to_V1[i] + c_ROXV1_to_V1[i] + c_RXV1_to_V1[i] + c_S_to_V1[i] + c_migration_V1[i]
c_V1_out[] <- c_V1_to_EAV1[i] + c_V1_to_EDV1[i] + c_V1_to_EOV1[i] + c_V1_to_EWV1[i] + c_V1_to_EXV1[i] + c_V1_to_S[i] + c_V1_to_V2[i] + c_death_V1[i]

update(V1[1]) <- V1[1] + c_V1_in[1] - c_V1_out[1] - c_aging_V1[1]
update(V1[2]) <- V1[2] + c_V1_in[2] - c_V1_out[2] + c_aging_V1[1] - c_aging_V1[2]
update(V1[3]) <- V1[3] + c_V1_in[3] - c_V1_out[3] + c_aging_V1[2] - c_aging_V1[3]
update(V1[4]) <- V1[4] + c_V1_in[4] - c_V1_out[4] + c_aging_V1[3] - c_aging_V1[4]
update(V1[5]) <- V1[5] + c_V1_in[5] - c_V1_out[5] + c_aging_V1[4] - c_aging_V1[5]
update(V1[6]) <- V1[6] + c_V1_in[6] - c_V1_out[6] + c_aging_V1[5] - c_aging_V1[6]
update(V1[7]) <- V1[7] + c_V1_in[7] - c_V1_out[7] + c_aging_V1[6] - c_aging_V1[7]
update(V1[8]) <- V1[8] + c_V1_in[8] - c_V1_out[8] + c_aging_V1[7] - c_aging_V1[8]
update(V1[9]) <- V1[9] + c_V1_in[9] - c_V1_out[9] + c_aging_V1[8] - c_aging_V1[9]
update(V1[10]) <- V1[10] + c_V1_in[10] - c_V1_out[10] + c_aging_V1[9] - c_aging_V1[10]
update(V1[11]) <- V1[11] + c_V1_in[11] - c_V1_out[11] + c_aging_V1[10] - c_aging_V1[11]
update(V1[12]) <- V1[12] + c_V1_in[12] - c_V1_out[12] + c_aging_V1[11] - c_aging_V1[12]
update(V1[13]) <- V1[13] + c_V1_in[13] - c_V1_out[13] + c_aging_V1[12] - c_aging_V1[13]
update(V1[14]) <- V1[14] + c_V1_in[14] - c_V1_out[14] + c_aging_V1[13] - c_aging_V1[14]
update(V1[15]) <- V1[15] + c_V1_in[15] - c_V1_out[15] + c_aging_V1[14] - c_aging_V1[15]
update(V1[16]) <- V1[16] + c_V1_in[16] - c_V1_out[16] + c_aging_V1[15]

# ----- 2.5. V2 Compartment -----------------------------------------------------

c_V2_in[] <- c_RAV2_to_V2[i] + c_RDV2_to_V2[i] + c_ROV2_to_V2[i] + c_RWV2_to_V2[i] + c_ROXV2_to_V2[i] + c_RXV2_to_V2[i] + c_V1_to_V2[i] + c_migration_V2[i]
c_V2_out[] <- c_V2_to_EAV2[i] + c_V2_to_EDV2[i] + c_V2_to_EOV2[i] + c_V2_to_EWV2[i] + c_V2_to_EXV2[i] + c_V2_to_S[i] + c_death_V2[i]

update(V2[1]) <- V2[1] + c_V2_in[1] - c_V2_out[1] - c_aging_V2[1]
update(V2[2]) <- V2[2] + c_V2_in[2] - c_V2_out[2] + c_aging_V2[1] - c_aging_V2[2]
update(V2[3]) <- V2[3] + c_V2_in[3] - c_V2_out[3] + c_aging_V2[2] - c_aging_V2[3]
update(V2[4]) <- V2[4] + c_V2_in[4] - c_V2_out[4] + c_aging_V2[3] - c_aging_V2[4]
update(V2[5]) <- V2[5] + c_V2_in[5] - c_V2_out[5] + c_aging_V2[4] - c_aging_V2[5]
update(V2[6]) <- V2[6] + c_V2_in[6] - c_V2_out[6] + c_aging_V2[5] - c_aging_V2[6]
update(V2[7]) <- V2[7] + c_V2_in[7] - c_V2_out[7] + c_aging_V2[6] - c_aging_V2[7]
update(V2[8]) <- V2[8] + c_V2_in[8] - c_V2_out[8] + c_aging_V2[7] - c_aging_V2[8]
update(V2[9]) <- V2[9] + c_V2_in[9] - c_V2_out[9] + c_aging_V2[8] - c_aging_V2[9]
update(V2[10]) <- V2[10] + c_V2_in[10] - c_V2_out[10] + c_aging_V2[9] - c_aging_V2[10]
update(V2[11]) <- V2[11] + c_V2_in[11] - c_V2_out[11] + c_aging_V2[10] - c_aging_V2[11]
update(V2[12]) <- V2[12] + c_V2_in[12] - c_V2_out[12] + c_aging_V2[11] - c_aging_V2[12]
update(V2[13]) <- V2[13] + c_V2_in[13] - c_V2_out[13] + c_aging_V2[12] - c_aging_V2[13]
update(V2[14]) <- V2[14] + c_V2_in[14] - c_V2_out[14] + c_aging_V2[13] - c_aging_V2[14]
update(V2[15]) <- V2[15] + c_V2_in[15] - c_V2_out[15] + c_aging_V2[14] - c_aging_V2[15]
update(V2[16]) <- V2[16] + c_V2_in[16] - c_V2_out[16] + c_aging_V2[15]

# ----- 2.7. EWS Compartment --------------------------------------------------

c_EWS_in[] <- c_EWV1_to_EWS[i] + c_EWV2_to_EWS[i] + c_S_to_EWS[i] + c_migration_EWS[i]  + c_S_exp_W[i]
c_EWS_out[] <- c_EWS_to_EWV1[i] + c_EWS_to_IWS[i] + c_death_EWS[i]

update(EWS[1]) <- EWS[1] + c_EWS_in[1] - c_EWS_out[1] - c_aging_EWS[1]
update(EWS[2]) <- EWS[2] + c_EWS_in[2] - c_EWS_out[2] + c_aging_EWS[1] - c_aging_EWS[2]
update(EWS[3]) <- EWS[3] + c_EWS_in[3] - c_EWS_out[3] + c_aging_EWS[2] - c_aging_EWS[3]
update(EWS[4]) <- EWS[4] + c_EWS_in[4] - c_EWS_out[4] + c_aging_EWS[3] - c_aging_EWS[4]
update(EWS[5]) <- EWS[5] + c_EWS_in[5] - c_EWS_out[5] + c_aging_EWS[4] - c_aging_EWS[5]
update(EWS[6]) <- EWS[6] + c_EWS_in[6] - c_EWS_out[6] + c_aging_EWS[5] - c_aging_EWS[6]
update(EWS[7]) <- EWS[7] + c_EWS_in[7] - c_EWS_out[7] + c_aging_EWS[6] - c_aging_EWS[7]
update(EWS[8]) <- EWS[8] + c_EWS_in[8] - c_EWS_out[8] + c_aging_EWS[7] - c_aging_EWS[8]
update(EWS[9]) <- EWS[9] + c_EWS_in[9] - c_EWS_out[9] + c_aging_EWS[8] - c_aging_EWS[9]
update(EWS[10]) <- EWS[10] + c_EWS_in[10] - c_EWS_out[10] + c_aging_EWS[9] - c_aging_EWS[10]
update(EWS[11]) <- EWS[11] + c_EWS_in[11] - c_EWS_out[11] + c_aging_EWS[10] - c_aging_EWS[11]
update(EWS[12]) <- EWS[12] + c_EWS_in[12] - c_EWS_out[12] + c_aging_EWS[11] - c_aging_EWS[12]
update(EWS[13]) <- EWS[13] + c_EWS_in[13] - c_EWS_out[13] + c_aging_EWS[12] - c_aging_EWS[13]
update(EWS[14]) <- EWS[14] + c_EWS_in[14] - c_EWS_out[14] + c_aging_EWS[13] - c_aging_EWS[14]
update(EWS[15]) <- EWS[15] + c_EWS_in[15] - c_EWS_out[15] + c_aging_EWS[14] - c_aging_EWS[15]
update(EWS[16]) <- EWS[16] + c_EWS_in[16] - c_EWS_out[16] + c_aging_EWS[15]

# ----- 2.8. EWV1 Compartment --------------------------------------------------

c_EWV1_in[] <- c_EWS_to_EWV1[i] + c_V1_to_EWV1[i] + c_migration_EWV1[i]
c_EWV1_out[] <- c_EWV1_to_EWS[i] + c_EWV1_to_EWV2[i] + c_EWV1_to_IWV1[i] + c_death_EWV1[i]

update(EWV1[1]) <- EWV1[1] + c_EWV1_in[1] - c_EWV1_out[1] - c_aging_EWV1[1]
update(EWV1[2]) <- EWV1[2] + c_EWV1_in[2] - c_EWV1_out[2] + c_aging_EWV1[1] - c_aging_EWV1[2]
update(EWV1[3]) <- EWV1[3] + c_EWV1_in[3] - c_EWV1_out[3] + c_aging_EWV1[2] - c_aging_EWV1[3]
update(EWV1[4]) <- EWV1[4] + c_EWV1_in[4] - c_EWV1_out[4] + c_aging_EWV1[3] - c_aging_EWV1[4]
update(EWV1[5]) <- EWV1[5] + c_EWV1_in[5] - c_EWV1_out[5] + c_aging_EWV1[4] - c_aging_EWV1[5]
update(EWV1[6]) <- EWV1[6] + c_EWV1_in[6] - c_EWV1_out[6] + c_aging_EWV1[5] - c_aging_EWV1[6]
update(EWV1[7]) <- EWV1[7] + c_EWV1_in[7] - c_EWV1_out[7] + c_aging_EWV1[6] - c_aging_EWV1[7]
update(EWV1[8]) <- EWV1[8] + c_EWV1_in[8] - c_EWV1_out[8] + c_aging_EWV1[7] - c_aging_EWV1[8]
update(EWV1[9]) <- EWV1[9] + c_EWV1_in[9] - c_EWV1_out[9] + c_aging_EWV1[8] - c_aging_EWV1[9]
update(EWV1[10]) <- EWV1[10] + c_EWV1_in[10] - c_EWV1_out[10] + c_aging_EWV1[9] - c_aging_EWV1[10]
update(EWV1[11]) <- EWV1[11] + c_EWV1_in[11] - c_EWV1_out[11] + c_aging_EWV1[10] - c_aging_EWV1[11]
update(EWV1[12]) <- EWV1[12] + c_EWV1_in[12] - c_EWV1_out[12] + c_aging_EWV1[11] - c_aging_EWV1[12]
update(EWV1[13]) <- EWV1[13] + c_EWV1_in[13] - c_EWV1_out[13] + c_aging_EWV1[12] - c_aging_EWV1[13]
update(EWV1[14]) <- EWV1[14] + c_EWV1_in[14] - c_EWV1_out[14] + c_aging_EWV1[13] - c_aging_EWV1[14]
update(EWV1[15]) <- EWV1[15] + c_EWV1_in[15] - c_EWV1_out[15] + c_aging_EWV1[14] - c_aging_EWV1[15]
update(EWV1[16]) <- EWV1[16] + c_EWV1_in[16] - c_EWV1_out[16] + c_aging_EWV1[15]

# ----- 2.9. EWV2 Compartment --------------------------------------------------

c_EWV2_in[] <- c_EWV1_to_EWV2[i] + c_V2_to_EWV2[i] + c_migration_EWV2[i]
c_EWV2_out[] <- c_EWV2_to_EWS[i] + c_EWV2_to_IWV2[i] + c_death_EWV2[i]

update(EWV2[1]) <- EWV2[1] + c_EWV2_in[1] - c_EWV2_out[1] - c_aging_EWV2[1]
update(EWV2[2]) <- EWV2[2] + c_EWV2_in[2] - c_EWV2_out[2] + c_aging_EWV2[1] - c_aging_EWV2[2]
update(EWV2[3]) <- EWV2[3] + c_EWV2_in[3] - c_EWV2_out[3] + c_aging_EWV2[2] - c_aging_EWV2[3]
update(EWV2[4]) <- EWV2[4] + c_EWV2_in[4] - c_EWV2_out[4] + c_aging_EWV2[3] - c_aging_EWV2[4]
update(EWV2[5]) <- EWV2[5] + c_EWV2_in[5] - c_EWV2_out[5] + c_aging_EWV2[4] - c_aging_EWV2[5]
update(EWV2[6]) <- EWV2[6] + c_EWV2_in[6] - c_EWV2_out[6] + c_aging_EWV2[5] - c_aging_EWV2[6]
update(EWV2[7]) <- EWV2[7] + c_EWV2_in[7] - c_EWV2_out[7] + c_aging_EWV2[6] - c_aging_EWV2[7]
update(EWV2[8]) <- EWV2[8] + c_EWV2_in[8] - c_EWV2_out[8] + c_aging_EWV2[7] - c_aging_EWV2[8]
update(EWV2[9]) <- EWV2[9] + c_EWV2_in[9] - c_EWV2_out[9] + c_aging_EWV2[8] - c_aging_EWV2[9]
update(EWV2[10]) <- EWV2[10] + c_EWV2_in[10] - c_EWV2_out[10] + c_aging_EWV2[9] - c_aging_EWV2[10]
update(EWV2[11]) <- EWV2[11] + c_EWV2_in[11] - c_EWV2_out[11] + c_aging_EWV2[10] - c_aging_EWV2[11]
update(EWV2[12]) <- EWV2[12] + c_EWV2_in[12] - c_EWV2_out[12] + c_aging_EWV2[11] - c_aging_EWV2[12]
update(EWV2[13]) <- EWV2[13] + c_EWV2_in[13] - c_EWV2_out[13] + c_aging_EWV2[12] - c_aging_EWV2[13]
update(EWV2[14]) <- EWV2[14] + c_EWV2_in[14] - c_EWV2_out[14] + c_aging_EWV2[13] - c_aging_EWV2[14]
update(EWV2[15]) <- EWV2[15] + c_EWV2_in[15] - c_EWV2_out[15] + c_aging_EWV2[14] - c_aging_EWV2[15]
update(EWV2[16]) <- EWV2[16] + c_EWV2_in[16] - c_EWV2_out[16] + c_aging_EWV2[15]

# ----- 2.11. EAS Compartment --------------------------------------------------

c_EAS_in[] <- c_EAV1_to_EAS[i] + c_EAV2_to_EAS[i] + c_S_to_EAS[i] + c_migration_EAS[i]
c_EAS_out[] <- c_EAS_to_EAV1[i] + c_EAS_to_IAS[i] + c_death_EAS[i]

update(EAS[1]) <- EAS[1] + c_EAS_in[1] - c_EAS_out[1] - c_aging_EAS[1]
update(EAS[2]) <- EAS[2] + c_EAS_in[2] - c_EAS_out[2] + c_aging_EAS[1] - c_aging_EAS[2]
update(EAS[3]) <- EAS[3] + c_EAS_in[3] - c_EAS_out[3] + c_aging_EAS[2] - c_aging_EAS[3]
update(EAS[4]) <- EAS[4] + c_EAS_in[4] - c_EAS_out[4] + c_aging_EAS[3] - c_aging_EAS[4]
update(EAS[5]) <- EAS[5] + c_EAS_in[5] - c_EAS_out[5] + c_aging_EAS[4] - c_aging_EAS[5]
update(EAS[6]) <- EAS[6] + c_EAS_in[6] - c_EAS_out[6] + c_aging_EAS[5] - c_aging_EAS[6]
update(EAS[7]) <- EAS[7] + c_EAS_in[7] - c_EAS_out[7] + c_aging_EAS[6] - c_aging_EAS[7]
update(EAS[8]) <- EAS[8] + c_EAS_in[8] - c_EAS_out[8] + c_aging_EAS[7] - c_aging_EAS[8]
update(EAS[9]) <- EAS[9] + c_EAS_in[9] - c_EAS_out[9] + c_aging_EAS[8] - c_aging_EAS[9]
update(EAS[10]) <- EAS[10] + c_EAS_in[10] - c_EAS_out[10] + c_aging_EAS[9] - c_aging_EAS[10]
update(EAS[11]) <- EAS[11] + c_EAS_in[11] - c_EAS_out[11] + c_aging_EAS[10] - c_aging_EAS[11]
update(EAS[12]) <- EAS[12] + c_EAS_in[12] - c_EAS_out[12] + c_aging_EAS[11] - c_aging_EAS[12]
update(EAS[13]) <- EAS[13] + c_EAS_in[13] - c_EAS_out[13] + c_aging_EAS[12] - c_aging_EAS[13]
update(EAS[14]) <- EAS[14] + c_EAS_in[14] - c_EAS_out[14] + c_aging_EAS[13] - c_aging_EAS[14]
update(EAS[15]) <- EAS[15] + c_EAS_in[15] - c_EAS_out[15] + c_aging_EAS[14] - c_aging_EAS[15]
update(EAS[16]) <- EAS[16] + c_EAS_in[16] - c_EAS_out[16] + c_aging_EAS[15]

# ----- 2.12. EAV1 Compartment --------------------------------------------------

c_EAV1_in[] <- c_EAS_to_EAV1[i] + c_V1_to_EAV1[i] + c_migration_EAV1[i]
c_EAV1_out[] <- c_EAV1_to_EAS[i] + c_EAV1_to_EAV2[i] + c_EAV1_to_IAV1[i] + c_death_EAV1[i]

update(EAV1[1]) <- EAV1[1] + c_EAV1_in[1] - c_EAV1_out[1] - c_aging_EAV1[1]
update(EAV1[2]) <- EAV1[2] + c_EAV1_in[2] - c_EAV1_out[2] + c_aging_EAV1[1] - c_aging_EAV1[2]
update(EAV1[3]) <- EAV1[3] + c_EAV1_in[3] - c_EAV1_out[3] + c_aging_EAV1[2] - c_aging_EAV1[3]
update(EAV1[4]) <- EAV1[4] + c_EAV1_in[4] - c_EAV1_out[4] + c_aging_EAV1[3] - c_aging_EAV1[4]
update(EAV1[5]) <- EAV1[5] + c_EAV1_in[5] - c_EAV1_out[5] + c_aging_EAV1[4] - c_aging_EAV1[5]
update(EAV1[6]) <- EAV1[6] + c_EAV1_in[6] - c_EAV1_out[6] + c_aging_EAV1[5] - c_aging_EAV1[6]
update(EAV1[7]) <- EAV1[7] + c_EAV1_in[7] - c_EAV1_out[7] + c_aging_EAV1[6] - c_aging_EAV1[7]
update(EAV1[8]) <- EAV1[8] + c_EAV1_in[8] - c_EAV1_out[8] + c_aging_EAV1[7] - c_aging_EAV1[8]
update(EAV1[9]) <- EAV1[9] + c_EAV1_in[9] - c_EAV1_out[9] + c_aging_EAV1[8] - c_aging_EAV1[9]
update(EAV1[10]) <- EAV1[10] + c_EAV1_in[10] - c_EAV1_out[10] + c_aging_EAV1[9] - c_aging_EAV1[10]
update(EAV1[11]) <- EAV1[11] + c_EAV1_in[11] - c_EAV1_out[11] + c_aging_EAV1[10] - c_aging_EAV1[11]
update(EAV1[12]) <- EAV1[12] + c_EAV1_in[12] - c_EAV1_out[12] + c_aging_EAV1[11] - c_aging_EAV1[12]
update(EAV1[13]) <- EAV1[13] + c_EAV1_in[13] - c_EAV1_out[13] + c_aging_EAV1[12] - c_aging_EAV1[13]
update(EAV1[14]) <- EAV1[14] + c_EAV1_in[14] - c_EAV1_out[14] + c_aging_EAV1[13] - c_aging_EAV1[14]
update(EAV1[15]) <- EAV1[15] + c_EAV1_in[15] - c_EAV1_out[15] + c_aging_EAV1[14] - c_aging_EAV1[15]
update(EAV1[16]) <- EAV1[16] + c_EAV1_in[16] - c_EAV1_out[16] + c_aging_EAV1[15]

# ----- 2.13. EAV2 Compartment --------------------------------------------------

c_EAV2_in[] <- c_EAV1_to_EAV2[i] + c_V2_to_EAV2[i] + c_migration_EAV2[i]
c_EAV2_out[] <- c_EAV2_to_EAS[i] + c_EAV2_to_IAV2[i] + c_death_EAV2[i]

update(EAV2[1]) <- EAV2[1] + c_EAV2_in[1] - c_EAV2_out[1] - c_aging_EAV2[1]
update(EAV2[2]) <- EAV2[2] + c_EAV2_in[2] - c_EAV2_out[2] + c_aging_EAV2[1] - c_aging_EAV2[2]
update(EAV2[3]) <- EAV2[3] + c_EAV2_in[3] - c_EAV2_out[3] + c_aging_EAV2[2] - c_aging_EAV2[3]
update(EAV2[4]) <- EAV2[4] + c_EAV2_in[4] - c_EAV2_out[4] + c_aging_EAV2[3] - c_aging_EAV2[4]
update(EAV2[5]) <- EAV2[5] + c_EAV2_in[5] - c_EAV2_out[5] + c_aging_EAV2[4] - c_aging_EAV2[5]
update(EAV2[6]) <- EAV2[6] + c_EAV2_in[6] - c_EAV2_out[6] + c_aging_EAV2[5] - c_aging_EAV2[6]
update(EAV2[7]) <- EAV2[7] + c_EAV2_in[7] - c_EAV2_out[7] + c_aging_EAV2[6] - c_aging_EAV2[7]
update(EAV2[8]) <- EAV2[8] + c_EAV2_in[8] - c_EAV2_out[8] + c_aging_EAV2[7] - c_aging_EAV2[8]
update(EAV2[9]) <- EAV2[9] + c_EAV2_in[9] - c_EAV2_out[9] + c_aging_EAV2[8] - c_aging_EAV2[9]
update(EAV2[10]) <- EAV2[10] + c_EAV2_in[10] - c_EAV2_out[10] + c_aging_EAV2[9] - c_aging_EAV2[10]
update(EAV2[11]) <- EAV2[11] + c_EAV2_in[11] - c_EAV2_out[11] + c_aging_EAV2[10] - c_aging_EAV2[11]
update(EAV2[12]) <- EAV2[12] + c_EAV2_in[12] - c_EAV2_out[12] + c_aging_EAV2[11] - c_aging_EAV2[12]
update(EAV2[13]) <- EAV2[13] + c_EAV2_in[13] - c_EAV2_out[13] + c_aging_EAV2[12] - c_aging_EAV2[13]
update(EAV2[14]) <- EAV2[14] + c_EAV2_in[14] - c_EAV2_out[14] + c_aging_EAV2[13] - c_aging_EAV2[14]
update(EAV2[15]) <- EAV2[15] + c_EAV2_in[15] - c_EAV2_out[15] + c_aging_EAV2[14] - c_aging_EAV2[15]
update(EAV2[16]) <- EAV2[16] + c_EAV2_in[16] - c_EAV2_out[16] + c_aging_EAV2[15]

# ----- 2.15. EDS Compartment --------------------------------------------------

c_EDS_in[] <- c_EDV1_to_EDS[i] + c_EDV2_to_EDS[i] + c_S_to_EDS[i] + c_migration_EDS[i]
c_EDS_out[] <- c_EDS_to_EDV1[i] + c_EDS_to_IDS[i] + c_death_EDS[i]

update(EDS[1]) <- EDS[1] + c_EDS_in[1] - c_EDS_out[1] - c_aging_EDS[1]
update(EDS[2]) <- EDS[2] + c_EDS_in[2] - c_EDS_out[2] + c_aging_EDS[1] - c_aging_EDS[2]
update(EDS[3]) <- EDS[3] + c_EDS_in[3] - c_EDS_out[3] + c_aging_EDS[2] - c_aging_EDS[3]
update(EDS[4]) <- EDS[4] + c_EDS_in[4] - c_EDS_out[4] + c_aging_EDS[3] - c_aging_EDS[4]
update(EDS[5]) <- EDS[5] + c_EDS_in[5] - c_EDS_out[5] + c_aging_EDS[4] - c_aging_EDS[5]
update(EDS[6]) <- EDS[6] + c_EDS_in[6] - c_EDS_out[6] + c_aging_EDS[5] - c_aging_EDS[6]
update(EDS[7]) <- EDS[7] + c_EDS_in[7] - c_EDS_out[7] + c_aging_EDS[6] - c_aging_EDS[7]
update(EDS[8]) <- EDS[8] + c_EDS_in[8] - c_EDS_out[8] + c_aging_EDS[7] - c_aging_EDS[8]
update(EDS[9]) <- EDS[9] + c_EDS_in[9] - c_EDS_out[9] + c_aging_EDS[8] - c_aging_EDS[9]
update(EDS[10]) <- EDS[10] + c_EDS_in[10] - c_EDS_out[10] + c_aging_EDS[9] - c_aging_EDS[10]
update(EDS[11]) <- EDS[11] + c_EDS_in[11] - c_EDS_out[11] + c_aging_EDS[10] - c_aging_EDS[11]
update(EDS[12]) <- EDS[12] + c_EDS_in[12] - c_EDS_out[12] + c_aging_EDS[11] - c_aging_EDS[12]
update(EDS[13]) <- EDS[13] + c_EDS_in[13] - c_EDS_out[13] + c_aging_EDS[12] - c_aging_EDS[13]
update(EDS[14]) <- EDS[14] + c_EDS_in[14] - c_EDS_out[14] + c_aging_EDS[13] - c_aging_EDS[14]
update(EDS[15]) <- EDS[15] + c_EDS_in[15] - c_EDS_out[15] + c_aging_EDS[14] - c_aging_EDS[15]
update(EDS[16]) <- EDS[16] + c_EDS_in[16] - c_EDS_out[16] + c_aging_EDS[15]

# ----- 2.16. EDV1 Compartment --------------------------------------------------

c_EDV1_in[] <- c_EDS_to_EDV1[i] + c_V1_to_EDV1[i] + c_migration_EDV1[i]
c_EDV1_out[] <- c_EDV1_to_EDS[i] + c_EDV1_to_EDV2[i] + c_EDV1_to_IDV1[i] + c_death_EDV1[i]

update(EDV1[1]) <- EDV1[1] + c_EDV1_in[1] - c_EDV1_out[1] - c_aging_EDV1[1]
update(EDV1[2]) <- EDV1[2] + c_EDV1_in[2] - c_EDV1_out[2] + c_aging_EDV1[1] - c_aging_EDV1[2]
update(EDV1[3]) <- EDV1[3] + c_EDV1_in[3] - c_EDV1_out[3] + c_aging_EDV1[2] - c_aging_EDV1[3]
update(EDV1[4]) <- EDV1[4] + c_EDV1_in[4] - c_EDV1_out[4] + c_aging_EDV1[3] - c_aging_EDV1[4]
update(EDV1[5]) <- EDV1[5] + c_EDV1_in[5] - c_EDV1_out[5] + c_aging_EDV1[4] - c_aging_EDV1[5]
update(EDV1[6]) <- EDV1[6] + c_EDV1_in[6] - c_EDV1_out[6] + c_aging_EDV1[5] - c_aging_EDV1[6]
update(EDV1[7]) <- EDV1[7] + c_EDV1_in[7] - c_EDV1_out[7] + c_aging_EDV1[6] - c_aging_EDV1[7]
update(EDV1[8]) <- EDV1[8] + c_EDV1_in[8] - c_EDV1_out[8] + c_aging_EDV1[7] - c_aging_EDV1[8]
update(EDV1[9]) <- EDV1[9] + c_EDV1_in[9] - c_EDV1_out[9] + c_aging_EDV1[8] - c_aging_EDV1[9]
update(EDV1[10]) <- EDV1[10] + c_EDV1_in[10] - c_EDV1_out[10] + c_aging_EDV1[9] - c_aging_EDV1[10]
update(EDV1[11]) <- EDV1[11] + c_EDV1_in[11] - c_EDV1_out[11] + c_aging_EDV1[10] - c_aging_EDV1[11]
update(EDV1[12]) <- EDV1[12] + c_EDV1_in[12] - c_EDV1_out[12] + c_aging_EDV1[11] - c_aging_EDV1[12]
update(EDV1[13]) <- EDV1[13] + c_EDV1_in[13] - c_EDV1_out[13] + c_aging_EDV1[12] - c_aging_EDV1[13]
update(EDV1[14]) <- EDV1[14] + c_EDV1_in[14] - c_EDV1_out[14] + c_aging_EDV1[13] - c_aging_EDV1[14]
update(EDV1[15]) <- EDV1[15] + c_EDV1_in[15] - c_EDV1_out[15] + c_aging_EDV1[14] - c_aging_EDV1[15]
update(EDV1[16]) <- EDV1[16] + c_EDV1_in[16] - c_EDV1_out[16] + c_aging_EDV1[15]

# ----- 2.17. EDV2 Compartment --------------------------------------------------

c_EDV2_in[] <- c_EDV1_to_EDV2[i] + c_V2_to_EDV2[i] + c_migration_EDV2[i]
c_EDV2_out[] <- c_EDV2_to_EDS[i] + c_EDV2_to_IDV2[i] + c_death_EDV2[i]

update(EDV2[1]) <- EDV2[1] + c_EDV2_in[1] - c_EDV2_out[1] - c_aging_EDV2[1]
update(EDV2[2]) <- EDV2[2] + c_EDV2_in[2] - c_EDV2_out[2] + c_aging_EDV2[1] - c_aging_EDV2[2]
update(EDV2[3]) <- EDV2[3] + c_EDV2_in[3] - c_EDV2_out[3] + c_aging_EDV2[2] - c_aging_EDV2[3]
update(EDV2[4]) <- EDV2[4] + c_EDV2_in[4] - c_EDV2_out[4] + c_aging_EDV2[3] - c_aging_EDV2[4]
update(EDV2[5]) <- EDV2[5] + c_EDV2_in[5] - c_EDV2_out[5] + c_aging_EDV2[4] - c_aging_EDV2[5]
update(EDV2[6]) <- EDV2[6] + c_EDV2_in[6] - c_EDV2_out[6] + c_aging_EDV2[5] - c_aging_EDV2[6]
update(EDV2[7]) <- EDV2[7] + c_EDV2_in[7] - c_EDV2_out[7] + c_aging_EDV2[6] - c_aging_EDV2[7]
update(EDV2[8]) <- EDV2[8] + c_EDV2_in[8] - c_EDV2_out[8] + c_aging_EDV2[7] - c_aging_EDV2[8]
update(EDV2[9]) <- EDV2[9] + c_EDV2_in[9] - c_EDV2_out[9] + c_aging_EDV2[8] - c_aging_EDV2[9]
update(EDV2[10]) <- EDV2[10] + c_EDV2_in[10] - c_EDV2_out[10] + c_aging_EDV2[9] - c_aging_EDV2[10]
update(EDV2[11]) <- EDV2[11] + c_EDV2_in[11] - c_EDV2_out[11] + c_aging_EDV2[10] - c_aging_EDV2[11]
update(EDV2[12]) <- EDV2[12] + c_EDV2_in[12] - c_EDV2_out[12] + c_aging_EDV2[11] - c_aging_EDV2[12]
update(EDV2[13]) <- EDV2[13] + c_EDV2_in[13] - c_EDV2_out[13] + c_aging_EDV2[12] - c_aging_EDV2[13]
update(EDV2[14]) <- EDV2[14] + c_EDV2_in[14] - c_EDV2_out[14] + c_aging_EDV2[13] - c_aging_EDV2[14]
update(EDV2[15]) <- EDV2[15] + c_EDV2_in[15] - c_EDV2_out[15] + c_aging_EDV2[14] - c_aging_EDV2[15]
update(EDV2[16]) <- EDV2[16] + c_EDV2_in[16] - c_EDV2_out[16] + c_aging_EDV2[15]

# ----- 2.19. EOS Compartment --------------------------------------------------

c_EOS_in[] <- c_EOV1_to_EOS[i] + c_EOV2_to_EOS[i] + c_S_to_EOS[i] + c_migration_EOS[i]
c_EOS_out[] <- c_EOS_to_EOV1[i] + c_EOS_to_IOS[i] + c_death_EOS[i]

update(EOS[1]) <- EOS[1] + c_EOS_in[1] - c_EOS_out[1] - c_aging_EOS[1]
update(EOS[2]) <- EOS[2] + c_EOS_in[2] - c_EOS_out[2] + c_aging_EOS[1] - c_aging_EOS[2]
update(EOS[3]) <- EOS[3] + c_EOS_in[3] - c_EOS_out[3] + c_aging_EOS[2] - c_aging_EOS[3]
update(EOS[4]) <- EOS[4] + c_EOS_in[4] - c_EOS_out[4] + c_aging_EOS[3] - c_aging_EOS[4]
update(EOS[5]) <- EOS[5] + c_EOS_in[5] - c_EOS_out[5] + c_aging_EOS[4] - c_aging_EOS[5]
update(EOS[6]) <- EOS[6] + c_EOS_in[6] - c_EOS_out[6] + c_aging_EOS[5] - c_aging_EOS[6]
update(EOS[7]) <- EOS[7] + c_EOS_in[7] - c_EOS_out[7] + c_aging_EOS[6] - c_aging_EOS[7]
update(EOS[8]) <- EOS[8] + c_EOS_in[8] - c_EOS_out[8] + c_aging_EOS[7] - c_aging_EOS[8]
update(EOS[9]) <- EOS[9] + c_EOS_in[9] - c_EOS_out[9] + c_aging_EOS[8] - c_aging_EOS[9]
update(EOS[10]) <- EOS[10] + c_EOS_in[10] - c_EOS_out[10] + c_aging_EOS[9] - c_aging_EOS[10]
update(EOS[11]) <- EOS[11] + c_EOS_in[11] - c_EOS_out[11] + c_aging_EOS[10] - c_aging_EOS[11]
update(EOS[12]) <- EOS[12] + c_EOS_in[12] - c_EOS_out[12] + c_aging_EOS[11] - c_aging_EOS[12]
update(EOS[13]) <- EOS[13] + c_EOS_in[13] - c_EOS_out[13] + c_aging_EOS[12] - c_aging_EOS[13]
update(EOS[14]) <- EOS[14] + c_EOS_in[14] - c_EOS_out[14] + c_aging_EOS[13] - c_aging_EOS[14]
update(EOS[15]) <- EOS[15] + c_EOS_in[15] - c_EOS_out[15] + c_aging_EOS[14] - c_aging_EOS[15]
update(EOS[16]) <- EOS[16] + c_EOS_in[16] - c_EOS_out[16] + c_aging_EOS[15]

# ----- 2.20. EOV1 Compartment --------------------------------------------------

c_EOV1_in[] <- c_EOS_to_EOV1[i] + c_V1_to_EOV1[i] + c_migration_EOV1[i]
c_EOV1_out[] <- c_EOV1_to_EOS[i] + c_EOV1_to_EOV2[i] + c_EOV1_to_IOV1[i] + c_death_EOV1[i]

update(EOV1[1]) <- EOV1[1] + c_EOV1_in[1] - c_EOV1_out[1] - c_aging_EOV1[1]
update(EOV1[2]) <- EOV1[2] + c_EOV1_in[2] - c_EOV1_out[2] + c_aging_EOV1[1] - c_aging_EOV1[2]
update(EOV1[3]) <- EOV1[3] + c_EOV1_in[3] - c_EOV1_out[3] + c_aging_EOV1[2] - c_aging_EOV1[3]
update(EOV1[4]) <- EOV1[4] + c_EOV1_in[4] - c_EOV1_out[4] + c_aging_EOV1[3] - c_aging_EOV1[4]
update(EOV1[5]) <- EOV1[5] + c_EOV1_in[5] - c_EOV1_out[5] + c_aging_EOV1[4] - c_aging_EOV1[5]
update(EOV1[6]) <- EOV1[6] + c_EOV1_in[6] - c_EOV1_out[6] + c_aging_EOV1[5] - c_aging_EOV1[6]
update(EOV1[7]) <- EOV1[7] + c_EOV1_in[7] - c_EOV1_out[7] + c_aging_EOV1[6] - c_aging_EOV1[7]
update(EOV1[8]) <- EOV1[8] + c_EOV1_in[8] - c_EOV1_out[8] + c_aging_EOV1[7] - c_aging_EOV1[8]
update(EOV1[9]) <- EOV1[9] + c_EOV1_in[9] - c_EOV1_out[9] + c_aging_EOV1[8] - c_aging_EOV1[9]
update(EOV1[10]) <- EOV1[10] + c_EOV1_in[10] - c_EOV1_out[10] + c_aging_EOV1[9] - c_aging_EOV1[10]
update(EOV1[11]) <- EOV1[11] + c_EOV1_in[11] - c_EOV1_out[11] + c_aging_EOV1[10] - c_aging_EOV1[11]
update(EOV1[12]) <- EOV1[12] + c_EOV1_in[12] - c_EOV1_out[12] + c_aging_EOV1[11] - c_aging_EOV1[12]
update(EOV1[13]) <- EOV1[13] + c_EOV1_in[13] - c_EOV1_out[13] + c_aging_EOV1[12] - c_aging_EOV1[13]
update(EOV1[14]) <- EOV1[14] + c_EOV1_in[14] - c_EOV1_out[14] + c_aging_EOV1[13] - c_aging_EOV1[14]
update(EOV1[15]) <- EOV1[15] + c_EOV1_in[15] - c_EOV1_out[15] + c_aging_EOV1[14] - c_aging_EOV1[15]
update(EOV1[16]) <- EOV1[16] + c_EOV1_in[16] - c_EOV1_out[16] + c_aging_EOV1[15]

# ----- 2.21. EOV2 Compartment --------------------------------------------------

c_EOV2_in[] <- c_EOV1_to_EOV2[i] + c_V2_to_EOV2[i] + c_migration_EOV2[i]
c_EOV2_out[] <- c_EOV2_to_EOS[i] + c_EOV2_to_IOV2[i] + c_death_EOV2[i]

update(EOV2[1]) <- EOV2[1] + c_EOV2_in[1] - c_EOV2_out[1] - c_aging_EOV2[1]
update(EOV2[2]) <- EOV2[2] + c_EOV2_in[2] - c_EOV2_out[2] + c_aging_EOV2[1] - c_aging_EOV2[2]
update(EOV2[3]) <- EOV2[3] + c_EOV2_in[3] - c_EOV2_out[3] + c_aging_EOV2[2] - c_aging_EOV2[3]
update(EOV2[4]) <- EOV2[4] + c_EOV2_in[4] - c_EOV2_out[4] + c_aging_EOV2[3] - c_aging_EOV2[4]
update(EOV2[5]) <- EOV2[5] + c_EOV2_in[5] - c_EOV2_out[5] + c_aging_EOV2[4] - c_aging_EOV2[5]
update(EOV2[6]) <- EOV2[6] + c_EOV2_in[6] - c_EOV2_out[6] + c_aging_EOV2[5] - c_aging_EOV2[6]
update(EOV2[7]) <- EOV2[7] + c_EOV2_in[7] - c_EOV2_out[7] + c_aging_EOV2[6] - c_aging_EOV2[7]
update(EOV2[8]) <- EOV2[8] + c_EOV2_in[8] - c_EOV2_out[8] + c_aging_EOV2[7] - c_aging_EOV2[8]
update(EOV2[9]) <- EOV2[9] + c_EOV2_in[9] - c_EOV2_out[9] + c_aging_EOV2[8] - c_aging_EOV2[9]
update(EOV2[10]) <- EOV2[10] + c_EOV2_in[10] - c_EOV2_out[10] + c_aging_EOV2[9] - c_aging_EOV2[10]
update(EOV2[11]) <- EOV2[11] + c_EOV2_in[11] - c_EOV2_out[11] + c_aging_EOV2[10] - c_aging_EOV2[11]
update(EOV2[12]) <- EOV2[12] + c_EOV2_in[12] - c_EOV2_out[12] + c_aging_EOV2[11] - c_aging_EOV2[12]
update(EOV2[13]) <- EOV2[13] + c_EOV2_in[13] - c_EOV2_out[13] + c_aging_EOV2[12] - c_aging_EOV2[13]
update(EOV2[14]) <- EOV2[14] + c_EOV2_in[14] - c_EOV2_out[14] + c_aging_EOV2[13] - c_aging_EOV2[14]
update(EOV2[15]) <- EOV2[15] + c_EOV2_in[15] - c_EOV2_out[15] + c_aging_EOV2[14] - c_aging_EOV2[15]
update(EOV2[16]) <- EOV2[16] + c_EOV2_in[16] - c_EOV2_out[16] + c_aging_EOV2[15]

# ----- 2.23. EXS Compartment --------------------------------------------------

c_EXS_in[] <- c_EXV1_to_EXS[i] + c_EXV2_to_EXS[i] + c_S_to_EXS[i] + c_migration_EXS[i]  + c_S_exp_X[i]
c_EXS_out[] <- c_EXS_to_EXV1[i] + c_EXS_to_IXS[i] + c_death_EXS[i]

update(EXS[1]) <- EXS[1] + c_EXS_in[1] - c_EXS_out[1] - c_aging_EXS[1]
update(EXS[2]) <- EXS[2] + c_EXS_in[2] - c_EXS_out[2] + c_aging_EXS[1] - c_aging_EXS[2]
update(EXS[3]) <- EXS[3] + c_EXS_in[3] - c_EXS_out[3] + c_aging_EXS[2] - c_aging_EXS[3]
update(EXS[4]) <- EXS[4] + c_EXS_in[4] - c_EXS_out[4] + c_aging_EXS[3] - c_aging_EXS[4]
update(EXS[5]) <- EXS[5] + c_EXS_in[5] - c_EXS_out[5] + c_aging_EXS[4] - c_aging_EXS[5]
update(EXS[6]) <- EXS[6] + c_EXS_in[6] - c_EXS_out[6] + c_aging_EXS[5] - c_aging_EXS[6]
update(EXS[7]) <- EXS[7] + c_EXS_in[7] - c_EXS_out[7] + c_aging_EXS[6] - c_aging_EXS[7]
update(EXS[8]) <- EXS[8] + c_EXS_in[8] - c_EXS_out[8] + c_aging_EXS[7] - c_aging_EXS[8]
update(EXS[9]) <- EXS[9] + c_EXS_in[9] - c_EXS_out[9] + c_aging_EXS[8] - c_aging_EXS[9]
update(EXS[10]) <- EXS[10] + c_EXS_in[10] - c_EXS_out[10] + c_aging_EXS[9] - c_aging_EXS[10]
update(EXS[11]) <- EXS[11] + c_EXS_in[11] - c_EXS_out[11] + c_aging_EXS[10] - c_aging_EXS[11]
update(EXS[12]) <- EXS[12] + c_EXS_in[12] - c_EXS_out[12] + c_aging_EXS[11] - c_aging_EXS[12]
update(EXS[13]) <- EXS[13] + c_EXS_in[13] - c_EXS_out[13] + c_aging_EXS[12] - c_aging_EXS[13]
update(EXS[14]) <- EXS[14] + c_EXS_in[14] - c_EXS_out[14] + c_aging_EXS[13] - c_aging_EXS[14]
update(EXS[15]) <- EXS[15] + c_EXS_in[15] - c_EXS_out[15] + c_aging_EXS[14] - c_aging_EXS[15]
update(EXS[16]) <- EXS[16] + c_EXS_in[16] - c_EXS_out[16] + c_aging_EXS[15]

# ----- 2.24. EXV1 Compartment --------------------------------------------------

c_EXV1_in[] <- c_EXS_to_EXV1[i] + c_V1_to_EXV1[i] + c_migration_EXV1[i]
c_EXV1_out[] <- c_EXV1_to_EXS[i] + c_EXV1_to_EXV2[i] + c_EXV1_to_IXV1[i] + c_death_EXV1[i]

update(EXV1[1]) <- EXV1[1] + c_EXV1_in[1] - c_EXV1_out[1] - c_aging_EXV1[1]
update(EXV1[2]) <- EXV1[2] + c_EXV1_in[2] - c_EXV1_out[2] + c_aging_EXV1[1] - c_aging_EXV1[2]
update(EXV1[3]) <- EXV1[3] + c_EXV1_in[3] - c_EXV1_out[3] + c_aging_EXV1[2] - c_aging_EXV1[3]
update(EXV1[4]) <- EXV1[4] + c_EXV1_in[4] - c_EXV1_out[4] + c_aging_EXV1[3] - c_aging_EXV1[4]
update(EXV1[5]) <- EXV1[5] + c_EXV1_in[5] - c_EXV1_out[5] + c_aging_EXV1[4] - c_aging_EXV1[5]
update(EXV1[6]) <- EXV1[6] + c_EXV1_in[6] - c_EXV1_out[6] + c_aging_EXV1[5] - c_aging_EXV1[6]
update(EXV1[7]) <- EXV1[7] + c_EXV1_in[7] - c_EXV1_out[7] + c_aging_EXV1[6] - c_aging_EXV1[7]
update(EXV1[8]) <- EXV1[8] + c_EXV1_in[8] - c_EXV1_out[8] + c_aging_EXV1[7] - c_aging_EXV1[8]
update(EXV1[9]) <- EXV1[9] + c_EXV1_in[9] - c_EXV1_out[9] + c_aging_EXV1[8] - c_aging_EXV1[9]
update(EXV1[10]) <- EXV1[10] + c_EXV1_in[10] - c_EXV1_out[10] + c_aging_EXV1[9] - c_aging_EXV1[10]
update(EXV1[11]) <- EXV1[11] + c_EXV1_in[11] - c_EXV1_out[11] + c_aging_EXV1[10] - c_aging_EXV1[11]
update(EXV1[12]) <- EXV1[12] + c_EXV1_in[12] - c_EXV1_out[12] + c_aging_EXV1[11] - c_aging_EXV1[12]
update(EXV1[13]) <- EXV1[13] + c_EXV1_in[13] - c_EXV1_out[13] + c_aging_EXV1[12] - c_aging_EXV1[13]
update(EXV1[14]) <- EXV1[14] + c_EXV1_in[14] - c_EXV1_out[14] + c_aging_EXV1[13] - c_aging_EXV1[14]
update(EXV1[15]) <- EXV1[15] + c_EXV1_in[15] - c_EXV1_out[15] + c_aging_EXV1[14] - c_aging_EXV1[15]
update(EXV1[16]) <- EXV1[16] + c_EXV1_in[16] - c_EXV1_out[16] + c_aging_EXV1[15]

# ----- 2.25. EXV2 Compartment --------------------------------------------------

c_EXV2_in[] <- c_EXV1_to_EXV2[i] + c_V2_to_EXV2[i] + c_migration_EXV2[i]
c_EXV2_out[] <- c_EXV2_to_EXS[i] + c_EXV2_to_IXV2[i] + c_death_EXV2[i]

update(EXV2[1]) <- EXV2[1] + c_EXV2_in[1] - c_EXV2_out[1] - c_aging_EXV2[1]
update(EXV2[2]) <- EXV2[2] + c_EXV2_in[2] - c_EXV2_out[2] + c_aging_EXV2[1] - c_aging_EXV2[2]
update(EXV2[3]) <- EXV2[3] + c_EXV2_in[3] - c_EXV2_out[3] + c_aging_EXV2[2] - c_aging_EXV2[3]
update(EXV2[4]) <- EXV2[4] + c_EXV2_in[4] - c_EXV2_out[4] + c_aging_EXV2[3] - c_aging_EXV2[4]
update(EXV2[5]) <- EXV2[5] + c_EXV2_in[5] - c_EXV2_out[5] + c_aging_EXV2[4] - c_aging_EXV2[5]
update(EXV2[6]) <- EXV2[6] + c_EXV2_in[6] - c_EXV2_out[6] + c_aging_EXV2[5] - c_aging_EXV2[6]
update(EXV2[7]) <- EXV2[7] + c_EXV2_in[7] - c_EXV2_out[7] + c_aging_EXV2[6] - c_aging_EXV2[7]
update(EXV2[8]) <- EXV2[8] + c_EXV2_in[8] - c_EXV2_out[8] + c_aging_EXV2[7] - c_aging_EXV2[8]
update(EXV2[9]) <- EXV2[9] + c_EXV2_in[9] - c_EXV2_out[9] + c_aging_EXV2[8] - c_aging_EXV2[9]
update(EXV2[10]) <- EXV2[10] + c_EXV2_in[10] - c_EXV2_out[10] + c_aging_EXV2[9] - c_aging_EXV2[10]
update(EXV2[11]) <- EXV2[11] + c_EXV2_in[11] - c_EXV2_out[11] + c_aging_EXV2[10] - c_aging_EXV2[11]
update(EXV2[12]) <- EXV2[12] + c_EXV2_in[12] - c_EXV2_out[12] + c_aging_EXV2[11] - c_aging_EXV2[12]
update(EXV2[13]) <- EXV2[13] + c_EXV2_in[13] - c_EXV2_out[13] + c_aging_EXV2[12] - c_aging_EXV2[13]
update(EXV2[14]) <- EXV2[14] + c_EXV2_in[14] - c_EXV2_out[14] + c_aging_EXV2[13] - c_aging_EXV2[14]
update(EXV2[15]) <- EXV2[15] + c_EXV2_in[15] - c_EXV2_out[15] + c_aging_EXV2[14] - c_aging_EXV2[15]
update(EXV2[16]) <- EXV2[16] + c_EXV2_in[16] - c_EXV2_out[16] + c_aging_EXV2[15]

# ----- 2.27. IWS Compartment --------------------------------------------------

c_IWS_in[] <- c_EWS_to_IWS[i] + c_IWV1_to_IWS[i] + c_IWV2_to_IWS[i] + c_migration_IWS[i]
c_IWS_out[] <- c_IWS_to_IWV1[i] + c_IWS_to_RWS[i] + c_death_inf_IWS[i] + c_death_IWS[i]

update(IWS[1]) <- IWS[1] + c_IWS_in[1] - c_IWS_out[1] - c_aging_IWS[1]
update(IWS[2]) <- IWS[2] + c_IWS_in[2] - c_IWS_out[2] + c_aging_IWS[1] - c_aging_IWS[2]
update(IWS[3]) <- IWS[3] + c_IWS_in[3] - c_IWS_out[3] + c_aging_IWS[2] - c_aging_IWS[3]
update(IWS[4]) <- IWS[4] + c_IWS_in[4] - c_IWS_out[4] + c_aging_IWS[3] - c_aging_IWS[4]
update(IWS[5]) <- IWS[5] + c_IWS_in[5] - c_IWS_out[5] + c_aging_IWS[4] - c_aging_IWS[5]
update(IWS[6]) <- IWS[6] + c_IWS_in[6] - c_IWS_out[6] + c_aging_IWS[5] - c_aging_IWS[6]
update(IWS[7]) <- IWS[7] + c_IWS_in[7] - c_IWS_out[7] + c_aging_IWS[6] - c_aging_IWS[7]
update(IWS[8]) <- IWS[8] + c_IWS_in[8] - c_IWS_out[8] + c_aging_IWS[7] - c_aging_IWS[8]
update(IWS[9]) <- IWS[9] + c_IWS_in[9] - c_IWS_out[9] + c_aging_IWS[8] - c_aging_IWS[9]
update(IWS[10]) <- IWS[10] + c_IWS_in[10] - c_IWS_out[10] + c_aging_IWS[9] - c_aging_IWS[10]
update(IWS[11]) <- IWS[11] + c_IWS_in[11] - c_IWS_out[11] + c_aging_IWS[10] - c_aging_IWS[11]
update(IWS[12]) <- IWS[12] + c_IWS_in[12] - c_IWS_out[12] + c_aging_IWS[11] - c_aging_IWS[12]
update(IWS[13]) <- IWS[13] + c_IWS_in[13] - c_IWS_out[13] + c_aging_IWS[12] - c_aging_IWS[13]
update(IWS[14]) <- IWS[14] + c_IWS_in[14] - c_IWS_out[14] + c_aging_IWS[13] - c_aging_IWS[14]
update(IWS[15]) <- IWS[15] + c_IWS_in[15] - c_IWS_out[15] + c_aging_IWS[14] - c_aging_IWS[15]
update(IWS[16]) <- IWS[16] + c_IWS_in[16] - c_IWS_out[16] + c_aging_IWS[15]

# ----- 2.28. IWV1 Compartment --------------------------------------------------

c_IWV1_in[] <- c_EWV1_to_IWV1[i] + c_IWS_to_IWV1[i] + c_migration_IWV1[i]
c_IWV1_out[] <- c_IWV1_to_IWS[i] + c_IWV1_to_IWV2[i] + c_IWV1_to_RWV1[i] + c_death_inf_IWV1[i] + c_death_IWV1[i]

update(IWV1[1]) <- IWV1[1] + c_IWV1_in[1] - c_IWV1_out[1] - c_aging_IWV1[1]
update(IWV1[2]) <- IWV1[2] + c_IWV1_in[2] - c_IWV1_out[2] + c_aging_IWV1[1] - c_aging_IWV1[2]
update(IWV1[3]) <- IWV1[3] + c_IWV1_in[3] - c_IWV1_out[3] + c_aging_IWV1[2] - c_aging_IWV1[3]
update(IWV1[4]) <- IWV1[4] + c_IWV1_in[4] - c_IWV1_out[4] + c_aging_IWV1[3] - c_aging_IWV1[4]
update(IWV1[5]) <- IWV1[5] + c_IWV1_in[5] - c_IWV1_out[5] + c_aging_IWV1[4] - c_aging_IWV1[5]
update(IWV1[6]) <- IWV1[6] + c_IWV1_in[6] - c_IWV1_out[6] + c_aging_IWV1[5] - c_aging_IWV1[6]
update(IWV1[7]) <- IWV1[7] + c_IWV1_in[7] - c_IWV1_out[7] + c_aging_IWV1[6] - c_aging_IWV1[7]
update(IWV1[8]) <- IWV1[8] + c_IWV1_in[8] - c_IWV1_out[8] + c_aging_IWV1[7] - c_aging_IWV1[8]
update(IWV1[9]) <- IWV1[9] + c_IWV1_in[9] - c_IWV1_out[9] + c_aging_IWV1[8] - c_aging_IWV1[9]
update(IWV1[10]) <- IWV1[10] + c_IWV1_in[10] - c_IWV1_out[10] + c_aging_IWV1[9] - c_aging_IWV1[10]
update(IWV1[11]) <- IWV1[11] + c_IWV1_in[11] - c_IWV1_out[11] + c_aging_IWV1[10] - c_aging_IWV1[11]
update(IWV1[12]) <- IWV1[12] + c_IWV1_in[12] - c_IWV1_out[12] + c_aging_IWV1[11] - c_aging_IWV1[12]
update(IWV1[13]) <- IWV1[13] + c_IWV1_in[13] - c_IWV1_out[13] + c_aging_IWV1[12] - c_aging_IWV1[13]
update(IWV1[14]) <- IWV1[14] + c_IWV1_in[14] - c_IWV1_out[14] + c_aging_IWV1[13] - c_aging_IWV1[14]
update(IWV1[15]) <- IWV1[15] + c_IWV1_in[15] - c_IWV1_out[15] + c_aging_IWV1[14] - c_aging_IWV1[15]
update(IWV1[16]) <- IWV1[16] + c_IWV1_in[16] - c_IWV1_out[16] + c_aging_IWV1[15]

# ----- 2.29. IWV2 Compartment --------------------------------------------------

c_IWV2_in[] <- c_EWV2_to_IWV2[i] + c_IWV1_to_IWV2[i] + c_migration_IWV2[i]
c_IWV2_out[] <- c_IWV2_to_IWS[i] + c_IWV2_to_RWV2[i] + c_death_inf_IWV2[i] + c_death_IWV2[i]

update(IWV2[1]) <- IWV2[1] + c_IWV2_in[1] - c_IWV2_out[1] - c_aging_IWV2[1]
update(IWV2[2]) <- IWV2[2] + c_IWV2_in[2] - c_IWV2_out[2] + c_aging_IWV2[1] - c_aging_IWV2[2]
update(IWV2[3]) <- IWV2[3] + c_IWV2_in[3] - c_IWV2_out[3] + c_aging_IWV2[2] - c_aging_IWV2[3]
update(IWV2[4]) <- IWV2[4] + c_IWV2_in[4] - c_IWV2_out[4] + c_aging_IWV2[3] - c_aging_IWV2[4]
update(IWV2[5]) <- IWV2[5] + c_IWV2_in[5] - c_IWV2_out[5] + c_aging_IWV2[4] - c_aging_IWV2[5]
update(IWV2[6]) <- IWV2[6] + c_IWV2_in[6] - c_IWV2_out[6] + c_aging_IWV2[5] - c_aging_IWV2[6]
update(IWV2[7]) <- IWV2[7] + c_IWV2_in[7] - c_IWV2_out[7] + c_aging_IWV2[6] - c_aging_IWV2[7]
update(IWV2[8]) <- IWV2[8] + c_IWV2_in[8] - c_IWV2_out[8] + c_aging_IWV2[7] - c_aging_IWV2[8]
update(IWV2[9]) <- IWV2[9] + c_IWV2_in[9] - c_IWV2_out[9] + c_aging_IWV2[8] - c_aging_IWV2[9]
update(IWV2[10]) <- IWV2[10] + c_IWV2_in[10] - c_IWV2_out[10] + c_aging_IWV2[9] - c_aging_IWV2[10]
update(IWV2[11]) <- IWV2[11] + c_IWV2_in[11] - c_IWV2_out[11] + c_aging_IWV2[10] - c_aging_IWV2[11]
update(IWV2[12]) <- IWV2[12] + c_IWV2_in[12] - c_IWV2_out[12] + c_aging_IWV2[11] - c_aging_IWV2[12]
update(IWV2[13]) <- IWV2[13] + c_IWV2_in[13] - c_IWV2_out[13] + c_aging_IWV2[12] - c_aging_IWV2[13]
update(IWV2[14]) <- IWV2[14] + c_IWV2_in[14] - c_IWV2_out[14] + c_aging_IWV2[13] - c_aging_IWV2[14]
update(IWV2[15]) <- IWV2[15] + c_IWV2_in[15] - c_IWV2_out[15] + c_aging_IWV2[14] - c_aging_IWV2[15]
update(IWV2[16]) <- IWV2[16] + c_IWV2_in[16] - c_IWV2_out[16] + c_aging_IWV2[15]

# ----- 2.31. IAS Compartment --------------------------------------------------

c_IAS_in[] <- c_EAS_to_IAS[i] + c_IAV1_to_IAS[i] + c_IAV2_to_IAS[i] + c_migration_IAS[i]
c_IAS_out[] <- c_IAS_to_IAV1[i] + c_IAS_to_RAS[i] + c_death_inf_IAS[i] + c_death_IAS[i]

update(IAS[1]) <- IAS[1] + c_IAS_in[1] - c_IAS_out[1] - c_aging_IAS[1]
update(IAS[2]) <- IAS[2] + c_IAS_in[2] - c_IAS_out[2] + c_aging_IAS[1] - c_aging_IAS[2]
update(IAS[3]) <- IAS[3] + c_IAS_in[3] - c_IAS_out[3] + c_aging_IAS[2] - c_aging_IAS[3]
update(IAS[4]) <- IAS[4] + c_IAS_in[4] - c_IAS_out[4] + c_aging_IAS[3] - c_aging_IAS[4]
update(IAS[5]) <- IAS[5] + c_IAS_in[5] - c_IAS_out[5] + c_aging_IAS[4] - c_aging_IAS[5]
update(IAS[6]) <- IAS[6] + c_IAS_in[6] - c_IAS_out[6] + c_aging_IAS[5] - c_aging_IAS[6]
update(IAS[7]) <- IAS[7] + c_IAS_in[7] - c_IAS_out[7] + c_aging_IAS[6] - c_aging_IAS[7]
update(IAS[8]) <- IAS[8] + c_IAS_in[8] - c_IAS_out[8] + c_aging_IAS[7] - c_aging_IAS[8]
update(IAS[9]) <- IAS[9] + c_IAS_in[9] - c_IAS_out[9] + c_aging_IAS[8] - c_aging_IAS[9]
update(IAS[10]) <- IAS[10] + c_IAS_in[10] - c_IAS_out[10] + c_aging_IAS[9] - c_aging_IAS[10]
update(IAS[11]) <- IAS[11] + c_IAS_in[11] - c_IAS_out[11] + c_aging_IAS[10] - c_aging_IAS[11]
update(IAS[12]) <- IAS[12] + c_IAS_in[12] - c_IAS_out[12] + c_aging_IAS[11] - c_aging_IAS[12]
update(IAS[13]) <- IAS[13] + c_IAS_in[13] - c_IAS_out[13] + c_aging_IAS[12] - c_aging_IAS[13]
update(IAS[14]) <- IAS[14] + c_IAS_in[14] - c_IAS_out[14] + c_aging_IAS[13] - c_aging_IAS[14]
update(IAS[15]) <- IAS[15] + c_IAS_in[15] - c_IAS_out[15] + c_aging_IAS[14] - c_aging_IAS[15]
update(IAS[16]) <- IAS[16] + c_IAS_in[16] - c_IAS_out[16] + c_aging_IAS[15]

# ----- 2.32. IAV1 Compartment --------------------------------------------------

c_IAV1_in[] <- c_EAV1_to_IAV1[i] + c_IAS_to_IAV1[i] + c_migration_IAV1[i]
c_IAV1_out[] <- c_IAV1_to_IAS[i] + c_IAV1_to_IAV2[i] + c_IAV1_to_RAV1[i] + c_death_inf_IAV1[i] + c_death_IAV1[i]

update(IAV1[1]) <- IAV1[1] + c_IAV1_in[1] - c_IAV1_out[1] - c_aging_IAV1[1]
update(IAV1[2]) <- IAV1[2] + c_IAV1_in[2] - c_IAV1_out[2] + c_aging_IAV1[1] - c_aging_IAV1[2]
update(IAV1[3]) <- IAV1[3] + c_IAV1_in[3] - c_IAV1_out[3] + c_aging_IAV1[2] - c_aging_IAV1[3]
update(IAV1[4]) <- IAV1[4] + c_IAV1_in[4] - c_IAV1_out[4] + c_aging_IAV1[3] - c_aging_IAV1[4]
update(IAV1[5]) <- IAV1[5] + c_IAV1_in[5] - c_IAV1_out[5] + c_aging_IAV1[4] - c_aging_IAV1[5]
update(IAV1[6]) <- IAV1[6] + c_IAV1_in[6] - c_IAV1_out[6] + c_aging_IAV1[5] - c_aging_IAV1[6]
update(IAV1[7]) <- IAV1[7] + c_IAV1_in[7] - c_IAV1_out[7] + c_aging_IAV1[6] - c_aging_IAV1[7]
update(IAV1[8]) <- IAV1[8] + c_IAV1_in[8] - c_IAV1_out[8] + c_aging_IAV1[7] - c_aging_IAV1[8]
update(IAV1[9]) <- IAV1[9] + c_IAV1_in[9] - c_IAV1_out[9] + c_aging_IAV1[8] - c_aging_IAV1[9]
update(IAV1[10]) <- IAV1[10] + c_IAV1_in[10] - c_IAV1_out[10] + c_aging_IAV1[9] - c_aging_IAV1[10]
update(IAV1[11]) <- IAV1[11] + c_IAV1_in[11] - c_IAV1_out[11] + c_aging_IAV1[10] - c_aging_IAV1[11]
update(IAV1[12]) <- IAV1[12] + c_IAV1_in[12] - c_IAV1_out[12] + c_aging_IAV1[11] - c_aging_IAV1[12]
update(IAV1[13]) <- IAV1[13] + c_IAV1_in[13] - c_IAV1_out[13] + c_aging_IAV1[12] - c_aging_IAV1[13]
update(IAV1[14]) <- IAV1[14] + c_IAV1_in[14] - c_IAV1_out[14] + c_aging_IAV1[13] - c_aging_IAV1[14]
update(IAV1[15]) <- IAV1[15] + c_IAV1_in[15] - c_IAV1_out[15] + c_aging_IAV1[14] - c_aging_IAV1[15]
update(IAV1[16]) <- IAV1[16] + c_IAV1_in[16] - c_IAV1_out[16] + c_aging_IAV1[15]

# ----- 2.33. IAV2 Compartment --------------------------------------------------

c_IAV2_in[] <- c_EAV2_to_IAV2[i] + c_IAV1_to_IAV2[i] + c_migration_IAV2[i]
c_IAV2_out[] <- c_IAV2_to_IAS[i] + c_IAV2_to_RAV2[i] + c_death_inf_IAV2[i] + c_death_IAV2[i]

update(IAV2[1]) <- IAV2[1] + c_IAV2_in[1] - c_IAV2_out[1] - c_aging_IAV2[1]
update(IAV2[2]) <- IAV2[2] + c_IAV2_in[2] - c_IAV2_out[2] + c_aging_IAV2[1] - c_aging_IAV2[2]
update(IAV2[3]) <- IAV2[3] + c_IAV2_in[3] - c_IAV2_out[3] + c_aging_IAV2[2] - c_aging_IAV2[3]
update(IAV2[4]) <- IAV2[4] + c_IAV2_in[4] - c_IAV2_out[4] + c_aging_IAV2[3] - c_aging_IAV2[4]
update(IAV2[5]) <- IAV2[5] + c_IAV2_in[5] - c_IAV2_out[5] + c_aging_IAV2[4] - c_aging_IAV2[5]
update(IAV2[6]) <- IAV2[6] + c_IAV2_in[6] - c_IAV2_out[6] + c_aging_IAV2[5] - c_aging_IAV2[6]
update(IAV2[7]) <- IAV2[7] + c_IAV2_in[7] - c_IAV2_out[7] + c_aging_IAV2[6] - c_aging_IAV2[7]
update(IAV2[8]) <- IAV2[8] + c_IAV2_in[8] - c_IAV2_out[8] + c_aging_IAV2[7] - c_aging_IAV2[8]
update(IAV2[9]) <- IAV2[9] + c_IAV2_in[9] - c_IAV2_out[9] + c_aging_IAV2[8] - c_aging_IAV2[9]
update(IAV2[10]) <- IAV2[10] + c_IAV2_in[10] - c_IAV2_out[10] + c_aging_IAV2[9] - c_aging_IAV2[10]
update(IAV2[11]) <- IAV2[11] + c_IAV2_in[11] - c_IAV2_out[11] + c_aging_IAV2[10] - c_aging_IAV2[11]
update(IAV2[12]) <- IAV2[12] + c_IAV2_in[12] - c_IAV2_out[12] + c_aging_IAV2[11] - c_aging_IAV2[12]
update(IAV2[13]) <- IAV2[13] + c_IAV2_in[13] - c_IAV2_out[13] + c_aging_IAV2[12] - c_aging_IAV2[13]
update(IAV2[14]) <- IAV2[14] + c_IAV2_in[14] - c_IAV2_out[14] + c_aging_IAV2[13] - c_aging_IAV2[14]
update(IAV2[15]) <- IAV2[15] + c_IAV2_in[15] - c_IAV2_out[15] + c_aging_IAV2[14] - c_aging_IAV2[15]
update(IAV2[16]) <- IAV2[16] + c_IAV2_in[16] - c_IAV2_out[16] + c_aging_IAV2[15]

# ----- 2.35. IDS Compartment --------------------------------------------------

c_IDS_in[] <- c_EDS_to_IDS[i] + c_IDV1_to_IDS[i] + c_IDV2_to_IDS[i] + c_migration_IDS[i]
c_IDS_out[] <- c_IDS_to_IDV1[i] + c_IDS_to_RDS[i] + c_death_inf_IDS[i] + c_death_IDS[i]

update(IDS[1]) <- IDS[1] + c_IDS_in[1] - c_IDS_out[1] - c_aging_IDS[1]
update(IDS[2]) <- IDS[2] + c_IDS_in[2] - c_IDS_out[2] + c_aging_IDS[1] - c_aging_IDS[2]
update(IDS[3]) <- IDS[3] + c_IDS_in[3] - c_IDS_out[3] + c_aging_IDS[2] - c_aging_IDS[3]
update(IDS[4]) <- IDS[4] + c_IDS_in[4] - c_IDS_out[4] + c_aging_IDS[3] - c_aging_IDS[4]
update(IDS[5]) <- IDS[5] + c_IDS_in[5] - c_IDS_out[5] + c_aging_IDS[4] - c_aging_IDS[5]
update(IDS[6]) <- IDS[6] + c_IDS_in[6] - c_IDS_out[6] + c_aging_IDS[5] - c_aging_IDS[6]
update(IDS[7]) <- IDS[7] + c_IDS_in[7] - c_IDS_out[7] + c_aging_IDS[6] - c_aging_IDS[7]
update(IDS[8]) <- IDS[8] + c_IDS_in[8] - c_IDS_out[8] + c_aging_IDS[7] - c_aging_IDS[8]
update(IDS[9]) <- IDS[9] + c_IDS_in[9] - c_IDS_out[9] + c_aging_IDS[8] - c_aging_IDS[9]
update(IDS[10]) <- IDS[10] + c_IDS_in[10] - c_IDS_out[10] + c_aging_IDS[9] - c_aging_IDS[10]
update(IDS[11]) <- IDS[11] + c_IDS_in[11] - c_IDS_out[11] + c_aging_IDS[10] - c_aging_IDS[11]
update(IDS[12]) <- IDS[12] + c_IDS_in[12] - c_IDS_out[12] + c_aging_IDS[11] - c_aging_IDS[12]
update(IDS[13]) <- IDS[13] + c_IDS_in[13] - c_IDS_out[13] + c_aging_IDS[12] - c_aging_IDS[13]
update(IDS[14]) <- IDS[14] + c_IDS_in[14] - c_IDS_out[14] + c_aging_IDS[13] - c_aging_IDS[14]
update(IDS[15]) <- IDS[15] + c_IDS_in[15] - c_IDS_out[15] + c_aging_IDS[14] - c_aging_IDS[15]
update(IDS[16]) <- IDS[16] + c_IDS_in[16] - c_IDS_out[16] + c_aging_IDS[15]

# ----- 2.36. IDV1 Compartment --------------------------------------------------

c_IDV1_in[] <- c_EDV1_to_IDV1[i] + c_IDS_to_IDV1[i] + c_migration_IDV1[i]
c_IDV1_out[] <- c_IDV1_to_IDS[i] + c_IDV1_to_IDV2[i] + c_IDV1_to_RDV1[i] + c_death_inf_IDV1[i] + c_death_IDV1[i]

update(IDV1[1]) <- IDV1[1] + c_IDV1_in[1] - c_IDV1_out[1] - c_aging_IDV1[1]
update(IDV1[2]) <- IDV1[2] + c_IDV1_in[2] - c_IDV1_out[2] + c_aging_IDV1[1] - c_aging_IDV1[2]
update(IDV1[3]) <- IDV1[3] + c_IDV1_in[3] - c_IDV1_out[3] + c_aging_IDV1[2] - c_aging_IDV1[3]
update(IDV1[4]) <- IDV1[4] + c_IDV1_in[4] - c_IDV1_out[4] + c_aging_IDV1[3] - c_aging_IDV1[4]
update(IDV1[5]) <- IDV1[5] + c_IDV1_in[5] - c_IDV1_out[5] + c_aging_IDV1[4] - c_aging_IDV1[5]
update(IDV1[6]) <- IDV1[6] + c_IDV1_in[6] - c_IDV1_out[6] + c_aging_IDV1[5] - c_aging_IDV1[6]
update(IDV1[7]) <- IDV1[7] + c_IDV1_in[7] - c_IDV1_out[7] + c_aging_IDV1[6] - c_aging_IDV1[7]
update(IDV1[8]) <- IDV1[8] + c_IDV1_in[8] - c_IDV1_out[8] + c_aging_IDV1[7] - c_aging_IDV1[8]
update(IDV1[9]) <- IDV1[9] + c_IDV1_in[9] - c_IDV1_out[9] + c_aging_IDV1[8] - c_aging_IDV1[9]
update(IDV1[10]) <- IDV1[10] + c_IDV1_in[10] - c_IDV1_out[10] + c_aging_IDV1[9] - c_aging_IDV1[10]
update(IDV1[11]) <- IDV1[11] + c_IDV1_in[11] - c_IDV1_out[11] + c_aging_IDV1[10] - c_aging_IDV1[11]
update(IDV1[12]) <- IDV1[12] + c_IDV1_in[12] - c_IDV1_out[12] + c_aging_IDV1[11] - c_aging_IDV1[12]
update(IDV1[13]) <- IDV1[13] + c_IDV1_in[13] - c_IDV1_out[13] + c_aging_IDV1[12] - c_aging_IDV1[13]
update(IDV1[14]) <- IDV1[14] + c_IDV1_in[14] - c_IDV1_out[14] + c_aging_IDV1[13] - c_aging_IDV1[14]
update(IDV1[15]) <- IDV1[15] + c_IDV1_in[15] - c_IDV1_out[15] + c_aging_IDV1[14] - c_aging_IDV1[15]
update(IDV1[16]) <- IDV1[16] + c_IDV1_in[16] - c_IDV1_out[16] + c_aging_IDV1[15]

# ----- 2.37. IDV2 Compartment --------------------------------------------------

c_IDV2_in[] <- c_EDV2_to_IDV2[i] + c_IDV1_to_IDV2[i] + c_migration_IDV2[i]
c_IDV2_out[] <- c_IDV2_to_IDS[i] + c_IDV2_to_RDV2[i] + c_death_inf_IDV2[i] + c_death_IDV2[i]

update(IDV2[1]) <- IDV2[1] + c_IDV2_in[1] - c_IDV2_out[1] - c_aging_IDV2[1]
update(IDV2[2]) <- IDV2[2] + c_IDV2_in[2] - c_IDV2_out[2] + c_aging_IDV2[1] - c_aging_IDV2[2]
update(IDV2[3]) <- IDV2[3] + c_IDV2_in[3] - c_IDV2_out[3] + c_aging_IDV2[2] - c_aging_IDV2[3]
update(IDV2[4]) <- IDV2[4] + c_IDV2_in[4] - c_IDV2_out[4] + c_aging_IDV2[3] - c_aging_IDV2[4]
update(IDV2[5]) <- IDV2[5] + c_IDV2_in[5] - c_IDV2_out[5] + c_aging_IDV2[4] - c_aging_IDV2[5]
update(IDV2[6]) <- IDV2[6] + c_IDV2_in[6] - c_IDV2_out[6] + c_aging_IDV2[5] - c_aging_IDV2[6]
update(IDV2[7]) <- IDV2[7] + c_IDV2_in[7] - c_IDV2_out[7] + c_aging_IDV2[6] - c_aging_IDV2[7]
update(IDV2[8]) <- IDV2[8] + c_IDV2_in[8] - c_IDV2_out[8] + c_aging_IDV2[7] - c_aging_IDV2[8]
update(IDV2[9]) <- IDV2[9] + c_IDV2_in[9] - c_IDV2_out[9] + c_aging_IDV2[8] - c_aging_IDV2[9]
update(IDV2[10]) <- IDV2[10] + c_IDV2_in[10] - c_IDV2_out[10] + c_aging_IDV2[9] - c_aging_IDV2[10]
update(IDV2[11]) <- IDV2[11] + c_IDV2_in[11] - c_IDV2_out[11] + c_aging_IDV2[10] - c_aging_IDV2[11]
update(IDV2[12]) <- IDV2[12] + c_IDV2_in[12] - c_IDV2_out[12] + c_aging_IDV2[11] - c_aging_IDV2[12]
update(IDV2[13]) <- IDV2[13] + c_IDV2_in[13] - c_IDV2_out[13] + c_aging_IDV2[12] - c_aging_IDV2[13]
update(IDV2[14]) <- IDV2[14] + c_IDV2_in[14] - c_IDV2_out[14] + c_aging_IDV2[13] - c_aging_IDV2[14]
update(IDV2[15]) <- IDV2[15] + c_IDV2_in[15] - c_IDV2_out[15] + c_aging_IDV2[14] - c_aging_IDV2[15]
update(IDV2[16]) <- IDV2[16] + c_IDV2_in[16] - c_IDV2_out[16] + c_aging_IDV2[15]

# ----- 2.39. IOS Compartment --------------------------------------------------

c_IOS_in[] <- c_EOS_to_IOS[i] + c_IOV1_to_IOS[i] + c_IOV2_to_IOS[i] + c_migration_IOS[i]
c_IOS_out[] <- c_IOS_to_IOV1[i] + c_IOS_to_ROS[i] + c_death_inf_IOS[i] + c_death_IOS[i]

update(IOS[1]) <- IOS[1] + c_IOS_in[1] - c_IOS_out[1] - c_aging_IOS[1]
update(IOS[2]) <- IOS[2] + c_IOS_in[2] - c_IOS_out[2] + c_aging_IOS[1] - c_aging_IOS[2]
update(IOS[3]) <- IOS[3] + c_IOS_in[3] - c_IOS_out[3] + c_aging_IOS[2] - c_aging_IOS[3]
update(IOS[4]) <- IOS[4] + c_IOS_in[4] - c_IOS_out[4] + c_aging_IOS[3] - c_aging_IOS[4]
update(IOS[5]) <- IOS[5] + c_IOS_in[5] - c_IOS_out[5] + c_aging_IOS[4] - c_aging_IOS[5]
update(IOS[6]) <- IOS[6] + c_IOS_in[6] - c_IOS_out[6] + c_aging_IOS[5] - c_aging_IOS[6]
update(IOS[7]) <- IOS[7] + c_IOS_in[7] - c_IOS_out[7] + c_aging_IOS[6] - c_aging_IOS[7]
update(IOS[8]) <- IOS[8] + c_IOS_in[8] - c_IOS_out[8] + c_aging_IOS[7] - c_aging_IOS[8]
update(IOS[9]) <- IOS[9] + c_IOS_in[9] - c_IOS_out[9] + c_aging_IOS[8] - c_aging_IOS[9]
update(IOS[10]) <- IOS[10] + c_IOS_in[10] - c_IOS_out[10] + c_aging_IOS[9] - c_aging_IOS[10]
update(IOS[11]) <- IOS[11] + c_IOS_in[11] - c_IOS_out[11] + c_aging_IOS[10] - c_aging_IOS[11]
update(IOS[12]) <- IOS[12] + c_IOS_in[12] - c_IOS_out[12] + c_aging_IOS[11] - c_aging_IOS[12]
update(IOS[13]) <- IOS[13] + c_IOS_in[13] - c_IOS_out[13] + c_aging_IOS[12] - c_aging_IOS[13]
update(IOS[14]) <- IOS[14] + c_IOS_in[14] - c_IOS_out[14] + c_aging_IOS[13] - c_aging_IOS[14]
update(IOS[15]) <- IOS[15] + c_IOS_in[15] - c_IOS_out[15] + c_aging_IOS[14] - c_aging_IOS[15]
update(IOS[16]) <- IOS[16] + c_IOS_in[16] - c_IOS_out[16] + c_aging_IOS[15]

# ----- 2.40. IOV1 Compartment --------------------------------------------------

c_IOV1_in[] <- c_EOV1_to_IOV1[i] + c_IOS_to_IOV1[i] + c_migration_IOV1[i]
c_IOV1_out[] <- c_IOV1_to_IOS[i] + c_IOV1_to_IOV2[i] + c_IOV1_to_ROV1[i] + c_death_inf_IOV1[i] + c_death_IOV1[i]

update(IOV1[1]) <- IOV1[1] + c_IOV1_in[1] - c_IOV1_out[1] - c_aging_IOV1[1]
update(IOV1[2]) <- IOV1[2] + c_IOV1_in[2] - c_IOV1_out[2] + c_aging_IOV1[1] - c_aging_IOV1[2]
update(IOV1[3]) <- IOV1[3] + c_IOV1_in[3] - c_IOV1_out[3] + c_aging_IOV1[2] - c_aging_IOV1[3]
update(IOV1[4]) <- IOV1[4] + c_IOV1_in[4] - c_IOV1_out[4] + c_aging_IOV1[3] - c_aging_IOV1[4]
update(IOV1[5]) <- IOV1[5] + c_IOV1_in[5] - c_IOV1_out[5] + c_aging_IOV1[4] - c_aging_IOV1[5]
update(IOV1[6]) <- IOV1[6] + c_IOV1_in[6] - c_IOV1_out[6] + c_aging_IOV1[5] - c_aging_IOV1[6]
update(IOV1[7]) <- IOV1[7] + c_IOV1_in[7] - c_IOV1_out[7] + c_aging_IOV1[6] - c_aging_IOV1[7]
update(IOV1[8]) <- IOV1[8] + c_IOV1_in[8] - c_IOV1_out[8] + c_aging_IOV1[7] - c_aging_IOV1[8]
update(IOV1[9]) <- IOV1[9] + c_IOV1_in[9] - c_IOV1_out[9] + c_aging_IOV1[8] - c_aging_IOV1[9]
update(IOV1[10]) <- IOV1[10] + c_IOV1_in[10] - c_IOV1_out[10] + c_aging_IOV1[9] - c_aging_IOV1[10]
update(IOV1[11]) <- IOV1[11] + c_IOV1_in[11] - c_IOV1_out[11] + c_aging_IOV1[10] - c_aging_IOV1[11]
update(IOV1[12]) <- IOV1[12] + c_IOV1_in[12] - c_IOV1_out[12] + c_aging_IOV1[11] - c_aging_IOV1[12]
update(IOV1[13]) <- IOV1[13] + c_IOV1_in[13] - c_IOV1_out[13] + c_aging_IOV1[12] - c_aging_IOV1[13]
update(IOV1[14]) <- IOV1[14] + c_IOV1_in[14] - c_IOV1_out[14] + c_aging_IOV1[13] - c_aging_IOV1[14]
update(IOV1[15]) <- IOV1[15] + c_IOV1_in[15] - c_IOV1_out[15] + c_aging_IOV1[14] - c_aging_IOV1[15]
update(IOV1[16]) <- IOV1[16] + c_IOV1_in[16] - c_IOV1_out[16] + c_aging_IOV1[15]

# ----- 2.41. IOV2 Compartment --------------------------------------------------

c_IOV2_in[] <- c_EOV2_to_IOV2[i] + c_IOV1_to_IOV2[i] + c_migration_IOV2[i]
c_IOV2_out[] <- c_IOV2_to_IOS[i] + c_IOV2_to_ROV2[i] + c_death_inf_IOV2[i] + c_death_IOV2[i]

update(IOV2[1]) <- IOV2[1] + c_IOV2_in[1] - c_IOV2_out[1] - c_aging_IOV2[1]
update(IOV2[2]) <- IOV2[2] + c_IOV2_in[2] - c_IOV2_out[2] + c_aging_IOV2[1] - c_aging_IOV2[2]
update(IOV2[3]) <- IOV2[3] + c_IOV2_in[3] - c_IOV2_out[3] + c_aging_IOV2[2] - c_aging_IOV2[3]
update(IOV2[4]) <- IOV2[4] + c_IOV2_in[4] - c_IOV2_out[4] + c_aging_IOV2[3] - c_aging_IOV2[4]
update(IOV2[5]) <- IOV2[5] + c_IOV2_in[5] - c_IOV2_out[5] + c_aging_IOV2[4] - c_aging_IOV2[5]
update(IOV2[6]) <- IOV2[6] + c_IOV2_in[6] - c_IOV2_out[6] + c_aging_IOV2[5] - c_aging_IOV2[6]
update(IOV2[7]) <- IOV2[7] + c_IOV2_in[7] - c_IOV2_out[7] + c_aging_IOV2[6] - c_aging_IOV2[7]
update(IOV2[8]) <- IOV2[8] + c_IOV2_in[8] - c_IOV2_out[8] + c_aging_IOV2[7] - c_aging_IOV2[8]
update(IOV2[9]) <- IOV2[9] + c_IOV2_in[9] - c_IOV2_out[9] + c_aging_IOV2[8] - c_aging_IOV2[9]
update(IOV2[10]) <- IOV2[10] + c_IOV2_in[10] - c_IOV2_out[10] + c_aging_IOV2[9] - c_aging_IOV2[10]
update(IOV2[11]) <- IOV2[11] + c_IOV2_in[11] - c_IOV2_out[11] + c_aging_IOV2[10] - c_aging_IOV2[11]
update(IOV2[12]) <- IOV2[12] + c_IOV2_in[12] - c_IOV2_out[12] + c_aging_IOV2[11] - c_aging_IOV2[12]
update(IOV2[13]) <- IOV2[13] + c_IOV2_in[13] - c_IOV2_out[13] + c_aging_IOV2[12] - c_aging_IOV2[13]
update(IOV2[14]) <- IOV2[14] + c_IOV2_in[14] - c_IOV2_out[14] + c_aging_IOV2[13] - c_aging_IOV2[14]
update(IOV2[15]) <- IOV2[15] + c_IOV2_in[15] - c_IOV2_out[15] + c_aging_IOV2[14] - c_aging_IOV2[15]
update(IOV2[16]) <- IOV2[16] + c_IOV2_in[16] - c_IOV2_out[16] + c_aging_IOV2[15]

# ----- 2.43. IXS Compartment --------------------------------------------------

c_IXS_in[] <- c_EXS_to_IXS[i] + c_IXV1_to_IXS[i] + c_IXV2_to_IXS[i] + c_migration_IXS[i]
c_IXS_out[] <- c_IXS_to_IXV1[i] + c_IXS_to_RXS[i] + c_death_inf_IXS[i] + c_death_IXS[i]

update(IXS[1]) <- IXS[1] + c_IXS_in[1] - c_IXS_out[1] - c_aging_IXS[1]
update(IXS[2]) <- IXS[2] + c_IXS_in[2] - c_IXS_out[2] + c_aging_IXS[1] - c_aging_IXS[2]
update(IXS[3]) <- IXS[3] + c_IXS_in[3] - c_IXS_out[3] + c_aging_IXS[2] - c_aging_IXS[3]
update(IXS[4]) <- IXS[4] + c_IXS_in[4] - c_IXS_out[4] + c_aging_IXS[3] - c_aging_IXS[4]
update(IXS[5]) <- IXS[5] + c_IXS_in[5] - c_IXS_out[5] + c_aging_IXS[4] - c_aging_IXS[5]
update(IXS[6]) <- IXS[6] + c_IXS_in[6] - c_IXS_out[6] + c_aging_IXS[5] - c_aging_IXS[6]
update(IXS[7]) <- IXS[7] + c_IXS_in[7] - c_IXS_out[7] + c_aging_IXS[6] - c_aging_IXS[7]
update(IXS[8]) <- IXS[8] + c_IXS_in[8] - c_IXS_out[8] + c_aging_IXS[7] - c_aging_IXS[8]
update(IXS[9]) <- IXS[9] + c_IXS_in[9] - c_IXS_out[9] + c_aging_IXS[8] - c_aging_IXS[9]
update(IXS[10]) <- IXS[10] + c_IXS_in[10] - c_IXS_out[10] + c_aging_IXS[9] - c_aging_IXS[10]
update(IXS[11]) <- IXS[11] + c_IXS_in[11] - c_IXS_out[11] + c_aging_IXS[10] - c_aging_IXS[11]
update(IXS[12]) <- IXS[12] + c_IXS_in[12] - c_IXS_out[12] + c_aging_IXS[11] - c_aging_IXS[12]
update(IXS[13]) <- IXS[13] + c_IXS_in[13] - c_IXS_out[13] + c_aging_IXS[12] - c_aging_IXS[13]
update(IXS[14]) <- IXS[14] + c_IXS_in[14] - c_IXS_out[14] + c_aging_IXS[13] - c_aging_IXS[14]
update(IXS[15]) <- IXS[15] + c_IXS_in[15] - c_IXS_out[15] + c_aging_IXS[14] - c_aging_IXS[15]
update(IXS[16]) <- IXS[16] + c_IXS_in[16] - c_IXS_out[16] + c_aging_IXS[15]

# ----- 2.44. IXV1 Compartment --------------------------------------------------

c_IXV1_in[] <- c_EXV1_to_IXV1[i] + c_IXS_to_IXV1[i] + c_migration_IXV1[i]
c_IXV1_out[] <- c_IXV1_to_IXS[i] + c_IXV1_to_IXV2[i] + c_IXV1_to_RXV1[i] + c_death_inf_IXV1[i] + c_death_IXV1[i]

update(IXV1[1]) <- IXV1[1] + c_IXV1_in[1] - c_IXV1_out[1] - c_aging_IXV1[1]
update(IXV1[2]) <- IXV1[2] + c_IXV1_in[2] - c_IXV1_out[2] + c_aging_IXV1[1] - c_aging_IXV1[2]
update(IXV1[3]) <- IXV1[3] + c_IXV1_in[3] - c_IXV1_out[3] + c_aging_IXV1[2] - c_aging_IXV1[3]
update(IXV1[4]) <- IXV1[4] + c_IXV1_in[4] - c_IXV1_out[4] + c_aging_IXV1[3] - c_aging_IXV1[4]
update(IXV1[5]) <- IXV1[5] + c_IXV1_in[5] - c_IXV1_out[5] + c_aging_IXV1[4] - c_aging_IXV1[5]
update(IXV1[6]) <- IXV1[6] + c_IXV1_in[6] - c_IXV1_out[6] + c_aging_IXV1[5] - c_aging_IXV1[6]
update(IXV1[7]) <- IXV1[7] + c_IXV1_in[7] - c_IXV1_out[7] + c_aging_IXV1[6] - c_aging_IXV1[7]
update(IXV1[8]) <- IXV1[8] + c_IXV1_in[8] - c_IXV1_out[8] + c_aging_IXV1[7] - c_aging_IXV1[8]
update(IXV1[9]) <- IXV1[9] + c_IXV1_in[9] - c_IXV1_out[9] + c_aging_IXV1[8] - c_aging_IXV1[9]
update(IXV1[10]) <- IXV1[10] + c_IXV1_in[10] - c_IXV1_out[10] + c_aging_IXV1[9] - c_aging_IXV1[10]
update(IXV1[11]) <- IXV1[11] + c_IXV1_in[11] - c_IXV1_out[11] + c_aging_IXV1[10] - c_aging_IXV1[11]
update(IXV1[12]) <- IXV1[12] + c_IXV1_in[12] - c_IXV1_out[12] + c_aging_IXV1[11] - c_aging_IXV1[12]
update(IXV1[13]) <- IXV1[13] + c_IXV1_in[13] - c_IXV1_out[13] + c_aging_IXV1[12] - c_aging_IXV1[13]
update(IXV1[14]) <- IXV1[14] + c_IXV1_in[14] - c_IXV1_out[14] + c_aging_IXV1[13] - c_aging_IXV1[14]
update(IXV1[15]) <- IXV1[15] + c_IXV1_in[15] - c_IXV1_out[15] + c_aging_IXV1[14] - c_aging_IXV1[15]
update(IXV1[16]) <- IXV1[16] + c_IXV1_in[16] - c_IXV1_out[16] + c_aging_IXV1[15]

# ----- 2.45. IXV2 Compartment --------------------------------------------------

c_IXV2_in[] <- c_EXV2_to_IXV2[i] + c_IXV1_to_IXV2[i] + c_migration_IXV2[i]
c_IXV2_out[] <- c_IXV2_to_IXS[i] + c_IXV2_to_RXV2[i] + c_death_inf_IXV2[i] + c_death_IXV2[i]

update(IXV2[1]) <- IXV2[1] + c_IXV2_in[1] - c_IXV2_out[1] - c_aging_IXV2[1]
update(IXV2[2]) <- IXV2[2] + c_IXV2_in[2] - c_IXV2_out[2] + c_aging_IXV2[1] - c_aging_IXV2[2]
update(IXV2[3]) <- IXV2[3] + c_IXV2_in[3] - c_IXV2_out[3] + c_aging_IXV2[2] - c_aging_IXV2[3]
update(IXV2[4]) <- IXV2[4] + c_IXV2_in[4] - c_IXV2_out[4] + c_aging_IXV2[3] - c_aging_IXV2[4]
update(IXV2[5]) <- IXV2[5] + c_IXV2_in[5] - c_IXV2_out[5] + c_aging_IXV2[4] - c_aging_IXV2[5]
update(IXV2[6]) <- IXV2[6] + c_IXV2_in[6] - c_IXV2_out[6] + c_aging_IXV2[5] - c_aging_IXV2[6]
update(IXV2[7]) <- IXV2[7] + c_IXV2_in[7] - c_IXV2_out[7] + c_aging_IXV2[6] - c_aging_IXV2[7]
update(IXV2[8]) <- IXV2[8] + c_IXV2_in[8] - c_IXV2_out[8] + c_aging_IXV2[7] - c_aging_IXV2[8]
update(IXV2[9]) <- IXV2[9] + c_IXV2_in[9] - c_IXV2_out[9] + c_aging_IXV2[8] - c_aging_IXV2[9]
update(IXV2[10]) <- IXV2[10] + c_IXV2_in[10] - c_IXV2_out[10] + c_aging_IXV2[9] - c_aging_IXV2[10]
update(IXV2[11]) <- IXV2[11] + c_IXV2_in[11] - c_IXV2_out[11] + c_aging_IXV2[10] - c_aging_IXV2[11]
update(IXV2[12]) <- IXV2[12] + c_IXV2_in[12] - c_IXV2_out[12] + c_aging_IXV2[11] - c_aging_IXV2[12]
update(IXV2[13]) <- IXV2[13] + c_IXV2_in[13] - c_IXV2_out[13] + c_aging_IXV2[12] - c_aging_IXV2[13]
update(IXV2[14]) <- IXV2[14] + c_IXV2_in[14] - c_IXV2_out[14] + c_aging_IXV2[13] - c_aging_IXV2[14]
update(IXV2[15]) <- IXV2[15] + c_IXV2_in[15] - c_IXV2_out[15] + c_aging_IXV2[14] - c_aging_IXV2[15]
update(IXV2[16]) <- IXV2[16] + c_IXV2_in[16] - c_IXV2_out[16] + c_aging_IXV2[15]

# ----- 2.47. RWS Compartment --------------------------------------------------

c_RWS_in[] <- c_IWS_to_RWS[i] + c_RWV1_to_RWS[i] + c_RWV2_to_RWS[i] + c_migration_RWS[i]
c_RWS_out[] <- c_RWS_to_EARWS[i] + c_RWS_to_EDRWS[i] + c_RWS_to_EORWS[i] + c_RWS_to_RWV1[i] + c_RWS_to_S[i] + c_death_RWS[i]

update(RWS[1]) <- RWS[1] + c_RWS_in[1] - c_RWS_out[1] - c_aging_RWS[1]
update(RWS[2]) <- RWS[2] + c_RWS_in[2] - c_RWS_out[2] + c_aging_RWS[1] - c_aging_RWS[2]
update(RWS[3]) <- RWS[3] + c_RWS_in[3] - c_RWS_out[3] + c_aging_RWS[2] - c_aging_RWS[3]
update(RWS[4]) <- RWS[4] + c_RWS_in[4] - c_RWS_out[4] + c_aging_RWS[3] - c_aging_RWS[4]
update(RWS[5]) <- RWS[5] + c_RWS_in[5] - c_RWS_out[5] + c_aging_RWS[4] - c_aging_RWS[5]
update(RWS[6]) <- RWS[6] + c_RWS_in[6] - c_RWS_out[6] + c_aging_RWS[5] - c_aging_RWS[6]
update(RWS[7]) <- RWS[7] + c_RWS_in[7] - c_RWS_out[7] + c_aging_RWS[6] - c_aging_RWS[7]
update(RWS[8]) <- RWS[8] + c_RWS_in[8] - c_RWS_out[8] + c_aging_RWS[7] - c_aging_RWS[8]
update(RWS[9]) <- RWS[9] + c_RWS_in[9] - c_RWS_out[9] + c_aging_RWS[8] - c_aging_RWS[9]
update(RWS[10]) <- RWS[10] + c_RWS_in[10] - c_RWS_out[10] + c_aging_RWS[9] - c_aging_RWS[10]
update(RWS[11]) <- RWS[11] + c_RWS_in[11] - c_RWS_out[11] + c_aging_RWS[10] - c_aging_RWS[11]
update(RWS[12]) <- RWS[12] + c_RWS_in[12] - c_RWS_out[12] + c_aging_RWS[11] - c_aging_RWS[12]
update(RWS[13]) <- RWS[13] + c_RWS_in[13] - c_RWS_out[13] + c_aging_RWS[12] - c_aging_RWS[13]
update(RWS[14]) <- RWS[14] + c_RWS_in[14] - c_RWS_out[14] + c_aging_RWS[13] - c_aging_RWS[14]
update(RWS[15]) <- RWS[15] + c_RWS_in[15] - c_RWS_out[15] + c_aging_RWS[14] - c_aging_RWS[15]
update(RWS[16]) <- RWS[16] + c_RWS_in[16] - c_RWS_out[16] + c_aging_RWS[15]

# ----- 2.48. RWV1 Compartment --------------------------------------------------

c_RWV1_in[] <- c_IWV1_to_RWV1[i] + c_RWS_to_RWV1[i] + c_migration_RWV1[i]
c_RWV1_out[] <- c_RWV1_to_EARWV1[i] + c_RWV1_to_EDRWV1[i] + c_RWV1_to_EORWV1[i] + c_RWV1_to_RWS[i] + c_RWV1_to_RWV2[i] + c_RWV1_to_V1[i] + c_death_RWV1[i]

update(RWV1[1]) <- RWV1[1] + c_RWV1_in[1] - c_RWV1_out[1] - c_aging_RWV1[1]
update(RWV1[2]) <- RWV1[2] + c_RWV1_in[2] - c_RWV1_out[2] + c_aging_RWV1[1] - c_aging_RWV1[2]
update(RWV1[3]) <- RWV1[3] + c_RWV1_in[3] - c_RWV1_out[3] + c_aging_RWV1[2] - c_aging_RWV1[3]
update(RWV1[4]) <- RWV1[4] + c_RWV1_in[4] - c_RWV1_out[4] + c_aging_RWV1[3] - c_aging_RWV1[4]
update(RWV1[5]) <- RWV1[5] + c_RWV1_in[5] - c_RWV1_out[5] + c_aging_RWV1[4] - c_aging_RWV1[5]
update(RWV1[6]) <- RWV1[6] + c_RWV1_in[6] - c_RWV1_out[6] + c_aging_RWV1[5] - c_aging_RWV1[6]
update(RWV1[7]) <- RWV1[7] + c_RWV1_in[7] - c_RWV1_out[7] + c_aging_RWV1[6] - c_aging_RWV1[7]
update(RWV1[8]) <- RWV1[8] + c_RWV1_in[8] - c_RWV1_out[8] + c_aging_RWV1[7] - c_aging_RWV1[8]
update(RWV1[9]) <- RWV1[9] + c_RWV1_in[9] - c_RWV1_out[9] + c_aging_RWV1[8] - c_aging_RWV1[9]
update(RWV1[10]) <- RWV1[10] + c_RWV1_in[10] - c_RWV1_out[10] + c_aging_RWV1[9] - c_aging_RWV1[10]
update(RWV1[11]) <- RWV1[11] + c_RWV1_in[11] - c_RWV1_out[11] + c_aging_RWV1[10] - c_aging_RWV1[11]
update(RWV1[12]) <- RWV1[12] + c_RWV1_in[12] - c_RWV1_out[12] + c_aging_RWV1[11] - c_aging_RWV1[12]
update(RWV1[13]) <- RWV1[13] + c_RWV1_in[13] - c_RWV1_out[13] + c_aging_RWV1[12] - c_aging_RWV1[13]
update(RWV1[14]) <- RWV1[14] + c_RWV1_in[14] - c_RWV1_out[14] + c_aging_RWV1[13] - c_aging_RWV1[14]
update(RWV1[15]) <- RWV1[15] + c_RWV1_in[15] - c_RWV1_out[15] + c_aging_RWV1[14] - c_aging_RWV1[15]
update(RWV1[16]) <- RWV1[16] + c_RWV1_in[16] - c_RWV1_out[16] + c_aging_RWV1[15]

# ----- 2.49. RWV2 Compartment --------------------------------------------------

c_RWV2_in[] <- c_IWV2_to_RWV2[i] + c_RWV1_to_RWV2[i] + c_migration_RWV2[i]
c_RWV2_out[] <- c_RWV2_to_EARWV2[i] + c_RWV2_to_EDRWV2[i] + c_RWV2_to_EORWV2[i] + c_RWV2_to_RWS[i] + c_RWV2_to_V2[i] + c_death_RWV2[i]

update(RWV2[1]) <- RWV2[1] + c_RWV2_in[1] - c_RWV2_out[1] - c_aging_RWV2[1]
update(RWV2[2]) <- RWV2[2] + c_RWV2_in[2] - c_RWV2_out[2] + c_aging_RWV2[1] - c_aging_RWV2[2]
update(RWV2[3]) <- RWV2[3] + c_RWV2_in[3] - c_RWV2_out[3] + c_aging_RWV2[2] - c_aging_RWV2[3]
update(RWV2[4]) <- RWV2[4] + c_RWV2_in[4] - c_RWV2_out[4] + c_aging_RWV2[3] - c_aging_RWV2[4]
update(RWV2[5]) <- RWV2[5] + c_RWV2_in[5] - c_RWV2_out[5] + c_aging_RWV2[4] - c_aging_RWV2[5]
update(RWV2[6]) <- RWV2[6] + c_RWV2_in[6] - c_RWV2_out[6] + c_aging_RWV2[5] - c_aging_RWV2[6]
update(RWV2[7]) <- RWV2[7] + c_RWV2_in[7] - c_RWV2_out[7] + c_aging_RWV2[6] - c_aging_RWV2[7]
update(RWV2[8]) <- RWV2[8] + c_RWV2_in[8] - c_RWV2_out[8] + c_aging_RWV2[7] - c_aging_RWV2[8]
update(RWV2[9]) <- RWV2[9] + c_RWV2_in[9] - c_RWV2_out[9] + c_aging_RWV2[8] - c_aging_RWV2[9]
update(RWV2[10]) <- RWV2[10] + c_RWV2_in[10] - c_RWV2_out[10] + c_aging_RWV2[9] - c_aging_RWV2[10]
update(RWV2[11]) <- RWV2[11] + c_RWV2_in[11] - c_RWV2_out[11] + c_aging_RWV2[10] - c_aging_RWV2[11]
update(RWV2[12]) <- RWV2[12] + c_RWV2_in[12] - c_RWV2_out[12] + c_aging_RWV2[11] - c_aging_RWV2[12]
update(RWV2[13]) <- RWV2[13] + c_RWV2_in[13] - c_RWV2_out[13] + c_aging_RWV2[12] - c_aging_RWV2[13]
update(RWV2[14]) <- RWV2[14] + c_RWV2_in[14] - c_RWV2_out[14] + c_aging_RWV2[13] - c_aging_RWV2[14]
update(RWV2[15]) <- RWV2[15] + c_RWV2_in[15] - c_RWV2_out[15] + c_aging_RWV2[14] - c_aging_RWV2[15]
update(RWV2[16]) <- RWV2[16] + c_RWV2_in[16] - c_RWV2_out[16] + c_aging_RWV2[15]

# ----- 2.51. RAS Compartment --------------------------------------------------

c_RAS_in[] <- c_IARWS_to_RAS[i] + c_IAS_to_RAS[i] + c_RAV1_to_RAS[i] + c_RAV2_to_RAS[i] + c_migration_RAS[i]
c_RAS_out[] <- c_RAS_to_EDRAS[i] + c_RAS_to_EORAS[i] + c_RAS_to_RAV1[i] + c_RAS_to_S[i] + c_death_RAS[i]

update(RAS[1]) <- RAS[1] + c_RAS_in[1] - c_RAS_out[1] - c_aging_RAS[1]
update(RAS[2]) <- RAS[2] + c_RAS_in[2] - c_RAS_out[2] + c_aging_RAS[1] - c_aging_RAS[2]
update(RAS[3]) <- RAS[3] + c_RAS_in[3] - c_RAS_out[3] + c_aging_RAS[2] - c_aging_RAS[3]
update(RAS[4]) <- RAS[4] + c_RAS_in[4] - c_RAS_out[4] + c_aging_RAS[3] - c_aging_RAS[4]
update(RAS[5]) <- RAS[5] + c_RAS_in[5] - c_RAS_out[5] + c_aging_RAS[4] - c_aging_RAS[5]
update(RAS[6]) <- RAS[6] + c_RAS_in[6] - c_RAS_out[6] + c_aging_RAS[5] - c_aging_RAS[6]
update(RAS[7]) <- RAS[7] + c_RAS_in[7] - c_RAS_out[7] + c_aging_RAS[6] - c_aging_RAS[7]
update(RAS[8]) <- RAS[8] + c_RAS_in[8] - c_RAS_out[8] + c_aging_RAS[7] - c_aging_RAS[8]
update(RAS[9]) <- RAS[9] + c_RAS_in[9] - c_RAS_out[9] + c_aging_RAS[8] - c_aging_RAS[9]
update(RAS[10]) <- RAS[10] + c_RAS_in[10] - c_RAS_out[10] + c_aging_RAS[9] - c_aging_RAS[10]
update(RAS[11]) <- RAS[11] + c_RAS_in[11] - c_RAS_out[11] + c_aging_RAS[10] - c_aging_RAS[11]
update(RAS[12]) <- RAS[12] + c_RAS_in[12] - c_RAS_out[12] + c_aging_RAS[11] - c_aging_RAS[12]
update(RAS[13]) <- RAS[13] + c_RAS_in[13] - c_RAS_out[13] + c_aging_RAS[12] - c_aging_RAS[13]
update(RAS[14]) <- RAS[14] + c_RAS_in[14] - c_RAS_out[14] + c_aging_RAS[13] - c_aging_RAS[14]
update(RAS[15]) <- RAS[15] + c_RAS_in[15] - c_RAS_out[15] + c_aging_RAS[14] - c_aging_RAS[15]
update(RAS[16]) <- RAS[16] + c_RAS_in[16] - c_RAS_out[16] + c_aging_RAS[15]

# ----- 2.52. RAV1 Compartment --------------------------------------------------

c_RAV1_in[] <- c_IARWV1_to_RAV1[i] + c_IAV1_to_RAV1[i] + c_RAS_to_RAV1[i] + c_migration_RAV1[i]
c_RAV1_out[] <- c_RAV1_to_EDRAV1[i] + c_RAV1_to_EORAV1[i] + c_RAV1_to_RAS[i] + c_RAV1_to_RAV2[i] + c_RAV1_to_V1[i] + c_death_RAV1[i]

update(RAV1[1]) <- RAV1[1] + c_RAV1_in[1] - c_RAV1_out[1] - c_aging_RAV1[1]
update(RAV1[2]) <- RAV1[2] + c_RAV1_in[2] - c_RAV1_out[2] + c_aging_RAV1[1] - c_aging_RAV1[2]
update(RAV1[3]) <- RAV1[3] + c_RAV1_in[3] - c_RAV1_out[3] + c_aging_RAV1[2] - c_aging_RAV1[3]
update(RAV1[4]) <- RAV1[4] + c_RAV1_in[4] - c_RAV1_out[4] + c_aging_RAV1[3] - c_aging_RAV1[4]
update(RAV1[5]) <- RAV1[5] + c_RAV1_in[5] - c_RAV1_out[5] + c_aging_RAV1[4] - c_aging_RAV1[5]
update(RAV1[6]) <- RAV1[6] + c_RAV1_in[6] - c_RAV1_out[6] + c_aging_RAV1[5] - c_aging_RAV1[6]
update(RAV1[7]) <- RAV1[7] + c_RAV1_in[7] - c_RAV1_out[7] + c_aging_RAV1[6] - c_aging_RAV1[7]
update(RAV1[8]) <- RAV1[8] + c_RAV1_in[8] - c_RAV1_out[8] + c_aging_RAV1[7] - c_aging_RAV1[8]
update(RAV1[9]) <- RAV1[9] + c_RAV1_in[9] - c_RAV1_out[9] + c_aging_RAV1[8] - c_aging_RAV1[9]
update(RAV1[10]) <- RAV1[10] + c_RAV1_in[10] - c_RAV1_out[10] + c_aging_RAV1[9] - c_aging_RAV1[10]
update(RAV1[11]) <- RAV1[11] + c_RAV1_in[11] - c_RAV1_out[11] + c_aging_RAV1[10] - c_aging_RAV1[11]
update(RAV1[12]) <- RAV1[12] + c_RAV1_in[12] - c_RAV1_out[12] + c_aging_RAV1[11] - c_aging_RAV1[12]
update(RAV1[13]) <- RAV1[13] + c_RAV1_in[13] - c_RAV1_out[13] + c_aging_RAV1[12] - c_aging_RAV1[13]
update(RAV1[14]) <- RAV1[14] + c_RAV1_in[14] - c_RAV1_out[14] + c_aging_RAV1[13] - c_aging_RAV1[14]
update(RAV1[15]) <- RAV1[15] + c_RAV1_in[15] - c_RAV1_out[15] + c_aging_RAV1[14] - c_aging_RAV1[15]
update(RAV1[16]) <- RAV1[16] + c_RAV1_in[16] - c_RAV1_out[16] + c_aging_RAV1[15]

# ----- 2.53. RAV2 Compartment --------------------------------------------------

c_RAV2_in[] <- c_IARWV2_to_RAV2[i] + c_IAV2_to_RAV2[i] + c_RAV1_to_RAV2[i] + c_migration_RAV2[i]
c_RAV2_out[] <- c_RAV2_to_EDRAV2[i] + c_RAV2_to_EORAV2[i] + c_RAV2_to_RAS[i] + c_RAV2_to_V2[i] + c_death_RAV2[i]

update(RAV2[1]) <- RAV2[1] + c_RAV2_in[1] - c_RAV2_out[1] - c_aging_RAV2[1]
update(RAV2[2]) <- RAV2[2] + c_RAV2_in[2] - c_RAV2_out[2] + c_aging_RAV2[1] - c_aging_RAV2[2]
update(RAV2[3]) <- RAV2[3] + c_RAV2_in[3] - c_RAV2_out[3] + c_aging_RAV2[2] - c_aging_RAV2[3]
update(RAV2[4]) <- RAV2[4] + c_RAV2_in[4] - c_RAV2_out[4] + c_aging_RAV2[3] - c_aging_RAV2[4]
update(RAV2[5]) <- RAV2[5] + c_RAV2_in[5] - c_RAV2_out[5] + c_aging_RAV2[4] - c_aging_RAV2[5]
update(RAV2[6]) <- RAV2[6] + c_RAV2_in[6] - c_RAV2_out[6] + c_aging_RAV2[5] - c_aging_RAV2[6]
update(RAV2[7]) <- RAV2[7] + c_RAV2_in[7] - c_RAV2_out[7] + c_aging_RAV2[6] - c_aging_RAV2[7]
update(RAV2[8]) <- RAV2[8] + c_RAV2_in[8] - c_RAV2_out[8] + c_aging_RAV2[7] - c_aging_RAV2[8]
update(RAV2[9]) <- RAV2[9] + c_RAV2_in[9] - c_RAV2_out[9] + c_aging_RAV2[8] - c_aging_RAV2[9]
update(RAV2[10]) <- RAV2[10] + c_RAV2_in[10] - c_RAV2_out[10] + c_aging_RAV2[9] - c_aging_RAV2[10]
update(RAV2[11]) <- RAV2[11] + c_RAV2_in[11] - c_RAV2_out[11] + c_aging_RAV2[10] - c_aging_RAV2[11]
update(RAV2[12]) <- RAV2[12] + c_RAV2_in[12] - c_RAV2_out[12] + c_aging_RAV2[11] - c_aging_RAV2[12]
update(RAV2[13]) <- RAV2[13] + c_RAV2_in[13] - c_RAV2_out[13] + c_aging_RAV2[12] - c_aging_RAV2[13]
update(RAV2[14]) <- RAV2[14] + c_RAV2_in[14] - c_RAV2_out[14] + c_aging_RAV2[13] - c_aging_RAV2[14]
update(RAV2[15]) <- RAV2[15] + c_RAV2_in[15] - c_RAV2_out[15] + c_aging_RAV2[14] - c_aging_RAV2[15]
update(RAV2[16]) <- RAV2[16] + c_RAV2_in[16] - c_RAV2_out[16] + c_aging_RAV2[15]

# ----- 2.54. RAV2 Compartment --------------------------------------------------

c_RAV2_in[] <- c_IARWV2_to_RAV2[i] + c_IAV2_to_RAV2[i] + c_RAV1_to_RAV2[i] + c_migration_RAV2[i]
c_RAV2_out[] <- c_RAV2_to_EDRAV2[i] + c_RAV2_to_EORAV2[i] + c_RAV2_to_RAS[i] + c_RAV2_to_V2[i] + c_death_RAV2[i]

update(RAV2[1]) <- RAV2[1] + c_RAV2_in[1] - c_RAV2_out[1] - c_aging_RAV2[1]
update(RAV2[2]) <- RAV2[2] + c_RAV2_in[2] - c_RAV2_out[2] + c_aging_RAV2[1] - c_aging_RAV2[2]
update(RAV2[3]) <- RAV2[3] + c_RAV2_in[3] - c_RAV2_out[3] + c_aging_RAV2[2] - c_aging_RAV2[3]
update(RAV2[4]) <- RAV2[4] + c_RAV2_in[4] - c_RAV2_out[4] + c_aging_RAV2[3] - c_aging_RAV2[4]
update(RAV2[5]) <- RAV2[5] + c_RAV2_in[5] - c_RAV2_out[5] + c_aging_RAV2[4] - c_aging_RAV2[5]
update(RAV2[6]) <- RAV2[6] + c_RAV2_in[6] - c_RAV2_out[6] + c_aging_RAV2[5] - c_aging_RAV2[6]
update(RAV2[7]) <- RAV2[7] + c_RAV2_in[7] - c_RAV2_out[7] + c_aging_RAV2[6] - c_aging_RAV2[7]
update(RAV2[8]) <- RAV2[8] + c_RAV2_in[8] - c_RAV2_out[8] + c_aging_RAV2[7] - c_aging_RAV2[8]
update(RAV2[9]) <- RAV2[9] + c_RAV2_in[9] - c_RAV2_out[9] + c_aging_RAV2[8] - c_aging_RAV2[9]
update(RAV2[10]) <- RAV2[10] + c_RAV2_in[10] - c_RAV2_out[10] + c_aging_RAV2[9] - c_aging_RAV2[10]
update(RAV2[11]) <- RAV2[11] + c_RAV2_in[11] - c_RAV2_out[11] + c_aging_RAV2[10] - c_aging_RAV2[11]
update(RAV2[12]) <- RAV2[12] + c_RAV2_in[12] - c_RAV2_out[12] + c_aging_RAV2[11] - c_aging_RAV2[12]
update(RAV2[13]) <- RAV2[13] + c_RAV2_in[13] - c_RAV2_out[13] + c_aging_RAV2[12] - c_aging_RAV2[13]
update(RAV2[14]) <- RAV2[14] + c_RAV2_in[14] - c_RAV2_out[14] + c_aging_RAV2[13] - c_aging_RAV2[14]
update(RAV2[15]) <- RAV2[15] + c_RAV2_in[15] - c_RAV2_out[15] + c_aging_RAV2[14] - c_aging_RAV2[15]
update(RAV2[16]) <- RAV2[16] + c_RAV2_in[16] - c_RAV2_out[16] + c_aging_RAV2[15]

# ----- 2.56. RDS Compartment --------------------------------------------------

c_RDS_in[] <- c_IDRAS_to_RDS[i] + c_IDRWS_to_RDS[i] + c_IDS_to_RDS[i] + c_RDV1_to_RDS[i] + c_RDV2_to_RDS[i] + c_migration_RDS[i]
c_RDS_out[] <- c_RDS_to_EORDS[i] + c_RDS_to_RDV1[i] + c_RDS_to_S[i] + c_death_RDS[i]

update(RDS[1]) <- RDS[1] + c_RDS_in[1] - c_RDS_out[1] - c_aging_RDS[1]
update(RDS[2]) <- RDS[2] + c_RDS_in[2] - c_RDS_out[2] + c_aging_RDS[1] - c_aging_RDS[2]
update(RDS[3]) <- RDS[3] + c_RDS_in[3] - c_RDS_out[3] + c_aging_RDS[2] - c_aging_RDS[3]
update(RDS[4]) <- RDS[4] + c_RDS_in[4] - c_RDS_out[4] + c_aging_RDS[3] - c_aging_RDS[4]
update(RDS[5]) <- RDS[5] + c_RDS_in[5] - c_RDS_out[5] + c_aging_RDS[4] - c_aging_RDS[5]
update(RDS[6]) <- RDS[6] + c_RDS_in[6] - c_RDS_out[6] + c_aging_RDS[5] - c_aging_RDS[6]
update(RDS[7]) <- RDS[7] + c_RDS_in[7] - c_RDS_out[7] + c_aging_RDS[6] - c_aging_RDS[7]
update(RDS[8]) <- RDS[8] + c_RDS_in[8] - c_RDS_out[8] + c_aging_RDS[7] - c_aging_RDS[8]
update(RDS[9]) <- RDS[9] + c_RDS_in[9] - c_RDS_out[9] + c_aging_RDS[8] - c_aging_RDS[9]
update(RDS[10]) <- RDS[10] + c_RDS_in[10] - c_RDS_out[10] + c_aging_RDS[9] - c_aging_RDS[10]
update(RDS[11]) <- RDS[11] + c_RDS_in[11] - c_RDS_out[11] + c_aging_RDS[10] - c_aging_RDS[11]
update(RDS[12]) <- RDS[12] + c_RDS_in[12] - c_RDS_out[12] + c_aging_RDS[11] - c_aging_RDS[12]
update(RDS[13]) <- RDS[13] + c_RDS_in[13] - c_RDS_out[13] + c_aging_RDS[12] - c_aging_RDS[13]
update(RDS[14]) <- RDS[14] + c_RDS_in[14] - c_RDS_out[14] + c_aging_RDS[13] - c_aging_RDS[14]
update(RDS[15]) <- RDS[15] + c_RDS_in[15] - c_RDS_out[15] + c_aging_RDS[14] - c_aging_RDS[15]
update(RDS[16]) <- RDS[16] + c_RDS_in[16] - c_RDS_out[16] + c_aging_RDS[15]

# ----- 2.57. RDV1 Compartment --------------------------------------------------

c_RDV1_in[] <- c_IDRAV1_to_RDV1[i] + c_IDRWV1_to_RDV1[i] + c_IDV1_to_RDV1[i] + c_RDS_to_RDV1[i] + c_migration_RDV1[i]
c_RDV1_out[] <- c_RDV1_to_EORDV1[i] + c_RDV1_to_RDS[i] + c_RDV1_to_RDV2[i] + c_RDV1_to_V1[i] + c_death_RDV1[i]

update(RDV1[1]) <- RDV1[1] + c_RDV1_in[1] - c_RDV1_out[1] - c_aging_RDV1[1]
update(RDV1[2]) <- RDV1[2] + c_RDV1_in[2] - c_RDV1_out[2] + c_aging_RDV1[1] - c_aging_RDV1[2]
update(RDV1[3]) <- RDV1[3] + c_RDV1_in[3] - c_RDV1_out[3] + c_aging_RDV1[2] - c_aging_RDV1[3]
update(RDV1[4]) <- RDV1[4] + c_RDV1_in[4] - c_RDV1_out[4] + c_aging_RDV1[3] - c_aging_RDV1[4]
update(RDV1[5]) <- RDV1[5] + c_RDV1_in[5] - c_RDV1_out[5] + c_aging_RDV1[4] - c_aging_RDV1[5]
update(RDV1[6]) <- RDV1[6] + c_RDV1_in[6] - c_RDV1_out[6] + c_aging_RDV1[5] - c_aging_RDV1[6]
update(RDV1[7]) <- RDV1[7] + c_RDV1_in[7] - c_RDV1_out[7] + c_aging_RDV1[6] - c_aging_RDV1[7]
update(RDV1[8]) <- RDV1[8] + c_RDV1_in[8] - c_RDV1_out[8] + c_aging_RDV1[7] - c_aging_RDV1[8]
update(RDV1[9]) <- RDV1[9] + c_RDV1_in[9] - c_RDV1_out[9] + c_aging_RDV1[8] - c_aging_RDV1[9]
update(RDV1[10]) <- RDV1[10] + c_RDV1_in[10] - c_RDV1_out[10] + c_aging_RDV1[9] - c_aging_RDV1[10]
update(RDV1[11]) <- RDV1[11] + c_RDV1_in[11] - c_RDV1_out[11] + c_aging_RDV1[10] - c_aging_RDV1[11]
update(RDV1[12]) <- RDV1[12] + c_RDV1_in[12] - c_RDV1_out[12] + c_aging_RDV1[11] - c_aging_RDV1[12]
update(RDV1[13]) <- RDV1[13] + c_RDV1_in[13] - c_RDV1_out[13] + c_aging_RDV1[12] - c_aging_RDV1[13]
update(RDV1[14]) <- RDV1[14] + c_RDV1_in[14] - c_RDV1_out[14] + c_aging_RDV1[13] - c_aging_RDV1[14]
update(RDV1[15]) <- RDV1[15] + c_RDV1_in[15] - c_RDV1_out[15] + c_aging_RDV1[14] - c_aging_RDV1[15]
update(RDV1[16]) <- RDV1[16] + c_RDV1_in[16] - c_RDV1_out[16] + c_aging_RDV1[15]

# ----- 2.58. RDV2 Compartment --------------------------------------------------

c_RDV2_in[] <- c_IDRAV2_to_RDV2[i] + c_IDRWV2_to_RDV2[i] + c_IDV2_to_RDV2[i] + c_RDV1_to_RDV2[i] + c_migration_RDV2[i]
c_RDV2_out[] <- c_RDV2_to_EORDV2[i] + c_RDV2_to_RDS[i] + c_RDV2_to_V2[i] + c_death_RDV2[i]

update(RDV2[1]) <- RDV2[1] + c_RDV2_in[1] - c_RDV2_out[1] - c_aging_RDV2[1]
update(RDV2[2]) <- RDV2[2] + c_RDV2_in[2] - c_RDV2_out[2] + c_aging_RDV2[1] - c_aging_RDV2[2]
update(RDV2[3]) <- RDV2[3] + c_RDV2_in[3] - c_RDV2_out[3] + c_aging_RDV2[2] - c_aging_RDV2[3]
update(RDV2[4]) <- RDV2[4] + c_RDV2_in[4] - c_RDV2_out[4] + c_aging_RDV2[3] - c_aging_RDV2[4]
update(RDV2[5]) <- RDV2[5] + c_RDV2_in[5] - c_RDV2_out[5] + c_aging_RDV2[4] - c_aging_RDV2[5]
update(RDV2[6]) <- RDV2[6] + c_RDV2_in[6] - c_RDV2_out[6] + c_aging_RDV2[5] - c_aging_RDV2[6]
update(RDV2[7]) <- RDV2[7] + c_RDV2_in[7] - c_RDV2_out[7] + c_aging_RDV2[6] - c_aging_RDV2[7]
update(RDV2[8]) <- RDV2[8] + c_RDV2_in[8] - c_RDV2_out[8] + c_aging_RDV2[7] - c_aging_RDV2[8]
update(RDV2[9]) <- RDV2[9] + c_RDV2_in[9] - c_RDV2_out[9] + c_aging_RDV2[8] - c_aging_RDV2[9]
update(RDV2[10]) <- RDV2[10] + c_RDV2_in[10] - c_RDV2_out[10] + c_aging_RDV2[9] - c_aging_RDV2[10]
update(RDV2[11]) <- RDV2[11] + c_RDV2_in[11] - c_RDV2_out[11] + c_aging_RDV2[10] - c_aging_RDV2[11]
update(RDV2[12]) <- RDV2[12] + c_RDV2_in[12] - c_RDV2_out[12] + c_aging_RDV2[11] - c_aging_RDV2[12]
update(RDV2[13]) <- RDV2[13] + c_RDV2_in[13] - c_RDV2_out[13] + c_aging_RDV2[12] - c_aging_RDV2[13]
update(RDV2[14]) <- RDV2[14] + c_RDV2_in[14] - c_RDV2_out[14] + c_aging_RDV2[13] - c_aging_RDV2[14]
update(RDV2[15]) <- RDV2[15] + c_RDV2_in[15] - c_RDV2_out[15] + c_aging_RDV2[14] - c_aging_RDV2[15]
update(RDV2[16]) <- RDV2[16] + c_RDV2_in[16] - c_RDV2_out[16] + c_aging_RDV2[15]

# ----- 2.60. ROS Compartment --------------------------------------------------

c_ROS_in[] <- c_IORAS_to_ROS[i] + c_IORDS_to_ROS[i] + c_IORWS_to_ROS[i] + c_IOS_to_ROS[i] + c_ROV1_to_ROS[i] + c_ROV2_to_ROS[i] + c_migration_ROS[i]
c_ROS_out[] <- c_ROS_to_EXROS[i] + c_ROS_to_ROV1[i] + c_ROS_to_S[i] + c_death_ROS[i]

update(ROS[1]) <- ROS[1] + c_ROS_in[1] - c_ROS_out[1] - c_aging_ROS[1]
update(ROS[2]) <- ROS[2] + c_ROS_in[2] - c_ROS_out[2] + c_aging_ROS[1] - c_aging_ROS[2]
update(ROS[3]) <- ROS[3] + c_ROS_in[3] - c_ROS_out[3] + c_aging_ROS[2] - c_aging_ROS[3]
update(ROS[4]) <- ROS[4] + c_ROS_in[4] - c_ROS_out[4] + c_aging_ROS[3] - c_aging_ROS[4]
update(ROS[5]) <- ROS[5] + c_ROS_in[5] - c_ROS_out[5] + c_aging_ROS[4] - c_aging_ROS[5]
update(ROS[6]) <- ROS[6] + c_ROS_in[6] - c_ROS_out[6] + c_aging_ROS[5] - c_aging_ROS[6]
update(ROS[7]) <- ROS[7] + c_ROS_in[7] - c_ROS_out[7] + c_aging_ROS[6] - c_aging_ROS[7]
update(ROS[8]) <- ROS[8] + c_ROS_in[8] - c_ROS_out[8] + c_aging_ROS[7] - c_aging_ROS[8]
update(ROS[9]) <- ROS[9] + c_ROS_in[9] - c_ROS_out[9] + c_aging_ROS[8] - c_aging_ROS[9]
update(ROS[10]) <- ROS[10] + c_ROS_in[10] - c_ROS_out[10] + c_aging_ROS[9] - c_aging_ROS[10]
update(ROS[11]) <- ROS[11] + c_ROS_in[11] - c_ROS_out[11] + c_aging_ROS[10] - c_aging_ROS[11]
update(ROS[12]) <- ROS[12] + c_ROS_in[12] - c_ROS_out[12] + c_aging_ROS[11] - c_aging_ROS[12]
update(ROS[13]) <- ROS[13] + c_ROS_in[13] - c_ROS_out[13] + c_aging_ROS[12] - c_aging_ROS[13]
update(ROS[14]) <- ROS[14] + c_ROS_in[14] - c_ROS_out[14] + c_aging_ROS[13] - c_aging_ROS[14]
update(ROS[15]) <- ROS[15] + c_ROS_in[15] - c_ROS_out[15] + c_aging_ROS[14] - c_aging_ROS[15]
update(ROS[16]) <- ROS[16] + c_ROS_in[16] - c_ROS_out[16] + c_aging_ROS[15]

# ----- 2.61. ROV1 Compartment --------------------------------------------------

c_ROV1_in[] <- c_IORAV1_to_ROV1[i] + c_IORDV1_to_ROV1[i] + c_IORWV1_to_ROV1[i] + c_IOV1_to_ROV1[i] + c_ROS_to_ROV1[i] + c_migration_ROV1[i]
c_ROV1_out[] <- c_ROV1_to_EXROV1[i] + c_ROV1_to_ROS[i] + c_ROV1_to_ROV2[i] + c_ROV1_to_V1[i] + c_death_ROV1[i]

update(ROV1[1]) <- ROV1[1] + c_ROV1_in[1] - c_ROV1_out[1] - c_aging_ROV1[1]
update(ROV1[2]) <- ROV1[2] + c_ROV1_in[2] - c_ROV1_out[2] + c_aging_ROV1[1] - c_aging_ROV1[2]
update(ROV1[3]) <- ROV1[3] + c_ROV1_in[3] - c_ROV1_out[3] + c_aging_ROV1[2] - c_aging_ROV1[3]
update(ROV1[4]) <- ROV1[4] + c_ROV1_in[4] - c_ROV1_out[4] + c_aging_ROV1[3] - c_aging_ROV1[4]
update(ROV1[5]) <- ROV1[5] + c_ROV1_in[5] - c_ROV1_out[5] + c_aging_ROV1[4] - c_aging_ROV1[5]
update(ROV1[6]) <- ROV1[6] + c_ROV1_in[6] - c_ROV1_out[6] + c_aging_ROV1[5] - c_aging_ROV1[6]
update(ROV1[7]) <- ROV1[7] + c_ROV1_in[7] - c_ROV1_out[7] + c_aging_ROV1[6] - c_aging_ROV1[7]
update(ROV1[8]) <- ROV1[8] + c_ROV1_in[8] - c_ROV1_out[8] + c_aging_ROV1[7] - c_aging_ROV1[8]
update(ROV1[9]) <- ROV1[9] + c_ROV1_in[9] - c_ROV1_out[9] + c_aging_ROV1[8] - c_aging_ROV1[9]
update(ROV1[10]) <- ROV1[10] + c_ROV1_in[10] - c_ROV1_out[10] + c_aging_ROV1[9] - c_aging_ROV1[10]
update(ROV1[11]) <- ROV1[11] + c_ROV1_in[11] - c_ROV1_out[11] + c_aging_ROV1[10] - c_aging_ROV1[11]
update(ROV1[12]) <- ROV1[12] + c_ROV1_in[12] - c_ROV1_out[12] + c_aging_ROV1[11] - c_aging_ROV1[12]
update(ROV1[13]) <- ROV1[13] + c_ROV1_in[13] - c_ROV1_out[13] + c_aging_ROV1[12] - c_aging_ROV1[13]
update(ROV1[14]) <- ROV1[14] + c_ROV1_in[14] - c_ROV1_out[14] + c_aging_ROV1[13] - c_aging_ROV1[14]
update(ROV1[15]) <- ROV1[15] + c_ROV1_in[15] - c_ROV1_out[15] + c_aging_ROV1[14] - c_aging_ROV1[15]
update(ROV1[16]) <- ROV1[16] + c_ROV1_in[16] - c_ROV1_out[16] + c_aging_ROV1[15]

# ----- 2.62. ROV2 Compartment --------------------------------------------------

c_ROV2_in[] <- c_IORAV2_to_ROV2[i] + c_IORDV2_to_ROV2[i] + c_IORWV2_to_ROV2[i] + c_IOV2_to_ROV2[i] + c_ROV1_to_ROV2[i] + c_migration_ROV2[i]
c_ROV2_out[] <- c_ROV2_to_EXROV2[i] + c_ROV2_to_ROS[i] + c_ROV2_to_V2[i] + c_death_ROV2[i]

update(ROV2[1]) <- ROV2[1] + c_ROV2_in[1] - c_ROV2_out[1] - c_aging_ROV2[1]
update(ROV2[2]) <- ROV2[2] + c_ROV2_in[2] - c_ROV2_out[2] + c_aging_ROV2[1] - c_aging_ROV2[2]
update(ROV2[3]) <- ROV2[3] + c_ROV2_in[3] - c_ROV2_out[3] + c_aging_ROV2[2] - c_aging_ROV2[3]
update(ROV2[4]) <- ROV2[4] + c_ROV2_in[4] - c_ROV2_out[4] + c_aging_ROV2[3] - c_aging_ROV2[4]
update(ROV2[5]) <- ROV2[5] + c_ROV2_in[5] - c_ROV2_out[5] + c_aging_ROV2[4] - c_aging_ROV2[5]
update(ROV2[6]) <- ROV2[6] + c_ROV2_in[6] - c_ROV2_out[6] + c_aging_ROV2[5] - c_aging_ROV2[6]
update(ROV2[7]) <- ROV2[7] + c_ROV2_in[7] - c_ROV2_out[7] + c_aging_ROV2[6] - c_aging_ROV2[7]
update(ROV2[8]) <- ROV2[8] + c_ROV2_in[8] - c_ROV2_out[8] + c_aging_ROV2[7] - c_aging_ROV2[8]
update(ROV2[9]) <- ROV2[9] + c_ROV2_in[9] - c_ROV2_out[9] + c_aging_ROV2[8] - c_aging_ROV2[9]
update(ROV2[10]) <- ROV2[10] + c_ROV2_in[10] - c_ROV2_out[10] + c_aging_ROV2[9] - c_aging_ROV2[10]
update(ROV2[11]) <- ROV2[11] + c_ROV2_in[11] - c_ROV2_out[11] + c_aging_ROV2[10] - c_aging_ROV2[11]
update(ROV2[12]) <- ROV2[12] + c_ROV2_in[12] - c_ROV2_out[12] + c_aging_ROV2[11] - c_aging_ROV2[12]
update(ROV2[13]) <- ROV2[13] + c_ROV2_in[13] - c_ROV2_out[13] + c_aging_ROV2[12] - c_aging_ROV2[13]
update(ROV2[14]) <- ROV2[14] + c_ROV2_in[14] - c_ROV2_out[14] + c_aging_ROV2[13] - c_aging_ROV2[14]
update(ROV2[15]) <- ROV2[15] + c_ROV2_in[15] - c_ROV2_out[15] + c_aging_ROV2[14] - c_aging_ROV2[15]
update(ROV2[16]) <- ROV2[16] + c_ROV2_in[16] - c_ROV2_out[16] + c_aging_ROV2[15]

# ----- 2.64. RXS Compartment --------------------------------------------------

c_RXS_in[] <- c_IXS_to_RXS[i] + c_RXV1_to_RXS[i] + c_RXV2_to_RXS[i] + c_migration_RXS[i]
c_RXS_out[] <- c_RXS_to_EORXS[i] + c_RXS_to_RXV1[i] + c_RXS_to_S[i] + c_death_RXS[i]

update(RXS[1]) <- RXS[1] + c_RXS_in[1] - c_RXS_out[1] - c_aging_RXS[1]
update(RXS[2]) <- RXS[2] + c_RXS_in[2] - c_RXS_out[2] + c_aging_RXS[1] - c_aging_RXS[2]
update(RXS[3]) <- RXS[3] + c_RXS_in[3] - c_RXS_out[3] + c_aging_RXS[2] - c_aging_RXS[3]
update(RXS[4]) <- RXS[4] + c_RXS_in[4] - c_RXS_out[4] + c_aging_RXS[3] - c_aging_RXS[4]
update(RXS[5]) <- RXS[5] + c_RXS_in[5] - c_RXS_out[5] + c_aging_RXS[4] - c_aging_RXS[5]
update(RXS[6]) <- RXS[6] + c_RXS_in[6] - c_RXS_out[6] + c_aging_RXS[5] - c_aging_RXS[6]
update(RXS[7]) <- RXS[7] + c_RXS_in[7] - c_RXS_out[7] + c_aging_RXS[6] - c_aging_RXS[7]
update(RXS[8]) <- RXS[8] + c_RXS_in[8] - c_RXS_out[8] + c_aging_RXS[7] - c_aging_RXS[8]
update(RXS[9]) <- RXS[9] + c_RXS_in[9] - c_RXS_out[9] + c_aging_RXS[8] - c_aging_RXS[9]
update(RXS[10]) <- RXS[10] + c_RXS_in[10] - c_RXS_out[10] + c_aging_RXS[9] - c_aging_RXS[10]
update(RXS[11]) <- RXS[11] + c_RXS_in[11] - c_RXS_out[11] + c_aging_RXS[10] - c_aging_RXS[11]
update(RXS[12]) <- RXS[12] + c_RXS_in[12] - c_RXS_out[12] + c_aging_RXS[11] - c_aging_RXS[12]
update(RXS[13]) <- RXS[13] + c_RXS_in[13] - c_RXS_out[13] + c_aging_RXS[12] - c_aging_RXS[13]
update(RXS[14]) <- RXS[14] + c_RXS_in[14] - c_RXS_out[14] + c_aging_RXS[13] - c_aging_RXS[14]
update(RXS[15]) <- RXS[15] + c_RXS_in[15] - c_RXS_out[15] + c_aging_RXS[14] - c_aging_RXS[15]
update(RXS[16]) <- RXS[16] + c_RXS_in[16] - c_RXS_out[16] + c_aging_RXS[15]

# ----- 2.65. RXV1 Compartment --------------------------------------------------

c_RXV1_in[] <- c_IXV1_to_RXV1[i] + c_RXS_to_RXV1[i] + c_migration_RXV1[i]
c_RXV1_out[] <- c_RXV1_to_EORXV1[i] + c_RXV1_to_RXS[i] + c_RXV1_to_RXV2[i] + c_RXV1_to_V1[i] + c_death_RXV1[i]

update(RXV1[1]) <- RXV1[1] + c_RXV1_in[1] - c_RXV1_out[1] - c_aging_RXV1[1]
update(RXV1[2]) <- RXV1[2] + c_RXV1_in[2] - c_RXV1_out[2] + c_aging_RXV1[1] - c_aging_RXV1[2]
update(RXV1[3]) <- RXV1[3] + c_RXV1_in[3] - c_RXV1_out[3] + c_aging_RXV1[2] - c_aging_RXV1[3]
update(RXV1[4]) <- RXV1[4] + c_RXV1_in[4] - c_RXV1_out[4] + c_aging_RXV1[3] - c_aging_RXV1[4]
update(RXV1[5]) <- RXV1[5] + c_RXV1_in[5] - c_RXV1_out[5] + c_aging_RXV1[4] - c_aging_RXV1[5]
update(RXV1[6]) <- RXV1[6] + c_RXV1_in[6] - c_RXV1_out[6] + c_aging_RXV1[5] - c_aging_RXV1[6]
update(RXV1[7]) <- RXV1[7] + c_RXV1_in[7] - c_RXV1_out[7] + c_aging_RXV1[6] - c_aging_RXV1[7]
update(RXV1[8]) <- RXV1[8] + c_RXV1_in[8] - c_RXV1_out[8] + c_aging_RXV1[7] - c_aging_RXV1[8]
update(RXV1[9]) <- RXV1[9] + c_RXV1_in[9] - c_RXV1_out[9] + c_aging_RXV1[8] - c_aging_RXV1[9]
update(RXV1[10]) <- RXV1[10] + c_RXV1_in[10] - c_RXV1_out[10] + c_aging_RXV1[9] - c_aging_RXV1[10]
update(RXV1[11]) <- RXV1[11] + c_RXV1_in[11] - c_RXV1_out[11] + c_aging_RXV1[10] - c_aging_RXV1[11]
update(RXV1[12]) <- RXV1[12] + c_RXV1_in[12] - c_RXV1_out[12] + c_aging_RXV1[11] - c_aging_RXV1[12]
update(RXV1[13]) <- RXV1[13] + c_RXV1_in[13] - c_RXV1_out[13] + c_aging_RXV1[12] - c_aging_RXV1[13]
update(RXV1[14]) <- RXV1[14] + c_RXV1_in[14] - c_RXV1_out[14] + c_aging_RXV1[13] - c_aging_RXV1[14]
update(RXV1[15]) <- RXV1[15] + c_RXV1_in[15] - c_RXV1_out[15] + c_aging_RXV1[14] - c_aging_RXV1[15]
update(RXV1[16]) <- RXV1[16] + c_RXV1_in[16] - c_RXV1_out[16] + c_aging_RXV1[15]

# ----- 2.66. RXV2 Compartment --------------------------------------------------

c_RXV2_in[] <- c_IXV2_to_RXV2[i] + c_RXV1_to_RXV2[i] + c_migration_RXV2[i]
c_RXV2_out[] <- c_RXV2_to_EORXV2[i] + c_RXV2_to_RXS[i] + c_RXV2_to_V2[i] + c_death_RXV2[i]

update(RXV2[1]) <- RXV2[1] + c_RXV2_in[1] - c_RXV2_out[1] - c_aging_RXV2[1]
update(RXV2[2]) <- RXV2[2] + c_RXV2_in[2] - c_RXV2_out[2] + c_aging_RXV2[1] - c_aging_RXV2[2]
update(RXV2[3]) <- RXV2[3] + c_RXV2_in[3] - c_RXV2_out[3] + c_aging_RXV2[2] - c_aging_RXV2[3]
update(RXV2[4]) <- RXV2[4] + c_RXV2_in[4] - c_RXV2_out[4] + c_aging_RXV2[3] - c_aging_RXV2[4]
update(RXV2[5]) <- RXV2[5] + c_RXV2_in[5] - c_RXV2_out[5] + c_aging_RXV2[4] - c_aging_RXV2[5]
update(RXV2[6]) <- RXV2[6] + c_RXV2_in[6] - c_RXV2_out[6] + c_aging_RXV2[5] - c_aging_RXV2[6]
update(RXV2[7]) <- RXV2[7] + c_RXV2_in[7] - c_RXV2_out[7] + c_aging_RXV2[6] - c_aging_RXV2[7]
update(RXV2[8]) <- RXV2[8] + c_RXV2_in[8] - c_RXV2_out[8] + c_aging_RXV2[7] - c_aging_RXV2[8]
update(RXV2[9]) <- RXV2[9] + c_RXV2_in[9] - c_RXV2_out[9] + c_aging_RXV2[8] - c_aging_RXV2[9]
update(RXV2[10]) <- RXV2[10] + c_RXV2_in[10] - c_RXV2_out[10] + c_aging_RXV2[9] - c_aging_RXV2[10]
update(RXV2[11]) <- RXV2[11] + c_RXV2_in[11] - c_RXV2_out[11] + c_aging_RXV2[10] - c_aging_RXV2[11]
update(RXV2[12]) <- RXV2[12] + c_RXV2_in[12] - c_RXV2_out[12] + c_aging_RXV2[11] - c_aging_RXV2[12]
update(RXV2[13]) <- RXV2[13] + c_RXV2_in[13] - c_RXV2_out[13] + c_aging_RXV2[12] - c_aging_RXV2[13]
update(RXV2[14]) <- RXV2[14] + c_RXV2_in[14] - c_RXV2_out[14] + c_aging_RXV2[13] - c_aging_RXV2[14]
update(RXV2[15]) <- RXV2[15] + c_RXV2_in[15] - c_RXV2_out[15] + c_aging_RXV2[14] - c_aging_RXV2[15]
update(RXV2[16]) <- RXV2[16] + c_RXV2_in[16] - c_RXV2_out[16] + c_aging_RXV2[15]

# ----- 2.68. EARWS Compartment -----------------------------------------------

c_EARWS_in[] <- c_EARWV1_to_EARWS[i] + c_EARWV2_to_EARWS[i] + c_RWS_to_EARWS[i] + c_migration_EARWS[i]
c_EARWS_out[] <- c_EARWS_to_EARWV1[i] + c_EARWS_to_IARWS[i] + c_death_EARWS[i]

update(EARWS[1]) <- EARWS[1] + c_EARWS_in[1] - c_EARWS_out[1] - c_aging_EARWS[1]
update(EARWS[2]) <- EARWS[2] + c_EARWS_in[2] - c_EARWS_out[2] + c_aging_EARWS[1] - c_aging_EARWS[2]
update(EARWS[3]) <- EARWS[3] + c_EARWS_in[3] - c_EARWS_out[3] + c_aging_EARWS[2] - c_aging_EARWS[3]
update(EARWS[4]) <- EARWS[4] + c_EARWS_in[4] - c_EARWS_out[4] + c_aging_EARWS[3] - c_aging_EARWS[4]
update(EARWS[5]) <- EARWS[5] + c_EARWS_in[5] - c_EARWS_out[5] + c_aging_EARWS[4] - c_aging_EARWS[5]
update(EARWS[6]) <- EARWS[6] + c_EARWS_in[6] - c_EARWS_out[6] + c_aging_EARWS[5] - c_aging_EARWS[6]
update(EARWS[7]) <- EARWS[7] + c_EARWS_in[7] - c_EARWS_out[7] + c_aging_EARWS[6] - c_aging_EARWS[7]
update(EARWS[8]) <- EARWS[8] + c_EARWS_in[8] - c_EARWS_out[8] + c_aging_EARWS[7] - c_aging_EARWS[8]
update(EARWS[9]) <- EARWS[9] + c_EARWS_in[9] - c_EARWS_out[9] + c_aging_EARWS[8] - c_aging_EARWS[9]
update(EARWS[10]) <- EARWS[10] + c_EARWS_in[10] - c_EARWS_out[10] + c_aging_EARWS[9] - c_aging_EARWS[10]
update(EARWS[11]) <- EARWS[11] + c_EARWS_in[11] - c_EARWS_out[11] + c_aging_EARWS[10] - c_aging_EARWS[11]
update(EARWS[12]) <- EARWS[12] + c_EARWS_in[12] - c_EARWS_out[12] + c_aging_EARWS[11] - c_aging_EARWS[12]
update(EARWS[13]) <- EARWS[13] + c_EARWS_in[13] - c_EARWS_out[13] + c_aging_EARWS[12] - c_aging_EARWS[13]
update(EARWS[14]) <- EARWS[14] + c_EARWS_in[14] - c_EARWS_out[14] + c_aging_EARWS[13] - c_aging_EARWS[14]
update(EARWS[15]) <- EARWS[15] + c_EARWS_in[15] - c_EARWS_out[15] + c_aging_EARWS[14] - c_aging_EARWS[15]
update(EARWS[16]) <- EARWS[16] + c_EARWS_in[16] - c_EARWS_out[16] + c_aging_EARWS[15]

# ----- 2.69. EARWV1 Compartment -----------------------------------------------

c_EARWV1_in[] <- c_EARWS_to_EARWV1[i] + c_RWV1_to_EARWV1[i] + c_migration_EARWV1[i]
c_EARWV1_out[] <- c_EARWV1_to_EARWS[i] + c_EARWV1_to_EARWV2[i] + c_EARWV1_to_IARWV1[i] + c_death_EARWV1[i]

update(EARWV1[1]) <- EARWV1[1] + c_EARWV1_in[1] - c_EARWV1_out[1] - c_aging_EARWV1[1]
update(EARWV1[2]) <- EARWV1[2] + c_EARWV1_in[2] - c_EARWV1_out[2] + c_aging_EARWV1[1] - c_aging_EARWV1[2]
update(EARWV1[3]) <- EARWV1[3] + c_EARWV1_in[3] - c_EARWV1_out[3] + c_aging_EARWV1[2] - c_aging_EARWV1[3]
update(EARWV1[4]) <- EARWV1[4] + c_EARWV1_in[4] - c_EARWV1_out[4] + c_aging_EARWV1[3] - c_aging_EARWV1[4]
update(EARWV1[5]) <- EARWV1[5] + c_EARWV1_in[5] - c_EARWV1_out[5] + c_aging_EARWV1[4] - c_aging_EARWV1[5]
update(EARWV1[6]) <- EARWV1[6] + c_EARWV1_in[6] - c_EARWV1_out[6] + c_aging_EARWV1[5] - c_aging_EARWV1[6]
update(EARWV1[7]) <- EARWV1[7] + c_EARWV1_in[7] - c_EARWV1_out[7] + c_aging_EARWV1[6] - c_aging_EARWV1[7]
update(EARWV1[8]) <- EARWV1[8] + c_EARWV1_in[8] - c_EARWV1_out[8] + c_aging_EARWV1[7] - c_aging_EARWV1[8]
update(EARWV1[9]) <- EARWV1[9] + c_EARWV1_in[9] - c_EARWV1_out[9] + c_aging_EARWV1[8] - c_aging_EARWV1[9]
update(EARWV1[10]) <- EARWV1[10] + c_EARWV1_in[10] - c_EARWV1_out[10] + c_aging_EARWV1[9] - c_aging_EARWV1[10]
update(EARWV1[11]) <- EARWV1[11] + c_EARWV1_in[11] - c_EARWV1_out[11] + c_aging_EARWV1[10] - c_aging_EARWV1[11]
update(EARWV1[12]) <- EARWV1[12] + c_EARWV1_in[12] - c_EARWV1_out[12] + c_aging_EARWV1[11] - c_aging_EARWV1[12]
update(EARWV1[13]) <- EARWV1[13] + c_EARWV1_in[13] - c_EARWV1_out[13] + c_aging_EARWV1[12] - c_aging_EARWV1[13]
update(EARWV1[14]) <- EARWV1[14] + c_EARWV1_in[14] - c_EARWV1_out[14] + c_aging_EARWV1[13] - c_aging_EARWV1[14]
update(EARWV1[15]) <- EARWV1[15] + c_EARWV1_in[15] - c_EARWV1_out[15] + c_aging_EARWV1[14] - c_aging_EARWV1[15]
update(EARWV1[16]) <- EARWV1[16] + c_EARWV1_in[16] - c_EARWV1_out[16] + c_aging_EARWV1[15]

# ----- 2.70. EARWV2 Compartment -----------------------------------------------

c_EARWV2_in[] <- c_EARWV1_to_EARWV2[i] + c_RWV2_to_EARWV2[i] + c_migration_EARWV2[i]
c_EARWV2_out[] <- c_EARWV2_to_EARWS[i] + c_EARWV2_to_IARWV2[i] + c_death_EARWV2[i]

update(EARWV2[1]) <- EARWV2[1] + c_EARWV2_in[1] - c_EARWV2_out[1] - c_aging_EARWV2[1]
update(EARWV2[2]) <- EARWV2[2] + c_EARWV2_in[2] - c_EARWV2_out[2] + c_aging_EARWV2[1] - c_aging_EARWV2[2]
update(EARWV2[3]) <- EARWV2[3] + c_EARWV2_in[3] - c_EARWV2_out[3] + c_aging_EARWV2[2] - c_aging_EARWV2[3]
update(EARWV2[4]) <- EARWV2[4] + c_EARWV2_in[4] - c_EARWV2_out[4] + c_aging_EARWV2[3] - c_aging_EARWV2[4]
update(EARWV2[5]) <- EARWV2[5] + c_EARWV2_in[5] - c_EARWV2_out[5] + c_aging_EARWV2[4] - c_aging_EARWV2[5]
update(EARWV2[6]) <- EARWV2[6] + c_EARWV2_in[6] - c_EARWV2_out[6] + c_aging_EARWV2[5] - c_aging_EARWV2[6]
update(EARWV2[7]) <- EARWV2[7] + c_EARWV2_in[7] - c_EARWV2_out[7] + c_aging_EARWV2[6] - c_aging_EARWV2[7]
update(EARWV2[8]) <- EARWV2[8] + c_EARWV2_in[8] - c_EARWV2_out[8] + c_aging_EARWV2[7] - c_aging_EARWV2[8]
update(EARWV2[9]) <- EARWV2[9] + c_EARWV2_in[9] - c_EARWV2_out[9] + c_aging_EARWV2[8] - c_aging_EARWV2[9]
update(EARWV2[10]) <- EARWV2[10] + c_EARWV2_in[10] - c_EARWV2_out[10] + c_aging_EARWV2[9] - c_aging_EARWV2[10]
update(EARWV2[11]) <- EARWV2[11] + c_EARWV2_in[11] - c_EARWV2_out[11] + c_aging_EARWV2[10] - c_aging_EARWV2[11]
update(EARWV2[12]) <- EARWV2[12] + c_EARWV2_in[12] - c_EARWV2_out[12] + c_aging_EARWV2[11] - c_aging_EARWV2[12]
update(EARWV2[13]) <- EARWV2[13] + c_EARWV2_in[13] - c_EARWV2_out[13] + c_aging_EARWV2[12] - c_aging_EARWV2[13]
update(EARWV2[14]) <- EARWV2[14] + c_EARWV2_in[14] - c_EARWV2_out[14] + c_aging_EARWV2[13] - c_aging_EARWV2[14]
update(EARWV2[15]) <- EARWV2[15] + c_EARWV2_in[15] - c_EARWV2_out[15] + c_aging_EARWV2[14] - c_aging_EARWV2[15]
update(EARWV2[16]) <- EARWV2[16] + c_EARWV2_in[16] - c_EARWV2_out[16] + c_aging_EARWV2[15]

# ----- 2.72. EDRWS Compartment -----------------------------------------------

c_EDRWS_in[] <- c_EDRWV1_to_EDRWS[i] + c_EDRWV2_to_EDRWS[i] + c_RWS_to_EDRWS[i] + c_migration_EDRWS[i]
c_EDRWS_out[] <- c_EDRWS_to_EDRWV1[i] + c_EDRWS_to_IDRWS[i] + c_death_EDRWS[i]

update(EDRWS[1]) <- EDRWS[1] + c_EDRWS_in[1] - c_EDRWS_out[1] - c_aging_EDRWS[1]
update(EDRWS[2]) <- EDRWS[2] + c_EDRWS_in[2] - c_EDRWS_out[2] + c_aging_EDRWS[1] - c_aging_EDRWS[2]
update(EDRWS[3]) <- EDRWS[3] + c_EDRWS_in[3] - c_EDRWS_out[3] + c_aging_EDRWS[2] - c_aging_EDRWS[3]
update(EDRWS[4]) <- EDRWS[4] + c_EDRWS_in[4] - c_EDRWS_out[4] + c_aging_EDRWS[3] - c_aging_EDRWS[4]
update(EDRWS[5]) <- EDRWS[5] + c_EDRWS_in[5] - c_EDRWS_out[5] + c_aging_EDRWS[4] - c_aging_EDRWS[5]
update(EDRWS[6]) <- EDRWS[6] + c_EDRWS_in[6] - c_EDRWS_out[6] + c_aging_EDRWS[5] - c_aging_EDRWS[6]
update(EDRWS[7]) <- EDRWS[7] + c_EDRWS_in[7] - c_EDRWS_out[7] + c_aging_EDRWS[6] - c_aging_EDRWS[7]
update(EDRWS[8]) <- EDRWS[8] + c_EDRWS_in[8] - c_EDRWS_out[8] + c_aging_EDRWS[7] - c_aging_EDRWS[8]
update(EDRWS[9]) <- EDRWS[9] + c_EDRWS_in[9] - c_EDRWS_out[9] + c_aging_EDRWS[8] - c_aging_EDRWS[9]
update(EDRWS[10]) <- EDRWS[10] + c_EDRWS_in[10] - c_EDRWS_out[10] + c_aging_EDRWS[9] - c_aging_EDRWS[10]
update(EDRWS[11]) <- EDRWS[11] + c_EDRWS_in[11] - c_EDRWS_out[11] + c_aging_EDRWS[10] - c_aging_EDRWS[11]
update(EDRWS[12]) <- EDRWS[12] + c_EDRWS_in[12] - c_EDRWS_out[12] + c_aging_EDRWS[11] - c_aging_EDRWS[12]
update(EDRWS[13]) <- EDRWS[13] + c_EDRWS_in[13] - c_EDRWS_out[13] + c_aging_EDRWS[12] - c_aging_EDRWS[13]
update(EDRWS[14]) <- EDRWS[14] + c_EDRWS_in[14] - c_EDRWS_out[14] + c_aging_EDRWS[13] - c_aging_EDRWS[14]
update(EDRWS[15]) <- EDRWS[15] + c_EDRWS_in[15] - c_EDRWS_out[15] + c_aging_EDRWS[14] - c_aging_EDRWS[15]
update(EDRWS[16]) <- EDRWS[16] + c_EDRWS_in[16] - c_EDRWS_out[16] + c_aging_EDRWS[15]

# ----- 2.73. EDRWV1 Compartment -----------------------------------------------

c_EDRWV1_in[] <- c_EDRWS_to_EDRWV1[i] + c_RWV1_to_EDRWV1[i] + c_migration_EDRWV1[i]
c_EDRWV1_out[] <- c_EDRWV1_to_EDRWS[i] + c_EDRWV1_to_EDRWV2[i] + c_EDRWV1_to_IDRWV1[i] + c_death_EDRWV1[i]

update(EDRWV1[1]) <- EDRWV1[1] + c_EDRWV1_in[1] - c_EDRWV1_out[1] - c_aging_EDRWV1[1]
update(EDRWV1[2]) <- EDRWV1[2] + c_EDRWV1_in[2] - c_EDRWV1_out[2] + c_aging_EDRWV1[1] - c_aging_EDRWV1[2]
update(EDRWV1[3]) <- EDRWV1[3] + c_EDRWV1_in[3] - c_EDRWV1_out[3] + c_aging_EDRWV1[2] - c_aging_EDRWV1[3]
update(EDRWV1[4]) <- EDRWV1[4] + c_EDRWV1_in[4] - c_EDRWV1_out[4] + c_aging_EDRWV1[3] - c_aging_EDRWV1[4]
update(EDRWV1[5]) <- EDRWV1[5] + c_EDRWV1_in[5] - c_EDRWV1_out[5] + c_aging_EDRWV1[4] - c_aging_EDRWV1[5]
update(EDRWV1[6]) <- EDRWV1[6] + c_EDRWV1_in[6] - c_EDRWV1_out[6] + c_aging_EDRWV1[5] - c_aging_EDRWV1[6]
update(EDRWV1[7]) <- EDRWV1[7] + c_EDRWV1_in[7] - c_EDRWV1_out[7] + c_aging_EDRWV1[6] - c_aging_EDRWV1[7]
update(EDRWV1[8]) <- EDRWV1[8] + c_EDRWV1_in[8] - c_EDRWV1_out[8] + c_aging_EDRWV1[7] - c_aging_EDRWV1[8]
update(EDRWV1[9]) <- EDRWV1[9] + c_EDRWV1_in[9] - c_EDRWV1_out[9] + c_aging_EDRWV1[8] - c_aging_EDRWV1[9]
update(EDRWV1[10]) <- EDRWV1[10] + c_EDRWV1_in[10] - c_EDRWV1_out[10] + c_aging_EDRWV1[9] - c_aging_EDRWV1[10]
update(EDRWV1[11]) <- EDRWV1[11] + c_EDRWV1_in[11] - c_EDRWV1_out[11] + c_aging_EDRWV1[10] - c_aging_EDRWV1[11]
update(EDRWV1[12]) <- EDRWV1[12] + c_EDRWV1_in[12] - c_EDRWV1_out[12] + c_aging_EDRWV1[11] - c_aging_EDRWV1[12]
update(EDRWV1[13]) <- EDRWV1[13] + c_EDRWV1_in[13] - c_EDRWV1_out[13] + c_aging_EDRWV1[12] - c_aging_EDRWV1[13]
update(EDRWV1[14]) <- EDRWV1[14] + c_EDRWV1_in[14] - c_EDRWV1_out[14] + c_aging_EDRWV1[13] - c_aging_EDRWV1[14]
update(EDRWV1[15]) <- EDRWV1[15] + c_EDRWV1_in[15] - c_EDRWV1_out[15] + c_aging_EDRWV1[14] - c_aging_EDRWV1[15]
update(EDRWV1[16]) <- EDRWV1[16] + c_EDRWV1_in[16] - c_EDRWV1_out[16] + c_aging_EDRWV1[15]

# ----- 2.74. EDRWV2 Compartment -----------------------------------------------

c_EDRWV2_in[] <- c_EDRWV1_to_EDRWV2[i] + c_RWV2_to_EDRWV2[i] + c_migration_EDRWV2[i]
c_EDRWV2_out[] <- c_EDRWV2_to_EDRWS[i] + c_EDRWV2_to_IDRWV2[i] + c_death_EDRWV2[i]

update(EDRWV2[1]) <- EDRWV2[1] + c_EDRWV2_in[1] - c_EDRWV2_out[1] - c_aging_EDRWV2[1]
update(EDRWV2[2]) <- EDRWV2[2] + c_EDRWV2_in[2] - c_EDRWV2_out[2] + c_aging_EDRWV2[1] - c_aging_EDRWV2[2]
update(EDRWV2[3]) <- EDRWV2[3] + c_EDRWV2_in[3] - c_EDRWV2_out[3] + c_aging_EDRWV2[2] - c_aging_EDRWV2[3]
update(EDRWV2[4]) <- EDRWV2[4] + c_EDRWV2_in[4] - c_EDRWV2_out[4] + c_aging_EDRWV2[3] - c_aging_EDRWV2[4]
update(EDRWV2[5]) <- EDRWV2[5] + c_EDRWV2_in[5] - c_EDRWV2_out[5] + c_aging_EDRWV2[4] - c_aging_EDRWV2[5]
update(EDRWV2[6]) <- EDRWV2[6] + c_EDRWV2_in[6] - c_EDRWV2_out[6] + c_aging_EDRWV2[5] - c_aging_EDRWV2[6]
update(EDRWV2[7]) <- EDRWV2[7] + c_EDRWV2_in[7] - c_EDRWV2_out[7] + c_aging_EDRWV2[6] - c_aging_EDRWV2[7]
update(EDRWV2[8]) <- EDRWV2[8] + c_EDRWV2_in[8] - c_EDRWV2_out[8] + c_aging_EDRWV2[7] - c_aging_EDRWV2[8]
update(EDRWV2[9]) <- EDRWV2[9] + c_EDRWV2_in[9] - c_EDRWV2_out[9] + c_aging_EDRWV2[8] - c_aging_EDRWV2[9]
update(EDRWV2[10]) <- EDRWV2[10] + c_EDRWV2_in[10] - c_EDRWV2_out[10] + c_aging_EDRWV2[9] - c_aging_EDRWV2[10]
update(EDRWV2[11]) <- EDRWV2[11] + c_EDRWV2_in[11] - c_EDRWV2_out[11] + c_aging_EDRWV2[10] - c_aging_EDRWV2[11]
update(EDRWV2[12]) <- EDRWV2[12] + c_EDRWV2_in[12] - c_EDRWV2_out[12] + c_aging_EDRWV2[11] - c_aging_EDRWV2[12]
update(EDRWV2[13]) <- EDRWV2[13] + c_EDRWV2_in[13] - c_EDRWV2_out[13] + c_aging_EDRWV2[12] - c_aging_EDRWV2[13]
update(EDRWV2[14]) <- EDRWV2[14] + c_EDRWV2_in[14] - c_EDRWV2_out[14] + c_aging_EDRWV2[13] - c_aging_EDRWV2[14]
update(EDRWV2[15]) <- EDRWV2[15] + c_EDRWV2_in[15] - c_EDRWV2_out[15] + c_aging_EDRWV2[14] - c_aging_EDRWV2[15]
update(EDRWV2[16]) <- EDRWV2[16] + c_EDRWV2_in[16] - c_EDRWV2_out[16] + c_aging_EDRWV2[15]

# ----- 2.76. EDRAS Compartment -----------------------------------------------

c_EDRAS_in[] <- c_EDRAV1_to_EDRAS[i] + c_EDRAV2_to_EDRAS[i] + c_RAS_to_EDRAS[i] + c_migration_EDRAS[i]
c_EDRAS_out[] <- c_EDRAS_to_EDRAV1[i] + c_EDRAS_to_IDRAS[i] + c_death_EDRAS[i]

update(EDRAS[1]) <- EDRAS[1] + c_EDRAS_in[1] - c_EDRAS_out[1] - c_aging_EDRAS[1]
update(EDRAS[2]) <- EDRAS[2] + c_EDRAS_in[2] - c_EDRAS_out[2] + c_aging_EDRAS[1] - c_aging_EDRAS[2]
update(EDRAS[3]) <- EDRAS[3] + c_EDRAS_in[3] - c_EDRAS_out[3] + c_aging_EDRAS[2] - c_aging_EDRAS[3]
update(EDRAS[4]) <- EDRAS[4] + c_EDRAS_in[4] - c_EDRAS_out[4] + c_aging_EDRAS[3] - c_aging_EDRAS[4]
update(EDRAS[5]) <- EDRAS[5] + c_EDRAS_in[5] - c_EDRAS_out[5] + c_aging_EDRAS[4] - c_aging_EDRAS[5]
update(EDRAS[6]) <- EDRAS[6] + c_EDRAS_in[6] - c_EDRAS_out[6] + c_aging_EDRAS[5] - c_aging_EDRAS[6]
update(EDRAS[7]) <- EDRAS[7] + c_EDRAS_in[7] - c_EDRAS_out[7] + c_aging_EDRAS[6] - c_aging_EDRAS[7]
update(EDRAS[8]) <- EDRAS[8] + c_EDRAS_in[8] - c_EDRAS_out[8] + c_aging_EDRAS[7] - c_aging_EDRAS[8]
update(EDRAS[9]) <- EDRAS[9] + c_EDRAS_in[9] - c_EDRAS_out[9] + c_aging_EDRAS[8] - c_aging_EDRAS[9]
update(EDRAS[10]) <- EDRAS[10] + c_EDRAS_in[10] - c_EDRAS_out[10] + c_aging_EDRAS[9] - c_aging_EDRAS[10]
update(EDRAS[11]) <- EDRAS[11] + c_EDRAS_in[11] - c_EDRAS_out[11] + c_aging_EDRAS[10] - c_aging_EDRAS[11]
update(EDRAS[12]) <- EDRAS[12] + c_EDRAS_in[12] - c_EDRAS_out[12] + c_aging_EDRAS[11] - c_aging_EDRAS[12]
update(EDRAS[13]) <- EDRAS[13] + c_EDRAS_in[13] - c_EDRAS_out[13] + c_aging_EDRAS[12] - c_aging_EDRAS[13]
update(EDRAS[14]) <- EDRAS[14] + c_EDRAS_in[14] - c_EDRAS_out[14] + c_aging_EDRAS[13] - c_aging_EDRAS[14]
update(EDRAS[15]) <- EDRAS[15] + c_EDRAS_in[15] - c_EDRAS_out[15] + c_aging_EDRAS[14] - c_aging_EDRAS[15]
update(EDRAS[16]) <- EDRAS[16] + c_EDRAS_in[16] - c_EDRAS_out[16] + c_aging_EDRAS[15]

# ----- 2.77. EDRAV1 Compartment -----------------------------------------------

c_EDRAV1_in[] <- c_EDRAS_to_EDRAV1[i] + c_RAV1_to_EDRAV1[i] + c_migration_EDRAV1[i]
c_EDRAV1_out[] <- c_EDRAV1_to_EDRAS[i] + c_EDRAV1_to_EDRAV2[i] + c_EDRAV1_to_IDRAV1[i] + c_death_EDRAV1[i]

update(EDRAV1[1]) <- EDRAV1[1] + c_EDRAV1_in[1] - c_EDRAV1_out[1] - c_aging_EDRAV1[1]
update(EDRAV1[2]) <- EDRAV1[2] + c_EDRAV1_in[2] - c_EDRAV1_out[2] + c_aging_EDRAV1[1] - c_aging_EDRAV1[2]
update(EDRAV1[3]) <- EDRAV1[3] + c_EDRAV1_in[3] - c_EDRAV1_out[3] + c_aging_EDRAV1[2] - c_aging_EDRAV1[3]
update(EDRAV1[4]) <- EDRAV1[4] + c_EDRAV1_in[4] - c_EDRAV1_out[4] + c_aging_EDRAV1[3] - c_aging_EDRAV1[4]
update(EDRAV1[5]) <- EDRAV1[5] + c_EDRAV1_in[5] - c_EDRAV1_out[5] + c_aging_EDRAV1[4] - c_aging_EDRAV1[5]
update(EDRAV1[6]) <- EDRAV1[6] + c_EDRAV1_in[6] - c_EDRAV1_out[6] + c_aging_EDRAV1[5] - c_aging_EDRAV1[6]
update(EDRAV1[7]) <- EDRAV1[7] + c_EDRAV1_in[7] - c_EDRAV1_out[7] + c_aging_EDRAV1[6] - c_aging_EDRAV1[7]
update(EDRAV1[8]) <- EDRAV1[8] + c_EDRAV1_in[8] - c_EDRAV1_out[8] + c_aging_EDRAV1[7] - c_aging_EDRAV1[8]
update(EDRAV1[9]) <- EDRAV1[9] + c_EDRAV1_in[9] - c_EDRAV1_out[9] + c_aging_EDRAV1[8] - c_aging_EDRAV1[9]
update(EDRAV1[10]) <- EDRAV1[10] + c_EDRAV1_in[10] - c_EDRAV1_out[10] + c_aging_EDRAV1[9] - c_aging_EDRAV1[10]
update(EDRAV1[11]) <- EDRAV1[11] + c_EDRAV1_in[11] - c_EDRAV1_out[11] + c_aging_EDRAV1[10] - c_aging_EDRAV1[11]
update(EDRAV1[12]) <- EDRAV1[12] + c_EDRAV1_in[12] - c_EDRAV1_out[12] + c_aging_EDRAV1[11] - c_aging_EDRAV1[12]
update(EDRAV1[13]) <- EDRAV1[13] + c_EDRAV1_in[13] - c_EDRAV1_out[13] + c_aging_EDRAV1[12] - c_aging_EDRAV1[13]
update(EDRAV1[14]) <- EDRAV1[14] + c_EDRAV1_in[14] - c_EDRAV1_out[14] + c_aging_EDRAV1[13] - c_aging_EDRAV1[14]
update(EDRAV1[15]) <- EDRAV1[15] + c_EDRAV1_in[15] - c_EDRAV1_out[15] + c_aging_EDRAV1[14] - c_aging_EDRAV1[15]
update(EDRAV1[16]) <- EDRAV1[16] + c_EDRAV1_in[16] - c_EDRAV1_out[16] + c_aging_EDRAV1[15]

# ----- 2.78. EDRAV2 Compartment -----------------------------------------------

c_EDRAV2_in[] <- c_EDRAV1_to_EDRAV2[i] + c_RAV2_to_EDRAV2[i] + c_migration_EDRAV2[i]
c_EDRAV2_out[] <- c_EDRAV2_to_EDRAS[i] + c_EDRAV2_to_IDRAV2[i] + c_death_EDRAV2[i]

update(EDRAV2[1]) <- EDRAV2[1] + c_EDRAV2_in[1] - c_EDRAV2_out[1] - c_aging_EDRAV2[1]
update(EDRAV2[2]) <- EDRAV2[2] + c_EDRAV2_in[2] - c_EDRAV2_out[2] + c_aging_EDRAV2[1] - c_aging_EDRAV2[2]
update(EDRAV2[3]) <- EDRAV2[3] + c_EDRAV2_in[3] - c_EDRAV2_out[3] + c_aging_EDRAV2[2] - c_aging_EDRAV2[3]
update(EDRAV2[4]) <- EDRAV2[4] + c_EDRAV2_in[4] - c_EDRAV2_out[4] + c_aging_EDRAV2[3] - c_aging_EDRAV2[4]
update(EDRAV2[5]) <- EDRAV2[5] + c_EDRAV2_in[5] - c_EDRAV2_out[5] + c_aging_EDRAV2[4] - c_aging_EDRAV2[5]
update(EDRAV2[6]) <- EDRAV2[6] + c_EDRAV2_in[6] - c_EDRAV2_out[6] + c_aging_EDRAV2[5] - c_aging_EDRAV2[6]
update(EDRAV2[7]) <- EDRAV2[7] + c_EDRAV2_in[7] - c_EDRAV2_out[7] + c_aging_EDRAV2[6] - c_aging_EDRAV2[7]
update(EDRAV2[8]) <- EDRAV2[8] + c_EDRAV2_in[8] - c_EDRAV2_out[8] + c_aging_EDRAV2[7] - c_aging_EDRAV2[8]
update(EDRAV2[9]) <- EDRAV2[9] + c_EDRAV2_in[9] - c_EDRAV2_out[9] + c_aging_EDRAV2[8] - c_aging_EDRAV2[9]
update(EDRAV2[10]) <- EDRAV2[10] + c_EDRAV2_in[10] - c_EDRAV2_out[10] + c_aging_EDRAV2[9] - c_aging_EDRAV2[10]
update(EDRAV2[11]) <- EDRAV2[11] + c_EDRAV2_in[11] - c_EDRAV2_out[11] + c_aging_EDRAV2[10] - c_aging_EDRAV2[11]
update(EDRAV2[12]) <- EDRAV2[12] + c_EDRAV2_in[12] - c_EDRAV2_out[12] + c_aging_EDRAV2[11] - c_aging_EDRAV2[12]
update(EDRAV2[13]) <- EDRAV2[13] + c_EDRAV2_in[13] - c_EDRAV2_out[13] + c_aging_EDRAV2[12] - c_aging_EDRAV2[13]
update(EDRAV2[14]) <- EDRAV2[14] + c_EDRAV2_in[14] - c_EDRAV2_out[14] + c_aging_EDRAV2[13] - c_aging_EDRAV2[14]
update(EDRAV2[15]) <- EDRAV2[15] + c_EDRAV2_in[15] - c_EDRAV2_out[15] + c_aging_EDRAV2[14] - c_aging_EDRAV2[15]
update(EDRAV2[16]) <- EDRAV2[16] + c_EDRAV2_in[16] - c_EDRAV2_out[16] + c_aging_EDRAV2[15]

# ----- 2.80. EORWS Compartment -----------------------------------------------

c_EORWS_in[] <- c_EORWV1_to_EORWS[i] + c_EORWV2_to_EORWS[i] + c_RWS_to_EORWS[i] + c_migration_EORWS[i]
c_EORWS_out[] <- c_EORWS_to_EORWV1[i] + c_EORWS_to_IORWS[i] + c_death_EORWS[i]

update(EORWS[1]) <- EORWS[1] + c_EORWS_in[1] - c_EORWS_out[1] - c_aging_EORWS[1]
update(EORWS[2]) <- EORWS[2] + c_EORWS_in[2] - c_EORWS_out[2] + c_aging_EORWS[1] - c_aging_EORWS[2]
update(EORWS[3]) <- EORWS[3] + c_EORWS_in[3] - c_EORWS_out[3] + c_aging_EORWS[2] - c_aging_EORWS[3]
update(EORWS[4]) <- EORWS[4] + c_EORWS_in[4] - c_EORWS_out[4] + c_aging_EORWS[3] - c_aging_EORWS[4]
update(EORWS[5]) <- EORWS[5] + c_EORWS_in[5] - c_EORWS_out[5] + c_aging_EORWS[4] - c_aging_EORWS[5]
update(EORWS[6]) <- EORWS[6] + c_EORWS_in[6] - c_EORWS_out[6] + c_aging_EORWS[5] - c_aging_EORWS[6]
update(EORWS[7]) <- EORWS[7] + c_EORWS_in[7] - c_EORWS_out[7] + c_aging_EORWS[6] - c_aging_EORWS[7]
update(EORWS[8]) <- EORWS[8] + c_EORWS_in[8] - c_EORWS_out[8] + c_aging_EORWS[7] - c_aging_EORWS[8]
update(EORWS[9]) <- EORWS[9] + c_EORWS_in[9] - c_EORWS_out[9] + c_aging_EORWS[8] - c_aging_EORWS[9]
update(EORWS[10]) <- EORWS[10] + c_EORWS_in[10] - c_EORWS_out[10] + c_aging_EORWS[9] - c_aging_EORWS[10]
update(EORWS[11]) <- EORWS[11] + c_EORWS_in[11] - c_EORWS_out[11] + c_aging_EORWS[10] - c_aging_EORWS[11]
update(EORWS[12]) <- EORWS[12] + c_EORWS_in[12] - c_EORWS_out[12] + c_aging_EORWS[11] - c_aging_EORWS[12]
update(EORWS[13]) <- EORWS[13] + c_EORWS_in[13] - c_EORWS_out[13] + c_aging_EORWS[12] - c_aging_EORWS[13]
update(EORWS[14]) <- EORWS[14] + c_EORWS_in[14] - c_EORWS_out[14] + c_aging_EORWS[13] - c_aging_EORWS[14]
update(EORWS[15]) <- EORWS[15] + c_EORWS_in[15] - c_EORWS_out[15] + c_aging_EORWS[14] - c_aging_EORWS[15]
update(EORWS[16]) <- EORWS[16] + c_EORWS_in[16] - c_EORWS_out[16] + c_aging_EORWS[15]

# ----- 2.81. EORWV1 Compartment -----------------------------------------------

c_EORWV1_in[] <- c_EORWS_to_EORWV1[i] + c_RWV1_to_EORWV1[i] + c_migration_EORWV1[i]
c_EORWV1_out[] <- c_EORWV1_to_EORWS[i] + c_EORWV1_to_EORWV2[i] + c_EORWV1_to_IORWV1[i] + c_death_EORWV1[i]

update(EORWV1[1]) <- EORWV1[1] + c_EORWV1_in[1] - c_EORWV1_out[1] - c_aging_EORWV1[1]
update(EORWV1[2]) <- EORWV1[2] + c_EORWV1_in[2] - c_EORWV1_out[2] + c_aging_EORWV1[1] - c_aging_EORWV1[2]
update(EORWV1[3]) <- EORWV1[3] + c_EORWV1_in[3] - c_EORWV1_out[3] + c_aging_EORWV1[2] - c_aging_EORWV1[3]
update(EORWV1[4]) <- EORWV1[4] + c_EORWV1_in[4] - c_EORWV1_out[4] + c_aging_EORWV1[3] - c_aging_EORWV1[4]
update(EORWV1[5]) <- EORWV1[5] + c_EORWV1_in[5] - c_EORWV1_out[5] + c_aging_EORWV1[4] - c_aging_EORWV1[5]
update(EORWV1[6]) <- EORWV1[6] + c_EORWV1_in[6] - c_EORWV1_out[6] + c_aging_EORWV1[5] - c_aging_EORWV1[6]
update(EORWV1[7]) <- EORWV1[7] + c_EORWV1_in[7] - c_EORWV1_out[7] + c_aging_EORWV1[6] - c_aging_EORWV1[7]
update(EORWV1[8]) <- EORWV1[8] + c_EORWV1_in[8] - c_EORWV1_out[8] + c_aging_EORWV1[7] - c_aging_EORWV1[8]
update(EORWV1[9]) <- EORWV1[9] + c_EORWV1_in[9] - c_EORWV1_out[9] + c_aging_EORWV1[8] - c_aging_EORWV1[9]
update(EORWV1[10]) <- EORWV1[10] + c_EORWV1_in[10] - c_EORWV1_out[10] + c_aging_EORWV1[9] - c_aging_EORWV1[10]
update(EORWV1[11]) <- EORWV1[11] + c_EORWV1_in[11] - c_EORWV1_out[11] + c_aging_EORWV1[10] - c_aging_EORWV1[11]
update(EORWV1[12]) <- EORWV1[12] + c_EORWV1_in[12] - c_EORWV1_out[12] + c_aging_EORWV1[11] - c_aging_EORWV1[12]
update(EORWV1[13]) <- EORWV1[13] + c_EORWV1_in[13] - c_EORWV1_out[13] + c_aging_EORWV1[12] - c_aging_EORWV1[13]
update(EORWV1[14]) <- EORWV1[14] + c_EORWV1_in[14] - c_EORWV1_out[14] + c_aging_EORWV1[13] - c_aging_EORWV1[14]
update(EORWV1[15]) <- EORWV1[15] + c_EORWV1_in[15] - c_EORWV1_out[15] + c_aging_EORWV1[14] - c_aging_EORWV1[15]
update(EORWV1[16]) <- EORWV1[16] + c_EORWV1_in[16] - c_EORWV1_out[16] + c_aging_EORWV1[15]

# ----- 2.82. EORWV2 Compartment -----------------------------------------------

c_EORWV2_in[] <- c_EORWV1_to_EORWV2[i] + c_RWV2_to_EORWV2[i] + c_migration_EORWV2[i]
c_EORWV2_out[] <- c_EORWV2_to_EORWS[i] + c_EORWV2_to_IORWV2[i] + c_death_EORWV2[i]

update(EORWV2[1]) <- EORWV2[1] + c_EORWV2_in[1] - c_EORWV2_out[1] - c_aging_EORWV2[1]
update(EORWV2[2]) <- EORWV2[2] + c_EORWV2_in[2] - c_EORWV2_out[2] + c_aging_EORWV2[1] - c_aging_EORWV2[2]
update(EORWV2[3]) <- EORWV2[3] + c_EORWV2_in[3] - c_EORWV2_out[3] + c_aging_EORWV2[2] - c_aging_EORWV2[3]
update(EORWV2[4]) <- EORWV2[4] + c_EORWV2_in[4] - c_EORWV2_out[4] + c_aging_EORWV2[3] - c_aging_EORWV2[4]
update(EORWV2[5]) <- EORWV2[5] + c_EORWV2_in[5] - c_EORWV2_out[5] + c_aging_EORWV2[4] - c_aging_EORWV2[5]
update(EORWV2[6]) <- EORWV2[6] + c_EORWV2_in[6] - c_EORWV2_out[6] + c_aging_EORWV2[5] - c_aging_EORWV2[6]
update(EORWV2[7]) <- EORWV2[7] + c_EORWV2_in[7] - c_EORWV2_out[7] + c_aging_EORWV2[6] - c_aging_EORWV2[7]
update(EORWV2[8]) <- EORWV2[8] + c_EORWV2_in[8] - c_EORWV2_out[8] + c_aging_EORWV2[7] - c_aging_EORWV2[8]
update(EORWV2[9]) <- EORWV2[9] + c_EORWV2_in[9] - c_EORWV2_out[9] + c_aging_EORWV2[8] - c_aging_EORWV2[9]
update(EORWV2[10]) <- EORWV2[10] + c_EORWV2_in[10] - c_EORWV2_out[10] + c_aging_EORWV2[9] - c_aging_EORWV2[10]
update(EORWV2[11]) <- EORWV2[11] + c_EORWV2_in[11] - c_EORWV2_out[11] + c_aging_EORWV2[10] - c_aging_EORWV2[11]
update(EORWV2[12]) <- EORWV2[12] + c_EORWV2_in[12] - c_EORWV2_out[12] + c_aging_EORWV2[11] - c_aging_EORWV2[12]
update(EORWV2[13]) <- EORWV2[13] + c_EORWV2_in[13] - c_EORWV2_out[13] + c_aging_EORWV2[12] - c_aging_EORWV2[13]
update(EORWV2[14]) <- EORWV2[14] + c_EORWV2_in[14] - c_EORWV2_out[14] + c_aging_EORWV2[13] - c_aging_EORWV2[14]
update(EORWV2[15]) <- EORWV2[15] + c_EORWV2_in[15] - c_EORWV2_out[15] + c_aging_EORWV2[14] - c_aging_EORWV2[15]
update(EORWV2[16]) <- EORWV2[16] + c_EORWV2_in[16] - c_EORWV2_out[16] + c_aging_EORWV2[15]

# ----- 2.84. EORAS Compartment -----------------------------------------------

c_EORAS_in[] <- c_EORAV1_to_EORAS[i] + c_EORAV2_to_EORAS[i] + c_RAS_to_EORAS[i] + c_migration_EORAS[i]
c_EORAS_out[] <- c_EORAS_to_EORAV1[i] + c_EORAS_to_IORAS[i] + c_death_EORAS[i]

update(EORAS[1]) <- EORAS[1] + c_EORAS_in[1] - c_EORAS_out[1] - c_aging_EORAS[1]
update(EORAS[2]) <- EORAS[2] + c_EORAS_in[2] - c_EORAS_out[2] + c_aging_EORAS[1] - c_aging_EORAS[2]
update(EORAS[3]) <- EORAS[3] + c_EORAS_in[3] - c_EORAS_out[3] + c_aging_EORAS[2] - c_aging_EORAS[3]
update(EORAS[4]) <- EORAS[4] + c_EORAS_in[4] - c_EORAS_out[4] + c_aging_EORAS[3] - c_aging_EORAS[4]
update(EORAS[5]) <- EORAS[5] + c_EORAS_in[5] - c_EORAS_out[5] + c_aging_EORAS[4] - c_aging_EORAS[5]
update(EORAS[6]) <- EORAS[6] + c_EORAS_in[6] - c_EORAS_out[6] + c_aging_EORAS[5] - c_aging_EORAS[6]
update(EORAS[7]) <- EORAS[7] + c_EORAS_in[7] - c_EORAS_out[7] + c_aging_EORAS[6] - c_aging_EORAS[7]
update(EORAS[8]) <- EORAS[8] + c_EORAS_in[8] - c_EORAS_out[8] + c_aging_EORAS[7] - c_aging_EORAS[8]
update(EORAS[9]) <- EORAS[9] + c_EORAS_in[9] - c_EORAS_out[9] + c_aging_EORAS[8] - c_aging_EORAS[9]
update(EORAS[10]) <- EORAS[10] + c_EORAS_in[10] - c_EORAS_out[10] + c_aging_EORAS[9] - c_aging_EORAS[10]
update(EORAS[11]) <- EORAS[11] + c_EORAS_in[11] - c_EORAS_out[11] + c_aging_EORAS[10] - c_aging_EORAS[11]
update(EORAS[12]) <- EORAS[12] + c_EORAS_in[12] - c_EORAS_out[12] + c_aging_EORAS[11] - c_aging_EORAS[12]
update(EORAS[13]) <- EORAS[13] + c_EORAS_in[13] - c_EORAS_out[13] + c_aging_EORAS[12] - c_aging_EORAS[13]
update(EORAS[14]) <- EORAS[14] + c_EORAS_in[14] - c_EORAS_out[14] + c_aging_EORAS[13] - c_aging_EORAS[14]
update(EORAS[15]) <- EORAS[15] + c_EORAS_in[15] - c_EORAS_out[15] + c_aging_EORAS[14] - c_aging_EORAS[15]
update(EORAS[16]) <- EORAS[16] + c_EORAS_in[16] - c_EORAS_out[16] + c_aging_EORAS[15]

# ----- 2.85. EORAV1 Compartment -----------------------------------------------

c_EORAV1_in[] <- c_EORAS_to_EORAV1[i] + c_RAV1_to_EORAV1[i] + c_migration_EORAV1[i]
c_EORAV1_out[] <- c_EORAV1_to_EORAS[i] + c_EORAV1_to_EORAV2[i] + c_EORAV1_to_IORAV1[i] + c_death_EORAV1[i]

update(EORAV1[1]) <- EORAV1[1] + c_EORAV1_in[1] - c_EORAV1_out[1] - c_aging_EORAV1[1]
update(EORAV1[2]) <- EORAV1[2] + c_EORAV1_in[2] - c_EORAV1_out[2] + c_aging_EORAV1[1] - c_aging_EORAV1[2]
update(EORAV1[3]) <- EORAV1[3] + c_EORAV1_in[3] - c_EORAV1_out[3] + c_aging_EORAV1[2] - c_aging_EORAV1[3]
update(EORAV1[4]) <- EORAV1[4] + c_EORAV1_in[4] - c_EORAV1_out[4] + c_aging_EORAV1[3] - c_aging_EORAV1[4]
update(EORAV1[5]) <- EORAV1[5] + c_EORAV1_in[5] - c_EORAV1_out[5] + c_aging_EORAV1[4] - c_aging_EORAV1[5]
update(EORAV1[6]) <- EORAV1[6] + c_EORAV1_in[6] - c_EORAV1_out[6] + c_aging_EORAV1[5] - c_aging_EORAV1[6]
update(EORAV1[7]) <- EORAV1[7] + c_EORAV1_in[7] - c_EORAV1_out[7] + c_aging_EORAV1[6] - c_aging_EORAV1[7]
update(EORAV1[8]) <- EORAV1[8] + c_EORAV1_in[8] - c_EORAV1_out[8] + c_aging_EORAV1[7] - c_aging_EORAV1[8]
update(EORAV1[9]) <- EORAV1[9] + c_EORAV1_in[9] - c_EORAV1_out[9] + c_aging_EORAV1[8] - c_aging_EORAV1[9]
update(EORAV1[10]) <- EORAV1[10] + c_EORAV1_in[10] - c_EORAV1_out[10] + c_aging_EORAV1[9] - c_aging_EORAV1[10]
update(EORAV1[11]) <- EORAV1[11] + c_EORAV1_in[11] - c_EORAV1_out[11] + c_aging_EORAV1[10] - c_aging_EORAV1[11]
update(EORAV1[12]) <- EORAV1[12] + c_EORAV1_in[12] - c_EORAV1_out[12] + c_aging_EORAV1[11] - c_aging_EORAV1[12]
update(EORAV1[13]) <- EORAV1[13] + c_EORAV1_in[13] - c_EORAV1_out[13] + c_aging_EORAV1[12] - c_aging_EORAV1[13]
update(EORAV1[14]) <- EORAV1[14] + c_EORAV1_in[14] - c_EORAV1_out[14] + c_aging_EORAV1[13] - c_aging_EORAV1[14]
update(EORAV1[15]) <- EORAV1[15] + c_EORAV1_in[15] - c_EORAV1_out[15] + c_aging_EORAV1[14] - c_aging_EORAV1[15]
update(EORAV1[16]) <- EORAV1[16] + c_EORAV1_in[16] - c_EORAV1_out[16] + c_aging_EORAV1[15]

# ----- 2.86. EORAV2 Compartment -----------------------------------------------

c_EORAV2_in[] <- c_EORAV1_to_EORAV2[i] + c_RAV2_to_EORAV2[i] + c_migration_EORAV2[i]
c_EORAV2_out[] <- c_EORAV2_to_EORAS[i] + c_EORAV2_to_IORAV2[i] + c_death_EORAV2[i]

update(EORAV2[1]) <- EORAV2[1] + c_EORAV2_in[1] - c_EORAV2_out[1] - c_aging_EORAV2[1]
update(EORAV2[2]) <- EORAV2[2] + c_EORAV2_in[2] - c_EORAV2_out[2] + c_aging_EORAV2[1] - c_aging_EORAV2[2]
update(EORAV2[3]) <- EORAV2[3] + c_EORAV2_in[3] - c_EORAV2_out[3] + c_aging_EORAV2[2] - c_aging_EORAV2[3]
update(EORAV2[4]) <- EORAV2[4] + c_EORAV2_in[4] - c_EORAV2_out[4] + c_aging_EORAV2[3] - c_aging_EORAV2[4]
update(EORAV2[5]) <- EORAV2[5] + c_EORAV2_in[5] - c_EORAV2_out[5] + c_aging_EORAV2[4] - c_aging_EORAV2[5]
update(EORAV2[6]) <- EORAV2[6] + c_EORAV2_in[6] - c_EORAV2_out[6] + c_aging_EORAV2[5] - c_aging_EORAV2[6]
update(EORAV2[7]) <- EORAV2[7] + c_EORAV2_in[7] - c_EORAV2_out[7] + c_aging_EORAV2[6] - c_aging_EORAV2[7]
update(EORAV2[8]) <- EORAV2[8] + c_EORAV2_in[8] - c_EORAV2_out[8] + c_aging_EORAV2[7] - c_aging_EORAV2[8]
update(EORAV2[9]) <- EORAV2[9] + c_EORAV2_in[9] - c_EORAV2_out[9] + c_aging_EORAV2[8] - c_aging_EORAV2[9]
update(EORAV2[10]) <- EORAV2[10] + c_EORAV2_in[10] - c_EORAV2_out[10] + c_aging_EORAV2[9] - c_aging_EORAV2[10]
update(EORAV2[11]) <- EORAV2[11] + c_EORAV2_in[11] - c_EORAV2_out[11] + c_aging_EORAV2[10] - c_aging_EORAV2[11]
update(EORAV2[12]) <- EORAV2[12] + c_EORAV2_in[12] - c_EORAV2_out[12] + c_aging_EORAV2[11] - c_aging_EORAV2[12]
update(EORAV2[13]) <- EORAV2[13] + c_EORAV2_in[13] - c_EORAV2_out[13] + c_aging_EORAV2[12] - c_aging_EORAV2[13]
update(EORAV2[14]) <- EORAV2[14] + c_EORAV2_in[14] - c_EORAV2_out[14] + c_aging_EORAV2[13] - c_aging_EORAV2[14]
update(EORAV2[15]) <- EORAV2[15] + c_EORAV2_in[15] - c_EORAV2_out[15] + c_aging_EORAV2[14] - c_aging_EORAV2[15]
update(EORAV2[16]) <- EORAV2[16] + c_EORAV2_in[16] - c_EORAV2_out[16] + c_aging_EORAV2[15]

# ----- 2.88. EORDS Compartment -----------------------------------------------

c_EORDS_in[] <- c_EORDV1_to_EORDS[i] + c_EORDV2_to_EORDS[i] + c_RDS_to_EORDS[i] + c_migration_EORDS[i]
c_EORDS_out[] <- c_EORDS_to_EORDV1[i] + c_EORDS_to_IORDS[i] + c_death_EORDS[i]

update(EORDS[1]) <- EORDS[1] + c_EORDS_in[1] - c_EORDS_out[1] - c_aging_EORDS[1]
update(EORDS[2]) <- EORDS[2] + c_EORDS_in[2] - c_EORDS_out[2] + c_aging_EORDS[1] - c_aging_EORDS[2]
update(EORDS[3]) <- EORDS[3] + c_EORDS_in[3] - c_EORDS_out[3] + c_aging_EORDS[2] - c_aging_EORDS[3]
update(EORDS[4]) <- EORDS[4] + c_EORDS_in[4] - c_EORDS_out[4] + c_aging_EORDS[3] - c_aging_EORDS[4]
update(EORDS[5]) <- EORDS[5] + c_EORDS_in[5] - c_EORDS_out[5] + c_aging_EORDS[4] - c_aging_EORDS[5]
update(EORDS[6]) <- EORDS[6] + c_EORDS_in[6] - c_EORDS_out[6] + c_aging_EORDS[5] - c_aging_EORDS[6]
update(EORDS[7]) <- EORDS[7] + c_EORDS_in[7] - c_EORDS_out[7] + c_aging_EORDS[6] - c_aging_EORDS[7]
update(EORDS[8]) <- EORDS[8] + c_EORDS_in[8] - c_EORDS_out[8] + c_aging_EORDS[7] - c_aging_EORDS[8]
update(EORDS[9]) <- EORDS[9] + c_EORDS_in[9] - c_EORDS_out[9] + c_aging_EORDS[8] - c_aging_EORDS[9]
update(EORDS[10]) <- EORDS[10] + c_EORDS_in[10] - c_EORDS_out[10] + c_aging_EORDS[9] - c_aging_EORDS[10]
update(EORDS[11]) <- EORDS[11] + c_EORDS_in[11] - c_EORDS_out[11] + c_aging_EORDS[10] - c_aging_EORDS[11]
update(EORDS[12]) <- EORDS[12] + c_EORDS_in[12] - c_EORDS_out[12] + c_aging_EORDS[11] - c_aging_EORDS[12]
update(EORDS[13]) <- EORDS[13] + c_EORDS_in[13] - c_EORDS_out[13] + c_aging_EORDS[12] - c_aging_EORDS[13]
update(EORDS[14]) <- EORDS[14] + c_EORDS_in[14] - c_EORDS_out[14] + c_aging_EORDS[13] - c_aging_EORDS[14]
update(EORDS[15]) <- EORDS[15] + c_EORDS_in[15] - c_EORDS_out[15] + c_aging_EORDS[14] - c_aging_EORDS[15]
update(EORDS[16]) <- EORDS[16] + c_EORDS_in[16] - c_EORDS_out[16] + c_aging_EORDS[15]

# ----- 2.89. EORDV1 Compartment -----------------------------------------------

c_EORDV1_in[] <- c_EORDS_to_EORDV1[i] + c_RDV1_to_EORDV1[i] + c_migration_EORDV1[i]
c_EORDV1_out[] <- c_EORDV1_to_EORDS[i] + c_EORDV1_to_EORDV2[i] + c_EORDV1_to_IORDV1[i] + c_death_EORDV1[i]

update(EORDV1[1]) <- EORDV1[1] + c_EORDV1_in[1] - c_EORDV1_out[1] - c_aging_EORDV1[1]
update(EORDV1[2]) <- EORDV1[2] + c_EORDV1_in[2] - c_EORDV1_out[2] + c_aging_EORDV1[1] - c_aging_EORDV1[2]
update(EORDV1[3]) <- EORDV1[3] + c_EORDV1_in[3] - c_EORDV1_out[3] + c_aging_EORDV1[2] - c_aging_EORDV1[3]
update(EORDV1[4]) <- EORDV1[4] + c_EORDV1_in[4] - c_EORDV1_out[4] + c_aging_EORDV1[3] - c_aging_EORDV1[4]
update(EORDV1[5]) <- EORDV1[5] + c_EORDV1_in[5] - c_EORDV1_out[5] + c_aging_EORDV1[4] - c_aging_EORDV1[5]
update(EORDV1[6]) <- EORDV1[6] + c_EORDV1_in[6] - c_EORDV1_out[6] + c_aging_EORDV1[5] - c_aging_EORDV1[6]
update(EORDV1[7]) <- EORDV1[7] + c_EORDV1_in[7] - c_EORDV1_out[7] + c_aging_EORDV1[6] - c_aging_EORDV1[7]
update(EORDV1[8]) <- EORDV1[8] + c_EORDV1_in[8] - c_EORDV1_out[8] + c_aging_EORDV1[7] - c_aging_EORDV1[8]
update(EORDV1[9]) <- EORDV1[9] + c_EORDV1_in[9] - c_EORDV1_out[9] + c_aging_EORDV1[8] - c_aging_EORDV1[9]
update(EORDV1[10]) <- EORDV1[10] + c_EORDV1_in[10] - c_EORDV1_out[10] + c_aging_EORDV1[9] - c_aging_EORDV1[10]
update(EORDV1[11]) <- EORDV1[11] + c_EORDV1_in[11] - c_EORDV1_out[11] + c_aging_EORDV1[10] - c_aging_EORDV1[11]
update(EORDV1[12]) <- EORDV1[12] + c_EORDV1_in[12] - c_EORDV1_out[12] + c_aging_EORDV1[11] - c_aging_EORDV1[12]
update(EORDV1[13]) <- EORDV1[13] + c_EORDV1_in[13] - c_EORDV1_out[13] + c_aging_EORDV1[12] - c_aging_EORDV1[13]
update(EORDV1[14]) <- EORDV1[14] + c_EORDV1_in[14] - c_EORDV1_out[14] + c_aging_EORDV1[13] - c_aging_EORDV1[14]
update(EORDV1[15]) <- EORDV1[15] + c_EORDV1_in[15] - c_EORDV1_out[15] + c_aging_EORDV1[14] - c_aging_EORDV1[15]
update(EORDV1[16]) <- EORDV1[16] + c_EORDV1_in[16] - c_EORDV1_out[16] + c_aging_EORDV1[15]

# ----- 2.90. EORDV2 Compartment -----------------------------------------------

c_EORDV2_in[] <- c_EORDV1_to_EORDV2[i] + c_RDV2_to_EORDV2[i] + c_migration_EORDV2[i]
c_EORDV2_out[] <- c_EORDV2_to_EORDS[i] + c_EORDV2_to_IORDV2[i] + c_death_EORDV2[i]

update(EORDV2[1]) <- EORDV2[1] + c_EORDV2_in[1] - c_EORDV2_out[1] - c_aging_EORDV2[1]
update(EORDV2[2]) <- EORDV2[2] + c_EORDV2_in[2] - c_EORDV2_out[2] + c_aging_EORDV2[1] - c_aging_EORDV2[2]
update(EORDV2[3]) <- EORDV2[3] + c_EORDV2_in[3] - c_EORDV2_out[3] + c_aging_EORDV2[2] - c_aging_EORDV2[3]
update(EORDV2[4]) <- EORDV2[4] + c_EORDV2_in[4] - c_EORDV2_out[4] + c_aging_EORDV2[3] - c_aging_EORDV2[4]
update(EORDV2[5]) <- EORDV2[5] + c_EORDV2_in[5] - c_EORDV2_out[5] + c_aging_EORDV2[4] - c_aging_EORDV2[5]
update(EORDV2[6]) <- EORDV2[6] + c_EORDV2_in[6] - c_EORDV2_out[6] + c_aging_EORDV2[5] - c_aging_EORDV2[6]
update(EORDV2[7]) <- EORDV2[7] + c_EORDV2_in[7] - c_EORDV2_out[7] + c_aging_EORDV2[6] - c_aging_EORDV2[7]
update(EORDV2[8]) <- EORDV2[8] + c_EORDV2_in[8] - c_EORDV2_out[8] + c_aging_EORDV2[7] - c_aging_EORDV2[8]
update(EORDV2[9]) <- EORDV2[9] + c_EORDV2_in[9] - c_EORDV2_out[9] + c_aging_EORDV2[8] - c_aging_EORDV2[9]
update(EORDV2[10]) <- EORDV2[10] + c_EORDV2_in[10] - c_EORDV2_out[10] + c_aging_EORDV2[9] - c_aging_EORDV2[10]
update(EORDV2[11]) <- EORDV2[11] + c_EORDV2_in[11] - c_EORDV2_out[11] + c_aging_EORDV2[10] - c_aging_EORDV2[11]
update(EORDV2[12]) <- EORDV2[12] + c_EORDV2_in[12] - c_EORDV2_out[12] + c_aging_EORDV2[11] - c_aging_EORDV2[12]
update(EORDV2[13]) <- EORDV2[13] + c_EORDV2_in[13] - c_EORDV2_out[13] + c_aging_EORDV2[12] - c_aging_EORDV2[13]
update(EORDV2[14]) <- EORDV2[14] + c_EORDV2_in[14] - c_EORDV2_out[14] + c_aging_EORDV2[13] - c_aging_EORDV2[14]
update(EORDV2[15]) <- EORDV2[15] + c_EORDV2_in[15] - c_EORDV2_out[15] + c_aging_EORDV2[14] - c_aging_EORDV2[15]
update(EORDV2[16]) <- EORDV2[16] + c_EORDV2_in[16] - c_EORDV2_out[16] + c_aging_EORDV2[15]

# ----- 2.92. IARWS Compartment -----------------------------------------------

c_IARWS_in[] <- c_EARWS_to_IARWS[i] + c_IARWV1_to_IARWS[i] + c_IARWV2_to_IARWS[i] + c_migration_IARWS[i]
c_IARWS_out[] <- c_IARWS_to_IARWV1[i] + c_IARWS_to_RAS[i] + c_death_IARWS[i]

update(IARWS[1]) <- IARWS[1] + c_IARWS_in[1] - c_IARWS_out[1] - c_aging_IARWS[1]
update(IARWS[2]) <- IARWS[2] + c_IARWS_in[2] - c_IARWS_out[2] + c_aging_IARWS[1] - c_aging_IARWS[2]
update(IARWS[3]) <- IARWS[3] + c_IARWS_in[3] - c_IARWS_out[3] + c_aging_IARWS[2] - c_aging_IARWS[3]
update(IARWS[4]) <- IARWS[4] + c_IARWS_in[4] - c_IARWS_out[4] + c_aging_IARWS[3] - c_aging_IARWS[4]
update(IARWS[5]) <- IARWS[5] + c_IARWS_in[5] - c_IARWS_out[5] + c_aging_IARWS[4] - c_aging_IARWS[5]
update(IARWS[6]) <- IARWS[6] + c_IARWS_in[6] - c_IARWS_out[6] + c_aging_IARWS[5] - c_aging_IARWS[6]
update(IARWS[7]) <- IARWS[7] + c_IARWS_in[7] - c_IARWS_out[7] + c_aging_IARWS[6] - c_aging_IARWS[7]
update(IARWS[8]) <- IARWS[8] + c_IARWS_in[8] - c_IARWS_out[8] + c_aging_IARWS[7] - c_aging_IARWS[8]
update(IARWS[9]) <- IARWS[9] + c_IARWS_in[9] - c_IARWS_out[9] + c_aging_IARWS[8] - c_aging_IARWS[9]
update(IARWS[10]) <- IARWS[10] + c_IARWS_in[10] - c_IARWS_out[10] + c_aging_IARWS[9] - c_aging_IARWS[10]
update(IARWS[11]) <- IARWS[11] + c_IARWS_in[11] - c_IARWS_out[11] + c_aging_IARWS[10] - c_aging_IARWS[11]
update(IARWS[12]) <- IARWS[12] + c_IARWS_in[12] - c_IARWS_out[12] + c_aging_IARWS[11] - c_aging_IARWS[12]
update(IARWS[13]) <- IARWS[13] + c_IARWS_in[13] - c_IARWS_out[13] + c_aging_IARWS[12] - c_aging_IARWS[13]
update(IARWS[14]) <- IARWS[14] + c_IARWS_in[14] - c_IARWS_out[14] + c_aging_IARWS[13] - c_aging_IARWS[14]
update(IARWS[15]) <- IARWS[15] + c_IARWS_in[15] - c_IARWS_out[15] + c_aging_IARWS[14] - c_aging_IARWS[15]
update(IARWS[16]) <- IARWS[16] + c_IARWS_in[16] - c_IARWS_out[16] + c_aging_IARWS[15]

# ----- 2.93. IARWV1 Compartment -----------------------------------------------

c_IARWV1_in[] <- c_EARWV1_to_IARWV1[i] + c_IARWS_to_IARWV1[i] + c_migration_IARWV1[i]
c_IARWV1_out[] <- c_IARWV1_to_IARWS[i] + c_IARWV1_to_IARWV2[i] + c_IARWV1_to_RAV1[i] + c_death_IARWV1[i]

update(IARWV1[1]) <- IARWV1[1] + c_IARWV1_in[1] - c_IARWV1_out[1] - c_aging_IARWV1[1]
update(IARWV1[2]) <- IARWV1[2] + c_IARWV1_in[2] - c_IARWV1_out[2] + c_aging_IARWV1[1] - c_aging_IARWV1[2]
update(IARWV1[3]) <- IARWV1[3] + c_IARWV1_in[3] - c_IARWV1_out[3] + c_aging_IARWV1[2] - c_aging_IARWV1[3]
update(IARWV1[4]) <- IARWV1[4] + c_IARWV1_in[4] - c_IARWV1_out[4] + c_aging_IARWV1[3] - c_aging_IARWV1[4]
update(IARWV1[5]) <- IARWV1[5] + c_IARWV1_in[5] - c_IARWV1_out[5] + c_aging_IARWV1[4] - c_aging_IARWV1[5]
update(IARWV1[6]) <- IARWV1[6] + c_IARWV1_in[6] - c_IARWV1_out[6] + c_aging_IARWV1[5] - c_aging_IARWV1[6]
update(IARWV1[7]) <- IARWV1[7] + c_IARWV1_in[7] - c_IARWV1_out[7] + c_aging_IARWV1[6] - c_aging_IARWV1[7]
update(IARWV1[8]) <- IARWV1[8] + c_IARWV1_in[8] - c_IARWV1_out[8] + c_aging_IARWV1[7] - c_aging_IARWV1[8]
update(IARWV1[9]) <- IARWV1[9] + c_IARWV1_in[9] - c_IARWV1_out[9] + c_aging_IARWV1[8] - c_aging_IARWV1[9]
update(IARWV1[10]) <- IARWV1[10] + c_IARWV1_in[10] - c_IARWV1_out[10] + c_aging_IARWV1[9] - c_aging_IARWV1[10]
update(IARWV1[11]) <- IARWV1[11] + c_IARWV1_in[11] - c_IARWV1_out[11] + c_aging_IARWV1[10] - c_aging_IARWV1[11]
update(IARWV1[12]) <- IARWV1[12] + c_IARWV1_in[12] - c_IARWV1_out[12] + c_aging_IARWV1[11] - c_aging_IARWV1[12]
update(IARWV1[13]) <- IARWV1[13] + c_IARWV1_in[13] - c_IARWV1_out[13] + c_aging_IARWV1[12] - c_aging_IARWV1[13]
update(IARWV1[14]) <- IARWV1[14] + c_IARWV1_in[14] - c_IARWV1_out[14] + c_aging_IARWV1[13] - c_aging_IARWV1[14]
update(IARWV1[15]) <- IARWV1[15] + c_IARWV1_in[15] - c_IARWV1_out[15] + c_aging_IARWV1[14] - c_aging_IARWV1[15]
update(IARWV1[16]) <- IARWV1[16] + c_IARWV1_in[16] - c_IARWV1_out[16] + c_aging_IARWV1[15]

# ----- 2.94. IARWV2 Compartment -----------------------------------------------

c_IARWV2_in[] <- c_EARWV2_to_IARWV2[i] + c_IARWV1_to_IARWV2[i] + c_migration_IARWV2[i]
c_IARWV2_out[] <- c_IARWV2_to_IARWS[i] + c_IARWV2_to_RAV2[i] + c_death_IARWV2[i]

update(IARWV2[1]) <- IARWV2[1] + c_IARWV2_in[1] - c_IARWV2_out[1] - c_aging_IARWV2[1]
update(IARWV2[2]) <- IARWV2[2] + c_IARWV2_in[2] - c_IARWV2_out[2] + c_aging_IARWV2[1] - c_aging_IARWV2[2]
update(IARWV2[3]) <- IARWV2[3] + c_IARWV2_in[3] - c_IARWV2_out[3] + c_aging_IARWV2[2] - c_aging_IARWV2[3]
update(IARWV2[4]) <- IARWV2[4] + c_IARWV2_in[4] - c_IARWV2_out[4] + c_aging_IARWV2[3] - c_aging_IARWV2[4]
update(IARWV2[5]) <- IARWV2[5] + c_IARWV2_in[5] - c_IARWV2_out[5] + c_aging_IARWV2[4] - c_aging_IARWV2[5]
update(IARWV2[6]) <- IARWV2[6] + c_IARWV2_in[6] - c_IARWV2_out[6] + c_aging_IARWV2[5] - c_aging_IARWV2[6]
update(IARWV2[7]) <- IARWV2[7] + c_IARWV2_in[7] - c_IARWV2_out[7] + c_aging_IARWV2[6] - c_aging_IARWV2[7]
update(IARWV2[8]) <- IARWV2[8] + c_IARWV2_in[8] - c_IARWV2_out[8] + c_aging_IARWV2[7] - c_aging_IARWV2[8]
update(IARWV2[9]) <- IARWV2[9] + c_IARWV2_in[9] - c_IARWV2_out[9] + c_aging_IARWV2[8] - c_aging_IARWV2[9]
update(IARWV2[10]) <- IARWV2[10] + c_IARWV2_in[10] - c_IARWV2_out[10] + c_aging_IARWV2[9] - c_aging_IARWV2[10]
update(IARWV2[11]) <- IARWV2[11] + c_IARWV2_in[11] - c_IARWV2_out[11] + c_aging_IARWV2[10] - c_aging_IARWV2[11]
update(IARWV2[12]) <- IARWV2[12] + c_IARWV2_in[12] - c_IARWV2_out[12] + c_aging_IARWV2[11] - c_aging_IARWV2[12]
update(IARWV2[13]) <- IARWV2[13] + c_IARWV2_in[13] - c_IARWV2_out[13] + c_aging_IARWV2[12] - c_aging_IARWV2[13]
update(IARWV2[14]) <- IARWV2[14] + c_IARWV2_in[14] - c_IARWV2_out[14] + c_aging_IARWV2[13] - c_aging_IARWV2[14]
update(IARWV2[15]) <- IARWV2[15] + c_IARWV2_in[15] - c_IARWV2_out[15] + c_aging_IARWV2[14] - c_aging_IARWV2[15]
update(IARWV2[16]) <- IARWV2[16] + c_IARWV2_in[16] - c_IARWV2_out[16] + c_aging_IARWV2[15]

# ----- 2.96. IDRWS Compartment -----------------------------------------------

c_IDRWS_in[] <- c_EDRWS_to_IDRWS[i] + c_IDRWV1_to_IDRWS[i] + c_IDRWV2_to_IDRWS[i] + c_migration_IDRWS[i]
c_IDRWS_out[] <- c_IDRWS_to_IDRWV1[i] + c_IDRWS_to_RDS[i] + c_death_IDRWS[i]

update(IDRWS[1]) <- IDRWS[1] + c_IDRWS_in[1] - c_IDRWS_out[1] - c_aging_IDRWS[1]
update(IDRWS[2]) <- IDRWS[2] + c_IDRWS_in[2] - c_IDRWS_out[2] + c_aging_IDRWS[1] - c_aging_IDRWS[2]
update(IDRWS[3]) <- IDRWS[3] + c_IDRWS_in[3] - c_IDRWS_out[3] + c_aging_IDRWS[2] - c_aging_IDRWS[3]
update(IDRWS[4]) <- IDRWS[4] + c_IDRWS_in[4] - c_IDRWS_out[4] + c_aging_IDRWS[3] - c_aging_IDRWS[4]
update(IDRWS[5]) <- IDRWS[5] + c_IDRWS_in[5] - c_IDRWS_out[5] + c_aging_IDRWS[4] - c_aging_IDRWS[5]
update(IDRWS[6]) <- IDRWS[6] + c_IDRWS_in[6] - c_IDRWS_out[6] + c_aging_IDRWS[5] - c_aging_IDRWS[6]
update(IDRWS[7]) <- IDRWS[7] + c_IDRWS_in[7] - c_IDRWS_out[7] + c_aging_IDRWS[6] - c_aging_IDRWS[7]
update(IDRWS[8]) <- IDRWS[8] + c_IDRWS_in[8] - c_IDRWS_out[8] + c_aging_IDRWS[7] - c_aging_IDRWS[8]
update(IDRWS[9]) <- IDRWS[9] + c_IDRWS_in[9] - c_IDRWS_out[9] + c_aging_IDRWS[8] - c_aging_IDRWS[9]
update(IDRWS[10]) <- IDRWS[10] + c_IDRWS_in[10] - c_IDRWS_out[10] + c_aging_IDRWS[9] - c_aging_IDRWS[10]
update(IDRWS[11]) <- IDRWS[11] + c_IDRWS_in[11] - c_IDRWS_out[11] + c_aging_IDRWS[10] - c_aging_IDRWS[11]
update(IDRWS[12]) <- IDRWS[12] + c_IDRWS_in[12] - c_IDRWS_out[12] + c_aging_IDRWS[11] - c_aging_IDRWS[12]
update(IDRWS[13]) <- IDRWS[13] + c_IDRWS_in[13] - c_IDRWS_out[13] + c_aging_IDRWS[12] - c_aging_IDRWS[13]
update(IDRWS[14]) <- IDRWS[14] + c_IDRWS_in[14] - c_IDRWS_out[14] + c_aging_IDRWS[13] - c_aging_IDRWS[14]
update(IDRWS[15]) <- IDRWS[15] + c_IDRWS_in[15] - c_IDRWS_out[15] + c_aging_IDRWS[14] - c_aging_IDRWS[15]
update(IDRWS[16]) <- IDRWS[16] + c_IDRWS_in[16] - c_IDRWS_out[16] + c_aging_IDRWS[15]

# ----- 2.97. IDRWV1 Compartment -----------------------------------------------

c_IDRWV1_in[] <- c_EDRWV1_to_IDRWV1[i] + c_IDRWS_to_IDRWV1[i] + c_migration_IDRWV1[i]
c_IDRWV1_out[] <- c_IDRWV1_to_IDRWS[i] + c_IDRWV1_to_IDRWV2[i] + c_IDRWV1_to_RDV1[i] + c_death_IDRWV1[i]

update(IDRWV1[1]) <- IDRWV1[1] + c_IDRWV1_in[1] - c_IDRWV1_out[1] - c_aging_IDRWV1[1]
update(IDRWV1[2]) <- IDRWV1[2] + c_IDRWV1_in[2] - c_IDRWV1_out[2] + c_aging_IDRWV1[1] - c_aging_IDRWV1[2]
update(IDRWV1[3]) <- IDRWV1[3] + c_IDRWV1_in[3] - c_IDRWV1_out[3] + c_aging_IDRWV1[2] - c_aging_IDRWV1[3]
update(IDRWV1[4]) <- IDRWV1[4] + c_IDRWV1_in[4] - c_IDRWV1_out[4] + c_aging_IDRWV1[3] - c_aging_IDRWV1[4]
update(IDRWV1[5]) <- IDRWV1[5] + c_IDRWV1_in[5] - c_IDRWV1_out[5] + c_aging_IDRWV1[4] - c_aging_IDRWV1[5]
update(IDRWV1[6]) <- IDRWV1[6] + c_IDRWV1_in[6] - c_IDRWV1_out[6] + c_aging_IDRWV1[5] - c_aging_IDRWV1[6]
update(IDRWV1[7]) <- IDRWV1[7] + c_IDRWV1_in[7] - c_IDRWV1_out[7] + c_aging_IDRWV1[6] - c_aging_IDRWV1[7]
update(IDRWV1[8]) <- IDRWV1[8] + c_IDRWV1_in[8] - c_IDRWV1_out[8] + c_aging_IDRWV1[7] - c_aging_IDRWV1[8]
update(IDRWV1[9]) <- IDRWV1[9] + c_IDRWV1_in[9] - c_IDRWV1_out[9] + c_aging_IDRWV1[8] - c_aging_IDRWV1[9]
update(IDRWV1[10]) <- IDRWV1[10] + c_IDRWV1_in[10] - c_IDRWV1_out[10] + c_aging_IDRWV1[9] - c_aging_IDRWV1[10]
update(IDRWV1[11]) <- IDRWV1[11] + c_IDRWV1_in[11] - c_IDRWV1_out[11] + c_aging_IDRWV1[10] - c_aging_IDRWV1[11]
update(IDRWV1[12]) <- IDRWV1[12] + c_IDRWV1_in[12] - c_IDRWV1_out[12] + c_aging_IDRWV1[11] - c_aging_IDRWV1[12]
update(IDRWV1[13]) <- IDRWV1[13] + c_IDRWV1_in[13] - c_IDRWV1_out[13] + c_aging_IDRWV1[12] - c_aging_IDRWV1[13]
update(IDRWV1[14]) <- IDRWV1[14] + c_IDRWV1_in[14] - c_IDRWV1_out[14] + c_aging_IDRWV1[13] - c_aging_IDRWV1[14]
update(IDRWV1[15]) <- IDRWV1[15] + c_IDRWV1_in[15] - c_IDRWV1_out[15] + c_aging_IDRWV1[14] - c_aging_IDRWV1[15]
update(IDRWV1[16]) <- IDRWV1[16] + c_IDRWV1_in[16] - c_IDRWV1_out[16] + c_aging_IDRWV1[15]

# ----- 2.98. IDRWV2 Compartment -----------------------------------------------

c_IDRWV2_in[] <- c_EDRWV2_to_IDRWV2[i] + c_IDRWV1_to_IDRWV2[i] + c_migration_IDRWV2[i]
c_IDRWV2_out[] <- c_IDRWV2_to_IDRWS[i] + c_IDRWV2_to_RDV2[i] + c_death_IDRWV2[i]

update(IDRWV2[1]) <- IDRWV2[1] + c_IDRWV2_in[1] - c_IDRWV2_out[1] - c_aging_IDRWV2[1]
update(IDRWV2[2]) <- IDRWV2[2] + c_IDRWV2_in[2] - c_IDRWV2_out[2] + c_aging_IDRWV2[1] - c_aging_IDRWV2[2]
update(IDRWV2[3]) <- IDRWV2[3] + c_IDRWV2_in[3] - c_IDRWV2_out[3] + c_aging_IDRWV2[2] - c_aging_IDRWV2[3]
update(IDRWV2[4]) <- IDRWV2[4] + c_IDRWV2_in[4] - c_IDRWV2_out[4] + c_aging_IDRWV2[3] - c_aging_IDRWV2[4]
update(IDRWV2[5]) <- IDRWV2[5] + c_IDRWV2_in[5] - c_IDRWV2_out[5] + c_aging_IDRWV2[4] - c_aging_IDRWV2[5]
update(IDRWV2[6]) <- IDRWV2[6] + c_IDRWV2_in[6] - c_IDRWV2_out[6] + c_aging_IDRWV2[5] - c_aging_IDRWV2[6]
update(IDRWV2[7]) <- IDRWV2[7] + c_IDRWV2_in[7] - c_IDRWV2_out[7] + c_aging_IDRWV2[6] - c_aging_IDRWV2[7]
update(IDRWV2[8]) <- IDRWV2[8] + c_IDRWV2_in[8] - c_IDRWV2_out[8] + c_aging_IDRWV2[7] - c_aging_IDRWV2[8]
update(IDRWV2[9]) <- IDRWV2[9] + c_IDRWV2_in[9] - c_IDRWV2_out[9] + c_aging_IDRWV2[8] - c_aging_IDRWV2[9]
update(IDRWV2[10]) <- IDRWV2[10] + c_IDRWV2_in[10] - c_IDRWV2_out[10] + c_aging_IDRWV2[9] - c_aging_IDRWV2[10]
update(IDRWV2[11]) <- IDRWV2[11] + c_IDRWV2_in[11] - c_IDRWV2_out[11] + c_aging_IDRWV2[10] - c_aging_IDRWV2[11]
update(IDRWV2[12]) <- IDRWV2[12] + c_IDRWV2_in[12] - c_IDRWV2_out[12] + c_aging_IDRWV2[11] - c_aging_IDRWV2[12]
update(IDRWV2[13]) <- IDRWV2[13] + c_IDRWV2_in[13] - c_IDRWV2_out[13] + c_aging_IDRWV2[12] - c_aging_IDRWV2[13]
update(IDRWV2[14]) <- IDRWV2[14] + c_IDRWV2_in[14] - c_IDRWV2_out[14] + c_aging_IDRWV2[13] - c_aging_IDRWV2[14]
update(IDRWV2[15]) <- IDRWV2[15] + c_IDRWV2_in[15] - c_IDRWV2_out[15] + c_aging_IDRWV2[14] - c_aging_IDRWV2[15]
update(IDRWV2[16]) <- IDRWV2[16] + c_IDRWV2_in[16] - c_IDRWV2_out[16] + c_aging_IDRWV2[15]

# ----- 2.100. IDRAS Compartment -----------------------------------------------

c_IDRAS_in[] <- c_EDRAS_to_IDRAS[i] + c_IDRAV1_to_IDRAS[i] + c_IDRAV2_to_IDRAS[i] + c_migration_IDRAS[i]
c_IDRAS_out[] <- c_IDRAS_to_IDRAV1[i] + c_IDRAS_to_RDS[i] + c_death_IDRAS[i]

update(IDRAS[1]) <- IDRAS[1] + c_IDRAS_in[1] - c_IDRAS_out[1] - c_aging_IDRAS[1]
update(IDRAS[2]) <- IDRAS[2] + c_IDRAS_in[2] - c_IDRAS_out[2] + c_aging_IDRAS[1] - c_aging_IDRAS[2]
update(IDRAS[3]) <- IDRAS[3] + c_IDRAS_in[3] - c_IDRAS_out[3] + c_aging_IDRAS[2] - c_aging_IDRAS[3]
update(IDRAS[4]) <- IDRAS[4] + c_IDRAS_in[4] - c_IDRAS_out[4] + c_aging_IDRAS[3] - c_aging_IDRAS[4]
update(IDRAS[5]) <- IDRAS[5] + c_IDRAS_in[5] - c_IDRAS_out[5] + c_aging_IDRAS[4] - c_aging_IDRAS[5]
update(IDRAS[6]) <- IDRAS[6] + c_IDRAS_in[6] - c_IDRAS_out[6] + c_aging_IDRAS[5] - c_aging_IDRAS[6]
update(IDRAS[7]) <- IDRAS[7] + c_IDRAS_in[7] - c_IDRAS_out[7] + c_aging_IDRAS[6] - c_aging_IDRAS[7]
update(IDRAS[8]) <- IDRAS[8] + c_IDRAS_in[8] - c_IDRAS_out[8] + c_aging_IDRAS[7] - c_aging_IDRAS[8]
update(IDRAS[9]) <- IDRAS[9] + c_IDRAS_in[9] - c_IDRAS_out[9] + c_aging_IDRAS[8] - c_aging_IDRAS[9]
update(IDRAS[10]) <- IDRAS[10] + c_IDRAS_in[10] - c_IDRAS_out[10] + c_aging_IDRAS[9] - c_aging_IDRAS[10]
update(IDRAS[11]) <- IDRAS[11] + c_IDRAS_in[11] - c_IDRAS_out[11] + c_aging_IDRAS[10] - c_aging_IDRAS[11]
update(IDRAS[12]) <- IDRAS[12] + c_IDRAS_in[12] - c_IDRAS_out[12] + c_aging_IDRAS[11] - c_aging_IDRAS[12]
update(IDRAS[13]) <- IDRAS[13] + c_IDRAS_in[13] - c_IDRAS_out[13] + c_aging_IDRAS[12] - c_aging_IDRAS[13]
update(IDRAS[14]) <- IDRAS[14] + c_IDRAS_in[14] - c_IDRAS_out[14] + c_aging_IDRAS[13] - c_aging_IDRAS[14]
update(IDRAS[15]) <- IDRAS[15] + c_IDRAS_in[15] - c_IDRAS_out[15] + c_aging_IDRAS[14] - c_aging_IDRAS[15]
update(IDRAS[16]) <- IDRAS[16] + c_IDRAS_in[16] - c_IDRAS_out[16] + c_aging_IDRAS[15]

# ----- 2.101. IDRAV1 Compartment -----------------------------------------------

c_IDRAV1_in[] <- c_EDRAV1_to_IDRAV1[i] + c_IDRAS_to_IDRAV1[i] + c_migration_IDRAV1[i]
c_IDRAV1_out[] <- c_IDRAV1_to_IDRAS[i] + c_IDRAV1_to_IDRAV2[i] + c_IDRAV1_to_RDV1[i] + c_death_IDRAV1[i]

update(IDRAV1[1]) <- IDRAV1[1] + c_IDRAV1_in[1] - c_IDRAV1_out[1] - c_aging_IDRAV1[1]
update(IDRAV1[2]) <- IDRAV1[2] + c_IDRAV1_in[2] - c_IDRAV1_out[2] + c_aging_IDRAV1[1] - c_aging_IDRAV1[2]
update(IDRAV1[3]) <- IDRAV1[3] + c_IDRAV1_in[3] - c_IDRAV1_out[3] + c_aging_IDRAV1[2] - c_aging_IDRAV1[3]
update(IDRAV1[4]) <- IDRAV1[4] + c_IDRAV1_in[4] - c_IDRAV1_out[4] + c_aging_IDRAV1[3] - c_aging_IDRAV1[4]
update(IDRAV1[5]) <- IDRAV1[5] + c_IDRAV1_in[5] - c_IDRAV1_out[5] + c_aging_IDRAV1[4] - c_aging_IDRAV1[5]
update(IDRAV1[6]) <- IDRAV1[6] + c_IDRAV1_in[6] - c_IDRAV1_out[6] + c_aging_IDRAV1[5] - c_aging_IDRAV1[6]
update(IDRAV1[7]) <- IDRAV1[7] + c_IDRAV1_in[7] - c_IDRAV1_out[7] + c_aging_IDRAV1[6] - c_aging_IDRAV1[7]
update(IDRAV1[8]) <- IDRAV1[8] + c_IDRAV1_in[8] - c_IDRAV1_out[8] + c_aging_IDRAV1[7] - c_aging_IDRAV1[8]
update(IDRAV1[9]) <- IDRAV1[9] + c_IDRAV1_in[9] - c_IDRAV1_out[9] + c_aging_IDRAV1[8] - c_aging_IDRAV1[9]
update(IDRAV1[10]) <- IDRAV1[10] + c_IDRAV1_in[10] - c_IDRAV1_out[10] + c_aging_IDRAV1[9] - c_aging_IDRAV1[10]
update(IDRAV1[11]) <- IDRAV1[11] + c_IDRAV1_in[11] - c_IDRAV1_out[11] + c_aging_IDRAV1[10] - c_aging_IDRAV1[11]
update(IDRAV1[12]) <- IDRAV1[12] + c_IDRAV1_in[12] - c_IDRAV1_out[12] + c_aging_IDRAV1[11] - c_aging_IDRAV1[12]
update(IDRAV1[13]) <- IDRAV1[13] + c_IDRAV1_in[13] - c_IDRAV1_out[13] + c_aging_IDRAV1[12] - c_aging_IDRAV1[13]
update(IDRAV1[14]) <- IDRAV1[14] + c_IDRAV1_in[14] - c_IDRAV1_out[14] + c_aging_IDRAV1[13] - c_aging_IDRAV1[14]
update(IDRAV1[15]) <- IDRAV1[15] + c_IDRAV1_in[15] - c_IDRAV1_out[15] + c_aging_IDRAV1[14] - c_aging_IDRAV1[15]
update(IDRAV1[16]) <- IDRAV1[16] + c_IDRAV1_in[16] - c_IDRAV1_out[16] + c_aging_IDRAV1[15]

# ----- 2.102. IDRAV2 Compartment -----------------------------------------------

c_IDRAV2_in[] <- c_EDRAV2_to_IDRAV2[i] + c_IDRAV1_to_IDRAV2[i] + c_migration_IDRAV2[i]
c_IDRAV2_out[] <- c_IDRAV2_to_IDRAS[i] + c_IDRAV2_to_RDV2[i] + c_death_IDRAV2[i]

update(IDRAV2[1]) <- IDRAV2[1] + c_IDRAV2_in[1] - c_IDRAV2_out[1] - c_aging_IDRAV2[1]
update(IDRAV2[2]) <- IDRAV2[2] + c_IDRAV2_in[2] - c_IDRAV2_out[2] + c_aging_IDRAV2[1] - c_aging_IDRAV2[2]
update(IDRAV2[3]) <- IDRAV2[3] + c_IDRAV2_in[3] - c_IDRAV2_out[3] + c_aging_IDRAV2[2] - c_aging_IDRAV2[3]
update(IDRAV2[4]) <- IDRAV2[4] + c_IDRAV2_in[4] - c_IDRAV2_out[4] + c_aging_IDRAV2[3] - c_aging_IDRAV2[4]
update(IDRAV2[5]) <- IDRAV2[5] + c_IDRAV2_in[5] - c_IDRAV2_out[5] + c_aging_IDRAV2[4] - c_aging_IDRAV2[5]
update(IDRAV2[6]) <- IDRAV2[6] + c_IDRAV2_in[6] - c_IDRAV2_out[6] + c_aging_IDRAV2[5] - c_aging_IDRAV2[6]
update(IDRAV2[7]) <- IDRAV2[7] + c_IDRAV2_in[7] - c_IDRAV2_out[7] + c_aging_IDRAV2[6] - c_aging_IDRAV2[7]
update(IDRAV2[8]) <- IDRAV2[8] + c_IDRAV2_in[8] - c_IDRAV2_out[8] + c_aging_IDRAV2[7] - c_aging_IDRAV2[8]
update(IDRAV2[9]) <- IDRAV2[9] + c_IDRAV2_in[9] - c_IDRAV2_out[9] + c_aging_IDRAV2[8] - c_aging_IDRAV2[9]
update(IDRAV2[10]) <- IDRAV2[10] + c_IDRAV2_in[10] - c_IDRAV2_out[10] + c_aging_IDRAV2[9] - c_aging_IDRAV2[10]
update(IDRAV2[11]) <- IDRAV2[11] + c_IDRAV2_in[11] - c_IDRAV2_out[11] + c_aging_IDRAV2[10] - c_aging_IDRAV2[11]
update(IDRAV2[12]) <- IDRAV2[12] + c_IDRAV2_in[12] - c_IDRAV2_out[12] + c_aging_IDRAV2[11] - c_aging_IDRAV2[12]
update(IDRAV2[13]) <- IDRAV2[13] + c_IDRAV2_in[13] - c_IDRAV2_out[13] + c_aging_IDRAV2[12] - c_aging_IDRAV2[13]
update(IDRAV2[14]) <- IDRAV2[14] + c_IDRAV2_in[14] - c_IDRAV2_out[14] + c_aging_IDRAV2[13] - c_aging_IDRAV2[14]
update(IDRAV2[15]) <- IDRAV2[15] + c_IDRAV2_in[15] - c_IDRAV2_out[15] + c_aging_IDRAV2[14] - c_aging_IDRAV2[15]
update(IDRAV2[16]) <- IDRAV2[16] + c_IDRAV2_in[16] - c_IDRAV2_out[16] + c_aging_IDRAV2[15]

# ----- 2.104. IORWS Compartment -----------------------------------------------

c_IORWS_in[] <- c_EORWS_to_IORWS[i] + c_IORWV1_to_IORWS[i] + c_IORWV2_to_IORWS[i] + c_migration_IORWS[i]
c_IORWS_out[] <- c_IORWS_to_IORWV1[i] + c_IORWS_to_ROS[i] + c_death_IORWS[i]

update(IORWS[1]) <- IORWS[1] + c_IORWS_in[1] - c_IORWS_out[1] - c_aging_IORWS[1]
update(IORWS[2]) <- IORWS[2] + c_IORWS_in[2] - c_IORWS_out[2] + c_aging_IORWS[1] - c_aging_IORWS[2]
update(IORWS[3]) <- IORWS[3] + c_IORWS_in[3] - c_IORWS_out[3] + c_aging_IORWS[2] - c_aging_IORWS[3]
update(IORWS[4]) <- IORWS[4] + c_IORWS_in[4] - c_IORWS_out[4] + c_aging_IORWS[3] - c_aging_IORWS[4]
update(IORWS[5]) <- IORWS[5] + c_IORWS_in[5] - c_IORWS_out[5] + c_aging_IORWS[4] - c_aging_IORWS[5]
update(IORWS[6]) <- IORWS[6] + c_IORWS_in[6] - c_IORWS_out[6] + c_aging_IORWS[5] - c_aging_IORWS[6]
update(IORWS[7]) <- IORWS[7] + c_IORWS_in[7] - c_IORWS_out[7] + c_aging_IORWS[6] - c_aging_IORWS[7]
update(IORWS[8]) <- IORWS[8] + c_IORWS_in[8] - c_IORWS_out[8] + c_aging_IORWS[7] - c_aging_IORWS[8]
update(IORWS[9]) <- IORWS[9] + c_IORWS_in[9] - c_IORWS_out[9] + c_aging_IORWS[8] - c_aging_IORWS[9]
update(IORWS[10]) <- IORWS[10] + c_IORWS_in[10] - c_IORWS_out[10] + c_aging_IORWS[9] - c_aging_IORWS[10]
update(IORWS[11]) <- IORWS[11] + c_IORWS_in[11] - c_IORWS_out[11] + c_aging_IORWS[10] - c_aging_IORWS[11]
update(IORWS[12]) <- IORWS[12] + c_IORWS_in[12] - c_IORWS_out[12] + c_aging_IORWS[11] - c_aging_IORWS[12]
update(IORWS[13]) <- IORWS[13] + c_IORWS_in[13] - c_IORWS_out[13] + c_aging_IORWS[12] - c_aging_IORWS[13]
update(IORWS[14]) <- IORWS[14] + c_IORWS_in[14] - c_IORWS_out[14] + c_aging_IORWS[13] - c_aging_IORWS[14]
update(IORWS[15]) <- IORWS[15] + c_IORWS_in[15] - c_IORWS_out[15] + c_aging_IORWS[14] - c_aging_IORWS[15]
update(IORWS[16]) <- IORWS[16] + c_IORWS_in[16] - c_IORWS_out[16] + c_aging_IORWS[15]

# ----- 2.105. IORWV1 Compartment -----------------------------------------------

c_IORWV1_in[] <- c_EORWV1_to_IORWV1[i] + c_IORWS_to_IORWV1[i] + c_migration_IORWV1[i]
c_IORWV1_out[] <- c_IORWV1_to_IORWS[i] + c_IORWV1_to_IORWV2[i] + c_IORWV1_to_ROV1[i] + c_death_IORWV1[i]

update(IORWV1[1]) <- IORWV1[1] + c_IORWV1_in[1] - c_IORWV1_out[1] - c_aging_IORWV1[1]
update(IORWV1[2]) <- IORWV1[2] + c_IORWV1_in[2] - c_IORWV1_out[2] + c_aging_IORWV1[1] - c_aging_IORWV1[2]
update(IORWV1[3]) <- IORWV1[3] + c_IORWV1_in[3] - c_IORWV1_out[3] + c_aging_IORWV1[2] - c_aging_IORWV1[3]
update(IORWV1[4]) <- IORWV1[4] + c_IORWV1_in[4] - c_IORWV1_out[4] + c_aging_IORWV1[3] - c_aging_IORWV1[4]
update(IORWV1[5]) <- IORWV1[5] + c_IORWV1_in[5] - c_IORWV1_out[5] + c_aging_IORWV1[4] - c_aging_IORWV1[5]
update(IORWV1[6]) <- IORWV1[6] + c_IORWV1_in[6] - c_IORWV1_out[6] + c_aging_IORWV1[5] - c_aging_IORWV1[6]
update(IORWV1[7]) <- IORWV1[7] + c_IORWV1_in[7] - c_IORWV1_out[7] + c_aging_IORWV1[6] - c_aging_IORWV1[7]
update(IORWV1[8]) <- IORWV1[8] + c_IORWV1_in[8] - c_IORWV1_out[8] + c_aging_IORWV1[7] - c_aging_IORWV1[8]
update(IORWV1[9]) <- IORWV1[9] + c_IORWV1_in[9] - c_IORWV1_out[9] + c_aging_IORWV1[8] - c_aging_IORWV1[9]
update(IORWV1[10]) <- IORWV1[10] + c_IORWV1_in[10] - c_IORWV1_out[10] + c_aging_IORWV1[9] - c_aging_IORWV1[10]
update(IORWV1[11]) <- IORWV1[11] + c_IORWV1_in[11] - c_IORWV1_out[11] + c_aging_IORWV1[10] - c_aging_IORWV1[11]
update(IORWV1[12]) <- IORWV1[12] + c_IORWV1_in[12] - c_IORWV1_out[12] + c_aging_IORWV1[11] - c_aging_IORWV1[12]
update(IORWV1[13]) <- IORWV1[13] + c_IORWV1_in[13] - c_IORWV1_out[13] + c_aging_IORWV1[12] - c_aging_IORWV1[13]
update(IORWV1[14]) <- IORWV1[14] + c_IORWV1_in[14] - c_IORWV1_out[14] + c_aging_IORWV1[13] - c_aging_IORWV1[14]
update(IORWV1[15]) <- IORWV1[15] + c_IORWV1_in[15] - c_IORWV1_out[15] + c_aging_IORWV1[14] - c_aging_IORWV1[15]
update(IORWV1[16]) <- IORWV1[16] + c_IORWV1_in[16] - c_IORWV1_out[16] + c_aging_IORWV1[15]

# ----- 2.106. IORWV2 Compartment -----------------------------------------------

c_IORWV2_in[] <- c_EORWV2_to_IORWV2[i] + c_IORWV1_to_IORWV2[i] + c_migration_IORWV2[i]
c_IORWV2_out[] <- c_IORWV2_to_IORWS[i] + c_IORWV2_to_ROV2[i] + c_death_IORWV2[i]

update(IORWV2[1]) <- IORWV2[1] + c_IORWV2_in[1] - c_IORWV2_out[1] - c_aging_IORWV2[1]
update(IORWV2[2]) <- IORWV2[2] + c_IORWV2_in[2] - c_IORWV2_out[2] + c_aging_IORWV2[1] - c_aging_IORWV2[2]
update(IORWV2[3]) <- IORWV2[3] + c_IORWV2_in[3] - c_IORWV2_out[3] + c_aging_IORWV2[2] - c_aging_IORWV2[3]
update(IORWV2[4]) <- IORWV2[4] + c_IORWV2_in[4] - c_IORWV2_out[4] + c_aging_IORWV2[3] - c_aging_IORWV2[4]
update(IORWV2[5]) <- IORWV2[5] + c_IORWV2_in[5] - c_IORWV2_out[5] + c_aging_IORWV2[4] - c_aging_IORWV2[5]
update(IORWV2[6]) <- IORWV2[6] + c_IORWV2_in[6] - c_IORWV2_out[6] + c_aging_IORWV2[5] - c_aging_IORWV2[6]
update(IORWV2[7]) <- IORWV2[7] + c_IORWV2_in[7] - c_IORWV2_out[7] + c_aging_IORWV2[6] - c_aging_IORWV2[7]
update(IORWV2[8]) <- IORWV2[8] + c_IORWV2_in[8] - c_IORWV2_out[8] + c_aging_IORWV2[7] - c_aging_IORWV2[8]
update(IORWV2[9]) <- IORWV2[9] + c_IORWV2_in[9] - c_IORWV2_out[9] + c_aging_IORWV2[8] - c_aging_IORWV2[9]
update(IORWV2[10]) <- IORWV2[10] + c_IORWV2_in[10] - c_IORWV2_out[10] + c_aging_IORWV2[9] - c_aging_IORWV2[10]
update(IORWV2[11]) <- IORWV2[11] + c_IORWV2_in[11] - c_IORWV2_out[11] + c_aging_IORWV2[10] - c_aging_IORWV2[11]
update(IORWV2[12]) <- IORWV2[12] + c_IORWV2_in[12] - c_IORWV2_out[12] + c_aging_IORWV2[11] - c_aging_IORWV2[12]
update(IORWV2[13]) <- IORWV2[13] + c_IORWV2_in[13] - c_IORWV2_out[13] + c_aging_IORWV2[12] - c_aging_IORWV2[13]
update(IORWV2[14]) <- IORWV2[14] + c_IORWV2_in[14] - c_IORWV2_out[14] + c_aging_IORWV2[13] - c_aging_IORWV2[14]
update(IORWV2[15]) <- IORWV2[15] + c_IORWV2_in[15] - c_IORWV2_out[15] + c_aging_IORWV2[14] - c_aging_IORWV2[15]
update(IORWV2[16]) <- IORWV2[16] + c_IORWV2_in[16] - c_IORWV2_out[16] + c_aging_IORWV2[15]

# ----- 2.108. IORAS Compartment -----------------------------------------------

c_IORAS_in[] <- c_EORAS_to_IORAS[i] + c_IORAV1_to_IORAS[i] + c_IORAV2_to_IORAS[i] + c_migration_IORAS[i]
c_IORAS_out[] <- c_IORAS_to_IORAV1[i] + c_IORAS_to_ROS[i] + c_death_IORAS[i]

update(IORAS[1]) <- IORAS[1] + c_IORAS_in[1] - c_IORAS_out[1] - c_aging_IORAS[1]
update(IORAS[2]) <- IORAS[2] + c_IORAS_in[2] - c_IORAS_out[2] + c_aging_IORAS[1] - c_aging_IORAS[2]
update(IORAS[3]) <- IORAS[3] + c_IORAS_in[3] - c_IORAS_out[3] + c_aging_IORAS[2] - c_aging_IORAS[3]
update(IORAS[4]) <- IORAS[4] + c_IORAS_in[4] - c_IORAS_out[4] + c_aging_IORAS[3] - c_aging_IORAS[4]
update(IORAS[5]) <- IORAS[5] + c_IORAS_in[5] - c_IORAS_out[5] + c_aging_IORAS[4] - c_aging_IORAS[5]
update(IORAS[6]) <- IORAS[6] + c_IORAS_in[6] - c_IORAS_out[6] + c_aging_IORAS[5] - c_aging_IORAS[6]
update(IORAS[7]) <- IORAS[7] + c_IORAS_in[7] - c_IORAS_out[7] + c_aging_IORAS[6] - c_aging_IORAS[7]
update(IORAS[8]) <- IORAS[8] + c_IORAS_in[8] - c_IORAS_out[8] + c_aging_IORAS[7] - c_aging_IORAS[8]
update(IORAS[9]) <- IORAS[9] + c_IORAS_in[9] - c_IORAS_out[9] + c_aging_IORAS[8] - c_aging_IORAS[9]
update(IORAS[10]) <- IORAS[10] + c_IORAS_in[10] - c_IORAS_out[10] + c_aging_IORAS[9] - c_aging_IORAS[10]
update(IORAS[11]) <- IORAS[11] + c_IORAS_in[11] - c_IORAS_out[11] + c_aging_IORAS[10] - c_aging_IORAS[11]
update(IORAS[12]) <- IORAS[12] + c_IORAS_in[12] - c_IORAS_out[12] + c_aging_IORAS[11] - c_aging_IORAS[12]
update(IORAS[13]) <- IORAS[13] + c_IORAS_in[13] - c_IORAS_out[13] + c_aging_IORAS[12] - c_aging_IORAS[13]
update(IORAS[14]) <- IORAS[14] + c_IORAS_in[14] - c_IORAS_out[14] + c_aging_IORAS[13] - c_aging_IORAS[14]
update(IORAS[15]) <- IORAS[15] + c_IORAS_in[15] - c_IORAS_out[15] + c_aging_IORAS[14] - c_aging_IORAS[15]
update(IORAS[16]) <- IORAS[16] + c_IORAS_in[16] - c_IORAS_out[16] + c_aging_IORAS[15]

# ----- 2.109. IORAV1 Compartment -----------------------------------------------

c_IORAV1_in[] <- c_EORAV1_to_IORAV1[i] + c_IORAS_to_IORAV1[i] + c_migration_IORAV1[i]
c_IORAV1_out[] <- c_IORAV1_to_IORAS[i] + c_IORAV1_to_IORAV2[i] + c_IORAV1_to_ROV1[i] + c_death_IORAV1[i]

update(IORAV1[1]) <- IORAV1[1] + c_IORAV1_in[1] - c_IORAV1_out[1] - c_aging_IORAV1[1]
update(IORAV1[2]) <- IORAV1[2] + c_IORAV1_in[2] - c_IORAV1_out[2] + c_aging_IORAV1[1] - c_aging_IORAV1[2]
update(IORAV1[3]) <- IORAV1[3] + c_IORAV1_in[3] - c_IORAV1_out[3] + c_aging_IORAV1[2] - c_aging_IORAV1[3]
update(IORAV1[4]) <- IORAV1[4] + c_IORAV1_in[4] - c_IORAV1_out[4] + c_aging_IORAV1[3] - c_aging_IORAV1[4]
update(IORAV1[5]) <- IORAV1[5] + c_IORAV1_in[5] - c_IORAV1_out[5] + c_aging_IORAV1[4] - c_aging_IORAV1[5]
update(IORAV1[6]) <- IORAV1[6] + c_IORAV1_in[6] - c_IORAV1_out[6] + c_aging_IORAV1[5] - c_aging_IORAV1[6]
update(IORAV1[7]) <- IORAV1[7] + c_IORAV1_in[7] - c_IORAV1_out[7] + c_aging_IORAV1[6] - c_aging_IORAV1[7]
update(IORAV1[8]) <- IORAV1[8] + c_IORAV1_in[8] - c_IORAV1_out[8] + c_aging_IORAV1[7] - c_aging_IORAV1[8]
update(IORAV1[9]) <- IORAV1[9] + c_IORAV1_in[9] - c_IORAV1_out[9] + c_aging_IORAV1[8] - c_aging_IORAV1[9]
update(IORAV1[10]) <- IORAV1[10] + c_IORAV1_in[10] - c_IORAV1_out[10] + c_aging_IORAV1[9] - c_aging_IORAV1[10]
update(IORAV1[11]) <- IORAV1[11] + c_IORAV1_in[11] - c_IORAV1_out[11] + c_aging_IORAV1[10] - c_aging_IORAV1[11]
update(IORAV1[12]) <- IORAV1[12] + c_IORAV1_in[12] - c_IORAV1_out[12] + c_aging_IORAV1[11] - c_aging_IORAV1[12]
update(IORAV1[13]) <- IORAV1[13] + c_IORAV1_in[13] - c_IORAV1_out[13] + c_aging_IORAV1[12] - c_aging_IORAV1[13]
update(IORAV1[14]) <- IORAV1[14] + c_IORAV1_in[14] - c_IORAV1_out[14] + c_aging_IORAV1[13] - c_aging_IORAV1[14]
update(IORAV1[15]) <- IORAV1[15] + c_IORAV1_in[15] - c_IORAV1_out[15] + c_aging_IORAV1[14] - c_aging_IORAV1[15]
update(IORAV1[16]) <- IORAV1[16] + c_IORAV1_in[16] - c_IORAV1_out[16] + c_aging_IORAV1[15]

# ----- 2.110. IORAV2 Compartment -----------------------------------------------

c_IORAV2_in[] <- c_EORAV2_to_IORAV2[i] + c_IORAV1_to_IORAV2[i] + c_migration_IORAV2[i]
c_IORAV2_out[] <- c_IORAV2_to_IORAS[i] + c_IORAV2_to_ROV2[i] + c_death_IORAV2[i]

update(IORAV2[1]) <- IORAV2[1] + c_IORAV2_in[1] - c_IORAV2_out[1] - c_aging_IORAV2[1]
update(IORAV2[2]) <- IORAV2[2] + c_IORAV2_in[2] - c_IORAV2_out[2] + c_aging_IORAV2[1] - c_aging_IORAV2[2]
update(IORAV2[3]) <- IORAV2[3] + c_IORAV2_in[3] - c_IORAV2_out[3] + c_aging_IORAV2[2] - c_aging_IORAV2[3]
update(IORAV2[4]) <- IORAV2[4] + c_IORAV2_in[4] - c_IORAV2_out[4] + c_aging_IORAV2[3] - c_aging_IORAV2[4]
update(IORAV2[5]) <- IORAV2[5] + c_IORAV2_in[5] - c_IORAV2_out[5] + c_aging_IORAV2[4] - c_aging_IORAV2[5]
update(IORAV2[6]) <- IORAV2[6] + c_IORAV2_in[6] - c_IORAV2_out[6] + c_aging_IORAV2[5] - c_aging_IORAV2[6]
update(IORAV2[7]) <- IORAV2[7] + c_IORAV2_in[7] - c_IORAV2_out[7] + c_aging_IORAV2[6] - c_aging_IORAV2[7]
update(IORAV2[8]) <- IORAV2[8] + c_IORAV2_in[8] - c_IORAV2_out[8] + c_aging_IORAV2[7] - c_aging_IORAV2[8]
update(IORAV2[9]) <- IORAV2[9] + c_IORAV2_in[9] - c_IORAV2_out[9] + c_aging_IORAV2[8] - c_aging_IORAV2[9]
update(IORAV2[10]) <- IORAV2[10] + c_IORAV2_in[10] - c_IORAV2_out[10] + c_aging_IORAV2[9] - c_aging_IORAV2[10]
update(IORAV2[11]) <- IORAV2[11] + c_IORAV2_in[11] - c_IORAV2_out[11] + c_aging_IORAV2[10] - c_aging_IORAV2[11]
update(IORAV2[12]) <- IORAV2[12] + c_IORAV2_in[12] - c_IORAV2_out[12] + c_aging_IORAV2[11] - c_aging_IORAV2[12]
update(IORAV2[13]) <- IORAV2[13] + c_IORAV2_in[13] - c_IORAV2_out[13] + c_aging_IORAV2[12] - c_aging_IORAV2[13]
update(IORAV2[14]) <- IORAV2[14] + c_IORAV2_in[14] - c_IORAV2_out[14] + c_aging_IORAV2[13] - c_aging_IORAV2[14]
update(IORAV2[15]) <- IORAV2[15] + c_IORAV2_in[15] - c_IORAV2_out[15] + c_aging_IORAV2[14] - c_aging_IORAV2[15]
update(IORAV2[16]) <- IORAV2[16] + c_IORAV2_in[16] - c_IORAV2_out[16] + c_aging_IORAV2[15]

# ----- 2.112. IORDS Compartment -----------------------------------------------

c_IORDS_in[] <- c_EORDS_to_IORDS[i] + c_IORDV1_to_IORDS[i] + c_IORDV2_to_IORDS[i] + c_migration_IORDS[i]
c_IORDS_out[] <- c_IORDS_to_IORDV1[i] + c_IORDS_to_ROS[i] + c_death_IORDS[i]

update(IORDS[1]) <- IORDS[1] + c_IORDS_in[1] - c_IORDS_out[1] - c_aging_IORDS[1]
update(IORDS[2]) <- IORDS[2] + c_IORDS_in[2] - c_IORDS_out[2] + c_aging_IORDS[1] - c_aging_IORDS[2]
update(IORDS[3]) <- IORDS[3] + c_IORDS_in[3] - c_IORDS_out[3] + c_aging_IORDS[2] - c_aging_IORDS[3]
update(IORDS[4]) <- IORDS[4] + c_IORDS_in[4] - c_IORDS_out[4] + c_aging_IORDS[3] - c_aging_IORDS[4]
update(IORDS[5]) <- IORDS[5] + c_IORDS_in[5] - c_IORDS_out[5] + c_aging_IORDS[4] - c_aging_IORDS[5]
update(IORDS[6]) <- IORDS[6] + c_IORDS_in[6] - c_IORDS_out[6] + c_aging_IORDS[5] - c_aging_IORDS[6]
update(IORDS[7]) <- IORDS[7] + c_IORDS_in[7] - c_IORDS_out[7] + c_aging_IORDS[6] - c_aging_IORDS[7]
update(IORDS[8]) <- IORDS[8] + c_IORDS_in[8] - c_IORDS_out[8] + c_aging_IORDS[7] - c_aging_IORDS[8]
update(IORDS[9]) <- IORDS[9] + c_IORDS_in[9] - c_IORDS_out[9] + c_aging_IORDS[8] - c_aging_IORDS[9]
update(IORDS[10]) <- IORDS[10] + c_IORDS_in[10] - c_IORDS_out[10] + c_aging_IORDS[9] - c_aging_IORDS[10]
update(IORDS[11]) <- IORDS[11] + c_IORDS_in[11] - c_IORDS_out[11] + c_aging_IORDS[10] - c_aging_IORDS[11]
update(IORDS[12]) <- IORDS[12] + c_IORDS_in[12] - c_IORDS_out[12] + c_aging_IORDS[11] - c_aging_IORDS[12]
update(IORDS[13]) <- IORDS[13] + c_IORDS_in[13] - c_IORDS_out[13] + c_aging_IORDS[12] - c_aging_IORDS[13]
update(IORDS[14]) <- IORDS[14] + c_IORDS_in[14] - c_IORDS_out[14] + c_aging_IORDS[13] - c_aging_IORDS[14]
update(IORDS[15]) <- IORDS[15] + c_IORDS_in[15] - c_IORDS_out[15] + c_aging_IORDS[14] - c_aging_IORDS[15]
update(IORDS[16]) <- IORDS[16] + c_IORDS_in[16] - c_IORDS_out[16] + c_aging_IORDS[15]

# ----- 2.113. IORDV1 Compartment -----------------------------------------------

c_IORDV1_in[] <- c_EORDV1_to_IORDV1[i] + c_IORDS_to_IORDV1[i] + c_migration_IORDV1[i]
c_IORDV1_out[] <- c_IORDV1_to_IORDS[i] + c_IORDV1_to_IORDV2[i] + c_IORDV1_to_ROV1[i] + c_death_IORDV1[i]

update(IORDV1[1]) <- IORDV1[1] + c_IORDV1_in[1] - c_IORDV1_out[1] - c_aging_IORDV1[1]
update(IORDV1[2]) <- IORDV1[2] + c_IORDV1_in[2] - c_IORDV1_out[2] + c_aging_IORDV1[1] - c_aging_IORDV1[2]
update(IORDV1[3]) <- IORDV1[3] + c_IORDV1_in[3] - c_IORDV1_out[3] + c_aging_IORDV1[2] - c_aging_IORDV1[3]
update(IORDV1[4]) <- IORDV1[4] + c_IORDV1_in[4] - c_IORDV1_out[4] + c_aging_IORDV1[3] - c_aging_IORDV1[4]
update(IORDV1[5]) <- IORDV1[5] + c_IORDV1_in[5] - c_IORDV1_out[5] + c_aging_IORDV1[4] - c_aging_IORDV1[5]
update(IORDV1[6]) <- IORDV1[6] + c_IORDV1_in[6] - c_IORDV1_out[6] + c_aging_IORDV1[5] - c_aging_IORDV1[6]
update(IORDV1[7]) <- IORDV1[7] + c_IORDV1_in[7] - c_IORDV1_out[7] + c_aging_IORDV1[6] - c_aging_IORDV1[7]
update(IORDV1[8]) <- IORDV1[8] + c_IORDV1_in[8] - c_IORDV1_out[8] + c_aging_IORDV1[7] - c_aging_IORDV1[8]
update(IORDV1[9]) <- IORDV1[9] + c_IORDV1_in[9] - c_IORDV1_out[9] + c_aging_IORDV1[8] - c_aging_IORDV1[9]
update(IORDV1[10]) <- IORDV1[10] + c_IORDV1_in[10] - c_IORDV1_out[10] + c_aging_IORDV1[9] - c_aging_IORDV1[10]
update(IORDV1[11]) <- IORDV1[11] + c_IORDV1_in[11] - c_IORDV1_out[11] + c_aging_IORDV1[10] - c_aging_IORDV1[11]
update(IORDV1[12]) <- IORDV1[12] + c_IORDV1_in[12] - c_IORDV1_out[12] + c_aging_IORDV1[11] - c_aging_IORDV1[12]
update(IORDV1[13]) <- IORDV1[13] + c_IORDV1_in[13] - c_IORDV1_out[13] + c_aging_IORDV1[12] - c_aging_IORDV1[13]
update(IORDV1[14]) <- IORDV1[14] + c_IORDV1_in[14] - c_IORDV1_out[14] + c_aging_IORDV1[13] - c_aging_IORDV1[14]
update(IORDV1[15]) <- IORDV1[15] + c_IORDV1_in[15] - c_IORDV1_out[15] + c_aging_IORDV1[14] - c_aging_IORDV1[15]
update(IORDV1[16]) <- IORDV1[16] + c_IORDV1_in[16] - c_IORDV1_out[16] + c_aging_IORDV1[15]

# ----- 2.114. IORDV2 Compartment -----------------------------------------------

c_IORDV2_in[] <- c_EORDV2_to_IORDV2[i] + c_IORDV1_to_IORDV2[i] + c_migration_IORDV2[i]
c_IORDV2_out[] <- c_IORDV2_to_IORDS[i] + c_IORDV2_to_ROV2[i] + c_death_IORDV2[i]

update(IORDV2[1]) <- IORDV2[1] + c_IORDV2_in[1] - c_IORDV2_out[1] - c_aging_IORDV2[1]
update(IORDV2[2]) <- IORDV2[2] + c_IORDV2_in[2] - c_IORDV2_out[2] + c_aging_IORDV2[1] - c_aging_IORDV2[2]
update(IORDV2[3]) <- IORDV2[3] + c_IORDV2_in[3] - c_IORDV2_out[3] + c_aging_IORDV2[2] - c_aging_IORDV2[3]
update(IORDV2[4]) <- IORDV2[4] + c_IORDV2_in[4] - c_IORDV2_out[4] + c_aging_IORDV2[3] - c_aging_IORDV2[4]
update(IORDV2[5]) <- IORDV2[5] + c_IORDV2_in[5] - c_IORDV2_out[5] + c_aging_IORDV2[4] - c_aging_IORDV2[5]
update(IORDV2[6]) <- IORDV2[6] + c_IORDV2_in[6] - c_IORDV2_out[6] + c_aging_IORDV2[5] - c_aging_IORDV2[6]
update(IORDV2[7]) <- IORDV2[7] + c_IORDV2_in[7] - c_IORDV2_out[7] + c_aging_IORDV2[6] - c_aging_IORDV2[7]
update(IORDV2[8]) <- IORDV2[8] + c_IORDV2_in[8] - c_IORDV2_out[8] + c_aging_IORDV2[7] - c_aging_IORDV2[8]
update(IORDV2[9]) <- IORDV2[9] + c_IORDV2_in[9] - c_IORDV2_out[9] + c_aging_IORDV2[8] - c_aging_IORDV2[9]
update(IORDV2[10]) <- IORDV2[10] + c_IORDV2_in[10] - c_IORDV2_out[10] + c_aging_IORDV2[9] - c_aging_IORDV2[10]
update(IORDV2[11]) <- IORDV2[11] + c_IORDV2_in[11] - c_IORDV2_out[11] + c_aging_IORDV2[10] - c_aging_IORDV2[11]
update(IORDV2[12]) <- IORDV2[12] + c_IORDV2_in[12] - c_IORDV2_out[12] + c_aging_IORDV2[11] - c_aging_IORDV2[12]
update(IORDV2[13]) <- IORDV2[13] + c_IORDV2_in[13] - c_IORDV2_out[13] + c_aging_IORDV2[12] - c_aging_IORDV2[13]
update(IORDV2[14]) <- IORDV2[14] + c_IORDV2_in[14] - c_IORDV2_out[14] + c_aging_IORDV2[13] - c_aging_IORDV2[14]
update(IORDV2[15]) <- IORDV2[15] + c_IORDV2_in[15] - c_IORDV2_out[15] + c_aging_IORDV2[14] - c_aging_IORDV2[15]
update(IORDV2[16]) <- IORDV2[16] + c_IORDV2_in[16] - c_IORDV2_out[16] + c_aging_IORDV2[15]

# ----- 2.116. EXROS Compartment -----------------------------------------------

c_EXROS_in[] <- c_EXROV1_to_EXROS[i] + c_EXROV2_to_EXROS[i] + c_ROS_to_EXROS[i] + c_migration_EXROS[i]
c_EXROS_out[] <- c_EXROS_to_EXROV1[i] + c_EXROS_to_IXROS[i] + c_death_EXROS[i]

update(EXROS[1]) <- EXROS[1] + c_EXROS_in[1] - c_EXROS_out[1] - c_aging_EXROS[1]
update(EXROS[2]) <- EXROS[2] + c_EXROS_in[2] - c_EXROS_out[2] + c_aging_EXROS[1] - c_aging_EXROS[2]
update(EXROS[3]) <- EXROS[3] + c_EXROS_in[3] - c_EXROS_out[3] + c_aging_EXROS[2] - c_aging_EXROS[3]
update(EXROS[4]) <- EXROS[4] + c_EXROS_in[4] - c_EXROS_out[4] + c_aging_EXROS[3] - c_aging_EXROS[4]
update(EXROS[5]) <- EXROS[5] + c_EXROS_in[5] - c_EXROS_out[5] + c_aging_EXROS[4] - c_aging_EXROS[5]
update(EXROS[6]) <- EXROS[6] + c_EXROS_in[6] - c_EXROS_out[6] + c_aging_EXROS[5] - c_aging_EXROS[6]
update(EXROS[7]) <- EXROS[7] + c_EXROS_in[7] - c_EXROS_out[7] + c_aging_EXROS[6] - c_aging_EXROS[7]
update(EXROS[8]) <- EXROS[8] + c_EXROS_in[8] - c_EXROS_out[8] + c_aging_EXROS[7] - c_aging_EXROS[8]
update(EXROS[9]) <- EXROS[9] + c_EXROS_in[9] - c_EXROS_out[9] + c_aging_EXROS[8] - c_aging_EXROS[9]
update(EXROS[10]) <- EXROS[10] + c_EXROS_in[10] - c_EXROS_out[10] + c_aging_EXROS[9] - c_aging_EXROS[10]
update(EXROS[11]) <- EXROS[11] + c_EXROS_in[11] - c_EXROS_out[11] + c_aging_EXROS[10] - c_aging_EXROS[11]
update(EXROS[12]) <- EXROS[12] + c_EXROS_in[12] - c_EXROS_out[12] + c_aging_EXROS[11] - c_aging_EXROS[12]
update(EXROS[13]) <- EXROS[13] + c_EXROS_in[13] - c_EXROS_out[13] + c_aging_EXROS[12] - c_aging_EXROS[13]
update(EXROS[14]) <- EXROS[14] + c_EXROS_in[14] - c_EXROS_out[14] + c_aging_EXROS[13] - c_aging_EXROS[14]
update(EXROS[15]) <- EXROS[15] + c_EXROS_in[15] - c_EXROS_out[15] + c_aging_EXROS[14] - c_aging_EXROS[15]
update(EXROS[16]) <- EXROS[16] + c_EXROS_in[16] - c_EXROS_out[16] + c_aging_EXROS[15]

# ----- 2.117. EXROV1 Compartment -----------------------------------------------

c_EXROV1_in[] <- c_EXROS_to_EXROV1[i] + c_ROV1_to_EXROV1[i] + c_migration_EXROV1[i]
c_EXROV1_out[] <- c_EXROV1_to_EXROS[i] + c_EXROV1_to_EXROV2[i] + c_EXROV1_to_IXROV1[i] + c_death_EXROV1[i]

update(EXROV1[1]) <- EXROV1[1] + c_EXROV1_in[1] - c_EXROV1_out[1] - c_aging_EXROV1[1]
update(EXROV1[2]) <- EXROV1[2] + c_EXROV1_in[2] - c_EXROV1_out[2] + c_aging_EXROV1[1] - c_aging_EXROV1[2]
update(EXROV1[3]) <- EXROV1[3] + c_EXROV1_in[3] - c_EXROV1_out[3] + c_aging_EXROV1[2] - c_aging_EXROV1[3]
update(EXROV1[4]) <- EXROV1[4] + c_EXROV1_in[4] - c_EXROV1_out[4] + c_aging_EXROV1[3] - c_aging_EXROV1[4]
update(EXROV1[5]) <- EXROV1[5] + c_EXROV1_in[5] - c_EXROV1_out[5] + c_aging_EXROV1[4] - c_aging_EXROV1[5]
update(EXROV1[6]) <- EXROV1[6] + c_EXROV1_in[6] - c_EXROV1_out[6] + c_aging_EXROV1[5] - c_aging_EXROV1[6]
update(EXROV1[7]) <- EXROV1[7] + c_EXROV1_in[7] - c_EXROV1_out[7] + c_aging_EXROV1[6] - c_aging_EXROV1[7]
update(EXROV1[8]) <- EXROV1[8] + c_EXROV1_in[8] - c_EXROV1_out[8] + c_aging_EXROV1[7] - c_aging_EXROV1[8]
update(EXROV1[9]) <- EXROV1[9] + c_EXROV1_in[9] - c_EXROV1_out[9] + c_aging_EXROV1[8] - c_aging_EXROV1[9]
update(EXROV1[10]) <- EXROV1[10] + c_EXROV1_in[10] - c_EXROV1_out[10] + c_aging_EXROV1[9] - c_aging_EXROV1[10]
update(EXROV1[11]) <- EXROV1[11] + c_EXROV1_in[11] - c_EXROV1_out[11] + c_aging_EXROV1[10] - c_aging_EXROV1[11]
update(EXROV1[12]) <- EXROV1[12] + c_EXROV1_in[12] - c_EXROV1_out[12] + c_aging_EXROV1[11] - c_aging_EXROV1[12]
update(EXROV1[13]) <- EXROV1[13] + c_EXROV1_in[13] - c_EXROV1_out[13] + c_aging_EXROV1[12] - c_aging_EXROV1[13]
update(EXROV1[14]) <- EXROV1[14] + c_EXROV1_in[14] - c_EXROV1_out[14] + c_aging_EXROV1[13] - c_aging_EXROV1[14]
update(EXROV1[15]) <- EXROV1[15] + c_EXROV1_in[15] - c_EXROV1_out[15] + c_aging_EXROV1[14] - c_aging_EXROV1[15]
update(EXROV1[16]) <- EXROV1[16] + c_EXROV1_in[16] - c_EXROV1_out[16] + c_aging_EXROV1[15]

# ----- 2.118. EXROV2 Compartment -----------------------------------------------

c_EXROV2_in[] <- c_EXROV1_to_EXROV2[i] + c_ROV2_to_EXROV2[i] + c_migration_EXROV2[i]
c_EXROV2_out[] <- c_EXROV2_to_EXROS[i] + c_EXROV2_to_IXROV2[i] + c_death_EXROV2[i]

update(EXROV2[1]) <- EXROV2[1] + c_EXROV2_in[1] - c_EXROV2_out[1] - c_aging_EXROV2[1]
update(EXROV2[2]) <- EXROV2[2] + c_EXROV2_in[2] - c_EXROV2_out[2] + c_aging_EXROV2[1] - c_aging_EXROV2[2]
update(EXROV2[3]) <- EXROV2[3] + c_EXROV2_in[3] - c_EXROV2_out[3] + c_aging_EXROV2[2] - c_aging_EXROV2[3]
update(EXROV2[4]) <- EXROV2[4] + c_EXROV2_in[4] - c_EXROV2_out[4] + c_aging_EXROV2[3] - c_aging_EXROV2[4]
update(EXROV2[5]) <- EXROV2[5] + c_EXROV2_in[5] - c_EXROV2_out[5] + c_aging_EXROV2[4] - c_aging_EXROV2[5]
update(EXROV2[6]) <- EXROV2[6] + c_EXROV2_in[6] - c_EXROV2_out[6] + c_aging_EXROV2[5] - c_aging_EXROV2[6]
update(EXROV2[7]) <- EXROV2[7] + c_EXROV2_in[7] - c_EXROV2_out[7] + c_aging_EXROV2[6] - c_aging_EXROV2[7]
update(EXROV2[8]) <- EXROV2[8] + c_EXROV2_in[8] - c_EXROV2_out[8] + c_aging_EXROV2[7] - c_aging_EXROV2[8]
update(EXROV2[9]) <- EXROV2[9] + c_EXROV2_in[9] - c_EXROV2_out[9] + c_aging_EXROV2[8] - c_aging_EXROV2[9]
update(EXROV2[10]) <- EXROV2[10] + c_EXROV2_in[10] - c_EXROV2_out[10] + c_aging_EXROV2[9] - c_aging_EXROV2[10]
update(EXROV2[11]) <- EXROV2[11] + c_EXROV2_in[11] - c_EXROV2_out[11] + c_aging_EXROV2[10] - c_aging_EXROV2[11]
update(EXROV2[12]) <- EXROV2[12] + c_EXROV2_in[12] - c_EXROV2_out[12] + c_aging_EXROV2[11] - c_aging_EXROV2[12]
update(EXROV2[13]) <- EXROV2[13] + c_EXROV2_in[13] - c_EXROV2_out[13] + c_aging_EXROV2[12] - c_aging_EXROV2[13]
update(EXROV2[14]) <- EXROV2[14] + c_EXROV2_in[14] - c_EXROV2_out[14] + c_aging_EXROV2[13] - c_aging_EXROV2[14]
update(EXROV2[15]) <- EXROV2[15] + c_EXROV2_in[15] - c_EXROV2_out[15] + c_aging_EXROV2[14] - c_aging_EXROV2[15]
update(EXROV2[16]) <- EXROV2[16] + c_EXROV2_in[16] - c_EXROV2_out[16] + c_aging_EXROV2[15]

# ----- 2.120. EORXS Compartment -----------------------------------------------

c_EORXS_in[] <- c_EORXV1_to_EORXS[i] + c_EORXV2_to_EORXS[i] + c_RXS_to_EORXS[i] + c_migration_EORXS[i]
c_EORXS_out[] <- c_EORXS_to_EORXV1[i] + c_EORXS_to_IORXS[i] + c_death_EORXS[i]

update(EORXS[1]) <- EORXS[1] + c_EORXS_in[1] - c_EORXS_out[1] - c_aging_EORXS[1]
update(EORXS[2]) <- EORXS[2] + c_EORXS_in[2] - c_EORXS_out[2] + c_aging_EORXS[1] - c_aging_EORXS[2]
update(EORXS[3]) <- EORXS[3] + c_EORXS_in[3] - c_EORXS_out[3] + c_aging_EORXS[2] - c_aging_EORXS[3]
update(EORXS[4]) <- EORXS[4] + c_EORXS_in[4] - c_EORXS_out[4] + c_aging_EORXS[3] - c_aging_EORXS[4]
update(EORXS[5]) <- EORXS[5] + c_EORXS_in[5] - c_EORXS_out[5] + c_aging_EORXS[4] - c_aging_EORXS[5]
update(EORXS[6]) <- EORXS[6] + c_EORXS_in[6] - c_EORXS_out[6] + c_aging_EORXS[5] - c_aging_EORXS[6]
update(EORXS[7]) <- EORXS[7] + c_EORXS_in[7] - c_EORXS_out[7] + c_aging_EORXS[6] - c_aging_EORXS[7]
update(EORXS[8]) <- EORXS[8] + c_EORXS_in[8] - c_EORXS_out[8] + c_aging_EORXS[7] - c_aging_EORXS[8]
update(EORXS[9]) <- EORXS[9] + c_EORXS_in[9] - c_EORXS_out[9] + c_aging_EORXS[8] - c_aging_EORXS[9]
update(EORXS[10]) <- EORXS[10] + c_EORXS_in[10] - c_EORXS_out[10] + c_aging_EORXS[9] - c_aging_EORXS[10]
update(EORXS[11]) <- EORXS[11] + c_EORXS_in[11] - c_EORXS_out[11] + c_aging_EORXS[10] - c_aging_EORXS[11]
update(EORXS[12]) <- EORXS[12] + c_EORXS_in[12] - c_EORXS_out[12] + c_aging_EORXS[11] - c_aging_EORXS[12]
update(EORXS[13]) <- EORXS[13] + c_EORXS_in[13] - c_EORXS_out[13] + c_aging_EORXS[12] - c_aging_EORXS[13]
update(EORXS[14]) <- EORXS[14] + c_EORXS_in[14] - c_EORXS_out[14] + c_aging_EORXS[13] - c_aging_EORXS[14]
update(EORXS[15]) <- EORXS[15] + c_EORXS_in[15] - c_EORXS_out[15] + c_aging_EORXS[14] - c_aging_EORXS[15]
update(EORXS[16]) <- EORXS[16] + c_EORXS_in[16] - c_EORXS_out[16] + c_aging_EORXS[15]

# ----- 2.121. EORXV1 Compartment -----------------------------------------------

c_EORXV1_in[] <- c_EORXS_to_EORXV1[i] + c_RXV1_to_EORXV1[i] + c_migration_EORXV1[i]
c_EORXV1_out[] <- c_EORXV1_to_EORXS[i] + c_EORXV1_to_EORXV2[i] + c_EORXV1_to_IORXV1[i] + c_death_EORXV1[i]

update(EORXV1[1]) <- EORXV1[1] + c_EORXV1_in[1] - c_EORXV1_out[1] - c_aging_EORXV1[1]
update(EORXV1[2]) <- EORXV1[2] + c_EORXV1_in[2] - c_EORXV1_out[2] + c_aging_EORXV1[1] - c_aging_EORXV1[2]
update(EORXV1[3]) <- EORXV1[3] + c_EORXV1_in[3] - c_EORXV1_out[3] + c_aging_EORXV1[2] - c_aging_EORXV1[3]
update(EORXV1[4]) <- EORXV1[4] + c_EORXV1_in[4] - c_EORXV1_out[4] + c_aging_EORXV1[3] - c_aging_EORXV1[4]
update(EORXV1[5]) <- EORXV1[5] + c_EORXV1_in[5] - c_EORXV1_out[5] + c_aging_EORXV1[4] - c_aging_EORXV1[5]
update(EORXV1[6]) <- EORXV1[6] + c_EORXV1_in[6] - c_EORXV1_out[6] + c_aging_EORXV1[5] - c_aging_EORXV1[6]
update(EORXV1[7]) <- EORXV1[7] + c_EORXV1_in[7] - c_EORXV1_out[7] + c_aging_EORXV1[6] - c_aging_EORXV1[7]
update(EORXV1[8]) <- EORXV1[8] + c_EORXV1_in[8] - c_EORXV1_out[8] + c_aging_EORXV1[7] - c_aging_EORXV1[8]
update(EORXV1[9]) <- EORXV1[9] + c_EORXV1_in[9] - c_EORXV1_out[9] + c_aging_EORXV1[8] - c_aging_EORXV1[9]
update(EORXV1[10]) <- EORXV1[10] + c_EORXV1_in[10] - c_EORXV1_out[10] + c_aging_EORXV1[9] - c_aging_EORXV1[10]
update(EORXV1[11]) <- EORXV1[11] + c_EORXV1_in[11] - c_EORXV1_out[11] + c_aging_EORXV1[10] - c_aging_EORXV1[11]
update(EORXV1[12]) <- EORXV1[12] + c_EORXV1_in[12] - c_EORXV1_out[12] + c_aging_EORXV1[11] - c_aging_EORXV1[12]
update(EORXV1[13]) <- EORXV1[13] + c_EORXV1_in[13] - c_EORXV1_out[13] + c_aging_EORXV1[12] - c_aging_EORXV1[13]
update(EORXV1[14]) <- EORXV1[14] + c_EORXV1_in[14] - c_EORXV1_out[14] + c_aging_EORXV1[13] - c_aging_EORXV1[14]
update(EORXV1[15]) <- EORXV1[15] + c_EORXV1_in[15] - c_EORXV1_out[15] + c_aging_EORXV1[14] - c_aging_EORXV1[15]
update(EORXV1[16]) <- EORXV1[16] + c_EORXV1_in[16] - c_EORXV1_out[16] + c_aging_EORXV1[15]

# ----- 2.122. EORXV2 Compartment -----------------------------------------------

c_EORXV2_in[] <- c_EORXV1_to_EORXV2[i] + c_RXV2_to_EORXV2[i] + c_migration_EORXV2[i]
c_EORXV2_out[] <- c_EORXV2_to_EORXS[i] + c_EORXV2_to_IORXV2[i] + c_death_EORXV2[i]

update(EORXV2[1]) <- EORXV2[1] + c_EORXV2_in[1] - c_EORXV2_out[1] - c_aging_EORXV2[1]
update(EORXV2[2]) <- EORXV2[2] + c_EORXV2_in[2] - c_EORXV2_out[2] + c_aging_EORXV2[1] - c_aging_EORXV2[2]
update(EORXV2[3]) <- EORXV2[3] + c_EORXV2_in[3] - c_EORXV2_out[3] + c_aging_EORXV2[2] - c_aging_EORXV2[3]
update(EORXV2[4]) <- EORXV2[4] + c_EORXV2_in[4] - c_EORXV2_out[4] + c_aging_EORXV2[3] - c_aging_EORXV2[4]
update(EORXV2[5]) <- EORXV2[5] + c_EORXV2_in[5] - c_EORXV2_out[5] + c_aging_EORXV2[4] - c_aging_EORXV2[5]
update(EORXV2[6]) <- EORXV2[6] + c_EORXV2_in[6] - c_EORXV2_out[6] + c_aging_EORXV2[5] - c_aging_EORXV2[6]
update(EORXV2[7]) <- EORXV2[7] + c_EORXV2_in[7] - c_EORXV2_out[7] + c_aging_EORXV2[6] - c_aging_EORXV2[7]
update(EORXV2[8]) <- EORXV2[8] + c_EORXV2_in[8] - c_EORXV2_out[8] + c_aging_EORXV2[7] - c_aging_EORXV2[8]
update(EORXV2[9]) <- EORXV2[9] + c_EORXV2_in[9] - c_EORXV2_out[9] + c_aging_EORXV2[8] - c_aging_EORXV2[9]
update(EORXV2[10]) <- EORXV2[10] + c_EORXV2_in[10] - c_EORXV2_out[10] + c_aging_EORXV2[9] - c_aging_EORXV2[10]
update(EORXV2[11]) <- EORXV2[11] + c_EORXV2_in[11] - c_EORXV2_out[11] + c_aging_EORXV2[10] - c_aging_EORXV2[11]
update(EORXV2[12]) <- EORXV2[12] + c_EORXV2_in[12] - c_EORXV2_out[12] + c_aging_EORXV2[11] - c_aging_EORXV2[12]
update(EORXV2[13]) <- EORXV2[13] + c_EORXV2_in[13] - c_EORXV2_out[13] + c_aging_EORXV2[12] - c_aging_EORXV2[13]
update(EORXV2[14]) <- EORXV2[14] + c_EORXV2_in[14] - c_EORXV2_out[14] + c_aging_EORXV2[13] - c_aging_EORXV2[14]
update(EORXV2[15]) <- EORXV2[15] + c_EORXV2_in[15] - c_EORXV2_out[15] + c_aging_EORXV2[14] - c_aging_EORXV2[15]
update(EORXV2[16]) <- EORXV2[16] + c_EORXV2_in[16] - c_EORXV2_out[16] + c_aging_EORXV2[15]

# ----- 2.125. IXROS Compartment -----------------------------------------------

c_IXROS_in[] <- c_EXROS_to_IXROS[i] + c_IXROV1_to_IXROS[i] + c_IXROV2_to_IXROS[i] + c_migration_IXROS[i]
c_IXROS_out[] <- c_IXROS_to_IXROV1[i] + c_IXROS_to_ROXS[i] + c_death_IXROS[i]

update(IXROS[1]) <- IXROS[1] + c_IXROS_in[1] - c_IXROS_out[1] - c_aging_IXROS[1]
update(IXROS[2]) <- IXROS[2] + c_IXROS_in[2] - c_IXROS_out[2] + c_aging_IXROS[1] - c_aging_IXROS[2]
update(IXROS[3]) <- IXROS[3] + c_IXROS_in[3] - c_IXROS_out[3] + c_aging_IXROS[2] - c_aging_IXROS[3]
update(IXROS[4]) <- IXROS[4] + c_IXROS_in[4] - c_IXROS_out[4] + c_aging_IXROS[3] - c_aging_IXROS[4]
update(IXROS[5]) <- IXROS[5] + c_IXROS_in[5] - c_IXROS_out[5] + c_aging_IXROS[4] - c_aging_IXROS[5]
update(IXROS[6]) <- IXROS[6] + c_IXROS_in[6] - c_IXROS_out[6] + c_aging_IXROS[5] - c_aging_IXROS[6]
update(IXROS[7]) <- IXROS[7] + c_IXROS_in[7] - c_IXROS_out[7] + c_aging_IXROS[6] - c_aging_IXROS[7]
update(IXROS[8]) <- IXROS[8] + c_IXROS_in[8] - c_IXROS_out[8] + c_aging_IXROS[7] - c_aging_IXROS[8]
update(IXROS[9]) <- IXROS[9] + c_IXROS_in[9] - c_IXROS_out[9] + c_aging_IXROS[8] - c_aging_IXROS[9]
update(IXROS[10]) <- IXROS[10] + c_IXROS_in[10] - c_IXROS_out[10] + c_aging_IXROS[9] - c_aging_IXROS[10]
update(IXROS[11]) <- IXROS[11] + c_IXROS_in[11] - c_IXROS_out[11] + c_aging_IXROS[10] - c_aging_IXROS[11]
update(IXROS[12]) <- IXROS[12] + c_IXROS_in[12] - c_IXROS_out[12] + c_aging_IXROS[11] - c_aging_IXROS[12]
update(IXROS[13]) <- IXROS[13] + c_IXROS_in[13] - c_IXROS_out[13] + c_aging_IXROS[12] - c_aging_IXROS[13]
update(IXROS[14]) <- IXROS[14] + c_IXROS_in[14] - c_IXROS_out[14] + c_aging_IXROS[13] - c_aging_IXROS[14]
update(IXROS[15]) <- IXROS[15] + c_IXROS_in[15] - c_IXROS_out[15] + c_aging_IXROS[14] - c_aging_IXROS[15]
update(IXROS[16]) <- IXROS[16] + c_IXROS_in[16] - c_IXROS_out[16] + c_aging_IXROS[15]

# ----- 2.126. IXROV1 Compartment -----------------------------------------------

c_IXROV1_in[] <- c_EXROV1_to_IXROV1[i] + c_IXROS_to_IXROV1[i] + c_migration_IXROV1[i]
c_IXROV1_out[] <- c_IXROV1_to_IXROS[i] + c_IXROV1_to_IXROV2[i] + c_IXROV1_to_ROXV1[i] + c_death_IXROV1[i]

update(IXROV1[1]) <- IXROV1[1] + c_IXROV1_in[1] - c_IXROV1_out[1] - c_aging_IXROV1[1]
update(IXROV1[2]) <- IXROV1[2] + c_IXROV1_in[2] - c_IXROV1_out[2] + c_aging_IXROV1[1] - c_aging_IXROV1[2]
update(IXROV1[3]) <- IXROV1[3] + c_IXROV1_in[3] - c_IXROV1_out[3] + c_aging_IXROV1[2] - c_aging_IXROV1[3]
update(IXROV1[4]) <- IXROV1[4] + c_IXROV1_in[4] - c_IXROV1_out[4] + c_aging_IXROV1[3] - c_aging_IXROV1[4]
update(IXROV1[5]) <- IXROV1[5] + c_IXROV1_in[5] - c_IXROV1_out[5] + c_aging_IXROV1[4] - c_aging_IXROV1[5]
update(IXROV1[6]) <- IXROV1[6] + c_IXROV1_in[6] - c_IXROV1_out[6] + c_aging_IXROV1[5] - c_aging_IXROV1[6]
update(IXROV1[7]) <- IXROV1[7] + c_IXROV1_in[7] - c_IXROV1_out[7] + c_aging_IXROV1[6] - c_aging_IXROV1[7]
update(IXROV1[8]) <- IXROV1[8] + c_IXROV1_in[8] - c_IXROV1_out[8] + c_aging_IXROV1[7] - c_aging_IXROV1[8]
update(IXROV1[9]) <- IXROV1[9] + c_IXROV1_in[9] - c_IXROV1_out[9] + c_aging_IXROV1[8] - c_aging_IXROV1[9]
update(IXROV1[10]) <- IXROV1[10] + c_IXROV1_in[10] - c_IXROV1_out[10] + c_aging_IXROV1[9] - c_aging_IXROV1[10]
update(IXROV1[11]) <- IXROV1[11] + c_IXROV1_in[11] - c_IXROV1_out[11] + c_aging_IXROV1[10] - c_aging_IXROV1[11]
update(IXROV1[12]) <- IXROV1[12] + c_IXROV1_in[12] - c_IXROV1_out[12] + c_aging_IXROV1[11] - c_aging_IXROV1[12]
update(IXROV1[13]) <- IXROV1[13] + c_IXROV1_in[13] - c_IXROV1_out[13] + c_aging_IXROV1[12] - c_aging_IXROV1[13]
update(IXROV1[14]) <- IXROV1[14] + c_IXROV1_in[14] - c_IXROV1_out[14] + c_aging_IXROV1[13] - c_aging_IXROV1[14]
update(IXROV1[15]) <- IXROV1[15] + c_IXROV1_in[15] - c_IXROV1_out[15] + c_aging_IXROV1[14] - c_aging_IXROV1[15]
update(IXROV1[16]) <- IXROV1[16] + c_IXROV1_in[16] - c_IXROV1_out[16] + c_aging_IXROV1[15]

# ----- 2.127. IXROV2 Compartment -----------------------------------------------

c_IXROV2_in[] <- c_EXROV2_to_IXROV2[i] + c_IXROV1_to_IXROV2[i] + c_migration_IXROV2[i]
c_IXROV2_out[] <- c_IXROV2_to_IXROS[i] + c_IXROV2_to_ROXV2[i] + c_death_IXROV2[i]

update(IXROV2[1]) <- IXROV2[1] + c_IXROV2_in[1] - c_IXROV2_out[1] - c_aging_IXROV2[1]
update(IXROV2[2]) <- IXROV2[2] + c_IXROV2_in[2] - c_IXROV2_out[2] + c_aging_IXROV2[1] - c_aging_IXROV2[2]
update(IXROV2[3]) <- IXROV2[3] + c_IXROV2_in[3] - c_IXROV2_out[3] + c_aging_IXROV2[2] - c_aging_IXROV2[3]
update(IXROV2[4]) <- IXROV2[4] + c_IXROV2_in[4] - c_IXROV2_out[4] + c_aging_IXROV2[3] - c_aging_IXROV2[4]
update(IXROV2[5]) <- IXROV2[5] + c_IXROV2_in[5] - c_IXROV2_out[5] + c_aging_IXROV2[4] - c_aging_IXROV2[5]
update(IXROV2[6]) <- IXROV2[6] + c_IXROV2_in[6] - c_IXROV2_out[6] + c_aging_IXROV2[5] - c_aging_IXROV2[6]
update(IXROV2[7]) <- IXROV2[7] + c_IXROV2_in[7] - c_IXROV2_out[7] + c_aging_IXROV2[6] - c_aging_IXROV2[7]
update(IXROV2[8]) <- IXROV2[8] + c_IXROV2_in[8] - c_IXROV2_out[8] + c_aging_IXROV2[7] - c_aging_IXROV2[8]
update(IXROV2[9]) <- IXROV2[9] + c_IXROV2_in[9] - c_IXROV2_out[9] + c_aging_IXROV2[8] - c_aging_IXROV2[9]
update(IXROV2[10]) <- IXROV2[10] + c_IXROV2_in[10] - c_IXROV2_out[10] + c_aging_IXROV2[9] - c_aging_IXROV2[10]
update(IXROV2[11]) <- IXROV2[11] + c_IXROV2_in[11] - c_IXROV2_out[11] + c_aging_IXROV2[10] - c_aging_IXROV2[11]
update(IXROV2[12]) <- IXROV2[12] + c_IXROV2_in[12] - c_IXROV2_out[12] + c_aging_IXROV2[11] - c_aging_IXROV2[12]
update(IXROV2[13]) <- IXROV2[13] + c_IXROV2_in[13] - c_IXROV2_out[13] + c_aging_IXROV2[12] - c_aging_IXROV2[13]
update(IXROV2[14]) <- IXROV2[14] + c_IXROV2_in[14] - c_IXROV2_out[14] + c_aging_IXROV2[13] - c_aging_IXROV2[14]
update(IXROV2[15]) <- IXROV2[15] + c_IXROV2_in[15] - c_IXROV2_out[15] + c_aging_IXROV2[14] - c_aging_IXROV2[15]
update(IXROV2[16]) <- IXROV2[16] + c_IXROV2_in[16] - c_IXROV2_out[16] + c_aging_IXROV2[15]

# ----- 2.129. IORXS Compartment -----------------------------------------------

c_IORXS_in[] <- c_EORXS_to_IORXS[i] + c_IORXV1_to_IORXS[i] + c_IORXV2_to_IORXS[i] + c_migration_IORXS[i]
c_IORXS_out[] <- c_IORXS_to_IORXV1[i] + c_IORXS_to_ROXS[i] + c_death_IORXS[i]

update(IORXS[1]) <- IORXS[1] + c_IORXS_in[1] - c_IORXS_out[1] - c_aging_IORXS[1]
update(IORXS[2]) <- IORXS[2] + c_IORXS_in[2] - c_IORXS_out[2] + c_aging_IORXS[1] - c_aging_IORXS[2]
update(IORXS[3]) <- IORXS[3] + c_IORXS_in[3] - c_IORXS_out[3] + c_aging_IORXS[2] - c_aging_IORXS[3]
update(IORXS[4]) <- IORXS[4] + c_IORXS_in[4] - c_IORXS_out[4] + c_aging_IORXS[3] - c_aging_IORXS[4]
update(IORXS[5]) <- IORXS[5] + c_IORXS_in[5] - c_IORXS_out[5] + c_aging_IORXS[4] - c_aging_IORXS[5]
update(IORXS[6]) <- IORXS[6] + c_IORXS_in[6] - c_IORXS_out[6] + c_aging_IORXS[5] - c_aging_IORXS[6]
update(IORXS[7]) <- IORXS[7] + c_IORXS_in[7] - c_IORXS_out[7] + c_aging_IORXS[6] - c_aging_IORXS[7]
update(IORXS[8]) <- IORXS[8] + c_IORXS_in[8] - c_IORXS_out[8] + c_aging_IORXS[7] - c_aging_IORXS[8]
update(IORXS[9]) <- IORXS[9] + c_IORXS_in[9] - c_IORXS_out[9] + c_aging_IORXS[8] - c_aging_IORXS[9]
update(IORXS[10]) <- IORXS[10] + c_IORXS_in[10] - c_IORXS_out[10] + c_aging_IORXS[9] - c_aging_IORXS[10]
update(IORXS[11]) <- IORXS[11] + c_IORXS_in[11] - c_IORXS_out[11] + c_aging_IORXS[10] - c_aging_IORXS[11]
update(IORXS[12]) <- IORXS[12] + c_IORXS_in[12] - c_IORXS_out[12] + c_aging_IORXS[11] - c_aging_IORXS[12]
update(IORXS[13]) <- IORXS[13] + c_IORXS_in[13] - c_IORXS_out[13] + c_aging_IORXS[12] - c_aging_IORXS[13]
update(IORXS[14]) <- IORXS[14] + c_IORXS_in[14] - c_IORXS_out[14] + c_aging_IORXS[13] - c_aging_IORXS[14]
update(IORXS[15]) <- IORXS[15] + c_IORXS_in[15] - c_IORXS_out[15] + c_aging_IORXS[14] - c_aging_IORXS[15]
update(IORXS[16]) <- IORXS[16] + c_IORXS_in[16] - c_IORXS_out[16] + c_aging_IORXS[15]

# ----- 2.130. IORXV1 Compartment -----------------------------------------------

c_IORXV1_in[] <- c_EORXV1_to_IORXV1[i] + c_IORXS_to_IORXV1[i] + c_migration_IORXV1[i]
c_IORXV1_out[] <- c_IORXV1_to_IORXS[i] + c_IORXV1_to_IORXV2[i] + c_IORXV1_to_ROXV1[i] + c_death_IORXV1[i]

update(IORXV1[1]) <- IORXV1[1] + c_IORXV1_in[1] - c_IORXV1_out[1] - c_aging_IORXV1[1]
update(IORXV1[2]) <- IORXV1[2] + c_IORXV1_in[2] - c_IORXV1_out[2] + c_aging_IORXV1[1] - c_aging_IORXV1[2]
update(IORXV1[3]) <- IORXV1[3] + c_IORXV1_in[3] - c_IORXV1_out[3] + c_aging_IORXV1[2] - c_aging_IORXV1[3]
update(IORXV1[4]) <- IORXV1[4] + c_IORXV1_in[4] - c_IORXV1_out[4] + c_aging_IORXV1[3] - c_aging_IORXV1[4]
update(IORXV1[5]) <- IORXV1[5] + c_IORXV1_in[5] - c_IORXV1_out[5] + c_aging_IORXV1[4] - c_aging_IORXV1[5]
update(IORXV1[6]) <- IORXV1[6] + c_IORXV1_in[6] - c_IORXV1_out[6] + c_aging_IORXV1[5] - c_aging_IORXV1[6]
update(IORXV1[7]) <- IORXV1[7] + c_IORXV1_in[7] - c_IORXV1_out[7] + c_aging_IORXV1[6] - c_aging_IORXV1[7]
update(IORXV1[8]) <- IORXV1[8] + c_IORXV1_in[8] - c_IORXV1_out[8] + c_aging_IORXV1[7] - c_aging_IORXV1[8]
update(IORXV1[9]) <- IORXV1[9] + c_IORXV1_in[9] - c_IORXV1_out[9] + c_aging_IORXV1[8] - c_aging_IORXV1[9]
update(IORXV1[10]) <- IORXV1[10] + c_IORXV1_in[10] - c_IORXV1_out[10] + c_aging_IORXV1[9] - c_aging_IORXV1[10]
update(IORXV1[11]) <- IORXV1[11] + c_IORXV1_in[11] - c_IORXV1_out[11] + c_aging_IORXV1[10] - c_aging_IORXV1[11]
update(IORXV1[12]) <- IORXV1[12] + c_IORXV1_in[12] - c_IORXV1_out[12] + c_aging_IORXV1[11] - c_aging_IORXV1[12]
update(IORXV1[13]) <- IORXV1[13] + c_IORXV1_in[13] - c_IORXV1_out[13] + c_aging_IORXV1[12] - c_aging_IORXV1[13]
update(IORXV1[14]) <- IORXV1[14] + c_IORXV1_in[14] - c_IORXV1_out[14] + c_aging_IORXV1[13] - c_aging_IORXV1[14]
update(IORXV1[15]) <- IORXV1[15] + c_IORXV1_in[15] - c_IORXV1_out[15] + c_aging_IORXV1[14] - c_aging_IORXV1[15]
update(IORXV1[16]) <- IORXV1[16] + c_IORXV1_in[16] - c_IORXV1_out[16] + c_aging_IORXV1[15]

# ----- 2.131. IORXV2 Compartment -----------------------------------------------

c_IORXV2_in[] <- c_EORXV2_to_IORXV2[i] + c_IORXV1_to_IORXV2[i] + c_migration_IORXV2[i]
c_IORXV2_out[] <- c_IORXV2_to_IORXS[i] + c_IORXV2_to_ROXV2[i] + c_death_IORXV2[i]

update(IORXV2[1]) <- IORXV2[1] + c_IORXV2_in[1] - c_IORXV2_out[1] - c_aging_IORXV2[1]
update(IORXV2[2]) <- IORXV2[2] + c_IORXV2_in[2] - c_IORXV2_out[2] + c_aging_IORXV2[1] - c_aging_IORXV2[2]
update(IORXV2[3]) <- IORXV2[3] + c_IORXV2_in[3] - c_IORXV2_out[3] + c_aging_IORXV2[2] - c_aging_IORXV2[3]
update(IORXV2[4]) <- IORXV2[4] + c_IORXV2_in[4] - c_IORXV2_out[4] + c_aging_IORXV2[3] - c_aging_IORXV2[4]
update(IORXV2[5]) <- IORXV2[5] + c_IORXV2_in[5] - c_IORXV2_out[5] + c_aging_IORXV2[4] - c_aging_IORXV2[5]
update(IORXV2[6]) <- IORXV2[6] + c_IORXV2_in[6] - c_IORXV2_out[6] + c_aging_IORXV2[5] - c_aging_IORXV2[6]
update(IORXV2[7]) <- IORXV2[7] + c_IORXV2_in[7] - c_IORXV2_out[7] + c_aging_IORXV2[6] - c_aging_IORXV2[7]
update(IORXV2[8]) <- IORXV2[8] + c_IORXV2_in[8] - c_IORXV2_out[8] + c_aging_IORXV2[7] - c_aging_IORXV2[8]
update(IORXV2[9]) <- IORXV2[9] + c_IORXV2_in[9] - c_IORXV2_out[9] + c_aging_IORXV2[8] - c_aging_IORXV2[9]
update(IORXV2[10]) <- IORXV2[10] + c_IORXV2_in[10] - c_IORXV2_out[10] + c_aging_IORXV2[9] - c_aging_IORXV2[10]
update(IORXV2[11]) <- IORXV2[11] + c_IORXV2_in[11] - c_IORXV2_out[11] + c_aging_IORXV2[10] - c_aging_IORXV2[11]
update(IORXV2[12]) <- IORXV2[12] + c_IORXV2_in[12] - c_IORXV2_out[12] + c_aging_IORXV2[11] - c_aging_IORXV2[12]
update(IORXV2[13]) <- IORXV2[13] + c_IORXV2_in[13] - c_IORXV2_out[13] + c_aging_IORXV2[12] - c_aging_IORXV2[13]
update(IORXV2[14]) <- IORXV2[14] + c_IORXV2_in[14] - c_IORXV2_out[14] + c_aging_IORXV2[13] - c_aging_IORXV2[14]
update(IORXV2[15]) <- IORXV2[15] + c_IORXV2_in[15] - c_IORXV2_out[15] + c_aging_IORXV2[14] - c_aging_IORXV2[15]
update(IORXV2[16]) <- IORXV2[16] + c_IORXV2_in[16] - c_IORXV2_out[16] + c_aging_IORXV2[15]

# ----- 2.133. ROXS Compartment -----------------------------------------------

c_ROXS_in[] <- c_IORXS_to_ROXS[i] + c_IXROS_to_ROXS[i] + c_ROXV1_to_ROXS[i] + c_ROXV2_to_ROXS[i] + c_migration_ROXS[i]
c_ROXS_out[] <- c_ROXS_to_ROXV1[i] + c_ROXS_to_S[i] + c_death_ROXS[i]

update(ROXS[1]) <- ROXS[1] + c_ROXS_in[1] - c_ROXS_out[1] - c_aging_ROXS[1]
update(ROXS[2]) <- ROXS[2] + c_ROXS_in[2] - c_ROXS_out[2] + c_aging_ROXS[1] - c_aging_ROXS[2]
update(ROXS[3]) <- ROXS[3] + c_ROXS_in[3] - c_ROXS_out[3] + c_aging_ROXS[2] - c_aging_ROXS[3]
update(ROXS[4]) <- ROXS[4] + c_ROXS_in[4] - c_ROXS_out[4] + c_aging_ROXS[3] - c_aging_ROXS[4]
update(ROXS[5]) <- ROXS[5] + c_ROXS_in[5] - c_ROXS_out[5] + c_aging_ROXS[4] - c_aging_ROXS[5]
update(ROXS[6]) <- ROXS[6] + c_ROXS_in[6] - c_ROXS_out[6] + c_aging_ROXS[5] - c_aging_ROXS[6]
update(ROXS[7]) <- ROXS[7] + c_ROXS_in[7] - c_ROXS_out[7] + c_aging_ROXS[6] - c_aging_ROXS[7]
update(ROXS[8]) <- ROXS[8] + c_ROXS_in[8] - c_ROXS_out[8] + c_aging_ROXS[7] - c_aging_ROXS[8]
update(ROXS[9]) <- ROXS[9] + c_ROXS_in[9] - c_ROXS_out[9] + c_aging_ROXS[8] - c_aging_ROXS[9]
update(ROXS[10]) <- ROXS[10] + c_ROXS_in[10] - c_ROXS_out[10] + c_aging_ROXS[9] - c_aging_ROXS[10]
update(ROXS[11]) <- ROXS[11] + c_ROXS_in[11] - c_ROXS_out[11] + c_aging_ROXS[10] - c_aging_ROXS[11]
update(ROXS[12]) <- ROXS[12] + c_ROXS_in[12] - c_ROXS_out[12] + c_aging_ROXS[11] - c_aging_ROXS[12]
update(ROXS[13]) <- ROXS[13] + c_ROXS_in[13] - c_ROXS_out[13] + c_aging_ROXS[12] - c_aging_ROXS[13]
update(ROXS[14]) <- ROXS[14] + c_ROXS_in[14] - c_ROXS_out[14] + c_aging_ROXS[13] - c_aging_ROXS[14]
update(ROXS[15]) <- ROXS[15] + c_ROXS_in[15] - c_ROXS_out[15] + c_aging_ROXS[14] - c_aging_ROXS[15]
update(ROXS[16]) <- ROXS[16] + c_ROXS_in[16] - c_ROXS_out[16] + c_aging_ROXS[15]

# ----- 2.134. ROXV1 Compartment -----------------------------------------------

c_ROXV1_in[] <- c_IORXV1_to_ROXV1[i] + c_IXROV1_to_ROXV1[i] + c_ROXS_to_ROXV1[i] + c_migration_ROXV1[i]
c_ROXV1_out[] <- c_ROXV1_to_ROXV2[i] + c_ROXV1_to_ROXS[i] + c_ROXV1_to_V1[i] + c_death_ROXV1[i]

update(ROXV1[1]) <- ROXV1[1] + c_ROXV1_in[1] - c_ROXV1_out[1] - c_aging_ROXV1[1]
update(ROXV1[2]) <- ROXV1[2] + c_ROXV1_in[2] - c_ROXV1_out[2] + c_aging_ROXV1[1] - c_aging_ROXV1[2]
update(ROXV1[3]) <- ROXV1[3] + c_ROXV1_in[3] - c_ROXV1_out[3] + c_aging_ROXV1[2] - c_aging_ROXV1[3]
update(ROXV1[4]) <- ROXV1[4] + c_ROXV1_in[4] - c_ROXV1_out[4] + c_aging_ROXV1[3] - c_aging_ROXV1[4]
update(ROXV1[5]) <- ROXV1[5] + c_ROXV1_in[5] - c_ROXV1_out[5] + c_aging_ROXV1[4] - c_aging_ROXV1[5]
update(ROXV1[6]) <- ROXV1[6] + c_ROXV1_in[6] - c_ROXV1_out[6] + c_aging_ROXV1[5] - c_aging_ROXV1[6]
update(ROXV1[7]) <- ROXV1[7] + c_ROXV1_in[7] - c_ROXV1_out[7] + c_aging_ROXV1[6] - c_aging_ROXV1[7]
update(ROXV1[8]) <- ROXV1[8] + c_ROXV1_in[8] - c_ROXV1_out[8] + c_aging_ROXV1[7] - c_aging_ROXV1[8]
update(ROXV1[9]) <- ROXV1[9] + c_ROXV1_in[9] - c_ROXV1_out[9] + c_aging_ROXV1[8] - c_aging_ROXV1[9]
update(ROXV1[10]) <- ROXV1[10] + c_ROXV1_in[10] - c_ROXV1_out[10] + c_aging_ROXV1[9] - c_aging_ROXV1[10]
update(ROXV1[11]) <- ROXV1[11] + c_ROXV1_in[11] - c_ROXV1_out[11] + c_aging_ROXV1[10] - c_aging_ROXV1[11]
update(ROXV1[12]) <- ROXV1[12] + c_ROXV1_in[12] - c_ROXV1_out[12] + c_aging_ROXV1[11] - c_aging_ROXV1[12]
update(ROXV1[13]) <- ROXV1[13] + c_ROXV1_in[13] - c_ROXV1_out[13] + c_aging_ROXV1[12] - c_aging_ROXV1[13]
update(ROXV1[14]) <- ROXV1[14] + c_ROXV1_in[14] - c_ROXV1_out[14] + c_aging_ROXV1[13] - c_aging_ROXV1[14]
update(ROXV1[15]) <- ROXV1[15] + c_ROXV1_in[15] - c_ROXV1_out[15] + c_aging_ROXV1[14] - c_aging_ROXV1[15]
update(ROXV1[16]) <- ROXV1[16] + c_ROXV1_in[16] - c_ROXV1_out[16] + c_aging_ROXV1[15]

# ----- 2.135. ROXV2 Compartment -----------------------------------------------

c_ROXV2_in[] <- c_IORXV2_to_ROXV2[i] + c_IXROV2_to_ROXV2[i] + c_ROXV1_to_ROXV2[i] + c_migration_ROXV2[i]
c_ROXV2_out[] <- c_ROXV2_to_ROXS[i] + c_ROXV2_to_V2[i] + c_death_ROXV2[i]

update(ROXV2[1]) <- ROXV2[1] + c_ROXV2_in[1] - c_ROXV2_out[1] - c_aging_ROXV2[1]
update(ROXV2[2]) <- ROXV2[2] + c_ROXV2_in[2] - c_ROXV2_out[2] + c_aging_ROXV2[1] - c_aging_ROXV2[2]
update(ROXV2[3]) <- ROXV2[3] + c_ROXV2_in[3] - c_ROXV2_out[3] + c_aging_ROXV2[2] - c_aging_ROXV2[3]
update(ROXV2[4]) <- ROXV2[4] + c_ROXV2_in[4] - c_ROXV2_out[4] + c_aging_ROXV2[3] - c_aging_ROXV2[4]
update(ROXV2[5]) <- ROXV2[5] + c_ROXV2_in[5] - c_ROXV2_out[5] + c_aging_ROXV2[4] - c_aging_ROXV2[5]
update(ROXV2[6]) <- ROXV2[6] + c_ROXV2_in[6] - c_ROXV2_out[6] + c_aging_ROXV2[5] - c_aging_ROXV2[6]
update(ROXV2[7]) <- ROXV2[7] + c_ROXV2_in[7] - c_ROXV2_out[7] + c_aging_ROXV2[6] - c_aging_ROXV2[7]
update(ROXV2[8]) <- ROXV2[8] + c_ROXV2_in[8] - c_ROXV2_out[8] + c_aging_ROXV2[7] - c_aging_ROXV2[8]
update(ROXV2[9]) <- ROXV2[9] + c_ROXV2_in[9] - c_ROXV2_out[9] + c_aging_ROXV2[8] - c_aging_ROXV2[9]
update(ROXV2[10]) <- ROXV2[10] + c_ROXV2_in[10] - c_ROXV2_out[10] + c_aging_ROXV2[9] - c_aging_ROXV2[10]
update(ROXV2[11]) <- ROXV2[11] + c_ROXV2_in[11] - c_ROXV2_out[11] + c_aging_ROXV2[10] - c_aging_ROXV2[11]
update(ROXV2[12]) <- ROXV2[12] + c_ROXV2_in[12] - c_ROXV2_out[12] + c_aging_ROXV2[11] - c_aging_ROXV2[12]
update(ROXV2[13]) <- ROXV2[13] + c_ROXV2_in[13] - c_ROXV2_out[13] + c_aging_ROXV2[12] - c_aging_ROXV2[13]
update(ROXV2[14]) <- ROXV2[14] + c_ROXV2_in[14] - c_ROXV2_out[14] + c_aging_ROXV2[13] - c_aging_ROXV2[14]
update(ROXV2[15]) <- ROXV2[15] + c_ROXV2_in[15] - c_ROXV2_out[15] + c_aging_ROXV2[14] - c_aging_ROXV2[15]
update(ROXV2[16]) <- ROXV2[16] + c_ROXV2_in[16] - c_ROXV2_out[16] + c_aging_ROXV2[15]

# ----- 2.135. Trackers --------------------------------------------------------

update(tracker_IX) <- sum(IX_tot)
update(tracker_I2) <- sum(I2_tot)

update(tracker_N) <- sum(N)

update(tracker_V) <- sum(V1) +
  sum(V2) +
  sum(EWV1) +
  sum(IWV1) +
  sum(RWV1) +
  sum(EWV2) +
  sum(IWV2) +
  sum(RWV2) +
  sum(EAV1) +
  sum(IAV1) +
  sum(RAV1) +
  sum(EARWV1) +
  sum(IARWV1) +
  sum(EAV2) +
  sum(IAV2) +
  sum(RAV2) +
  sum(EARWV2) +
  sum(IARWV2) +
  sum(EDV1) +
  sum(IDV1) +
  sum(RDV1) +
  sum(EDRWV1) +
  sum(IDRWV1) +
  sum(EDRAV1) +
  sum(IDRAV1) +
  sum(EDV2) +
  sum(IDV2) +
  sum(RDV2) +
  sum(EDRWV2) +
  sum(IDRWV2) +
  sum(EDRAV2) +
  sum(IDRAV2) +
  sum(EOV1) +
  sum(IOV1) +
  sum(ROV1) +
  sum(EORWV1) +
  sum(IORWV1) +
  sum(EORAV1) +
  sum(IORAV1) +
  sum(EORDV1) +
  sum(IORDV1) +
  sum(EOV2) +
  sum(IOV2) +
  sum(ROV2) +
  sum(EORWV2) +
  sum(IORWV2) +
  sum(EORAV2) +
  sum(IORAV2) +
  sum(EORDV2) +
  sum(IORDV2) +
  sum(EXV1) +
  sum(IXV1) +
  sum(RXV1) +
  sum(EXV2) +
  sum(IXV2) +
  sum(RXV2) +
  sum(EORXV1) +
  sum(IORXV1) +
  sum(EORXV2) +
  sum(IORXV2) +
  sum(EXROV1) +
  sum(IXROV1) +
  sum(EXROV2) +
  sum(IXROV2) +
  sum(ROXV1) +
  sum(ROXV2)


