sim <- new.env(parent = emptyenv())
# Population basis for simulation
sim$pop <- data.frame()
# Number of subjects in simulation
sim$N <- 0
# PK/PD simulation model
sim$pkpd <- NA
# Initializer factory for sim model
sim$inits_fac <- NA
# Default parameters
sim$params.default <- numeric()
