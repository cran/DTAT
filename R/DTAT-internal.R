sim <- new.env(parent = emptyenv())
# Population basis for simulation
sim$pop <- data.frame()
# Number of subjects in simulation
sim$N <- 0
# PK/PD simulation model
sim$pkpd <- NA
# Default parameters
sim$params.default <- numeric()
