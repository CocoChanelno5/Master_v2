############## COMMENTS ####################
### REGIMES
# 1 - EXPANSION
# 0 - RECESSION
set.seed(42)
setwd("~/Desktop/Master_v2")
source("MC_MG_HF_conv_functions.R")
load("~/Desktop/Master_v2/danedane.RData")

# the loop through all files makes the posterior estimation
posterior_a <- list()
# preparing matryx Y for GDP in USA
dftemp<-USA_GDP_ch
Y<-dftemp%>% select(c(ID,Period, Value)) %>% pivot_wider(names_from = ID,values_from = Value)
Y <- as.matrix(Y[,-1])
W<-W_USA
table(is.na(Y))

######################### PARAMETRY DLA USA GDP  ##########
N <- n_states
theta0 <- list(rho = 0.5,
               mu_1 = rep(4, N),
               mu_0 = rep(2.3, N),
               omega_d = rep(1, N), #VARIANCES (already squared)
               p_00 = rep(0.8, N),
               p_11 = rep(0.8, N))

hyperpar0 = list(alpha_prior = matrix(c(8, 2, 1, 9), nrow = 2, byrow = TRUE),
                 v_prior = 6,
                 delta_prior = 2,  ## changed from 0.4
                 m_prior = matrix(c(2.3,4), nrow = 2),
                 M_prior = diag(2))


start <- Sys.time()
posterior_a <- sample_posterior(initial = theta0, hyperpar = hyperpar0, S = 20000, S0 = 4000, S_rho = 10000, S0_rho = 3000, Y = Y, W = W)
end <- Sys.time()
print(end - start)
save.image(paste0("~/Desktop/Master_v2/post_simul/posterior_USA_GDP_", format(Sys.time(), "%b%d"), ".RData"))
