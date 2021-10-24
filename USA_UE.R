############## COMMENTS ####################
### REGIMES
# 1 - EXPANSION
# 0 - RECESSION
set.seed(42)
setwd("~/Desktop/Master_v2")
source("MC_MG_HF_conv_functions.R")
load("~/Desktop/Master_v2/danedane.RData")

## preparing matryx Y for unemployment rate in USA
dftemp<-USA_UE_ch
Y<-dftemp%>% select(c(Name,Date, Value)) %>% pivot_wider(names_from = Name,values_from = Value)
Y <- as.matrix(Y[,-1])
W<-W_USA
table(is.na(Y))

######################### PARAMETRY DLA USA STOPA BEZROBOCIA  ##########
N <- n_states
theta0 <- list(rho = 0.5,
               mu_1 = rep(-0.2, N),
               mu_0 = rep(0.3, N),
               omega_d = rep(1, N), #VARIANCES (already squared)
               p_00 = rep(0.8, N),
               p_11 = rep(0.8, N))

hyperpar0 = list(alpha_prior = matrix(c(8, 2, 1, 9), nrow = 2, byrow = TRUE),
                 v_prior = 6,
                 delta_prior = 8,
                 m_prior = matrix(c(0.5,-1), nrow = 2),
                 M_prior = diag(2))


start <- Sys.time()
posterior_a <- sample_posterior(initial = theta0, hyperpar = hyperpar0, S = 20000, S0 = 4000, S_rho = 5000, S0_rho = 1000, Y = Y, W = W)
end <- Sys.time()
print(end - start)
save.image(paste0("~/Desktop/Master_v2/post_simul/posterior_USA_UE_", format(Sys.time(), "%b%d"), ".RData"))

