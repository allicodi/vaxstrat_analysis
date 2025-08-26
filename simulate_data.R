# ---------------------------------------------------------------------------
# Simulate data based on PROVIDE rotavirus trial
# ---------------------------------------------------------------------------

# Parameters (eventually add to config)
# inflation <- 0        # if negative, protective effect in the doomed
# doomed_epsilon <- 1  # if != 1 violates hudgens assumptions
# epsilon <- 1          # where epsilon <= 1 ; if epsilon not 1, violates our assumption
# n <- 1e5              # Sample size to simulate

simulate_data <- function(seed = 12345,
                          effect_protect = TRUE,
                          inflation = 0,
                          doomed_epsilon = 1,
                          nat_inf_epsilon = 1,
                          n = 1e5){
  set.seed(seed)
  data <- data.frame(id = 1:n)
  
  # Covariates -----------------------------------------------------------------
  
  # Week 10 HAZ - N(mean=-0.97, sd=0.90)
  data$wk10_haz <- rnorm(n, mean = -0.97, sd = 0.90)
  
  # Gender - Bernoulli(0.5)
  data$gender <- rbinom(n, 1, 0.5)
  data$gender <- ifelse(data$gender == 0, "Female", "Male")
  
  # num_hh_sleep
  # NegativeBinomial(mu = 5.26, sigma = 2.5) - discrete, minimum of 1
  
  mu <- 5.26
  sigma <- 2.5
  size <- (mu^2) / (sigma^2 - mu)  # Note: variance = mu + mu^2 / size
  prob <- size / (size + mu)
  
  data$num_hh_sleep <- rnbinom(n, size = size, prob = prob)
  data$num_hh_sleep <- pmax(1, data$num_hh_sleep)  # enforce minimum 1
  
  # Principal Strata ------------------------------------------------------------
  
  # logit(P(Doomed | X)) = -2.16 + 0.81Gender + 0.18HAZ + 0.06Sleep 
  # data$p_doomed__x <- plogis(-2.16 + 
  #                            0.81*as.numeric(data$gender == "Male") +
  #                            0.18*data$wk10_haz +
  #                            0.06*data$num_hh_sleep)
  
  # logit(P(Immune | X)) = 1.29 + -0.30*Gender + 0.10HAZ + -0.09Sleep
  
  # data$p_immune__x <- plogis(1.29 +
  #                              -0.30*as.numeric(data$gender == "Male") +
  #                              0.10*data$wk10_haz +
  #                              -0.08*data$num_hh_sleep)
  
  # P(Protected | X) = 1 - P(Doomed | X) - P(Immune | X)
  # data$p_protected__x <- 1 - data$p_doomed__x - data$p_immune__x
  
  # ^ not all >0 :(
  # ~1000 / 1e6 < 0
  # male with large HAZ
  
  
  # chatgpt recommends softmax? to guarantee [0,1]
  
  # adjusted from -2.16 to -1.2 to try to get marginal probability closer to observed
  log_odds_doomed__x <- -1.2 + 0.81*as.numeric(data$gender == "Male") +
    0.18*data$wk10_haz +
    0.06*data$num_hh_sleep
  
  # adjusted from 1.29 to 1.5 to try to get marginal probability closer to observed
  log_odds_immune__x <- 1.5 - 0.30*as.numeric(data$gender == "Male") +
    0.10*data$wk10_haz -
    0.08*data$num_hh_sleep
  
  # Softmax transformation
  denom <- 1 + exp(log_odds_doomed__x) + exp(log_odds_immune__x)
  
  data$p_doomed__x <- exp(log_odds_doomed__x) / denom
  #^the doomed probability dist smaller mean than model in original data, no softmax
  data$p_immune__x <- exp(log_odds_immune__x) / denom
  data$p_protected__x <- 1 / denom
  
  ## Sample the strata
  probs <- cbind(data$p_doomed__x, data$p_immune__x, data$p_protected__x)
  strata <- c("Doomed", "Immune", "Protected")
  data$stratum <- apply(probs, 1, function(p) sample(strata, size = 1, prob = p))
  
  # true marginal probabilities = 17% doomed, 64% immune, 19% protected
  # marginally 7.3% doomed, 60.2% immune, 32.4% protected based on original models
  # tweaking intercepts --> 15% doomed, 59% immune, 26% protected 
  # close enough and maybe good to have more protected in simulation??
  
  # Outcome Probabilities --------------------------------------------------------
  
  # P(Y(1) = 1 | Doomed) does not have to be equal - controls size of estimand - increase this increases estimate
  
  # P(Y(0) = 1 | Doomed) = P(Y(1) = 1 | Doomed) = P(Y(0) = 1 | Protect) 
  # big ish 
  # Y ~ X | V = 1, S = 1
  
  # if negative, protective effect in the doomed
  # inflation <- 0
  # doomed_epsilon <- 1
  
  # violate hudgens:
  # P(Y(0) = 1 | Doomed)
  data$p_abx_0__doomed <-  plogis(-0.70 +
                                    0.78 * as.numeric(data$gender == "Male") +
                                    -1.44 * data$wk10_haz +
                                    0.49 * data$num_hh_sleep)
  
  # P(Y(0) = 1 | Protect)
  data$p_abx_0__protect <- data$p_abx_0__doomed * doomed_epsilon
  
  # P(Y(1) = 1 | Doomed)
  data$p_abx_1__doomed <- plogis(qlogis(data$p_abx_0__doomed) + inflation)
  
  # P(Y(0) = 1 | Immune) does not have to be equal & doesn't matter for size of ours or hudgens
  # equal implies no effect of intervention in the immune & should be true logically
  
  # P(Y(0) = 1 | Immune) = P(Y(1) = 1 | Immune) = P(Y(1) = 1 | Protect) 
  # Y ~ X | V = 0, S = 0
  
  # data$p_abx_01__immune <- plogis(-0.29 +
  #                                   0.41 * as.numeric(data$gender == "Male") +
  #                                   -0.10 * data$wk10_haz +
  #                                   0.13 * data$num_hh_sleep)
  
  # NEW / CHECK IF OK- flag to make protected effect = 0??
  if(effect_protect){
    data$p_abx_01__immune <- plogis(-0.29 +
                                      0.41 * as.numeric(data$gender == "Male") +
                                      -0.10 * data$wk10_haz +
                                      0.13 * data$num_hh_sleep)
  } else{
    # set probability in immune = to probability of protected without abx
    data$p_abx_01__immune <- data$p_abx_0__protect
  }
  
  data$p_abx_1__protect <- data$p_abx_01__immune * nat_inf_epsilon
  
  # Vaccine & Rotaepi ------------------------------------------------------------
  
  data$rotaarm <- rbinom(n, 1, 0.5)
  
  # Rotavirus + Abx Outcome ----------------------------------------------------------------------
  
  data$rotaepi <- NA
  data$any_abx_wk52 <- NA
  
  # Doomed, V = 1:
  is_doomed_v1 <- data$stratum == "Doomed" & data$rotaarm == 1
  data$rotaepi[is_doomed_v1] <- 1
  data$any_abx_wk52[is_doomed_v1] <- rbinom(sum(is_doomed_v1), 1, data$p_abx_1__doomed[is_doomed_v1])
  
  # Doomed, V = 0:
  is_doomed_v0 <- data$stratum == "Doomed" & data$rotaarm == 0
  data$rotaepi[is_doomed_v0] <- 1
  data$any_abx_wk52[is_doomed_v0] <- rbinom(sum(is_doomed_v0), 1, data$p_abx_0__doomed[is_doomed_v0])
  
  # Immune: 
  is_immune <- data$stratum == "Immune"
  data$rotaepi[is_immune] <- 0
  data$any_abx_wk52[is_immune] <- rbinom(sum(is_immune), 1, data$p_abx_01__immune[is_immune])
  
  # Protected, V=1:
  is_protected_v1 <- data$stratum == "Protected" & data$rotaarm == 1
  data$rotaepi[is_protected_v1] <- 0
  data$any_abx_wk52[is_protected_v1] <- rbinom(sum(is_protected_v1), 1, data$p_abx_1__protect[is_protected_v1])
  
  # Protected, V=0: 
  is_protected_v0 <- data$stratum == "Protected" & data$rotaarm == 0
  data$rotaepi[is_protected_v0] <- 1 
  data$any_abx_wk52[is_protected_v0] <- rbinom(sum(is_protected_v0), 1, data$p_abx_0__protect[is_protected_v0])
  
  return(data)
  #saveRDS(data, here::here(paste0("sim_data/sim_data_", n, ".Rds")))
  
}

# sims for different sample size, epsilons, inflation
# flexible models -- earth, gam, stepwise interaction, ... 

# marginal effect, doomed only (pt identification), 
# doomed only bound, naturally infected (pt identification), naturally infected bound

# figure out what bound we care abt / one sided hypothesis test
# null > 0, alt < 0 --> we care abt upper bound (least neg val); reject if upper CI of upper bound < 0

# Checks: --------------------------------------------------------------------

# Real data:
# P(abx) = 0.74
# P(abx | vax = 1) = 0.73
# P(abx | vax = 0) = 0.74
# P(abx | rotaepi = 1) = 0.91
# P(abx | rotaepi = 1 & vax = 1) = 0.96

# Sim data:
# P(abx) = 0.73
# P(abx | vax = 1) = 0.70
# P(abx | vax = 0) = 0.77
# P(abx | rotaepi = 1) = 0.94
# P(abx | rotaepi = 1 & vax = 1) = 0.942

# P(Y(1) = 1 | Protected) == P(Y(1) = 1 | Immune)
# mean(data$any_abx_wk52[data$stratum == "Protected" & data$rotaarm == 1])
# [1] 0.6705603
# > mean(data$any_abx_wk52[data$stratum == "Immune" & data$rotaarm == 1])
# [1] 0.6492297

# mean(data$wk52_haz[data$stratum == "Protected" & data$rotaarm == 1])
# -1.544
# mean(data$wk52_haz[data$stratum == "Immune" & data$rotaarm == 1])
# -1.5001

# P(Y(0) = 1 | Protected) == P(Y(0) = 1 | Doomed)
# > mean(data$any_abx_wk52[data$stratum == "Protected" & data$rotaarm == 0])
# [1] 0.9388416
# > mean(data$any_abx_wk52[data$stratum == "Doomed" & data$rotaarm == 0])
# [1] 0.9413367

# mean(data$wk52_haz[data$stratum == "Protected" & data$rotaarm == 0])
# -1.60
# mean(data$wk52_haz[data$stratum == "Doomed" & data$rotaarm == 0])
# -1.47


# ^ everything close enough?? for abx, not really for growth


# # to get truth
# data <- simulate_data(
#   seed = 12345,
#   inflation = 0,
#   doomed_epsilon = 1,
#   nat_inf_epsilon = 1,
#   n = 1e6
# )
# 
# naturally infected estimands
# # E[Y(1) | Protected or Doomed] = E[Y(1) | rotaarm = 1, Protected or Doomed]
# # E[ any_abx_wk52 | rotaarm = 1, Protected or Doomed ]
# mean(data$any_abx_wk52[
#   data$rotaarm == 1 & 
#   data$stratum %in% c("Protected", "Doomed")
# ])
# 
# mean(data$any_abx_wk52[
#   data$rotaarm == 0 & 
#     data$stratum %in% c("Protected", "Doomed")
# ])
# 
# # hudgens estimands
# mean(data$any_abx_wk52[
#   data$rotaarm == 1 & 
#     data$stratum %in% c("Doomed")
# ])
# 
# mean(data$any_abx_wk52[
#   data$rotaarm == 0 & 
#     data$stratum %in% c("Doomed")
# ])
# 
# # population estimands
# mean(data$any_abx_wk52[
#   data$rotaarm == 1
# ])
# 
# mean(data$any_abx_wk52[
#   data$rotaarm == 0
# ])
