# ---------------------------------------------------------------------------
# Simulate data based on PROVIDE rotavirus trial
# ---------------------------------------------------------------------------

# Parameters (eventually add to config)
# inflation <- 0        # if negative, protective effect in the doomed
# doomed_epsilon <- 1  # if != 1 violates hudgens assumptions
# epsilon <- 1          # where epsilon <= 1 ; if epsilon not 1, violates our assumption
# n <- 1e5              # Sample size to simulate

#' Function to simulate data to mimic PROVIDE rotavirus trial
simulate_data_provide <- function(seed = 12345,
                          effect_protect = TRUE,
                          doomed_inflation = 0,
                          protected_inflation = 0,
                          doomed_epsilon = 1,
                          protected_epsilon = 1,
                          immune_delta = 0,
                          protected_delta = 0,
                          n = 1e5){
  set.seed(seed)
  data <- data.frame(id = 1:n)
  
  # Covariates -----------------------------------------------------------------
  
  # Week 10 HAZ - N(mean=-0.97, sd=0.90)
  data$wk10_haz <- rnorm(n, mean = -0.97, sd = 0.90)
  
  # Gender - Bernoulli(0.5)
  # let 0 = female, 1 = male
  data$gender <- rbinom(n, 1, 0.5)
  #data$gender <- ifelse(data$gender == 0, "Female", "Male")
  
  # num_hh_sleep
  # NegativeBinomial(mu = 5.26, sigma = 2.5) - discrete, minimum of 1
  
  mu <- 5.26
  sigma <- 2.5
  size <- (mu^2) / (sigma^2 - mu)  # Note: variance = mu + mu^2 / size
  prob <- size / (size + mu)
  
  data$num_hh_sleep <- rnbinom(n, size = size, prob = prob)
  data$num_hh_sleep <- pmax(1, data$num_hh_sleep)  # enforce minimum 1
  
  # Principal Strata ------------------------------------------------------------
  
  # Softmax to guarantee [0,1]

  # adjusted from -2.16 fit in real data to -1.2 to try to get marginal probability closer to observed
  log_odds_doomed__x <- -1.2 + 0.81*as.numeric(data$gender == 1) +
    0.18*data$wk10_haz +
    0.06*data$num_hh_sleep 
  
  # adjusted from 1.29 fit in real data to 1.5 to try to get marginal probability closer to observed
  log_odds_immune__x <- 1.5 - 0.30*as.numeric(data$gender == 1) +
    0.10*data$wk10_haz -
    0.08*data$num_hh_sleep 
  
  # increase immune
  log_odds_immune__x <- log_odds_immune__x + immune_delta
  
  # increase protected (decreasing doomed and immune)
  log_odds_doomed__x <- log_odds_doomed__x + protected_delta
  log_odds_immune__x <- log_odds_immune__x + protected_delta
  
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
  
  #table(data$stratum) / (sum(table(data$stratum)))
  
  # Outcome Probabilities --------------------------------------------------------
  
  # P(Y(1) = 1 | Doomed, X)
  # approximated from data using regression of Y ~ X | Z = 1, S = 1
  data$p_abx_1__doomed <-  plogis(
    -0.70 + 0.78 * as.numeric(data$gender == 1) +
    -1.44 * data$wk10_haz + 0.49 * data$num_hh_sleep
  )
  
  # P(Y(0) = 1 | Doomed, X)
  # doomed_inflation controls size of effect in Doomed principal stratum
  # here we subtract doomed_inflation so that positive values of doomed_inflation lead
  # to larger positive effects in the Doomed
  data$p_abx_0__doomed <- plogis(qlogis(data$p_abx_1__doomed) - doomed_inflation)
  
  # P(Y(0) = 1 | Immune, X)
  data$p_abx_0__immune <- plogis(
    -0.29 + 0.41 * as.numeric(data$gender == 1) +
    -0.10 * data$wk10_haz + 0.13 * data$num_hh_sleep
  )
  # P(Y(1) = 1 | Immune, X)
  # should equal P(Y(0) = 1 | Immune, X) to satisfy stochastic exclusion restriction
  data$p_abx_1__immune <- data$p_abx_0__immune 

  # P(Y(1) = 1 | Immune, X)
  # should equal P(Y(1) = 1 | Immune, X) to satisfy partial principal ignorability
  # required for identification of the Naturally Infected Estimand
  # use protected_epsilon to violate this assumption
  data$p_abx_1__immune <- data$p_abx_1__protect * protected_epsilon

  # P(Y(0) = 1 | Protected, X)
  # protected_inflation controls size of effect in Protected principal stratum
  # NOTE: because this does not equal p_abx_0__doomed, we have a violation
  # of prin ignorability needed to id the doomed estimand. but we are ok with that.
  data$p_abx_0__protect <- plogis(qlogis(data$p_abx_1__protect) - protected_inflation)


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
  
  # Immune, V = 1: 
  is_immune_v1 <- data$stratum == "Immune" & data$rotaarm == 1
  data$rotaepi[is_immune_v1] <- 0
  data$any_abx_wk52[is_immune_v1] <- rbinom(sum(is_immune_v1), 1, data$p_abx_1__immune[is_immune_v1])
  
  # Immune, V = 0: 
  is_immune_v0 <- data$stratum == "Immune" & data$rotaarm == 0
  data$rotaepi[is_immune_v0] <- 0
  data$any_abx_wk52[is_immune_v0] <- rbinom(sum(is_immune_v0), 1, data$p_abx_0__immune[is_immune_v0])
  
  # Protected, V=1:
  is_protected_v1 <- data$stratum == "Protected" & data$rotaarm == 1
  data$rotaepi[is_protected_v1] <- 0
  data$any_abx_wk52[is_protected_v1] <- rbinom(sum(is_protected_v1), 1, data$p_abx_1__protect[is_protected_v1])
  
  # Protected, V=0: 
  is_protected_v0 <- data$stratum == "Protected" & data$rotaarm == 0
  data$rotaepi[is_protected_v0] <- 1 
  data$any_abx_wk52[is_protected_v0] <- rbinom(sum(is_protected_v0), 1, data$p_abx_0__protect[is_protected_v0])
  
  return(data)
  
}

simulate_data_contour <- function(seed = 12345,
                                  effect_protect = TRUE,
                                  doomed_inflation = 0,
                                  protected_inflation = 0,
                                  doomed_epsilon = 1,
                                  protected_epsilon = 1,
                                  n = 1e5){
  set.seed(seed)
  data <- data.frame(id = 1:n)
  
  # Covariates -----------------------------------------------------------------
  
  # make all binary or categorical
  
  # X1 - Binomial(3, 0.25)
  # X2 - Bernoulli(0.25)
  # X3 - Bernoulli(0.75)
  
  data$X1 <- rbinom(n, 3, 0.25)
  data$X2 <- rbinom(n, 1, 0.25)
  data$X3 <- rbinom(n, 1, 0.75)
  
  # Principal Strata ------------------------------------------------------------
  
  data$p_doomed__x <- plogis(-1 + 0.25*data$X1 - 0.4*data$X2 - 0.2*data$X3)
  data$p_immune__x <- plogis(-0.5 + -0.5*data$X1 + 0.25*data$X2 + 0.2*data$X3)
  data$p_protected__x <- 1 - data$p_doomed__x - data$p_immune__x
  
  # Sample the strata
  probs <- cbind(data$p_doomed__x, data$p_immune__x, data$p_protected__x)
  strata <- c("Doomed", "Immune", "Protected")
  data$stratum <- apply(probs, 1, function(p) sample(strata, size = 1, prob = p))
  
  round(prop.table(table(data$stratum)), 3)
  # Doomed    Immune Protected 
  # 0.242     0.343     0.415 
  
  # Outcome Probabilities --------------------------------------------------------
  
  # P(Y(1) = 1 | Doomed) does not have to be equal - controls size of estimand - increase this increases estimate
  
  # P(Y(0) = 1 | Doomed) = P(Y(1) = 1 | Doomed) = P(Y(0) = 1 | Protect) 
  # big ish 
  # Y ~ X | V = 1, S = 1
  
  # if negative, protective effect in the doomed
  # doomed_inflation <- 0
  # doomed_epsilon <- 1
  
  # violate hudgens:
  # P(Y(0) = 1 | Doomed)
  data$p_Y0__doomed <-  plogis(-1 +
                                 0.25 * data$X1 +
                                 -2 * data$X2 +
                                 0.5 * data$X3)
  
  # P(Y(0) = 1 | Protect)
  data$p_Y0__protect <- data$p_Y0__doomed * doomed_epsilon
  
  # P(Y(1) = 1 | Doomed)
  data$p_Y1__doomed <- plogis(qlogis(data$p_Y0__doomed) + doomed_inflation)
  
  # P(Y(0) = 1 | Immune) does not have to be equal & doesn't matter for size of ours or hudgens
  # equal implies no effect of intervention in the immune & should be true logically
  
  # P(Y(0) = 1 | Immune) = P(Y(1) = 1 | Immune) = P(Y(1) = 1 | Protect) 
  # Y ~ X | V = 0, S = 0
  
  # flag to make protected effect = 0
  if(effect_protect){
    # data$p_Y01__immune <- plogis(1 +
    #                                0.1 * data$X1 +
    #                                -1 * data$X2 +
    #                                0.25 * data$X3)
    data$p_Y01__immune <- plogis(-1 +
                                   0.25 * data$X1 +
                                   -2 * data$X2 +
                                   0.5 * data$X3)
  } else{
    # set probability in immune = to probability of protected without abx
    data$p_Y01__immune <- data$p_Y0__protect
  }
  
  data$p_Y1__protect <- data$p_Y01__immune * protected_epsilon
  
  data$p_Y1__protect <- plogis(qlogis(data$p_Y01__immune) + protected_inflation)
  
  # Vaccine (Z) & Infection (S) & Outcome (Y) ----------------------------------------
  
  # P(Z | X) ~ Bernoulli(plogis(X1 + X2 + X3))
  # This keeps it ~50% (IQR 0.465,0.5275); min = 0.44, max = 0.67
  p_Z__X <- plogis(-0.14 + 0.25*data$X1 + 0.1*data$X2 -0.1*data$X3)
  
  data$Z <- rbinom(n, 1, p_Z__X )
  data$S <- NA
  data$Y <- NA
  
  # Doomed, Z = 1:
  is_doomed_z1 <- data$stratum == "Doomed" & data$Z == 1
  data$S[is_doomed_z1] <- 1
  data$Y[is_doomed_z1] <- rbinom(sum(is_doomed_z1), 1, data$p_Y1__doomed[is_doomed_z1])
  
  # Doomed, Z = 0:
  is_doomed_z0 <- data$stratum == "Doomed" & data$Z == 0
  data$S[is_doomed_z0] <- 1
  data$Y[is_doomed_z0] <- rbinom(sum(is_doomed_z0), 1, data$p_Y0__doomed[is_doomed_z0])
  
  # Immune: 
  is_immune <- data$stratum == "Immune"
  data$S[is_immune] <- 0
  data$Y[is_immune] <- rbinom(sum(is_immune), 1, data$p_Y01__immune[is_immune])
  
  # Protected, Z=1:
  is_protected_z1 <- data$stratum == "Protected" & data$Z == 1
  data$S[is_protected_z1] <- 0
  data$Y[is_protected_z1] <- rbinom(sum(is_protected_z1), 1, data$p_Y1__protect[is_protected_z1])
  
  # Protected, Z=0: 
  is_protected_z0 <- data$stratum == "Protected" & data$Z == 0
  data$S[is_protected_z0] <- 1 
  data$Y[is_protected_z0] <- rbinom(sum(is_protected_z0), 1, data$p_Y0__protect[is_protected_z0])
  
  return(data)
  
}


# Function to simulate generic data 
simulate_data_generic <- function(seed = 12345,
                                  effect_protect = TRUE,
                                  doomed_inflation = 0,
                                  doomed_epsilon = 1,
                                  protected_epsilon = 1,
                                  immune_epsilon = 1,
                                  n = 1e5){
  set.seed(seed)
  data <- data.frame(id = 1:n)
  
  # Covariates -----------------------------------------------------------------
  
  # make all binary or categorical
  
  # X1 - Binomial(3, 0.25)
  # X2 - Bernoulli(0.25)
  # X3 - Bernoulli(0.75)
  
  data$X1 <- rbinom(n, 1, 0.5)
  data$X2 <- rbinom(n, 1, 0.5)
  data$X3 <- rbinom(n, 1, 0.5)
  
  # Principal Strata ------------------------------------------------------------
  
  data$p_doomed__x <- plogis(-1 + 0.5*data$X1 - 1*data$X1*data$X2 - 0.5*data$X3)
  data$p_immune__x <- plogis(-1 + 0.5*data$X1 - 1*data$X3*data$X1 - 0.5*data$X3)
  data$p_protected__x <- 1 - data$p_doomed__x - data$p_immune__x
  
  # Sample the strata
  probs <- cbind(data$p_doomed__x, data$p_immune__x, data$p_protected__x)
  strata <- c("Doomed", "Immune", "Protected")
  data$stratum <- apply(probs, 1, function(p) sample(strata, size = 1, prob = p))
  
  #round(prop.table(table(data$stratum)), 3)
  
  # Outcome Probabilities --------------------------------------------------------
  
  # P(Y(1) = 1 | Doomed) does not have to be equal - controls size of estimand - increase this increases estimate
  
  # P(Y(0) = 1 | Doomed) = P(Y(1) = 1 | Doomed) = P(Y(0) = 1 | Protect) 
  # big ish 
  # Y ~ X | V = 1, S = 1
  
  # if negative, protective effect in the doomed
  # doomed_inflation <- 0
  # doomed_epsilon <- 1
  
  # violate hudgens:
  # P(Y(0) = 1 | Doomed)
  data$p_Y0__doomed <-  plogis(-1 +
                               0.5 * data$X1 +                                                     
                              -1 * data$X2*data$X1 +
                               0.5 * data$X3)
  
  # P(Y(0) = 1 | Protect)
  data$p_Y0__protect <- data$p_Y0__doomed * doomed_epsilon
  
  # P(Y(1) = 1 | Doomed)
  data$p_Y1__doomed <- plogis(qlogis(data$p_Y0__doomed) + doomed_inflation)
  
  # P(Y(0) = 1 | Immune) does not have to be equal & doesn't matter for size of ours or hudgens
  # equal implies no effect of intervention in the immune & should be true logically
  
  # P(Y(0) = 1 | Immune) = P(Y(1) = 1 | Immune) = P(Y(1) = 1 | Protect) 
  # Y ~ X | V = 0, S = 0
  
  # flag to make protected effect = 0
  if(effect_protect){
    # ORIGINAL
    # data$p_Y01__immune <- plogis(-0.5 +                                                            
    #                              0.5 * data$X1 +                                                   
    #                              -1 * data$X3 * data$X1 +                                          
    #                              0.5 * data$X2)
    
    data$p_Y0__immune <- plogis(-0.5 +                                                            
                                   0.5 * data$X1 +                                                   
                                   -1 * data$X3 * data$X1 +                                          
                                   0.5 * data$X2)
    data$p_Y1__immune <- data$p_Y0__immune * immune_epsilon
    
  } else{
    # ORIGINAL 
    # set probability in immune = to probability of protected without abx
    # data$p_Y01__immune <- data$p_Y0__protect
    data$p_Y0__immune <- data$p_Y0__protect
    data$p_Y1__immune <- data$p_Y0__immune * immune_epsilon
  }
  
  # is this right? assume vaccinated protected == vaccinated immune? 
  # vaccinated immune may or may not equal unvaccinated immune now based on epsilon
  data$p_Y1__protect <- data$p_Y1__immune * protected_epsilon
  
  # Vaccine (Z) & Infection (S) & Outcome (Y) ----------------------------------------
  
  # P(Z | X) ~ Bernoulli(plogis(X1 + X2 + X3))
  p_Z__X <- plogis(-0.14 - 0.5*data$X1 + 1*data$X1*data$X2 - 1.2*data$X3)
  
  data$Z <- rbinom(n, 1, p_Z__X )
  data$S <- NA
  data$Y <- NA
  
  # Doomed, Z = 1:
  is_doomed <- data$stratum == "Doomed"
  data$S0[is_doomed] <- 1
  data$S1[is_doomed] <- 1
  data$Y1[is_doomed] <- rbinom(sum(is_doomed), 1, data$p_Y1__doomed[is_doomed])
  data$Y0[is_doomed] <- rbinom(sum(is_doomed), 1, data$p_Y0__doomed[is_doomed])
  
  data$S[is_doomed & data$Z == 1] <- data$S1[is_doomed & data$Z == 1]
  data$S[is_doomed & data$Z == 0] <- data$S0[is_doomed & data$Z == 0]
  data$Y[is_doomed & data$Z == 1] <- data$Y1[is_doomed & data$Z == 1]
  data$Y[is_doomed & data$Z == 0] <- data$Y0[is_doomed & data$Z == 0]
  
  # Immune: 
  is_immune <- data$stratum == "Immune"
  data$S0[is_immune] <- 0
  data$S1[is_immune] <- 0
  # data$Y1[is_immune] <- rbinom(sum(is_immune), 1, data$p_Y01__immune[is_immune])
  # data$Y0[is_immune] <- data$Y1[is_immune]
  data$Y1[is_immune] <- rbinom(sum(is_immune), 1, data$p_Y1__immune[is_immune])
  data$Y0[is_immune] <- rbinom(sum(is_immune), 1, data$p_Y0__immune[is_immune])
  
  data$S[is_immune & data$Z == 1] <- data$S1[is_immune & data$Z == 1]
  data$S[is_immune & data$Z == 0] <- data$S0[is_immune & data$Z == 0]
  data$Y[is_immune & data$Z == 1] <- data$Y1[is_immune & data$Z == 1]
  data$Y[is_immune & data$Z == 0] <- data$Y0[is_immune & data$Z == 0]
  
  # Protected, Z=1:
  is_protected <- data$stratum == "Protected" 
  data$S0[is_protected] <- 1
  data$S1[is_protected] <- 0
  data$Y1[is_protected] <- rbinom(sum(is_protected), 1, data$p_Y1__protect[is_protected])
  data$Y0[is_protected] <- rbinom(sum(is_protected), 1, data$p_Y0__protect[is_protected])
  
  data$S[is_protected & data$Z == 1] <- data$S1[is_protected & data$Z == 1]
  data$S[is_protected & data$Z == 0] <- data$S0[is_protected & data$Z == 0]
  data$Y[is_protected & data$Z == 1] <- data$Y1[is_protected & data$Z == 1]
  data$Y[is_protected & data$Z == 0] <- data$Y0[is_protected & data$Z == 0]
  
  return(data)
   
}

