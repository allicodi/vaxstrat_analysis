
#' Function to simulate PROVIDE-like data 
#' 
#' @param seed seed for replicability
#' @param effect_protect boolean for effect in protected stratum, default TRUE
#' @param doomed_inflation numeric modify effect size in doomed stratum, default 0
#' @param protected_inflation numeric modify effect size in protected stratum, default 0
#' @param doomed_epsilon numeric violate assumption protected w abx == immune (assumption for naturally infected), default 1 (no violation)
#' @param protected_epsilon numeric violate assumption effect protected w/o abx == immune (assumption for doomed estimand), default 1 (no violation)
#' @param immune_delta numeric modify size of immune stratum
#' @param protected_delta numeric modify size of protected stratum
#' @param n sample size to generate
#' 
#' @returns simulated dataset
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
  
  # P(Y(1) = 1 | Doomed) does not have to be equal - controls size of estimand - increase this increases estimate
  
  # P(Y(0) = 1 | Doomed) = P(Y(1) = 1 | Doomed) = P(Y(0) = 1 | Protect) 
  # big ish 
  # Y ~ X | V = 1, S = 1
  
  # if negative, protective effect in the doomed
  # doomed_inflation <- 0
  # doomed_epsilon <- 1
  
  # violate hudgens doomed assumption:
  # P(Y(0) = 1 | Doomed)
  data$p_abx_0__doomed <-  plogis(-0.70 +
                                    0.78 * as.numeric(data$gender == 1) +
                                    -1.44 * data$wk10_haz +
                                    0.49 * data$num_hh_sleep)
  
  # P(Y(0) = 1 | Protect)
  data$p_abx_0__protect <- data$p_abx_0__doomed * doomed_epsilon
  
  # P(Y(1) = 1 | Doomed)
  data$p_abx_1__doomed <- plogis(qlogis(data$p_abx_0__doomed) + doomed_inflation)
  
  # P(Y(0) = 1 | Immune) does not have to be equal & doesn't matter for size of ours or hudgens
  # equal implies no effect of intervention in the immune & should be true logically
  
  # P(Y(0) = 1 | Immune) = P(Y(1) = 1 | Immune) = P(Y(1) = 1 | Protect) 
  # Y ~ X | V = 0, S = 0
  
  
  # flag to make protected effect = 0 if false, default should be true
  if(effect_protect){
    data$p_abx_01__immune <- plogis(-0.29 +
                                      0.41 * as.numeric(data$gender == 1) +
                                      -0.10 * data$wk10_haz +
                                      0.13 * data$num_hh_sleep)
  } else{
    # set probability in immune = to probability of protected without abx
    data$p_abx_01__immune <- data$p_abx_0__protect
  }
  
  data$p_abx_1__protect <- data$p_abx_01__immune * protected_epsilon
  
  data$p_abx_1__protect <- plogis(qlogis(data$p_abx_1__protect) + protected_inflation)
  
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
  
}