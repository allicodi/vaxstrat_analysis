
seed <- 1
n <- 1e7

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
deltas <- expand.grid(
  immune_delta = seq(-2, 0, by = 0.01),
  protected_delta = seq(-2, 0, by = 0.01)
) 

deltas$doomed <- NA
deltas$immune <- NA
deltas$protected <- NA

# final combos
deltas <- data.frame(immune_delta = c(-1.18, -0.32, 0.6, 0.18, 0.6, 1.3),
                     protected_delta = c(-0.13, -0.11, -0.17, 0.53, -0.17, -1.09))

for(i in 1:nrow(deltas)){
  
  log_odds_doomed__x <- -1.2 + 0.81*as.numeric(data$gender == 1) +
    0.18*data$wk10_haz +
    0.06*data$num_hh_sleep +
    deltas[i, "protected_delta"]
  
  log_odds_immune__x <- 1.5 - 0.30*as.numeric(data$gender == 1) +
    0.10*data$wk10_haz -
    0.08*data$num_hh_sleep +
    deltas[i, "immune_delta"] +
    deltas[i, "protected_delta"]
  
  # Softmax probabilities
  denom <- 1 + exp(log_odds_doomed__x) + exp(log_odds_immune__x)
  p_doomed <- mean(exp(log_odds_doomed__x) / denom)
  p_immune <- mean(exp(log_odds_immune__x) / denom)
  p_protected <- 1 - p_doomed - p_immune
  
  deltas[i, 3:5] <- c(p_doomed, p_immune, p_protected)
}

best_deltas <- function(deltas, doomed, immune, protected){
  target <- c(doomed = doomed, immune = immune, protected = protected)
  deltas$distance <- (deltas$doomed - target["doomed"])^2 +
    (deltas$immune - target["immune"])^2 +
    (deltas$protected - target["protected"])^2
  
  best_deltas <- deltas %>% arrange(distance) %>% slice(1)
  best_deltas
}

best_deltas(deltas, doomed = 0.10, immune = 0.70, protected = 0.20)
best_deltas(deltas, doomed = 0.23, immune = 0.30, protected = 0.47)
best_deltas(deltas, doomed = 0.17, immune = 0.50, protected = 0.33)
best_deltas(deltas, doomed = 0.15, immune = 0.70, protected = 0.15)
best_deltas(deltas, doomed = 0.05, immune = 0.70, protected = 0.25)
