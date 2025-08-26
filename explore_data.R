# ------------------------------------------------------------------------------
# Examine real PROVIDE data to use as base for simulation data
# ------------------------------------------------------------------------------

here::i_am("simulations/explore_data.R")

# Read in per-protocol cleaned data
data <- readRDS(here::here("provide_data/per_protocol_data.Rds"))

# drop person missing wk10_haz
data <- data[-which(is.na(data$wk10_haz)),]

# Describe marginally -------------------------------------------------------

# 1:1 vaccine vs placebo

table(data$rotaarm,data$rotaepi)
# 0   1
# 0 171  93
# 1 215  45

# marginally
# 93 / (93 + 171) = 0.35 placebo arm
# 45 / (45 + 215) = 0.17 vax arm

# lasso w/ all covariates, rotaapi

# Initial model fit
fit <- glm(rotaepi ~ rotaarm + gender + wk10_haz + 
             wk10_ageday + num_hh_lt_5 + num_hh_sleep + 
             fedu_bin + medu_bin + 
             inco + elec + gas + tv + toil_bin + watr_bin + food_avail_bin, 
           data = data,
           family = stats::binomial())

# Outcome as binary numeric
y <- as.numeric(data$rotaepi)

# Matrix of predictors — exclude the outcome variable
x <- model.matrix(~ rotaarm + gender + wk10_haz + 
                    wk10_ageday + num_hh_lt_5 + num_hh_sleep + 
                    fedu_bin + medu_bin + 
                    inco + elec + gas + tv + toil_bin + watr_bin + food_avail_bin, 
                  data = data)[, -1]  # remove intercept

# LASSO for variable selection
fit_lasso <- glmnet(x, y, family = "binomial", alpha = 1)
cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 1)
coef(cv_fit, s = "lambda.min")

# vaccine arm, gender, num_hh_sleep
# also throw in wk10_HAZ? then we've got vaccine, binary var (gender), cont var (HAZ),  discrete continuous (num_hh_sleep)

fit2 <- glm(rotaepi ~ rotaarm + gender + wk10_haz + num_hh_sleep, 
            data = data, 
            family = stats::binomial())

# Intercept: -1.28
# rotaarm: -0.97
# gender: 0.49
# wk10_haz: 0.0003 ---> let's pretend this is -0.1?? 
# num_hh_sleep: 0.08

# P(S | Z, X) ^^ without accounting for the strata

# gender 1:1
hist(as.numeric(data$gender))

# HAZ - normal (mean = -0.97, sd = 0.90)
hist(data$wk10_haz)
mean(data$wk10_haz)
sd(data$wk10_haz)

# num_hh_sleep - discrete gamma ? with minimum of 1
summary(data$num_hh_sleep)

mu <- 5.26
sigma <- 2.5

shape <- (mu / sigma)^2
rate <- mu / sigma^2

# ------------------------------------------------------------------------
# Get strata probs -------------------------------------------------------
# ------------------------------------------------------------------------

## Label strata we know

# Infected people in the vaccine arm are doomed- we don't know who is doomed in placebo arm
data$doomed <- ifelse(data$rotaarm == 1 & data$rotaepi == 1, 1, 0)
data$doomed <- ifelse(data$rotaarm == 0, NA, data$doomed)

# Uninfected people in the placebo arm are immune- we don't know who is immune in vaccine arm
data$immune <- ifelse(data$rotaarm == 0 & data$rotaepi == 0, 1, 0)
data$immune <- ifelse(data$rotaarm == 1, NA, data$immune)

# P(Doomed) = 0.17
# P(Immune) = 0.64

## 1. P(Doomed | X) ------------------------------------------------------------

p_doomed__x_fit <- glm(doomed ~ gender + wk10_haz + num_hh_sleep, 
                       data = data[data$rotaarm == 1,],
                       family = stats::binomial())

# Intercept: -2.16
# Gender (female = 0, male =1 ): 0.81
# HAZ: 0.18 <-- we want this to be opposite direction clinically?
# num_hh_sleep: 0.06

data$p_doomed__x <- predict(p_doomed__x_fit, newdata = data, type = 'response')

# But we observe the people in the vaccine arm and know if they're doomed so set to 1 or 0
data$p_doomed__x <- ifelse(data$rotaarm == 1, data$doomed, data$p_doomed__x)
# and we know people in placebo arm who did not get rota are immune, so set to 0
data$p_doomed__x <- ifelse(data$rotaarm == 0 & data$rotaepi == 0, 0, data$p_doomed__x)

data$p_notdoomed__x <- abs(1 - data$p_doomed__x)

# but for the people in the vaccine arm probability needs to be 1 or 0?? we directly observe them
# data$p_doomed__x <- ifelse(data$rotaarm == 1, data$doomed, data$p_doomed__x)
# data$p_notdoomed__x <- ifelse(data$rotaarm  == 1, abs(1 - data$doomed), 1 - data$p_doomed__x)

# and for the people we directly observe to be immune we also know their probability of being doomed is 0
#data$p_doomed__x <- ifelse(data$rotaarm == 0 & data$immune == 1, 0, data$p_doomed__x)
#data$p_notdoomed__x <- ifelse(data$rotaarm == 0 & data$immune == 1, 1, data$p_doomed__x)

# ^ maybe not because p_doomed__x, not p_doomed__x_z? 

## 2. P(Immune | X) ------------------------------------------------------------

p_immune__x_fit <- glm(immune ~ gender + wk10_haz + num_hh_sleep,
                       data = data[data$rotaarm == 0,],
                       family = stats::binomial())
# Intercept: 1.29
# Gender: -0.30
# wk10_haz: 0.10
# num_hh_sleep: -0.08

## 3. P(Protected | X) ------------------------------------------------------------

# Protected = 1 - Immune - Doomed

# ------------------------------------------------------------------------------
# Outcome Simulation - Binary
# ------------------------------------------------------------------------------

# P(Y(0) = 1 | Doomed) = P(Y(1) = 1 | Doomed) = P(Y(0) = 1 | Protect) 
# big ish 
# Y ~ X | V = 1, S = 1

abx_x__v_1_s_1 <- glm(any_abx_wk52 ~ gender + wk10_haz + num_hh_sleep,
                      data = data[which(data$rotaarm == 1 & data$rotaepi == 1),],
                      family = stats::binomial())

# Intercept: -0.70
# Gender: 0.78
# HAZ: -1.44
# Sleep: 0.49

summary(predict(abx_x__v_1_s_1, type = 'response'))
# Probabilities all >0.75, 1stQ 0.94 (on training data)
summary(predict(abx_x__v_1_s_1, newdata = data, type = 'response'))
# Probabilities >0.4, 1stQ 0.92 (on full data)

# is that too high?? 

abx_x__v_0_s_0 <- glm(any_abx_wk52 ~ gender + wk10_haz + num_hh_sleep,
                      data = data[which(data$rotaarm == 0 & data$rotaepi == 0),],
                      family = stats::binomial())

# Intercept: -0.29
# Gender: 0.41
# HAZ: -0.10
# Sleep: 0.13

summary(predict(abx_x__v_0_s_0, type = 'response'))
# Probabilities all >0.5, 1stQ 0.6, 3rdQ 0.71 (on training data)
summary(predict(abx_x__v_0_s_0, newdata = data, type = 'response'))
# Probabilities >0.5, 1stQ 0.61, 3rdQ 0.72 (on full data)


# ------------------------------------------------------------------------------
# Outcome Simulation - Continuous
# ------------------------------------------------------------------------------

# P(Y(0) = 1 | Doomed) = P(Y(1) = 1 | Doomed) = P(Y(0) = 1 | Protect) 
# big ish 
# Y ~ X | V = 1, S = 1

growth_x__v_1_s_1 <- glm(wk52_haz ~ gender + wk10_haz + num_hh_sleep,
                      data = data[which(data$rotaarm == 1 & data$rotaepi == 1),],
                      family = stats::gaussian())

# Intercept: -0.62
# Gender: 0.13
# HAZ: 0.85
# Sleep: -0.03

summary(predict(growth_x__v_1_s_1, type = 'response'))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -2.7971 -1.9727 -1.4782 -1.4673 -1.0730 -0.2673 

summary(predict(growth_x__v_1_s_1, newdata = data, type = 'response'))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -4.3520 -2.0456 -1.4992 -1.5513 -1.0400  0.8071 

growth_x__v_0_s_0 <- glm(wk52_haz ~ gender + wk10_haz + num_hh_sleep,
                      data = data[which(data$rotaarm == 0 & data$rotaepi == 0),],
                      family = stats::gaussian())

# Intercept: -0.92
# Gender: 0.04
# HAZ: 0.78
# Sleep: 0.03

summary(predict(growth_x__v_0_s_0 , type = 'response'))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -4.0830 -2.0135 -1.3110 -1.4160 -0.9385  0.5809 

summary(predict(growth_x__v_0_s_0 , newdata = data, type = 'response'))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -4.0830 -1.9265 -1.4228 -1.4588 -1.0080  0.6529 

