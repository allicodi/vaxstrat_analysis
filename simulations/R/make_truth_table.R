# -------------------------------------------------------------------
# Script to make table with truth from simulations for supplement
# -------------------------------------------------------------------

here::i_am("R/make_truth_table.R")

library(dplyr)
library(knitr)
library(kableExtra)

pi_1_er_1 <- readRDS(here::here("truth/default_truth.Rds"))
pi_1_er_0 <- readRDS(here::here("truth/violate_er_truth.Rds"))
pi_0_er_1 <- readRDS(here::here("truth/violate_cw_truth.Rds"))
pi_0_er_0 <- readRDS(here::here("truth/violate_cw_er_truth.Rds"))

data.frame(pi = c(1, 1, 0, 0),
           er = c(1, 0, 1, 0),
           E_Y1 = c(
             pi_1_er_1$E_Y1__protected_or_doomed,
             pi_1_er_0$E_Y1__protected_or_doomed,
             pi_0_er_1$E_Y1__protected_or_doomed, 
             pi_0_er_0$E_Y1__protected_or_doomed
           ),
           E_Y0 = c(
             pi_1_er_1$E_Y0__protected_or_doomed,
             pi_1_er_0$E_Y0__protected_or_doomed,
             pi_0_er_1$E_Y0__protected_or_doomed, 
             pi_0_er_0$E_Y0__protected_or_doomed
           ),
           additive_effect = c(
             pi_1_er_1$effect_nat_inf,
             pi_1_er_0$effect_nat_inf,
             pi_0_er_1$effect_nat_inf, 
             pi_0_er_0$effect_nat_inf
           ),
           mult_effect = c(
             pi_1_er_1$effect_nat_inf_mult,
             pi_1_er_0$effect_nat_inf_mult,
             pi_0_er_1$effect_nat_inf_mult, 
             pi_0_er_0$effect_nat_inf_mult
           ))

# make latex table where rows say 
# PI and ER satisfied 
# PI satisfied, ER violated
# PI violated, ER satisfied
# PI and ER violated

# columns are 
# $E\{Y(1) \mid S(0) = 1, S(1) = 0 \}
# $E\{Y(0) \mid S(0) = 1, S(1) = 0 \}
# Additive effect
# multiplicative effect

truth_df <- data.frame(
  pi = c(1, 1, 0, 0),
  er = c(1, 0, 1, 0),
  E_Y1 = c(
    pi_1_er_1$E_Y1__protected_or_doomed,
    pi_1_er_0$E_Y1__protected_or_doomed,
    pi_0_er_1$E_Y1__protected_or_doomed, 
    pi_0_er_0$E_Y1__protected_or_doomed
  ),
  E_Y0 = c(
    pi_1_er_1$E_Y0__protected_or_doomed,
    pi_1_er_0$E_Y0__protected_or_doomed,
    pi_0_er_1$E_Y0__protected_or_doomed, 
    pi_0_er_0$E_Y0__protected_or_doomed
  ),
  additive_effect = c(
    pi_1_er_1$effect_nat_inf,
    pi_1_er_0$effect_nat_inf,
    pi_0_er_1$effect_nat_inf, 
    pi_0_er_0$effect_nat_inf
  ),
  mult_effect = c(
    pi_1_er_1$effect_nat_inf_mult,
    pi_1_er_0$effect_nat_inf_mult,
    pi_0_er_1$effect_nat_inf_mult, 
    pi_0_er_0$effect_nat_inf_mult
  )
) %>%
  mutate(
    scenario = case_when(
      pi == 1 & er == 1 ~ "PI and ER satisfied",
      pi == 1 & er == 0 ~ "PI satisfied, ER violated",
      pi == 0 & er == 1 ~ "PI violated, ER satisfied",
      pi == 0 & er == 0 ~ "PI and ER violated"
    )
  ) %>%
  select(
    Scenario = scenario,
    E_Y1, E_Y0, additive_effect, mult_effect
  )

# LaTeX table
kable(
  truth_df,
  format = "latex",
  booktabs = TRUE,
  digits = 3,
  col.names = c(
    "",
    "$E\\{Y(1) \\mid S(0)=1, S(1)=0\\}$",
    "$E\\{Y(0) \\mid S(0)=1, S(1)=0\\}$",
    "Additive effect",
    "Multiplicative effect"
  ),
  escape = FALSE
) %>%
  kable_styling(latex_options = c("hold_position"))