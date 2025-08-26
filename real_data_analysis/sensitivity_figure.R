# ----------------------------------------------------------------------
# Make sensitivity figure for PROVIDE abx prescribing analysis
# ----------------------------------------------------------------------

library(ggplot2)

here::i_am("real_data_analysis/sensitivity_figure.R")

results <- readRDS(here::here("real_data_analysis/results/results_fewer_cov_same_lib.Rds"))

sens_res <- results$nat_inf$sens
sens_res_df <- data.frame(epsilon = sens_res$pt_est$epsilon,
                          psi_1 = sens_res$pt_est$psi_1_epsilon,
                          psi_0 = sens_res$pt_est$psi_0_aipw,
                          additive_effect = sens_res$pt_est$additive_effect,
                          lower = sens_res$pt_est$additive_effect - 1.96*sens_res$pt_est$additive_se,
                          upper = sens_res$pt_est$additive_effect + 1.96*sens_res$pt_est$additive_se)

# get the point estimate and CI at epsilon = 1
est_at1 <- sens_res_df$additive_effect[sens_res_df$epsilon == 1]
lower_at1 <- sens_res_df$lower[sens_res_df$epsilon == 1]
upper_at1 <- sens_res_df$upper[sens_res_df$epsilon == 1]

sens_plot <- ggplot(sens_res_df, aes(x = epsilon, y = additive_effect)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "#42B540FF") +
  geom_line(size = 1.2, color = "#42B540FF") +
  geom_hline(yintercept = est_at1, linetype = "dashed", color = "gray2") +
  annotate(
    "text",
    x = 1.02,  
    y = est_at1 + 0.15,
    hjust = 0,
    label = bquote(psi[1] - psi[0] == .(round(est_at1, 2)) ~ 
                     " (" * .(round(lower_at1, 2)) * ", " * .(round(upper_at1, 2)) * ")"),
    size = 5,
    color = "gray2"
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray2") +
  labs(
    x = expression(epsilon),
    y = expression(psi[1] - psi[0]~"(95 \u0025 CI)")
  ) +
  theme_minimal(base_size = 16)  +
  scale_x_log10(
    breaks = pretty(sens_res_df$epsilon),  # or use c(0.1, 0.2, 0.5, 1, 2, 5, 10) for custom ticks
    labels = scales::label_number()          # or scales::label_math() if you want expressions
  )

ggsave(here::here("real_data_analysis/results/sens_plot_provide.jpg"), sens_plot, width = 12, height = 8, dpi = 300)

# delta method stuff
t <- sens_res_df$psi_1 - sens_res_df$psi_0
t_star <- 0.5*log((1 + t) / (1 - t)) 

d__dpsi_1 <- 1 / ((sens_res_df$psi_1 - sens_res_df$psi_0)^2 - 1)
d__dpsi_0 <- -1 / ((sens_res_df$psi_1 - sens_res_df$psi_0)^2 - 1)

se_delta <- vector("numeric", length = length(sens_res$pt_est$cov_matrices))

for(i in 1:length(sens_res$pt_est$cov_matrices)){
  gradient <- c(d__dpsi_1[i], d__dpsi_0[i])
  
  se_delta[i] <- t(gradient) %*% sens_res$pt_est$cov_matrices[[i]] %*% gradient
}

sens_res_df <- cbind(sens_res_df, data.frame(t_star = t_star, 
                                             lower_t_star = (exp(2*(t_star - 1.96*se_delta)) - 1) / (exp(2*(t_star - 1.96*se_delta)) + 1),
                                             upper_t_star = (exp(2*(t_star + 1.96*se_delta)) - 1) / (exp(2*(t_star + 1.96*se_delta)) + 1)))

sens_plot_delta <- ggplot(sens_res_df, aes(x = epsilon, y = t)) +
  geom_ribbon(aes(ymin = lower_t_star, ymax = upper_t_star), alpha = 0.2, fill = "#42B540FF") +
  geom_line(size = 1.2, color = "#42B540FF") +
  geom_hline(yintercept = est_at1, linetype = "dashed", color = "gray2") +
  annotate(
    "text",
    x = 1.02,  
    y = est_at1 + 0.15,
    hjust = 0,
    label = bquote(psi[1] - psi[0] == .(round(est_at1, 2)) ~ 
                     " (" * .(round(lower_at1, 2)) * ", " * .(round(upper_at1, 2)) * ")"),
    size = 5,
    color = "gray2"
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray2") +
  labs(
    x = expression(epsilon),
    y = expression(psi[1] - psi[0]~"(95 \u0025 CI)")
  ) +
  theme_minimal(base_size = 16)  +
  scale_x_log10(
    breaks = pretty(sens_res_df$epsilon),  # or use c(0.1, 0.2, 0.5, 1, 2, 5, 10) for custom ticks
    labels = scales::label_number()          # or scales::label_math() if you want expressions
  )

ggsave(here::here("real_data_analysis/results/sens_plot_provide_fisher_z.jpg"), sens_plot_delta, width = 12, height = 8, dpi = 300)

