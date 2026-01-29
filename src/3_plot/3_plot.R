library(rstan); library(patchwork); library(tidyverse); library(lubridate); library(zoo); library(viridis)

orderly::orderly_dependency("1_data_cleaning", "latest()",
                            c("QLD_model_1_data.rds",
                              "QLD_model_2_data.rds",
                              "dose_df.rds",
                              "dose_df_QLD.rds",
                              "dose_df_sum_age_months.rds",
                              "QLD_model_1_data_in.rds",
                              "QLD_model_2_data_in.rds",
                              "QLD_births.rds"))

orderly::orderly_dependency("2_model_fit", "latest()",
                            c("fit_1.rds",
                              "fit_2.rds"))

QLD_births <- readRDS(file = "QLD_births.rds")

QLD_model_1_data <- readRDS(file = "QLD_model_1_data.rds")
QLD_model_2_data <- readRDS(file = "QLD_model_2_data.rds")
QLD_model_2_data_in <- readRDS(file = "QLD_model_2_data_in.rds")

dose_df <- readRDS(file = "dose_df.rds")
dose_df_QLD <- readRDS(file = "dose_df_QLD.rds")

fit_1 <- readRDS(file = "fit_1.rds")
fit_2 <- readRDS(file = "fit_2.rds")

##########################
##### raw data plots #####
##########################

### doses distributed

# model 1
doses_model_1 <- dose_df |> mutate(age_group = ifelse(age_in_months <= 8, "<=8 months", ifelse(age_in_months <= 24, "8-24 months", ">24 months"))) |>
  filter(year <= max(QLD_model_1_data$year)) |> group_by(age_group) |>
  summarise(QLD = sum(QLD_F) + sum(QLD_M) + sum(QLD_Unknown)) |>
  mutate(prop = QLD / sum(QLD))

sum(doses_model_1$QLD)

#######################
##### QLD model 1 #####
#######################

irr_pos_dis <- apply(rstan::extract(fit_1, "irr_pos_dis")$irr_pos_dis, 2, quantile, probs = c(0.025, 0.5, 0.975)) |>
  t() |>
  as.data.frame() |>
  mutate(week = QLD_model_1_data$week,
         year = QLD_model_1_data$year,
         treatment = QLD_model_1_data$treatment,
         week_year = QLD_model_1_data$week_year,
         week_cont = QLD_model_1_data$week_cont)

colnames(irr_pos_dis) <- c("l", "m", "u", "week", "year", "treatment", "week_year", "week_cont")

QLD_model_1_data <- QLD_model_1_data |>
  cbind(
    apply(rstan::extract(fit_1, "y_t_pos_dis")$y_t_pos_dis, 2, quantile, probs = c(0.025, 0.5, 0.975)) |>
      t() |>
      as.data.frame() |> rename(y_t_pos_l = 1, y_t_pos_m = 2, y_t_pos_u = 3)
  )

mu_gq_df <- apply(rstan::extract(fit_1, "mu_treat_gq")$mu_treat_gq, 2, quantile, probs = c(0.025, 0.5, 0.975)) |>
  t() |>
  as.data.frame() |>
  mutate(mu_c = 1:max(QLD_model_1_data$inc_greater_8m),
         treatment = 1) |>
  rbind(apply(rstan::extract(fit_1, "mu_untreat_gq")$mu_untreat_gq, 2, quantile, probs = c(0.025, 0.5, 0.975)) |>
          t() |>
          as.data.frame()|>
          mutate(mu_c = 1:max(QLD_model_1_data$inc_greater_8m),
                 treatment = 0)
  )

colnames(mu_gq_df)[1:3] <- c("l", "m", "u")

coef_df <- cbind(data.frame("coef" = c("alpha", "beta_1", "beta_2", "d")),
                 rbind(quantile(rstan::extract(fit_1, "intercept")[[1]], probs = c(0.025, 0.5, 0.975)),
                       quantile(rstan::extract(fit_1, "coef_inc_c")[[1]], probs = c(0.025, 0.5, 0.975)),
                       quantile(rstan::extract(fit_1, "coef_treat")[[1]], probs = c(0.025, 0.5, 0.975)),
                       quantile(rstan::extract(fit_1, "sigma_week_year")[[1]], probs = c(0.025, 0.5, 0.975))) |>
                   as.data.frame() |> rename("l" = 1, "m" = 2, "u" = 3)
                 )

fit_plots <-
  ggplot(data = QLD_model_1_data, aes(x = week_cont, y = inc_less_8m, fill = factor(treatment))) +
  geom_bar(stat = "identity", col = "grey30") +
  geom_pointrange(aes(x = week_cont, y = y_t_pos_m, ymin = y_t_pos_l, ymax = y_t_pos_u,
                      fill = factor(treatment)),
                  inherit.aes = FALSE,
                  shape = 21, size = 0.7, col = "grey30", alpha = 0.7) +
  theme_bw() + theme(text = element_text(size = 11), legend.position = c(0.8, 0.8), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
  scale_fill_manual(values = c("grey70", "skyblue"), name = "Nirsevimab\ndistribution",
                    labels = c("No", "Yes")) +
  xlab("Week-Year") + ylab("Incidence in those\n8 months old and under") +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  scale_x_continuous(labels = QLD_model_1_data$week_year[seq(1, 105, 10)], breaks = seq(1, 105, 10)) +

ggplot(data = QLD_model_1_data, aes(x = week_cont, y = inc_greater_8m)) +
  geom_bar(stat = "identity", col = "grey30", fill = "grey70") +
  theme_bw() + theme(text = element_text(size = 11)) +
  xlab("Week-Year") + ylab("Incidence in those\nolder than 8-months") +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  scale_x_continuous(labels = QLD_model_1_data$week_year[seq(1, 105, 10)], breaks = seq(1, 105, 10)) +

ggplot(data = QLD_model_1_data,
         aes(x = week_cont, y = inc_less_8m / inc_greater_8m, col = factor(treatment),
             fill = factor(treatment))) +
  geom_bar(stat = "identity", position = position_dodge(), col = "grey30") +
  theme_bw() + theme(text = element_text(size = 11), legend.position = c(0.8, 0.8), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
  xlab("Week-Year") + ylab("Incidence rate ratio of [0, 8]\nmonth olds relative to older ages") +
  geom_pointrange(data = irr_pos_dis,
                  aes(x = week_cont, y = m, ymin = l, ymax = u, fill = factor(treatment)),
                  inherit.aes = FALSE,
                  shape = 21, col = "grey30",
                  size = 0.7, alpha = 0.7) +
  scale_fill_manual(values = c("grey70", "skyblue"), name = "Nirsevimab\ndistribution",
                    labels = c("No", "Yes")) +
  scale_x_continuous(labels = QLD_model_1_data$week_year[seq(1, 105, 10)], breaks = seq(1, 105, 10)) +

  ggplot(data = QLD_model_1_data, aes(x = inc_less_8m, y = y_t_pos_m, ymin = y_t_pos_l,
                                      ymax = y_t_pos_u, col = factor(treatment), shape = factor(year))) +
  geom_pointrange(size = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme_bw()  + theme(text = element_text(size = 11), legend.position = c(0.3, 0.8), legend.direction = "horizontal", legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
  xlab("Sampled incidence in those\n8 months old and under") +
  ylab("Predicted incidence in\nthose 8 months old and under") +
  scale_colour_manual(values = c("grey70", "skyblue"), name = "Nirsevimab\ndistribution",
                    labels = c("No", "Yes")) +
  scale_shape_manual(values = c(15, 16), name = "year") +

  ggplot(data = mu_gq_df, aes(x = mu_c, y = m, ymin = l, ymax = u, fill = factor(treatment))) +
  geom_ribbon(alpha = 0.7) +
  ylab("Weekly incidence in\nthose 8 months old and under") +
  geom_point(inherit.aes = FALSE,
             data = QLD_model_1_data,
             aes(x = inc_greater_8m, y = inc_less_8m, fill = factor(treatment)),
             size = 3, shape = 21,
             col = "grey30") +
  theme_bw() + theme(text = element_text(size = 11), legend.position = c(0.3, 0.75), legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
  geom_line(aes(col = factor(treatment)), linewidth = 1) +
  xlab("Weekly incidence in those\nolder than 8 months") +
  scale_colour_manual(values = c("grey70", "skyblue"), name = "Nirsevimab\ndistribution",
                      labels = c("No", "Yes")) +
  scale_fill_manual(values = c("grey70", "skyblue"), name = "Nirsevimab\ndistribution",
                    labels = c("No", "Yes")) +

  ggplot(data = coef_df, aes(x = coef, y = m, ymin = l, ymax = u)) +
  geom_pointrange(alpha = 0.5) +
  theme_bw() + theme(text = element_text(size = 11)) +
  ylab("Fitted parameter value") +
  xlab("Parameter") +
  scale_x_discrete(labels = c("alpha" = parse(text = "alpha"),
                              "beta_1" = parse(text = "beta[1]"),
                              "beta_2" = parse(text = "beta[2]"),
                              "d" = parse(text = "d"))) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(limits = c(-0.75, 1.75), breaks = seq(-0.5, 1.5, 0.5)) +

  plot_annotation(tag_levels = "a") +
  plot_layout(nrow = 2, ncol = 3)

ggsave(plot = fit_plots,
       filename = "model_1_plots.pdf",
       height = 20,
       width = 40,
       units = "cm",
       device = "pdf")

round(quantile(exp(rstan::extract(fit_1, "coef_inc_c")[[1]]), probs = c(0.025, 0.5, 0.975)), digits = 3)

round(quantile(exp(rstan::extract(fit_1, "coef_treat")[[1]]), probs = c(0.025, 0.5, 0.975)), digits = 2)

ggsave(
  traceplot(fit_1, pars = c("intercept", "coef_inc_c", "coef_treat", "sigma_week_year")),
  file = "model_1_traceplot.pdf",
  height = 10,
  width = 20,
  units = "cm",
  device = "pdf"
)

#######################
##### QLD model 2 #####
#######################

ggsave(
  traceplot(fit_2, pars = c("intercept", "coef_year_raw", "coef_doses", "coef_age_l", "coef_age_s", "amplitude", "phase_shift_cos_months", "sigma_cohort_births", "sigma_cohort_month_years")),
  file = "model_2_traceplot.pdf",
  height = 20,
  width = 30,
  units = "cm",
  device = "pdf"
)

pairs(fit_2, pars = c("intercept", "coef_year_raw", "coef_doses", "coef_age_l", "coef_age_s", "amplitude", "phase_shift_cos_months", "sigma_cohort_births", "sigma_cohort_month_years"))

quantile(exp(rstan::extract(fit_2, "coef_doses")[[1]]), probs = c(0.025, 0.5, 0.975))

QLD_model_2_data <- cbind(QLD_model_2_data,
                          apply(rstan::extract(fit_2, "y_pos_dis")[[1]], 2, quantile, probs = c(0.025, 0.5, 0.975)) |>
                            t() |>
                            as.data.frame() |>
                            rename("low" = 1, "med" = 2, "up" = 3))

actual_fitted_plot <- ggplot(data = QLD_model_2_data,
       aes(x = inc, y = med, col = factor(cohort_birth_months))) +
  geom_pointrange(alpha = 0.25, aes(, ymin = low, ymax = up)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.75) +
  xlab("Sampled monthly incidence per cohort") +
  ylab("Predicted monthly incidence per cohort") +
  theme_bw() + theme(text = element_text(size = 11), legend.position = "none") +
  scale_colour_manual(values = cols, name = "Cohort")

coef_df_model_2 <- cbind(data.frame("coef" = c("alpha", "beta_1", "beta_2", "beta_3", "beta_4", "beta_5", "beta_6", "sigma_1", "sigma_2")),
                         rbind(quantile(rstan::extract(fit_2, "intercept")[[1]], probs = c(0.025, 0.5, 0.975)),
                               quantile(as.vector(rstan::extract(fit_2, "coef_year_raw")[[1]]), probs = c(0.025, 0.5, 0.975)),
                               quantile(rstan::extract(fit_2, "coef_month_c")[[1]], probs = c(0.025, 0.5, 0.975)),
                               quantile(rstan::extract(fit_2, "coef_month_s")[[1]], probs = c(0.025, 0.5, 0.975)),
                               quantile(rstan::extract(fit_2, "coef_age_l")[[1]], probs = c(0.025, 0.5, 0.975)),
                               quantile(rstan::extract(fit_2, "coef_age_s")[[1]], probs = c(0.025, 0.5, 0.975)),
                               quantile(rstan::extract(fit_2, "coef_doses")[[1]], probs = c(0.025, 0.5, 0.975)),
                               quantile(rstan::extract(fit_2, "sigma_cohort_births")[[1]], probs = c(0.025, 0.5, 0.975)),
                               quantile(rstan::extract(fit_2, "sigma_cohort_month_years")[[1]], probs = c(0.025, 0.5, 0.975))) |>
                           as.data.frame() |> rename("l" = 1, "m" = 2, "u" = 3)
                         )

coef_plot_model_2 <- ggplot(data = coef_df_model_2, aes(x = coef, y = m, ymin = l, ymax = u)) +
  geom_pointrange(alpha = 0.5) +
  theme_bw() + theme(text = element_text(size = 11)) +
  ylab("Fitted parameter value") +
  xlab("Parameter") +
  scale_x_discrete(labels = c("alpha" = parse(text = "alpha"),
                              "beta_1" = parse(text = "beta[1]"),
                              "beta_2" = parse(text = "beta[2]"),
                              "beta_3" = parse(text = "beta[3]"),
                              "beta_4" = parse(text = "beta[4]"),
                              "beta_5" = parse(text = "beta[5]"),
                              "beta_6" = parse(text = "beta[6]"),
                              "sigma_1" = parse(text = "sigma[1]"),
                              "sigma_2" = parse(text = "sigma[2]")
                              )) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_y_continuous(limits = c(-9, 1), breaks = seq(-9, 1, 1))

### variables with uncertainty

cohorts <- sort(unique(QLD_model_2_data$cohort_birth_months))
cols <- viridis(length(cohorts), option = "viridis")
cols <- setNames(cols, cohorts)

dose_age_QLD <- dose_df_QLD |> mutate(age = interval(as.Date(cohort_birth_months), start_month) %/% months(1)) |>
  group_by(age, cohort_birth_months) |>
  summarise(doses = sum(QLD)) |>
  subset(doses > 0)

sum(dose_age_QLD$doses) == sum(dose_df_QLD$QLD)

dose_plot <-
  ggplot(dose_df_QLD, aes(x = start_month, y = QLD, fill = factor(cohort_birth_months))) +
  geom_bar(stat = "identity") +
  theme_bw() + theme(text = element_text(size = 11),
                     legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
  ylab("Nirsevimab doses per month") +
  xlab("Month") +
  scale_y_continuous(limits = c(0, 5500), breaks = seq(0, 5500, 1000)) +
  scale_fill_manual(name = "Cohort", values = cols) +
  guides(fill = guide_legend(
    keywidth = unit(0.425, "cm"),
    keyheight = unit(0.425, "cm"),
    override.aes = list(size = 0.725) # Makes the dots inside the legend smaller
  ))

dose_age_plot <- ggplot(data = dose_age_QLD,
                        aes(x = age, y = doses, fill = factor(cohort_birth_months))
                        ) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Cohort", values = cols) +
  scale_colour_manual(name = "Cohort", values = cols) +
  theme_bw() + theme(text = element_text(size = 11), legend.position = "none") +
  xlab("Age in months at time of immunisation") +
  ylab("Number of Nirsevimab\ndoses administered") +
  scale_y_sqrt()

# births
births_plot <- ggplot(data = QLD_births, aes(x = cohort_birth_months, y = births,
                                             col = factor(cohort_birth_months))) +
  geom_point() +
  geom_line(data = QLD_births[0, ]) +
  ylab("Registered births") + xlab("Cohort") +
  theme_bw() + theme(text = element_text(size = 11), legend.position = "none") +
  scale_y_continuous(limits = c(0, 6000)) +
  scale_colour_manual(name = "Cohort", values = cols)

ggsave(
  file = "uncertain_variables_plot.pdf",
  dose_plot + dose_age_plot + births_plot + guide_area() + plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(tag_levels = c("a")),
  height = 17.5,
  width = 27.5,
  units = "cm",
  device = "pdf"
)

### covariates

# estimating the coverage for each data point
dose_check <- rep(NA, nrow(QLD_model_2_data))

dose_check[QLD_model_2_data_in$doses_on_inds] <- as.vector(QLD_model_2_data_in$doses_mat %*% QLD_model_2_data_in$dose_data)
dose_check[is.na(dose_check)] <- 0

QLD_model_2_data$cov <- dose_check / QLD_births$births[QLD_model_2_data$ind_cohort_births]

n_iter <- length(rstan::extract(fit_2, "intercept")[[1]])

intercept <- rep(1, nrow(QLD_model_2_data)) %o% rstan::extract(fit_2, "intercept")[[1]]
year <- (QLD_model_2_data$ind_years - 1) %o% as.vector(rstan::extract(fit_2, "coef_year_raw")[[1]])
log_offset_months <- log(as.vector(QLD_model_2_data$offset_month) %o% rep(1, n_iter))
log_offset_births <- log(QLD_births$births[QLD_model_2_data$ind_cohort_births] %o% rep(1, n_iter))

months <- (sin(QLD_model_2_data$month * 2 * pi / 12) %o% rstan::extract(fit_2, "coef_month_s")[[1]]) +
  (cos(QLD_model_2_data$month * 2 * pi / 12) %o% rstan::extract(fit_2, "coef_month_c")[[1]])

age <- (QLD_model_2_data$age_rsv_months %o% rstan::extract(fit_2, "coef_age_l")[[1]]) +
  (QLD_model_2_data$age_rsv_months^2 %o% rstan::extract(fit_2, "coef_age_s")[[1]])

dose <- QLD_model_2_data$cov %o% rstan::extract(fit_2, "coef_doses")[[1]]

re_cohort <- t(rstan::extract(fit_2, "coef_cohort_births")[[1]][, QLD_model_2_data$ind_cohort_births])

re_cohort_month_years <- t(rstan::extract(fit_2, "coef_cohort_month_years")[[1]][,QLD_model_2_data$ind_cohort_month_years])

##### average effect plots

# age
e_inc_m_age <- intercept + log_offset_months + log_offset_births + months + dose + re_cohort + re_cohort_month_years

ages_in <- seq(0, max(QLD_model_2_data$age_rsv_months), 0.1)

n_ages <- length(ages_in)

ages_coef <- (ages_in %o% rstan::extract(fit_2, "coef_age_l")[[1]]) + (ages_in %o% rstan::extract(fit_2, "coef_age_s")[[1]])

avg_eff_age <- sapply(1:n_ages,
                      function(i){
                        quantile(colMeans(exp(e_inc_m_age + t(replicate(nrow(QLD_model_2_data), ages_coef[i,])))), probs = c(0.025, 0.5, 0.975))
                        }
                      )

avg_eff_age <- as.data.frame(t(avg_eff_age)) |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate(age_rsv_months = ages_in)

age_data_plot <-
  ggplot(data = QLD_model_2_data,
       aes(x = age_rsv_months, y = inc, group = cohort_birth_months, col = factor(cohort_birth_months))) +
  geom_ribbon(data = avg_eff_age, aes(x = age_rsv_months, ymin = l, ymax = u), inherit.aes = FALSE, alpha = 0.25) +
  geom_point(alpha = 0.75) +
  geom_line(data = avg_eff_age, aes(x = age_rsv_months, y = m), inherit.aes = FALSE, size = 1) +
  theme_bw() + theme(text = element_text(size = 11)) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Age in months") +
  scale_colour_manual(name = "Cohort", values = cols) +
  theme(legend.position = "none")

# seasonality
e_inc_m_season <- intercept + log_offset_months + log_offset_births + age + dose + re_cohort + re_cohort_month_years

months_in <- seq(1, 12, 0.1)

n_months <- length(months_in)

months_coef <- (sin(months_in * 2 * pi / 12) %o% rstan::extract(fit_2, "coef_month_s")[[1]]) +
  (cos(months_in * 2 * pi / 12) %o% rstan::extract(fit_2, "coef_month_c")[[1]])

avg_eff_season <- sapply(1:n_months,
                      function(i){
                        quantile(colMeans(exp(e_inc_m_season + t(replicate(nrow(QLD_model_2_data), months_coef[i,])))), probs = c(0.025, 0.5, 0.975))
                      }
)

avg_eff_season <- as.data.frame(t(avg_eff_season)) |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate(months = months_in)

season_data_plot <- ggplot(data = QLD_model_2_data,
                           aes(x = month, y = inc, col = factor(cohort_birth_months))) +
  geom_ribbon(data = avg_eff_season, aes(x = months, ymin = l, ymax = u), alpha = 0.25, inherit.aes = FALSE) +
  geom_point(alpha = 0.75) +
  geom_line(data = avg_eff_season, aes(x = months, y = m), inherit.aes = FALSE, size = 1) +
  theme_bw() + theme(text = element_text(size = 11)) +
  scale_colour_manual(name = "Cohort", values = cols) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Month") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(1, 12, 1), limits = c(1, 12))

# coverage

e_inc_m_cov <- intercept + log_offset_months + log_offset_births + age + months + re_cohort + re_cohort_month_years

cov_in <- seq(0, ceiling(max(QLD_model_2_data$cov)/0.1)*0.1, 0.025)

n_cov <- length(cov_in)

cov_coef <- cov_in %o% rstan::extract(fit_2, "coef_doses")[[1]]

avg_eff_cov <- sapply(1:n_cov,
                      function(i){
                        quantile(colMeans(exp(e_inc_m_cov + t(replicate(nrow(QLD_model_2_data), cov_coef[i,])))), probs = c(0.025, 0.5, 0.975))
                        }
)

avg_eff_cov <- as.data.frame(t(avg_eff_cov)) |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate(cov = cov_in)

cov_data_plot <- ggplot(data = QLD_model_2_data,
                           aes(x = cov, y = inc, col = factor(cohort_birth_months))) +
  geom_ribbon(data = avg_eff_cov, aes(x = cov, ymin = l, ymax = u), alpha = 0.25, inherit.aes = FALSE) +
  geom_point(alpha = 0.75) +
  geom_line(data = avg_eff_cov, aes(x = cov, y = m), inherit.aes = FALSE, size = 1) +
  theme_bw() + theme(text = element_text(size = 11)) +
  scale_colour_manual(name = "Cohort", values = cols) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Proportion of cohort immunised with Nirsevimab in previous 6-months") +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent) +
  guides(colour = guide_legend(
    keywidth = unit(0.4, "cm"),
    keyheight = unit(0.4, "cm"),
    override.aes = list(size = 0.7) # Makes the dots inside the legend smaller
  ))

### marginal effects

##### coefficient plots

ggsave(
  coef_plot_model_2 + actual_fitted_plot + age_data_plot +
  season_data_plot + cov_data_plot + guide_area() +
  plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = c("a")),
  file = "model_2_plots.pdf",
  device = "pdf",
  height = 20,
  width = 35,
  units = "cm"
)

