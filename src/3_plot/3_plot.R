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
QLD_model_1_data_in <- readRDS(file = "QLD_model_1_data_in.rds")
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

# mu_gq_df <- apply(rstan::extract(fit_1, "mu_treat_gq")$mu_treat_gq, 2, quantile, probs = c(0.025, 0.5, 0.975)) |>
#   t() |>
#   as.data.frame() |>
#   mutate(mu_c = 1:max(QLD_model_1_data$inc_greater_8m),
#          treatment = 1) |>
#   rbind(apply(rstan::extract(fit_1, "mu_untreat_gq")$mu_untreat_gq, 2, quantile, probs = c(0.025, 0.5, 0.975)) |>
#           t() |>
#           as.data.frame()|>
#           mutate(mu_c = 1:max(QLD_model_1_data$inc_greater_8m),
#                  treatment = 0)
#   )
#
# colnames(mu_gq_df)[1:3] <- c("l", "m", "u")

intercepts <- rstan::extract(fit_1, "intercept")[[1]]
re_week_year <- rstan::extract(fit_1, "coef_week_year")[[1]]
coef_inc_c <- rstan::extract(fit_1, "coef_inc_c")[[1]]
coef_treat <- rstan::extract(fit_1, "coef_treat")[[1]]

ind_week_year = sort(unique(QLD_model_1_data$week_cont))

mu_gq_df <- expand.grid(mu_c = seq(1, max(QLD_model_1_data$inc_greater_8m), 3),
                        treatment = c(0, 1))

pred <- lapply(1:nrow(mu_gq_df),
       function(i, mu_gq_df){
         out <- quantile(rowMeans(exp(t(log(7) + rep(1, ncol(re_week_year)) %o% (coef_inc_c * mu_gq_df[i, "mu_c"]/7 + coef_treat * mu_gq_df[i, "treatment"] + intercepts)) + re_week_year)), probs = c(0.025, 0.5, 0.975))
         return(data.frame("l" = out[1], "m" = out[2], "u" = out[3]))
         },
       mu_gq_df = mu_gq_df
       ) |> bind_rows() |> as.data.frame()

coef_df <- cbind(data.frame("coef" = c("alpha", "beta_1", "beta_2", "d")),
                 rbind(quantile(intercepts, probs = c(0.025, 0.5, 0.975)),
                       quantile(coef_inc_c, probs = c(0.025, 0.5, 0.975)),
                       quantile(coef_treat, probs = c(0.025, 0.5, 0.975)),
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
  scale_x_continuous(labels = arrange(subset(QLD_model_1_data, week_cont %in% seq(1, 105, 10)), week_cont) |> select(week_year) |> as.vector() |> unlist(), breaks = seq(1, 105, 10)) +

ggplot(data = QLD_model_1_data, aes(x = week_cont, y = inc_greater_8m)) +
  geom_bar(stat = "identity", col = "grey30", fill = "grey70") +
  theme_bw() + theme(text = element_text(size = 11)) +
  xlab("Week-Year") + ylab("Incidence in those\nolder than 8-months") +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  scale_x_continuous(labels = arrange(subset(QLD_model_1_data, week_cont %in% seq(1, 105, 10)), week_cont) |> select(week_year) |> as.vector() |> unlist(), breaks = seq(1, 105, 10)) +

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
  scale_x_continuous(labels = arrange(subset(QLD_model_1_data, week_cont %in% seq(1, 105, 10)), week_cont) |> select(week_year) |> as.vector() |> unlist(), breaks = seq(1, 105, 10)) +

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

  ggplot(data = cbind(mu_gq_df, pred), aes(x = mu_c, y = m, ymin = l, ymax = u, fill = factor(treatment))) +
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
       height = 25,
       width = 55,
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

round((1 - quantile(exp(rstan::extract(fit_2, "coef_doses")[[1]]), probs = c(0.025, 0.5, 0.975))) * 100, digits = 0)
round(quantile(exp(rstan::extract(fit_2, "coef_doses")[[1]]), probs = c(0.025, 0.5, 0.975)), digits = 2)

round(quantile(rstan::extract(fit_2, "phase_shift_cos_months")[[1]], probs = c(0.025, 0.5, 0.975)), digits = 2)
round(quantile(rstan::extract(fit_2, "amplitude")[[1]], probs = c(0.025, 0.5, 0.975)), digits = 2)
round(quantile(rstan::extract(fit_2, "coef_age_s")[[1]]/rstan::extract(fit_2, "coef_age_l")[[1]], probs = c(0.025, 0.5, 0.975)), digits = 2)

ggsave(
  traceplot(fit_2, pars = c("intercept", "coef_year_raw", "coef_doses", "coef_age_l", "coef_age_s", "amplitude", "phase_shift_cos_months", "sigma_cohort_births", "sigma_cohort_month_years")),
  file = "model_2_traceplot.pdf",
  height = 20,
  width = 30,
  units = "cm",
  device = "pdf"
)

pairs(fit_2, pars = c("intercept", "coef_year_raw", "coef_doses", "coef_age_l", "coef_age_s", "amplitude", "phase_shift_cos_months", "sigma_cohort_births", "sigma_cohort_month_years"))

QLD_model_2_data <- cbind(QLD_model_2_data,
                          apply(rstan::extract(fit_2, "y_pos_dis")[[1]], 2, quantile, probs = c(0.025, 0.5, 0.975)) |>
                            t() |>
                            as.data.frame() |>
                            rename("low" = 1, "med" = 2, "up" = 3))

cohorts <- sort(unique(QLD_model_2_data$cohort_birth_months))
cols <- viridis(length(cohorts), option = "viridis")
cols <- setNames(cols, cohorts)

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

dose_age_QLD <- dose_df_QLD |> mutate(age = interval(as.Date(cohort_birth_months), start_month) %/% months(1) + 0.01) |>
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
  ylab("Number of Nirsevimab\ndoses administered")

# births
births_plot <- ggplot(data = QLD_births, aes(x = cohort_birth_months, y = births,
                                             col = factor(cohort_birth_months))) +
  geom_point() +
  geom_line(data = QLD_births[0, ]) +
  ylab("Births") + xlab("Cohort") +
  theme_bw() + theme(text = element_text(size = 11), legend.position = "none") +
  scale_y_continuous(limits = c(0, 6000)) +
  scale_colour_manual(name = "Cohort", values = cols)

suppressWarnings(
  ggsave(
  file = "uncertain_variables_plot.pdf",
  dose_plot + dose_age_plot + births_plot + guide_area() + plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(tag_levels = c("a")),
  height = 17.5,
  width = 27.5,
  units = "cm",
  device = "pdf"
  )
)

######################
##### covariates #####
######################

# estimating the coverage for each data point
dose_check <- rep(NA, nrow(QLD_model_2_data))

dose_check[QLD_model_2_data_in$doses_on_inds] <- as.vector(QLD_model_2_data_in$doses_mat %*% QLD_model_2_data_in$dose_data)
dose_check[is.na(dose_check)] <- 0

QLD_model_2_data$cov <- dose_check / QLD_births$births[QLD_model_2_data$ind_cohort_births]

n_iter <- length(rstan::extract(fit_2, "intercept")[[1]])

# values to predict from GLM

glm_log_offset_months <- log(QLD_model_2_data$offset_month %o% rep(1, n_iter))
glm_log_offset_births <- log(QLD_births$births[QLD_model_2_data$ind_cohort_births] %o% rep(1, n_iter))

glm_intercept <- rep(1, nrow(QLD_model_2_data)) %o% rstan::extract(fit_2, "intercept")[[1]]

year_mat <- cbind(rep(0, n_iter), rstan::extract(fit_2, "coef_year_raw")[[1]])

glm_year <- t(year_mat[,QLD_model_2_data$ind_years])

glm_months <- (sin(QLD_model_2_data$month * 2 * pi / 12) %o% rstan::extract(fit_2, "coef_month_s")[[1]]) +
  (cos(QLD_model_2_data$month * 2 * pi / 12) %o% rstan::extract(fit_2, "coef_month_c")[[1]])

glm_age <- (-QLD_model_2_data$age_rsv_months %o% rstan::extract(fit_2, "coef_age_l")[[1]]) +
  (log(QLD_model_2_data$age_rsv_months) %o% rstan::extract(fit_2, "coef_age_s")[[1]])

glm_dose <- QLD_model_2_data$cov %o% rstan::extract(fit_2, "coef_doses")[[1]]

glm_re_cohort <- t(rstan::extract(fit_2, "coef_cohort_births")[[1]][, QLD_model_2_data$ind_cohort_births])

glm_re_cohort_month_years <- t(rstan::extract(fit_2, "coef_cohort_month_years")[[1]][,QLD_model_2_data$ind_cohort_month_years])

glm_overall_effect <- glm_intercept + glm_year + glm_log_offset_months + glm_log_offset_births + glm_age + glm_months + glm_dose + glm_re_cohort + glm_re_cohort_month_years

##### average effect plots

# age
e_inc_m_age <- glm_overall_effect - glm_age

ages_in <- seq(0.01, max(QLD_model_2_data$age_rsv_months), 0.1)

n_ages <- length(ages_in)

ages_coef <- (-ages_in %o% rstan::extract(fit_2, "coef_age_l")[[1]]) + (log(ages_in) %o% rstan::extract(fit_2, "coef_age_s")[[1]])

# ages_coef <- (ages_in %o% rstan::extract(fit_2, "coef_age_l")[[1]]) + (ages_in^2 %o% rstan::extract(fit_2, "coef_age_s")[[1]])

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
  geom_line(data = avg_eff_age, aes(x = age_rsv_months, y = m), inherit.aes = FALSE, linewidth = 1) +
  theme_bw() + theme(text = element_text(size = 11)) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Age in months") +
  scale_colour_manual(name = "Cohort", values = cols) +
  theme(legend.position = "none")

# seasonality
e_inc_m_season <- glm_overall_effect - glm_months

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
  geom_line(data = avg_eff_season, aes(x = months, y = m), inherit.aes = FALSE, linewidth = 1) +
  theme_bw() + theme(text = element_text(size = 11)) +
  scale_colour_manual(name = "Cohort", values = cols) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Month") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(1, 12, 1), limits = c(1, 12))

# coverage

e_inc_m_cov <- glm_overall_effect - glm_dose

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
  geom_line(data = avg_eff_cov, aes(x = cov, y = m), inherit.aes = FALSE, linewidth = 1) +
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

ggsave(
  coef_plot_model_2 + actual_fitted_plot + age_data_plot +
  season_data_plot + cov_data_plot + guide_area() + plot_layout(ncol = 3, guides = "collect") +
  plot_annotation(tag_levels = c("a")),
  file = "model_2_plots.pdf",
  device = "pdf",
  height = 20,
  width = 37.5,
  units = "cm"
)

##### counterfactual plots

QLD_model_2_data <- QLD_model_2_data |> mutate(cov_cf_on = ifelse(age_rsv_months <= 6 & month >= 4 & year >= 2024, 1, 0),
                                               cov_cf_off = 0) |>
  rowwise() |> mutate(cov_cf_on = max(cov_cf_on, cov))

cf_obs <- exp(glm_overall_effect)
cf_dose_on <- exp(glm_overall_effect - glm_dose + QLD_model_2_data$cov_cf_on %o% rstan::extract(fit_2, "coef_doses")[[1]])
cf_dose_off <- exp(glm_overall_effect - glm_dose + QLD_model_2_data$cov_cf_off %o% rstan::extract(fit_2, "coef_doses")[[1]])


QLD_year_month_all <- QLD_model_2_data$rsv_start_month

calc_cf <- function(cf, QLD_year_month){
  cf_summed_matrix <- rowsum(cf, group = QLD_year_month, reorder = FALSE)
  cf_cs_matrix <- apply(cf_summed_matrix, 2, cumsum)
  cf_cs <- apply(cf_cs_matrix, 1, quantile, probs = c(0.025, 0.5, 0.975)) |> t() |> as.data.frame() |>
  rename("l" = 1, "m" = 2, "u" = 3) |>
  mutate("rsv_start_month" = sort(unique(QLD_year_month)))
}

cf_cs <- rbind(calc_cf(cf_dose_on, QLD_year_month_all) |> mutate(model = "100% Nirsevimab coverage"),
               calc_cf(cf_dose_off, QLD_year_month_all) |> mutate(model = "0% Nirsevimab coverage"))

obs_cs <- QLD_model_2_data |> group_by(rsv_start_month) |> summarise(inc = sum(inc)) |>
  arrange(rsv_start_month) |>
  ungroup() |>
  mutate(c_inc = cumsum(inc))

cf_cs_obs <- calc_cf(cf_obs, QLD_year_month_all) |> mutate(model = "predicted")

cf_all_plot <- ggplot(data = cf_cs, aes(x = rsv_start_month, y = m, ymin = l, ymax = u, fill = model)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(aes(col = model)) +
  geom_ribbon(data = cf_cs_obs, aes(x = rsv_start_month, ymin = l, ymax = u), inherit.aes = FALSE, alpha = 0.25) +
  geom_line(data = cf_cs_obs, aes(x = rsv_start_month, y = m), inherit.aes = FALSE) +
  theme_bw() +
  geom_point(data = obs_cs, aes(x = rsv_start_month, y = c_inc), inherit.aes = FALSE) +
  ylab("Cumulative notifications in 0 - 60 month olds") +
  xlab("Month-Year") +
  scale_colour_manual(name = "Counterfactual", values = c("#CC79A7", "#56B4E9")) +
  scale_fill_manual(name = "Counterfactual", values = c("#CC79A7", "#56B4E9"))

# under 6-month olds only

age_index <- which(QLD_model_2_data$age_rsv_months <= 6)

cf_cs_6_month <- rbind(calc_cf(cf_dose_on[age_index,], QLD_year_month_all[age_index]) |> mutate(model = "100% Nirsevimab coverage"),
                       calc_cf(cf_dose_off[age_index,], QLD_year_month_all[age_index]) |> mutate(model = "0% Nirsevimab coverage"))

obs_cs_6_month <- QLD_model_2_data[age_index,] |> group_by(rsv_start_month) |> summarise(inc = sum(inc)) |>
  arrange(rsv_start_month) |>
  ungroup() |>
  mutate(c_inc = cumsum(inc))

cf_cs_obs_6_month <- calc_cf(cf_obs[age_index,], QLD_year_month_all[age_index]) |> mutate(model = "predicted")

cf_6m_plot <- ggplot(data = cf_cs_6_month, aes(x = rsv_start_month, y = m, ymin = l, ymax = u, fill = model)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(aes(col = model)) +
  geom_ribbon(data = cf_cs_obs_6_month, aes(x = rsv_start_month, ymin = l, ymax = u), inherit.aes = FALSE, alpha = 0.25) +
  geom_line(data = cf_cs_obs_6_month, aes(x = rsv_start_month, y = m), inherit.aes = FALSE) +
  theme_bw() +
  geom_point(data = obs_cs_6_month, aes(x = rsv_start_month, y = c_inc), inherit.aes = FALSE) +
  ylab("Cumulative notifications in 0 - 6 month olds") +
  xlab("Month-Year") +
  scale_colour_manual(name = "Counterfactual", values = c("#CC79A7", "#56B4E9")) +
  scale_fill_manual(name = "Counterfactual", values = c("#CC79A7", "#56B4E9"))

cf_all_plot + cf_6m_plot + plot_layout(guides = "collect") + plot_annotation(tag_levels = c("a"))

write.csv(rbind(cf_cs_6_month, cf_cs_obs_6_month, obs_cs_6_month |> rename("m" = c_inc) |> mutate(l = NA, u = NA, model = "observed") |> select(l, m, u, rsv_start_month, model)),
 file = "counterfactual_0_6_months.csv")

