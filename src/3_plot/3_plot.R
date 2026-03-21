library(cmdstanr); library(patchwork);
library(tidyverse); library(lubridate); library(zoo); library(viridis)
library(parallel); library(bayesplot)

theme_set(theme_bw() +
            theme(text = element_text(size = 11),
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 10)))

orderly::orderly_dependency("1_data_cleaning", "latest()",
                            c("QLD_model_1_data.rds",
                              "model_2_data_all.rds",
                              "model_2_data_ACT.rds",
                              "model_2_data_NSW.rds",
                              "model_2_data_QLD.rds",
                              "dose_df.rds",
                              "dose_df_sum_age_months_long.rds",
                              "QLD_model_1_data_in.rds",
                              "model_3_data_in.rds",
                              "model_2_data_ACT_in.rds",
                              "model_2_data_NSW_in.rds",
                              "model_2_data_QLD_in.rds",
                              "births_all.rds",
                              "ACT_births.rds",
                              "NSW_births.rds",
                              "QLD_births.rds"))

orderly::orderly_dependency("2_model_fit", "latest()",
                            c("fit_1.rds",
                              "fit_2_ACT.rds",
                              "fit_2_NSW.rds",
                              "fit_2_QLD.rds",
                              "fit_3.rds"))

births_all <- readRDS(file = "births_all.rds")
ACT_births <- readRDS(file = "ACT_births.rds")
NSW_births <- readRDS(file = "NSW_births.rds")
QLD_births <- readRDS(file = "QLD_births.rds")

QLD_model_1_data <- readRDS(file = "QLD_model_1_data.rds")
QLD_model_1_data_in <- readRDS(file = "QLD_model_1_data_in.rds")
model_2_data_ACT <- readRDS(file = "model_2_data_ACT.rds")
model_2_data_ACT_in <- readRDS(file = "model_2_data_ACT_in.rds")
model_2_data_NSW <- readRDS(file = "model_2_data_NSW.rds")
model_2_data_NSW_in <- readRDS(file = "model_2_data_NSW_in.rds")
model_2_data_QLD <- readRDS(file = "model_2_data_QLD.rds")
model_2_data_QLD_in <- readRDS(file = "model_2_data_QLD_in.rds")
model_3_data <- readRDS(file = "model_2_data_all.rds")
model_3_data_in <- readRDS(file = "model_3_data_in.rds")

dose_df <- readRDS(file = "dose_df.rds")
dose_df_sum_age_months_long <- readRDS(file = "dose_df_sum_age_months_long.rds")

fit_1 <- readRDS(file = "fit_1.rds")
fit_2_ACT <- readRDS(file = "fit_2_ACT.rds")
fit_2_NSW <- readRDS(file = "fit_2_NSW.rds")
fit_2_QLD <- readRDS(file = "fit_2_QLD.rds")
fit_3 <- readRDS(file = "fit_3.rds")

cohorts <- sort(unique(model_3_data$cohort_birth_months))
cols <- viridis(length(cohorts), option = "viridis")
cols <- setNames(cols, cohorts)
col_states <- c("ACT" = "#E69F00", "NSW" = "#CC79A7", "QLD" = "#56B4E9")
model_cols <- c("Pooled" = "#0072B2", "Independent" = "#009E73")

##########################
##### raw data plots #####
##########################

##### model 1
### doses distributed - all raw data

dose_sum_age_months_long <- dose_df_sum_age_months_long |>
  filter(year == 2024) |>
  mutate(cohort_birth_start_month = as.Date(cohort_birth_months),
         age_months = interval(cohort_birth_start_month, start_month) %/% months(1))

doses_model_1 <- dose_sum_age_months_long |>
  filter(state == "QLD") |> mutate(age_group = ifelse(age_months <= 8, "<=8 months", ifelse(age_months <= 24, "8-24 months", ">24 months"))) |>
  filter(year <= max(QLD_model_1_data$year)) |> group_by(age_group) |>
  summarise(QLD = sum(doses)) |>
  mutate(prop = QLD / sum(QLD))

sum(doses_model_1$QLD)

##### models 2 and 3
### variables with uncertainty

dose_df$cohort_birth_months <- factor(dose_df$cohort_birth_months, levels = cohorts)

# dose data actually used in the model (limited to cohorts modelled)
# maximum ages
dose_df |> subset(doses > 0) |> group_by(state) |> filter(age_months == max(age_months))

dose_df |> group_by(state) |> summarise(doses = sum(doses))

# age distribution
dose_df <- subset(dose_df, doses > 0)

dose_age <- dose_df |>
  summarise(doses = sum(doses), .by = c(age_months, state, cohort_birth_months)) |>
  subset(doses > 0)

sum(dose_age$doses) == sum(dose_df$doses)

sort(unique(dose_age$cohort_birth_months)) %in% sort(unique(dose_df$cohort_birth_months))

dose_plot <-
  ggplot(dose_df, aes(x = start_month, y = doses, fill = cohort_birth_months, group = state)) +
  geom_bar(stat = "identity")  +
  geom_bar(stat = "identity",
           data = dose_df |> summarise(doses = sum(doses), .by = c(start_month, state)),
           aes(x = start_month, y = doses), fill = NA, col = "grey30", inherit.aes = FALSE) +
  ylab("Nirsevimab doses per month") +
  xlab("Month") +
  scale_fill_manual(name = "Cohort", values = cols) +
  guides(fill = guide_legend(
    keywidth = unit(0.425, "cm"),
    keyheight = unit(0.425, "cm"),
    override.aes = list(size = 0.725))) +
  facet_wrap(~state, scales = "free_y") +
  scale_x_date(date_labels = "%b %Y",
               breaks = seq(min(dose_df$start_month), max(dose_df$start_month) - 30, by = "3 months"))

dose_age_plot <- ggplot(data = dose_age,
                        aes(x = age_months, y = doses,
                            fill = cohort_birth_months,
                            group = state)) +
  geom_bar(stat = "identity") +
  geom_bar(stat = "identity",
           data = dose_age |> summarise(doses = sum(doses), .by = c(age_months, state)),
           aes(x = age_months, y = doses), fill = NA, col = "grey30", inherit.aes = FALSE) +
  scale_fill_manual(name = "Cohort", values = cols) +
  xlab("Age in months at time of immunisation") +
  ylab("Number of Nirsevimab\ndoses administered") +
  facet_wrap(~state, scales = "free_y") +
  guides(fill = guide_legend(
    keywidth = unit(0.425, "cm"),
    keyheight = unit(0.425, "cm"),
    override.aes = list(size = 0.725))) +
  scale_x_continuous(breaks = seq(0, 50, 10))

# births
births_plot <- ggplot(data = births_all,
                      aes(x = cohort_birth_months, y = births,
                          fill = factor(cohort_birth_months), shape = state)) +
  geom_point(size = 1.5) +
  geom_line(data = births_all[0, ]) +
  ylab("Births") + xlab("Cohort") +
  theme_bw() +
  scale_fill_manual(name = "Cohort", values = cols) +
  scale_shape_manual(name = "State", values = c(24, 22, 21)) +
  scale_y_continuous(breaks = seq(0, 9000, 1000)) +
  guides(shape = guide_legend(order = 1),
         fill = "none") +
  theme(legend.position = c(0.925, 0.9),
        legend.background = element_rect(fill = "transparent"))

ggsave(
    file = "uncertain_variables_plot.pdf",
    (dose_plot + dose_age_plot + births_plot + guide_area()) +
     plot_layout(design = "
     AAAA
     BBBB
     CCDD", guides = "collect") +
      plot_annotation(tag_levels = c("a")) &
      theme(legend.box = "horizontal"),
    height = 25,
    width = 30,
    units = "cm",
    device = "pdf"
    )

#######################
##### QLD model 1 #####
#######################

probs_in <- c(0.025, 0.5, 0.975)

irr_pos_dis <- fit_1$summary("irr_pos_dis", extra_quantiles = ~posterior::quantile2(., probs = probs_in))[,2:4] |>
  as.data.frame() |>
  mutate(week = QLD_model_1_data$week,
         year = QLD_model_1_data$year,
         treatment = QLD_model_1_data$treatment,
         week_year = QLD_model_1_data$week_year,
         week_cont = QLD_model_1_data$week_cont) |>
  rename("l" = 1, "m" = 2, "u" = 3)

QLD_model_1_data <- QLD_model_1_data |>
  bind_cols(fit_1$summary("y_t_pos_dis", extra_quantiles = ~posterior::quantile2(., probs = probs_in))[,2:4] |>
              rename("y_t_pos_l" = 1, "y_t_pos_m" = 2, "y_t_pos_u" = 3))

intercepts <- fit_1$draws("intercept", format = "draws_matrix") |> as.vector()
re_week_year <- fit_1$draws("coef_week_year", format = "draws_matrix") |> as.matrix()
coef_inc_c <- fit_1$draws("coef_inc_c", format = "draws_matrix") |> as.vector()
coef_treat <- fit_1$draws("coef_treat", format = "draws_matrix") |> as.vector()

ind_week_year = sort(unique(QLD_model_1_data$week_cont))

mu_gq_df <- expand.grid(mu_c = seq(1, max(QLD_model_1_data$inc_greater_8m), 3),
                        treatment = c(0, 1))

pred <- lapply(1:nrow(mu_gq_df),
       function(i, mu_gq_df){
         out <- quantile(rowMeans(exp(t(log(7) +
                                          rep(1, ncol(re_week_year)) %o% (coef_inc_c * mu_gq_df[i, "mu_c"]/7 + coef_treat * mu_gq_df[i, "treatment"] + intercepts)
                                        ) +
                                        re_week_year)),
                         probs = probs_in)
         return(data.frame("l" = out[1], "m" = out[2], "u" = out[3]))
         },
       mu_gq_df = mu_gq_df
       ) |> bind_rows() |> as.data.frame()

rownames(pred) <- NULL

coef_df <- cbind(data.frame("coef" = c("alpha", "beta_1", "beta_2", "c")),
                 rbind(quantile(intercepts, probs = probs_in),
                       quantile(coef_inc_c, probs = probs_in),
                       quantile(coef_treat, probs = probs_in),
                       quantile(fit_1$draws("sigma_week_year", format = "draws_matrix"), probs = c(0.025, 0.5, 0.975))) |>
                   as.data.frame() |> rename("l" = 1, "m" = 2, "u" = 3)
                 )

fit_plots <-
  ggplot(data = QLD_model_1_data, aes(x = week_cont, y = inc_less_8m, fill = factor(treatment))) +
  geom_bar(stat = "identity", col = "grey30") +
  geom_pointrange(aes(x = week_cont, y = y_t_pos_m, ymin = y_t_pos_l, ymax = y_t_pos_u,
                      fill = factor(treatment)),
                  inherit.aes = FALSE,
                  shape = 21, size = 0.7, col = "grey30", alpha = 0.7) +
  theme(legend.position = c(0.8, 0.8)) +
  scale_fill_manual(values = c("grey70", "skyblue"), name = "Nirsevimab\ndistribution",
                    labels = c("No", "Yes")) +
  xlab("Week-Year") + ylab("Incidence in those\n8 months old and under") +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  scale_x_continuous(labels = arrange(subset(QLD_model_1_data, week_cont %in% seq(1, 105, 10)), week_cont) |> select(week_year) |> as.vector() |> unlist(), breaks = seq(1, 105, 10)) +

ggplot(data = QLD_model_1_data, aes(x = week_cont, y = inc_greater_8m)) +
  geom_bar(stat = "identity", col = "grey30", fill = "grey70") +
  xlab("Week-Year") + ylab("Incidence in those\nolder than 8-months") +
  scale_y_continuous(limits = c(0, 700), breaks = seq(0, 700, 100)) +
  scale_x_continuous(labels = arrange(subset(QLD_model_1_data, week_cont %in% seq(1, 105, 10)), week_cont) |> select(week_year) |> as.vector() |> unlist(), breaks = seq(1, 105, 10)) +

ggplot(data = QLD_model_1_data,
         aes(x = week_cont, y = inc_less_8m / inc_greater_8m, col = factor(treatment),
             fill = factor(treatment))) +
  geom_bar(stat = "identity", position = position_dodge(), col = "grey30") +
  theme(legend.position = c(0.8, 0.8)) +
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
  theme(legend.position = c(0.3, 0.8), legend.direction = "horizontal") +
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
  theme(legend.position = c(0.3, 0.75)) +
  geom_line(aes(col = factor(treatment)), linewidth = 1) +
  xlab("Weekly incidence in those\nolder than 8 months") +
  scale_colour_manual(values = c("grey70", "skyblue"), name = "Nirsevimab\ndistribution",
                      labels = c("No", "Yes")) +
  scale_fill_manual(values = c("grey70", "skyblue"), name = "Nirsevimab\ndistribution",
                    labels = c("No", "Yes")) +

  ggplot(data = coef_df, aes(x = coef, y = m, ymin = l, ymax = u)) +
  geom_pointrange(alpha = 0.5) +
  ylab("Fitted parameter value") +
  xlab("Parameter") +
  scale_x_discrete(labels = c("alpha" = parse(text = "alpha"),
                              "beta_1" = parse(text = "beta[1]"),
                              "beta_2" = parse(text = "beta[2]"),
                              "c" = parse(text = "c"))) +
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


round(quantile(exp(fit_1$draws("coef_inc_c", format = "draws_matrix")), probs = c(0.025, 0.5, 0.975)), digits = 3)

round(quantile(exp(fit_1$draws("coef_treat", format = "draws_matrix")), probs = c(0.025, 0.5, 0.975)), digits = 2)

ggsave(
  mcmc_trace(fit_1$draws(c("intercept", "coef_inc_c", "coef_treat", "sigma_week_year", "mu_c_prior_mean", "mu_c_prior_sd"))),
  file = "model_1_traceplot.pdf",
  height = 10,
  width = 20,
  units = "cm",
  device = "pdf"
)

#######################
##### cov model 2 #####
#######################

# posteriors
write.csv(fit_3$summary(c("phase_shift_cos_months", "amplitude"), extra_quantiles = ~posterior::quantile2(., probs = probs_in)),
          file = "model_3_phase_shift_amplitude.csv")

write.csv(
  rbind(fit_2_ACT$summary(c("phase_shift_cos_months", "amplitude"), extra_quantiles = ~posterior::quantile2(., probs = probs_in)) |>
  mutate(state = "ACT"),
  fit_2_NSW$summary(c("phase_shift_cos_months", "amplitude"), extra_quantiles = ~posterior::quantile2(., probs = probs_in)) |>
    mutate(state = "NSW"),
  fit_2_QLD$summary(c("phase_shift_cos_months", "amplitude"), extra_quantiles = ~posterior::quantile2(., probs = probs_in)) |>
    mutate(state = "QLD")
  ),
  file = "model_2_phase_shift_amplitude.csv")

write.csv(as.data.frame(
  rbind(round(quantile(fit_3$draws("coef_age_s", format = "draws_matrix")/fit_3$draws("coef_age_l", format = "draws_matrix"), probs = c(0.025, 0.5, 0.975)), digits = 2),
  round(quantile(fit_2_ACT$draws("coef_age_s", format = "draws_matrix")/fit_3$draws("coef_age_l", format = "draws_matrix"), probs = c(0.025, 0.5, 0.975)), digits = 2),
  round(quantile(fit_2_NSW$draws("coef_age_s", format = "draws_matrix")/fit_3$draws("coef_age_l", format = "draws_matrix"), probs = c(0.025, 0.5, 0.975)), digits = 2),
  round(quantile(fit_2_QLD$draws("coef_age_s", format = "draws_matrix")/fit_3$draws("coef_age_l", format = "draws_matrix"), probs = c(0.025, 0.5, 0.975)), digits = 2))
  ) |> mutate(state = c("pooled", "ACT", "NSW", "QLD"), model = c("pooled", rep("independent", 3))),
  file = "peak_age.csv")

# traceplots
ggsave(
  mcmc_trace(fit_2_ACT$draws(c("intercept", "coef_year_raw",
                               "coef_doses", "coef_age_l", "coef_age_s", "amplitude",
                               "phase_shift_cos_months", "sigma_cohorts", "sigma_observations",
                               "mu_doses_prior_mean", "mu_doses_prior_sd", "offset_cohort_births_mean", "offset_cohort_births_sd"))),
  file = "model_2_ACT_traceplot.pdf",
  height = 30,
  width = 40,
  units = "cm",
  device = "pdf"
)

ggsave(
  mcmc_trace(fit_2_NSW$draws(c("intercept", "coef_year_raw",
                               "coef_doses", "coef_age_l", "coef_age_s", "amplitude",
                               "phase_shift_cos_months", "sigma_cohorts", "sigma_observations",
                               "mu_doses_prior_mean", "mu_doses_prior_sd", "offset_cohort_births_mean", "offset_cohort_births_sd"))),
  file = "model_2_NSW_traceplot.pdf",
  height = 30,
  width = 40,
  units = "cm",
  device = "pdf"
)

ggsave(
  mcmc_trace(fit_2_QLD$draws(c("intercept", "coef_year_raw",
                               "coef_doses", "coef_age_l", "coef_age_s", "amplitude",
                               "phase_shift_cos_months", "sigma_cohorts", "sigma_observations",
                               "mu_doses_prior_mean", "mu_doses_prior_sd", "offset_cohort_births_mean", "offset_cohort_births_sd"))),
  file = "model_2_QLD_traceplot.pdf",
  height = 30,
  width = 40,
  units = "cm",
  device = "pdf"
)

ggsave(
  mcmc_trace(fit_3$draws(c("intercept", "coef_year_raw", "coef_state_raw",
                            "coef_doses", "coef_doses_diff", "coef_age_l", "coef_age_s", "amplitude",
                            "phase_shift_cos_months", "sigma_cohorts", "sigma_observations",
                           "mu_doses_prior_mean", "mu_doses_prior_sd", "offset_cohort_births_mean", "offset_cohort_births_sd"))
             ),
  file = "model_3_traceplot.pdf",
  height = 40,
  width = 50,
  units = "cm",
  device = "pdf"
)

##########################################################
##### predicted uncertainty in independent variables #####
##########################################################

extract_mu_doses <- function(data_in, fit, model_name, state_vec){

  data.frame("actual" = data_in$dose_data,
             "state" = state_vec[data_in$dose_data_state_index],
             "m" = apply(fit$draws("mu_doses", format = "draws_matrix"), 2, median),
             "l" = apply(fit$draws("mu_doses", format = "draws_matrix"), 2, quantile, probs = c(0.025)),
             "u" = apply(fit$draws("mu_doses", format = "draws_matrix"), 2, quantile, probs = c(0.975)),
             "model" = rep(model_name, length(data_in$dose_data))
             )
}

extract_offset_births <- function(data_in, fit, model_name, state_vec){

  if(length(state_vec) == 1){
    state <- rep(state_vec, length(data_in$sample_cohort_births))
  } else{
    state <- state_vec[data_in$cohort_births_state_index]
  }

  data.frame("actual" = data_in$sample_cohort_births,
             "state" = state,
             "m" = apply(fit$draws("offset_cohort_births", format = "draws_matrix"), 2, median),
             "l" = apply(fit$draws("offset_cohort_births", format = "draws_matrix"), 2, quantile, probs = c(0.025)),
             "u" = apply(fit$draws("offset_cohort_births", format = "draws_matrix"), 2, quantile, probs = c(0.975)),
             "model" = rep(model_name, length(data_in$sample_cohort_births))
             )
}

mu_doses_df <- rbind(extract_mu_doses(model_2_data_ACT_in, fit_2_ACT, "Independent", c("ACT")),
                     extract_mu_doses(model_2_data_NSW_in, fit_2_NSW, "Independent", c("NSW")),
                     extract_mu_doses(model_2_data_QLD_in, fit_2_QLD, "Independent", c("QLD")),
                     extract_mu_doses(model_3_data_in, fit_3, "Pooled", c("ACT", "NSW", "QLD")))

offset_births_df <- rbind(extract_offset_births(model_2_data_ACT_in, fit_2_ACT, "Independent", c("ACT")),
                     extract_offset_births(model_2_data_NSW_in, fit_2_NSW, "Independent", c("NSW")),
                     extract_offset_births(model_2_data_QLD_in, fit_2_QLD, "Independent", c("QLD")),
                     extract_offset_births(model_3_data_in, fit_3, "Pooled", c("ACT", "NSW", "QLD")))

# "#E69F00", "#CC79A7", "#009E73", "#0072B2", "#D55E00"

ggsave(
  ggplot(data = mu_doses_df,
       aes(x = actual, y = m, ymin = l, ymax = u, fill = model)) +
  geom_pointrange(size = 0.75, alpha = 0.75, shape = 21) +
  xlab("Sampled number of doses") + ylab("Modelled number of doses") +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme(text = element_text(size = 15)) +
  facet_wrap(~state, scales = "free") +
  scale_fill_manual(values = model_cols) +

  ggplot(data = offset_births_df,
       aes(x = actual, y = m, ymin = l, ymax = u, fill = model)) +
  geom_pointrange(size = 0.75, alpha = 0.75, shape = 21) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  theme(text = element_text(size = 15)) +
  facet_wrap(~ state, scales = "free") +
  xlab("Sampled number of births") + ylab("Modelled cohort size") +
  scale_fill_manual(values = model_cols) +
  plot_layout(nrow = 2, guides = "collect") +
  plot_annotation(tag_levels = c("a")),

  file = "measurement_error_plots.pdf",
  device = "pdf",
  height = 20,
  width = 37.5,
  units = "cm"
)

##################################
##### actual vs fitted plots #####
##################################

model_3_data <- cbind(model_3_data,
                          apply(fit_3$draws("y_pos_dis", format = "draws_matrix"), 2, quantile, probs = probs_in) |>
                            t() |>
                            as.data.frame() |>
                            rename("low" = 1, "med" = 2, "up" = 3))

model_2_data_ACT <- cbind(model_2_data_ACT,
                          apply(fit_2_ACT$draws("y_pos_dis", format = "draws_matrix"), 2, quantile, probs = probs_in) |>
                            t() |>
                            as.data.frame() |>
                            rename("low" = 1, "med" = 2, "up" = 3))

model_2_data_NSW <- cbind(model_2_data_NSW,
                          apply(fit_2_NSW$draws("y_pos_dis", format = "draws_matrix"), 2, quantile, probs = probs_in) |>
                            t() |>
                            as.data.frame() |>
                            rename("low" = 1, "med" = 2, "up" = 3))

model_2_data_QLD <- cbind(model_2_data_QLD,
                          apply(fit_2_QLD$draws("y_pos_dis", format = "draws_matrix"), 2, quantile, probs = probs_in) |>
                            t() |>
                            as.data.frame() |>
                            rename("low" = 1, "med" = 2, "up" = 3))

actual_fitted_plot_independent_models <- ggplot(
  data = rbind(model_2_data_ACT, model_2_data_NSW, model_2_data_QLD),
  aes(x = inc, y = med, col = factor(cohort_birth_months))) +
  geom_pointrange(alpha = 0.5, aes(ymin = low, ymax = up)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.75) +
  xlab("Sampled monthly incidence per cohort") +
  ylab("Predicted monthly incidence per cohort") +
  scale_colour_manual(values = cols, name = "Cohort") +
  facet_wrap(~state, scales = "free")

actual_fitted_plot_pooled_model <- ggplot(data = model_3_data,
                                          aes(x = inc, y = med, col = factor(cohort_birth_months))) +
  geom_pointrange(alpha = 0.5, aes(ymin = low, ymax = up)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.75) +
  xlab("Sampled monthly incidence per cohort") +
  ylab("Predicted monthly incidence per cohort") +
  scale_colour_manual(values = cols, name = "Cohort") +
  facet_wrap(~state, scales = "free")

######################
##### covariates #####
######################

n_iter <- length(fit_2_QLD$draws("intercept", format = "draws_matrix"))

# values to predict from GLM
glm_overall_effect_ACT <- t(fit_2_ACT$draws("lp_mu_gq", format = "draws_matrix"))
glm_overall_effect_NSW <- t(fit_2_NSW$draws("lp_mu_gq", format = "draws_matrix"))
glm_overall_effect_QLD <- t(fit_2_QLD$draws("lp_mu_gq", format = "draws_matrix"))
glm_overall_effect_pooled <- t(fit_3$draws("lp_mu_gq", format = "draws_matrix"))

##### average effect plots

# age
# model 2
gen_age_effect <- function(fit, model_data, glm_overall_effect){

  glm_age <- (-model_data$age_rsv_months %o% as.vector(fit$draws("coef_age_l", format = "draws_matrix"))) +
    (log(model_data$age_rsv_months) %o% as.vector(fit$draws("coef_age_s", format = "draws_matrix")))

  e_inc_m_age <- glm_overall_effect - glm_age

  ages_in <- seq(0.1, max(model_data$age_rsv_months), 0.1)
  n_ages <- length(ages_in)

  ages_coef <- (-ages_in %o% as.vector(fit$draws("coef_age_l", format = "draws_matrix"))) +
    (log(ages_in) %o% as.vector(fit$draws("coef_age_s", format = "draws_matrix")))

  avg_eff_age <- sapply(1:n_ages,
                      function(i){
                        quantile(colMeans(exp(e_inc_m_age + t(replicate(nrow(model_data), ages_coef[i,])))), probs = probs_in)
                        },
                      simplify = FALSE
                      ) |> bind_rows() |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate(age_rsv_months = ages_in)

  age_IRR_pooled <-
    apply(
      exp(ages_coef) /
        exp(rep(0.1, length(ages_in)) %o% as.vector(fit$draws("coef_age_l", format = "draws_matrix")) +
              (log(rep(0.1, length(ages_in))) %o% as.vector(fit$draws("coef_age_s", format = "draws_matrix")))),
      1, quantile, probs = probs_in) |>
    t() |> as.data.frame() |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate(age_rsv_months = ages_in)

  return(list("avg_eff_age" = avg_eff_age, "age_IRR_pooled" = age_IRR_pooled))
  }

age_eff_ACT <- gen_age_effect(fit_2_ACT, model_2_data_ACT, glm_overall_effect_ACT)
age_eff_NSW <- gen_age_effect(fit_2_NSW, model_2_data_NSW, glm_overall_effect_NSW)
age_eff_QLD <- gen_age_effect(fit_2_QLD, model_2_data_QLD, glm_overall_effect_QLD)
age_eff_pooled <- gen_age_effect(fit_3, model_3_data, glm_overall_effect_pooled)

avg_eff_age_ACT <- age_eff_ACT$avg_eff_age
avg_eff_age_NSW <- age_eff_NSW$avg_eff_age
avg_eff_age_QLD <- age_eff_QLD$avg_eff_age
avg_eff_age_pooled <- age_eff_pooled$avg_eff_age

avg_eff_age_df <- rbind(avg_eff_age_ACT |> mutate(state = "ACT", model = "independent"),
                        avg_eff_age_NSW |> mutate(state = "NSW", model = "independent"),
                        avg_eff_age_QLD |> mutate(state = "QLD", model = "independent"))

# pooled model
age_plot_pooled <- ggplot(data = model_3_data,
                          aes(x = age_rsv_months, y = inc, col = as.factor(cohort_birth_months))) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = avg_eff_age_pooled,
              aes(x = age_rsv_months, ymin = l, ymax = u), inherit.aes = FALSE, alpha = 0.25) +
  geom_line(data = avg_eff_age_pooled,
            aes(x = age_rsv_months, y = m), inherit.aes = FALSE, linewidth = 0.5) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Age in months") +
  scale_colour_manual(name = "Cohort", values = cols)

age_plot_independent <- ggplot(data = model_3_data,
                          aes(x = age_rsv_months, y = inc, col = as.factor(cohort_birth_months))) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = avg_eff_age_df,
              aes(x = age_rsv_months, ymin = l, ymax = u), inherit.aes = FALSE, alpha = 0.25) +
  geom_line(data = avg_eff_age_df,
            aes(x = age_rsv_months, y = m), inherit.aes = FALSE, linewidth = 0.5) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Age in months") +
  scale_colour_manual(name = "Cohort", values = cols) +
  facet_wrap(~state, scales = "free")

# incidence rate ratio
age_IRR_pooled_plot <- ggplot(data = age_eff_pooled$age_IRR_pooled,
                         aes(x = age_rsv_months, y = m, ymin = l, ymax = u)) +
  geom_ribbon(alpha = 0.25) + geom_line() +
  ylab("Incidence rate ratio relative to 0.1 month old infants") + xlab("Age in months")

age_IRR_independent_plot <- ggplot(data = rbind(age_eff_ACT$age_IRR_pooled |> mutate(state = "ACT"),
                                           age_eff_NSW$age_IRR_pooled |> mutate(state = "NSW"),
                                           age_eff_QLD$age_IRR_pooled |> mutate(state = "QLD")),
                              aes(x = age_rsv_months, y = m, ymin = l, ymax = u, fill = state)) +
  geom_ribbon(alpha = 0.25) + geom_line(aes(col = state), linewidth = 1) +
  ylab("Incidence rate ratio relative to 0.1 month old infants") + xlab("Age in months") +
  scale_colour_manual(values = col_states) +
  scale_fill_manual(values = col_states)

# seasonality

gen_season_effect <- function(fit, model_data, glm_overall_effect){

  glm_season <- (sin(model_data$month * 2 * pi / 12) %o% as.vector(fit$draws("coef_month_s", format = "draws_matrix"))) +
    (cos(model_data$month * 2 * pi / 12) %o% as.vector(fit$draws("coef_month_c", format = "draws_matrix")))

  e_inc_m_season <- glm_overall_effect - glm_season

  months_in <- seq(1, 12, 0.1)

  n_months <- length(months_in)

  months_coef <- ((sin(months_in * 2 * pi / 12) %o% as.vector(fit$draws("coef_month_s", format = "draws_matrix"))) +
                    (cos(months_in * 2 * pi / 12) %o% as.vector(fit$draws("coef_month_c", format = "draws_matrix"))))

  avg_eff_season <- sapply(1:n_months,
                        function(i){
                          quantile(colMeans(exp(e_inc_m_season + t(replicate(nrow(model_data), months_coef[i,])))), probs = probs_in)
                        },
                        simplify = FALSE
  ) |> bind_rows() |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate(months = months_in)

  IRR <- apply(exp(months_coef) / exp((sin(rep(1, n_months) * 2 * pi / 12) %o% as.vector(fit$draws("coef_month_s", format = "draws_matrix"))) +
                        (cos(rep(1, n_months) * 2 * pi / 12) %o% as.vector(fit$draws("coef_month_c", format = "draws_matrix")))),
        1, quantile, probs = probs_in) |>
    t() |> as.data.frame() |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate(months = months_in)

  return(list("avg_eff_season" = avg_eff_season, "IRR" = IRR))

}

eff_season_ACT <- gen_season_effect(fit_2_ACT, model_2_data_ACT, glm_overall_effect_ACT)
eff_season_NSW <- gen_season_effect(fit_2_NSW, model_2_data_NSW, glm_overall_effect_NSW)
eff_season_QLD <- gen_season_effect(fit_2_QLD, model_2_data_QLD, glm_overall_effect_QLD)

avg_eff_season_ACT <- eff_season_ACT$avg_eff_season
avg_eff_season_NSW <- eff_season_NSW$avg_eff_season
avg_eff_season_QLD <- eff_season_QLD$avg_eff_season

coef_month_s_mat <- t(fit_3$draws("coef_month_s", format = "draws_matrix")[,model_3_data$ind_states])
coef_month_l_mat <- t(fit_3$draws("coef_month_c", format = "draws_matrix")[,model_3_data$ind_states])

glm_months <- (sin(model_3_data$month * 2 * pi / 12) * coef_month_s_mat) + (cos(model_3_data$month * 2 * pi / 12) * coef_month_l_mat)

e_inc_m_season <- glm_overall_effect_pooled - glm_months

months_in <- seq(1, 12, 0.1)

n_months <- length(months_in)

state_group_counts <- as.vector(table(model_3_data$ind_states))

avg_eff_season_pooled <- sapply(1:n_months,
                      function(i){

                        months_coef <- t((sin(months_in[i] * 2 * pi / 12) * fit_3$draws("coef_month_s", format = "draws_matrix")[,model_3_data$ind_states]) +
                          (cos(months_in[i] * 2 * pi / 12) * fit_3$draws("coef_month_c", format = "draws_matrix")[,model_3_data$ind_states]))

                        mu <- rowsum(exp(e_inc_m_season + months_coef), group = model_3_data$ind_states) / state_group_counts

                        out <- as.data.frame(t(apply(mu, 1, quantile, probs = c(0.025, 0.5, 0.975)))) |> mutate(months = months_in[i])

                        out[,"state_index"] <- rownames(out)

                        rownames(out) <- NULL

                        return(out)
                      }, simplify = FALSE) |> bind_rows() |> rename("l" = 1, "m" = 2, "u" = 3)

avg_eff_season_pooled$state <- c("ACT", "NSW", "QLD")[as.numeric(avg_eff_season_pooled$state_index)]

IRR_season_pooled <- sapply(1:n_months,
                            function(i){

                              months_coef <- (sin(months_in[i] * 2 * pi / 12) * fit_3$draws("coef_month_s", format = "draws_matrix")) +
                                                 (cos(months_in[i] * 2 * pi / 12) * fit_3$draws("coef_month_c", format = "draws_matrix"))

                              months_coef_1 <- (sin(1 * 2 * pi / 12) * fit_3$draws("coef_month_s", format = "draws_matrix")) +
                                                   (cos(1 * 2 * pi / 12) * fit_3$draws("coef_month_c", format = "draws_matrix"))

                              out <- apply(exp(months_coef) / exp(months_coef_1), 2, quantile, probs = probs_in) |> t() |> as.data.frame() |>
                                mutate(months = months_in[i], state = c("ACT", "NSW", "QLD"))

                              rownames(out) <- NULL

                              return(out)
                            }, simplify = FALSE) |> bind_rows() |> rename("l" = 1, "m" = 2, "u" = 3)


avg_eff_season_df <- rbind(avg_eff_season_ACT |> mutate(state = "ACT", model = "independent"),
                           avg_eff_season_NSW |> mutate(state = "NSW", model = "independent"),
                           avg_eff_season_QLD |> mutate(state = "QLD", model = "independent"))

season_data_plot_independent <- ggplot(data = model_3_data,
                           aes(x = month, y = inc, col = as.character(cohort_birth_months))) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1)) +
  geom_ribbon(data = avg_eff_season_df, aes(x = months, ymin = l, ymax = u), alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = avg_eff_season_df, aes(x = months, y = m), inherit.aes = FALSE, linewidth = 0.5) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Month") +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = month.abb[seq(1, 12, 1)]) +
  scale_colour_manual(name = "Cohort", values = cols) +
  facet_wrap(~state, scales = "free_y")

season_data_plot_pooled <- ggplot(data = model_3_data,
                                  aes(x = month, y = inc, col = as.character(cohort_birth_months))) +
  geom_point(alpha = 0.5, position = position_jitter(width = 0.1)) +
  geom_ribbon(data = avg_eff_season_pooled, aes(x = months, ymin = l, ymax = u), alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = avg_eff_season_pooled, aes(x = months, y = m), inherit.aes = FALSE, linewidth = 0.5) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Month") +
  scale_x_continuous(breaks = seq(1, 12, 1), labels = month.abb[seq(1, 12, 1)]) +
  scale_colour_manual(name = "Cohort", values = cols) +
  facet_wrap(~state, scales = "free_y")

# incidence rate ratio
season_IRR_pooled_plot <- ggplot(data = IRR_season_pooled,
                         aes(x = months, y = m, ymin = l, ymax = u)) +
  geom_ribbon(alpha = 0.25, aes(fill = state)) + geom_line(aes(col = state), linewidth = 1) +
  ylab("Incidence rate ratio relative January") + xlab("Month") +
  scale_colour_manual(values = col_states, name = "State") +
  scale_fill_manual(values = col_states, name = "State") +
  scale_x_continuous(breaks = seq(1, 12), labels = month.abb)

season_IRR_independent_plot <- ggplot(data = rbind(eff_season_ACT$IRR |> mutate(state = "ACT"),
                                           eff_season_NSW$IRR |> mutate(state = "NSW"),
                                           eff_season_QLD$IRR |> mutate(state = "QLD")),
                              aes(x = months, y = m, ymin = l, ymax = u)) +
  geom_ribbon(alpha = 0.25, aes(fill = state)) + geom_line(aes(col = state), linewidth = 1) +
  ylab("Incidence rate ratio relative to January") + xlab("Month") +
  scale_colour_manual(values = col_states) +
  scale_fill_manual(values = col_states) +
  scale_x_continuous(breaks = seq(1, 12), labels = month.abb)

### coverage

# estimating the coverage for each data point
# mu_doses <- fit_3$draws("mu_doses", format = "draws_matrix")
# dose_check <- matrix(0, nrow = nrow(model_3_data), ncol = n_iter)
# dose_check[model_3_data_in$doses_on_inds, ] <- as.matrix(model_3_data_in$doses_mat %*% t(mu_doses))
# cov <- dose_check / t(fit_3$draws("offset_cohort_births", format = "draws_matrix")[,model_3_data$ind_cohorts])
# model_3_data$cov <- rowMeans(cov)
# glm_dose_pooled <- sweep(cov, 2, as.vector(fit_3$draws("coef_doses", format = "draws_matrix")), `*`)

dose_check <- rep(0, length = nrow(model_3_data))
dose_check[model_3_data_in$doses_on_inds] <- as.vector(model_3_data_in$doses_mat %*% model_3_data_in$dose_data)
model_3_data$cov <- dose_check / births_all[model_3_data$ind_cohorts, "births"]

glm_dose_state_specific <- model_3_data$cov * t(fit_3$draws("coef_doses_state", format = "draws_matrix")[, model_3_data$ind_states])

e_inc_m_cov_pooled <- glm_overall_effect_pooled - glm_dose_state_specific

cov_in <- seq(0, 1, 0.025)
n_cov <- length(cov_in)

avg_eff_cov_pooled <- sapply(1:n_cov,
                      function(i){
                        cov_coef <- t(cov_in[i] * fit_3$draws("coef_doses_state", format = "draws_matrix")[,model_3_data$ind_states])
                        mu <- rowsum(exp(e_inc_m_cov_pooled + cov_coef), group = model_3_data$ind_states) / state_group_counts
                        out <- as.data.frame(t(apply(mu, 1, quantile, probs = c(0.025, 0.5, 0.975)))) |> mutate(cov = cov_in[i])
                        out[,"state_index"] <- rownames(out)

                        cov_coef_pooled <- cov_in[i] %o% as.vector(fit_3$draws("coef_doses", format = "draws_matrix"))

                        out <- rbind(out,
                              quantile(colMeans(exp(e_inc_m_cov_pooled + t(replicate(nrow(model_3_data), as.vector(cov_coef_pooled))))), probs = probs_in) |>
                          t() |> as.data.frame() |> mutate(cov = cov_in[i], state_index = 4))

                        return(out)
                        }, simplify = FALSE) |> bind_rows() |> rename("l" = 1, "m" = 2, "u" = 3)

avg_eff_cov_pooled[,"state"] <- c("ACT", "NSW", "QLD", "pooled")[as.numeric(avg_eff_cov_pooled$state_index)]

gen_cov_effect <- function(fit, model_data, model_data_in, glm_overall_effect, births_in){

  dose_check <- rep(0, length = nrow(model_data))
  dose_check[model_data_in$doses_on_inds] <- as.vector(model_data_in$doses_mat %*% model_data_in$dose_data)
  model_data$cov <- dose_check / births_in[model_data$ind_cohorts, "births"]

  glm_cov <- model_data$cov %o% as.vector(fit$draws("coef_doses", format = "draws_matrix"))

  e_inc_m_cov <- glm_overall_effect - glm_cov

  cov_in <- seq(0, 1, 0.025)
  n_cov <- length(cov_in)

  cov_coef <- cov_in %o% as.vector(fit$draws("coef_doses", format = "draws_matrix"))

  avg_eff_season <- sapply(1:n_cov,
                           function(i){
                             quantile(colMeans(exp(e_inc_m_cov + t(replicate(nrow(model_data), cov_coef[i,])))), probs = probs_in)
                           },
                           simplify = FALSE
  ) |> bind_rows() |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate(cov = cov_in)

}

avg_eff_cov_ACT <- gen_cov_effect(fit_2_ACT, model_2_data_ACT, model_2_data_ACT_in, glm_overall_effect_ACT, ACT_births)
avg_eff_cov_NSW <- gen_cov_effect(fit_2_NSW, model_2_data_NSW, model_2_data_NSW_in, glm_overall_effect_NSW, NSW_births)
avg_eff_cov_QLD <- gen_cov_effect(fit_2_QLD, model_2_data_QLD, model_2_data_QLD_in, glm_overall_effect_QLD, QLD_births)

avg_eff_cov_df <- rbind(avg_eff_cov_ACT |> mutate(state = "ACT", model = "independent"),
                        avg_eff_cov_NSW |> mutate(state = "NSW", model = "independent"),
                        avg_eff_cov_QLD |> mutate(state = "QLD", model = "independent"))

cov_data_plot_pooled <- ggplot(data = model_3_data,
                           aes(x = cov, y = inc, col = as.character(cohort_birth_months))) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = subset(avg_eff_cov_pooled, state != "pooled"), aes(x = cov, ymin = l, ymax = u), alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = subset(avg_eff_cov_pooled, state != "pooled"), aes(x = cov, y = m), inherit.aes = FALSE, linewidth = 0.5) +
  scale_colour_manual(name = "Cohort", values = cols) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Proportion of cohort immunised with Nirsevimab in previous 6-months") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent) +
  facet_wrap(~state, scales = "free_y")

cov_data_plot_independent <- ggplot(data = model_3_data,
                               aes(x = cov, y = inc, col = as.character(cohort_birth_months))) +
  geom_point(alpha = 0.5) +
  geom_ribbon(data = avg_eff_cov_df, aes(x = cov, ymin = l, ymax = u), alpha = 0.25, inherit.aes = FALSE) +
  geom_line(data = avg_eff_cov_df, aes(x = cov, y = m), inherit.aes = FALSE, linewidth = 0.5) +
  scale_colour_manual(name = "Cohort", values = cols) +
  ylab("Monthly RSV notifications per cohort") +
  xlab("Proportion of cohort immunised with Nirsevimab in previous 6-months") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), labels = scales::percent) +
  facet_wrap(~state, scales = "free_y")

IRR_cov_pooled <- apply(exp(fit_3$draws("coef_doses_state", format = "draws_matrix")), 2, quantile, probs = probs_in) |>
  t() |> as.data.frame() |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate("state" = c("ACT", "NSW", "QLD")) |>
  rbind(quantile(exp(fit_3$draws("coef_doses", format = "draws_matrix")), probs = probs_in) |> t() |>
          as.data.frame() |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate("state" = c("Pooled")))

IRR_cov_pooled$state <- factor(IRR_cov_pooled$state, levels = c("ACT", "NSW", "QLD", "Pooled"))

IRR_cov_plot_pooled <- ggplot(data = IRR_cov_pooled,
                              aes(x = state, y = m, ymin = l, ymax = u, col = state)) +
  geom_pointrange() +
  xlab("State") + ylab("Incidence rate ratio of 100% Nirsevimab\ncoverage relative 0% coverage") +
  scale_colour_manual(values = c(col_states, "Pooled" = "black")) +
  theme(legend.position = "none")

IRR_cov_independent <- rbind(quantile(exp(fit_2_ACT$draws("coef_doses", format = "draws_matrix")), probs = probs_in) |> t() |>
                               as.data.frame() |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate("state" = c("ACT")),
                             quantile(exp(fit_2_NSW$draws("coef_doses", format = "draws_matrix")), probs = probs_in) |> t() |>
                               as.data.frame() |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate("state" = c("NSW")),
                             quantile(exp(fit_2_QLD$draws("coef_doses", format = "draws_matrix")), probs = probs_in) |> t() |>
                               as.data.frame() |> rename("l" = 1, "m" = 2, "u" = 3) |> mutate("state" = c("QLD"))
                             )

IRR_cov_independent$state <- factor(IRR_cov_independent$state, levels = c("ACT", "NSW", "QLD"))

IRR_cov_plot_independent <- ggplot(data = IRR_cov_independent,
                              aes(x = state, y = m, ymin = l, ymax = u, col = state)) +
  geom_pointrange() +
  xlab("State") + ylab("Incidence rate ratio of 100% Nirsevimab\ncoverage relative 0% coverage") +
  scale_colour_manual(values = col_states) +
  theme(legend.position = "none")

ggsave(
  (actual_fitted_plot_pooled_model) +
    ((age_plot_pooled + theme(legend.position = "none") + plot_spacer() + plot_layout(widths = c(1, 0.5)))) +
    (season_data_plot_pooled + theme(legend.position = "none")) +
    (cov_data_plot_pooled + theme(legend.position = "none")) +
    (age_IRR_pooled_plot + season_IRR_pooled_plot + IRR_cov_plot_pooled) +
    guide_area() +
    plot_layout(nrow = 3, guides = "collect") +
    plot_annotation(tag_levels = c("a")) & theme(legend.box = "horizontal"),
  file = "model_2_plots_pooled.pdf",
  device = "pdf",
  height = 37.5,
  width = 57.5,
  units = "cm"
)

ggsave(
  actual_fitted_plot_independent_models + (age_plot_independent + theme(legend.position = "none")) +
  (season_data_plot_independent + theme(legend.position = "none")) +
    (cov_data_plot_independent + theme(legend.position = "none")) +
    (age_IRR_independent_plot + season_IRR_independent_plot + IRR_cov_plot_independent) + guide_area() +
    plot_layout(nrow = 3, guides = "collect") +
  plot_annotation(tag_levels = c("a")) & theme(legend.box = "horizontal"),
  file = "model_2_plots_independent.pdf",
  device = "pdf",
  height = 37.5,
  width = 57.5,
  units = "cm"
)

################################
##### counterfactual plots #####
################################

# pooled model
model_3_data <- model_3_data |> mutate(cov_cf_on = ifelse(age_rsv_months <= 6 & year >= 2024, 1, 0),
                                       cov_cf_off = 0) |>
  rowwise() |> mutate(cov_cf_on = max(cov_cf_on, cov))

n_obs <- nrow(model_3_data)

cf_obs <- t(fit_3$draws("y_pos_dis", format = "draws_matrix"))

cf_dose_on <- matrix(rpois(n_obs * n_iter,
                           exp(glm_overall_effect_pooled - glm_dose_state_specific +
                                 (model_3_data$cov_cf_on * t(fit_3$draws("coef_doses_state", format = "draws_matrix")[,c(model_3_data_in$ind_states)])))
                           ),
                    nrow = n_obs)

cf_dose_on_pooled <- matrix(rpois(n_obs * n_iter,
                           exp(glm_overall_effect_pooled - glm_dose_state_specific + (model_3_data$cov_cf_on %o% as.vector(fit_3$draws("coef_doses", format = "draws_matrix"))))
                           ),
                           nrow = n_obs)

cf_dose_off <- matrix(rpois(n_obs * n_iter,
                            exp(glm_overall_effect_pooled - glm_dose_state_specific)
                            ),
                      nrow = n_obs)

QLD_inds <- which(model_3_data$state == "QLD" & model_3_data$rsv_start_month >= "2024-01-01")

year_month_all <- model_3_data$rsv_start_month[QLD_inds]

calc_cf <- function(cf, year_month){
  cf_summed_matrix <- rowsum(cf, group = year_month, reorder = TRUE)
  cf_cs_matrix <- apply(cf_summed_matrix, 2, cumsum)
  cf_cs <- apply(cf_cs_matrix, 1, quantile, probs = c(0.025, 0.5, 0.975)) |> t() |> as.data.frame() |>
  rename("l" = 1, "m" = 2, "u" = 3) |>
  mutate("rsv_start_month" = sort(unique(year_month)))
}

cf_cs <- rbind(calc_cf(cf_dose_on[QLD_inds,], year_month_all) |> mutate(model = "100% Nirsevimab coverage - QLD"),
               calc_cf(cf_dose_on_pooled[QLD_inds,], year_month_all) |> mutate(model = "100% Nirsevimab coverage - pooled"),
               calc_cf(cf_dose_off[QLD_inds], year_month_all) |> mutate(model = "0% Nirsevimab coverage"),
               calc_cf(cf_obs[QLD_inds], year_month_all) |> mutate(model = "Predicted"))

obs_cs <- model_3_data[QLD_inds,] |> filter(rsv_start_month >= "2024-01-01") |>
  ungroup() |> summarise(inc = sum(inc), .by = c(rsv_start_month)) |>
  arrange(rsv_start_month) |>
  ungroup() |>
  mutate(c_inc = cumsum(inc)) |>
  filter(rsv_start_month >= "2024-01-01")

cf_all_plot <- ggplot(data = cf_cs,
                      aes(x = rsv_start_month, y = m, ymin = l, ymax = u, fill = model)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(aes(col = model), linewidth = 1) +
  geom_point(data = obs_cs, aes(x = rsv_start_month, y = c_inc), inherit.aes = FALSE) +
  ylab("Cumulative notifications in 0 - 59 month olds") +
  xlab("Month-Year") +
  scale_colour_manual(name = "Counterfactual", values = c("#D55E00", "#C3D7A4", "#009E73", "#0072B2")) +
  scale_fill_manual(name = "Counterfactual", values = c("#D55E00", "#C3D7A4", "#009E73", "#0072B2")) +
  scale_y_continuous(limits = c(0, 21000), breaks = seq(0, 20000, 5000))

# under 6-month olds only

age_index <- which(model_3_data$age_rsv_months <= 6 & model_3_data$state == "QLD" & model_3_data$rsv_start_month >= "2024-01-01")

cf_cs_6_month <- rbind(calc_cf(cf_dose_on[age_index,], model_3_data$rsv_start_month[age_index]) |> mutate(model = "100% Nirsevimab coverage - QLD"),
                       calc_cf(cf_dose_off[age_index,], model_3_data$rsv_start_month[age_index]) |> mutate(model = "0% Nirsevimab coverage"),
                       calc_cf(cf_dose_on_pooled[age_index,], model_3_data$rsv_start_month[age_index]) |> mutate(model = "100% Nirsevimab coverage - pooled"),
                       calc_cf(cf_obs[age_index,], model_3_data$rsv_start_month[age_index]) |> mutate(model = "Predicted"))

obs_cs_6_month <- model_3_data[age_index,] |>
  filter(rsv_start_month >= "2024-01-01") |>
  ungroup() |> group_by(rsv_start_month) |> summarise(inc = sum(inc)) |>
  arrange(rsv_start_month) |>
  ungroup() |>
  mutate(c_inc = cumsum(inc))

cf_6m_plot <- ggplot(data = cf_cs_6_month, aes(x = rsv_start_month, y = m, ymin = l, ymax = u, fill = model)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(aes(col = model), linewidth = 1) +
  geom_point(data = obs_cs_6_month, aes(x = rsv_start_month, y = c_inc), inherit.aes = FALSE) +
  ylab("Cumulative notifications in 0 - 6 month olds") +
  xlab("Month-Year") +
  scale_colour_manual(name = "Counterfactual", values = c("#D55E00", "#C3D7A4", "#009E73", "#0072B2")) +
  scale_fill_manual(name = "Counterfactual", values = c("#D55E00", "#C3D7A4", "#009E73", "#0072B2")) +
  scale_y_continuous(limits = c(0, 3250), breaks = seq(0, 3000, 500))


ggsave(
cf_all_plot + cf_6m_plot + plot_layout(guides = "collect") + plot_annotation(tag_levels = c("a")),
file = "counterfactual_plots_pooled.pdf",
height = 12.5,
width = 35,
units = "cm",
device = "pdf"
)
# independent model

model_2_data_QLD <- model_2_data_QLD |> left_join(subset(model_3_data, state == "QLD")[, c("month", "year", "age_rsv_months", "cohort_birth_start_month", "cohort",
                                                                                           "min_date_nirsevimab", "cov")]) |>
  mutate(cov_cf_on = ifelse(age_rsv_months <= 6 & year >= 2024, 1, 0), cov_cf_off = 0) |>
  rowwise() |> mutate(cov_cf_on = max(cov_cf_on, cov))

n_obs <- nrow(model_2_data_QLD)

cf_obs_independent <- t(fit_2_QLD$draws("y_pos_dis", format = "draws_matrix"))

glm_dose_QLD <- model_2_data_QLD$cov %o% as.vector(fit_2_QLD$draws("coef_doses", format = "draws_matrix"))

cf_dose_on_independent <- matrix(rpois(n_obs * n_iter,
                                       exp(glm_overall_effect_QLD - glm_dose_QLD +
                                             (model_2_data_QLD$cov_cf_on %o% as.vector(fit_2_QLD$draws("coef_doses", format = "draws_matrix"))))
                                       ), nrow = n_obs)

cf_dose_off_independent <- matrix(rpois(n_obs * n_iter,
                                        exp(glm_overall_effect_QLD - glm_dose_QLD)
),
nrow = n_obs)

QLD_inds <- which(model_2_data_QLD$rsv_start_month >= "2024-01-01")
year_month_all <- model_2_data_QLD$rsv_start_month[QLD_inds]

cf_cs_independent <- rbind(calc_cf(cf_dose_on_independent[QLD_inds], year_month_all) |> mutate(model = "100% Nirsevimab coverage"),
                           calc_cf(cf_dose_off_independent[QLD_inds], year_month_all) |> mutate(model = "0% Nirsevimab coverage"),
                           calc_cf(cf_obs_independent[QLD_inds], year_month_all) |> mutate(model = "Predicted"))

obs_cs_independent <- model_2_data_QLD[QLD_inds,] |> filter(rsv_start_month >= "2024-01-01") |>
  ungroup() |> summarise(inc = sum(inc), .by = c(rsv_start_month)) |>
  arrange(rsv_start_month) |>
  ungroup() |>
  mutate(c_inc = cumsum(inc))

cf_all_plot_independent <- ggplot(data = cf_cs_independent,
                                  aes(x = rsv_start_month, y = m, ymin = l, ymax = u, fill = model)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(aes(col = model), linewidth = 1) +
  geom_point(data = obs_cs_independent, aes(x = rsv_start_month, y = c_inc), inherit.aes = FALSE) +
  ylab("Cumulative notifications in 0 - 59 month olds") +
  xlab("Month-Year") +
  scale_colour_manual(name = "Counterfactual", values = c("#D55E00", "#009E73", "#0072B2")) +
  scale_fill_manual(name = "Counterfactual", values = c("#D55E00", "#009E73", "#0072B2")) +
  scale_y_continuous(limits = c(0, 21000), breaks = seq(0, 20000, 5000))

# under 6-month olds only

age_index <- which(model_2_data_QLD$age_rsv_months <= 6 & model_2_data_QLD$rsv_start_month >= "2024-01-01")

cf_cs_6_month_independent <- rbind(calc_cf(cf_dose_on_independent[age_index,], model_2_data_QLD$rsv_start_month[age_index]) |> mutate(model = "100% Nirsevimab coverage"),
                                   calc_cf(cf_dose_off_independent[age_index,], model_2_data_QLD$rsv_start_month[age_index]) |> mutate(model = "0% Nirsevimab coverage"),
                                   calc_cf(cf_obs_independent[age_index,], model_2_data_QLD$rsv_start_month[age_index]) |> mutate(model = "Predicted"))

obs_cs_6_month_independent <- model_2_data_QLD[age_index, ] |>
  ungroup() |> group_by(rsv_start_month) |> summarise(inc = sum(inc)) |>
  arrange(rsv_start_month) |>
  ungroup() |>
  mutate(c_inc = cumsum(inc))

cf_6m_plot_independent <- ggplot(data = cf_cs_6_month_independent,
                                 aes(x = rsv_start_month, y = m, ymin = l, ymax = u, fill = model)) +
  geom_ribbon(alpha = 0.25) +
  geom_line(aes(col = model), linewidth = 1) +
  geom_point(data = obs_cs_6_month_independent, aes(x = rsv_start_month, y = c_inc), inherit.aes = FALSE) +
  ylab("Cumulative notifications in 0 - 6 month olds") +
  xlab("Month-Year") +
  scale_colour_manual(name = "Counterfactual", values = c("#D55E00", "#009E73", "#0072B2")) +
  scale_fill_manual(name = "Counterfactual", values = c("#D55E00", "#009E73", "#0072B2")) +
  scale_y_continuous(limits = c(0, 3250), breaks = seq(0, 3000, 500))

ggsave(
  cf_all_plot_independent + cf_6m_plot_independent + plot_layout(guides = "collect") + plot_annotation(tag_levels = c("a")),
  file = "counterfactual_plots_independent.pdf",
  height = 12.5,
  width = 35,
  units = "cm",
  device = "pdf"
)

##############################
##### model coefficients #####
##############################

coef_table_model_2_fun <- function(fit){

  coef_df_model_2 <- fit$summary(c("intercept", "coef_year_raw", "coef_month_c",
                                   "coef_month_s", "coef_age_l", "coef_age_s",
                                   "coef_doses", "sigma_cohorts", "sigma_observations",
                                   "mu_doses_prior_mean", "mu_doses_prior_sd", "offset_cohort_births_mean", "offset_cohort_births_sd"),
                                 extra_quantiles = ~posterior::quantile2(., probs = probs_in)) |>
  as.data.frame() |> rename("l" = 2, "m" = 3, "u" = 4) |>
  mutate(coef = c("$\\alpha$",
                  "$\\beta_{1}$", "$\\beta_{2}$", "$\\beta_{3}$",
                  "$\\beta_{4}$", "$\\beta_{5}$",
                  "$\\beta_{6}$", "$\\sigma_{1}$", "$\\sigma_{2}$",
                  "$a$", "$b$", "$c$", "$d$"
                  ),
         l = round(l, digits = 2),
         m = round(m, digits = 2),
         u = round(u, digits = 2)) |>
  relocate(c("variable","coef",  "m", "l", "u"))

}

write.csv(coef_table_model_2_fun(fit_2_ACT), file = "coef_model_2_ACT.csv")
write.csv(coef_table_model_2_fun(fit_2_NSW), file = "coef_model_2_NSW.csv")
write.csv(coef_table_model_2_fun(fit_2_QLD), file = "coef_model_2_QLD.csv")

coef_df_model_3 <- fit_3$summary(c("intercept", "coef_year_raw", "coef_state_raw", "coef_month_c",
                                   "coef_month_s", "coef_age_l", "coef_age_s",
                                   "coef_doses", "coef_doses_diff", "sigma_cohorts", "sigma_observations",
                                   "mu_doses_prior_mean", "mu_doses_prior_sd", "offset_cohort_births_mean", "offset_cohort_births_sd"),
                                 extra_quantiles = ~posterior::quantile2(., probs = probs_in)) |>
  as.data.frame() |> rename("l" = 2, "m" = 3, "u" = 4) |>
  mutate(coef = c("$\\alpha$",
                  "$\\beta_{1}$", "$\\beta_{7}$", "$\\beta_{8}$", "$\\beta^{ACT}_{2}$", "$\\beta^{NSW}_{2}$", "$\\beta^{QLD}_{2}$",
                  "$\\beta^{ACT}_{3}$", "$\\beta^{NSW}_{3}$", "$\\beta^{QLD}_{3}$",
                  "$\\beta_{4}$", "$\\beta_{5}$",
                  "$\\beta_{10}$", "$\\beta^{ACT}_{11}$", "$\\beta^{NSW}_{11}$", "$\\beta^{QLD}_{11}$", "$\\sigma_{1}$", "$\\sigma_{2}$",
                  "$a^{ACT}$", "$a^{NSW}$", "$a^{QLD}$", "$b^{ACT}$", "$b^{NSW}$", "$b^{QLD}$",
                  "$c^{ACT}$", "$c^{NSW}$", "$c^{QLD}$", "$d^{ACT}$", "$d^{NSW}$", "$d^{QLD}$"
  ),
  l = round(l, digits = 2),
  m = round(m, digits = 2),
  u = round(u, digits = 2)) |>
  relocate(c("variable","coef",  "m", "l", "u"))

write.csv(coef_df_model_3, file = "coef_model_3.csv")

write.csv(fit_1$summary(c("intercept", "coef_inc_c", "coef_treat", "sigma_week_year", "mu_c_prior_mean", "mu_c_prior_sd"),
              extra_quantiles = ~posterior::quantile2(., probs = probs_in)) |>
  as.data.frame() |> rename("l" = 2, "m" = 3, "u" = 4) |>
  mutate(coef = c("$\\alpha$", "$\\beta_{1}$", "$\\beta_{2}$", "$c$", "$a$", "$b$"),
         l = round(l, digits = 2),
         m = round(m, digits = 2),
         u = round(u, digits = 2)
         ) |>
  relocate(c("variable","coef",  "m", "l", "u")),
  file = "coef_model_1.csv")

####################
##### efficacy #####
####################

eff_model_pooled <- fit_3$summary(c("coef_doses", "coef_doses_state"), extra_quantiles = ~posterior::quantile2((1 - exp(.)), probs = probs_in)) |>
  rename("l" = 2, "m" = 3, "u" = 4) |>
  mutate(l = round(l, digits = 3) * 100,
         m = round(m, digits = 3) * 100,
         u = round(u, digits = 3) * 100,
         eff_model_pooled = paste0(m, " (", l, " - ", u, ")")) |>
  pull(eff_model_pooled)

IRR_model_pooled <- fit_3$summary(c("coef_doses", "coef_doses_state"), extra_quantiles = ~posterior::quantile2(exp(.), probs = probs_in)) |>
  rename("l" = 2, "m" = 3, "u" = 4) |>
  mutate(l = round(l, digits = 3),
         m = round(m, digits = 3),
         u = round(u, digits = 3),
         IRR_model_pooled = paste0(m, " (", l, " - ", u, ")")) |>
  pull(IRR_model_pooled)

eff_model_ind <- rbind(fit_2_ACT$summary("coef_doses", extra_quantiles = ~posterior::quantile2((1 - exp(.)), probs = probs_in)),
      fit_2_NSW$summary("coef_doses", extra_quantiles = ~posterior::quantile2((1 - exp(.)), probs = probs_in)),
      fit_2_QLD$summary("coef_doses", extra_quantiles = ~posterior::quantile2((1 - exp(.)), probs = probs_in))) |>
  rename("l" = 2, "m" = 3, "u" = 4) |>
  mutate(l = round(l, digits = 3) * 100,
         m = round(m, digits = 3) * 100,
         u = round(u, digits = 3) * 100,
         eff_model_ind = paste0(m, " (", l, " - ", u, ")")) |>
  pull(eff_model_ind)

IRR_model_ind <- rbind(fit_2_ACT$summary("coef_doses", extra_quantiles = ~posterior::quantile2(exp(.), probs = probs_in)),
                       fit_2_NSW$summary("coef_doses", extra_quantiles = ~posterior::quantile2(exp(.), probs = probs_in)),
                       fit_2_QLD$summary("coef_doses", extra_quantiles = ~posterior::quantile2(exp(.), probs = probs_in))) |>
  rename("l" = 2, "m" = 3, "u" = 4) |>
  mutate(l = round(l, digits = 3),
         m = round(m, digits = 3),
         u = round(u, digits = 3),
         IRR_model_ind = paste0(m, " (", l, " - ", u, ")")) |>
  pull(IRR_model_ind)

write.csv(data.frame("eff_model_ind" = c(NA, eff_model_ind),
                     "IRR_model_ind" = c(NA, IRR_model_ind),
                     "eff_model_pooled" = eff_model_pooled,
                     "IRR_model_pooled" = IRR_model_pooled),
          file = "efficacy.csv")
