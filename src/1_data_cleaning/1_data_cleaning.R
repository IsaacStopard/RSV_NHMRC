library(tidyverse)
library(readxl)
library(lubridate)
library(stringr)
library(patchwork)

#######################
##### births data #####
#######################

births <- read.csv(file = "births/ABS_BIRTHS_MONTH_OCCURRENCE.csv") |>
  select(Region, Month.of.occurence, TIME_PERIOD, OBS_VALUE) |> rename("state" = Region, "month" = Month.of.occurence, "year" = TIME_PERIOD, "births" = OBS_VALUE) |>
  mutate(state = case_when(state == "New South Wales" ~ "NSW",
                           state == "Australian Capital Territory" ~ "ACT",
                           state == "Queensland" ~ "QLD"),
         ind_states = case_when(state == "ACT" ~ 1,
                                 state == "NSW" ~ 2,
                                 state == "QLD" ~ 3),
         month = case_when(month == "January" ~ 1,
                           month == "February" ~ 2,
                           month == "March" ~ 3,
                           month == "April" ~ 4,
                           month == "May" ~ 5,
                           month == "June" ~ 6,
                           month == "July" ~ 7,
                           month == "August" ~ 8,
                           month == "September" ~ 9,
                           month == "October" ~ 10,
                           month == "November" ~ 11,
                           month == "December" ~ 12),
         cohort_birth_months = zoo::as.yearmon(make_date(year, month, 1)),
         cohort = paste0(state, " ", cohort_birth_months),
         births = ifelse(month == 12 & year == 2024 | month == 11 & year == 2024, NA, births) # birth data for these times are not correct
    ) |> na.omit() |>
  subset(cohort_birth_months > zoo::as.yearmon(make_date(2018, 01, 1)))

# fill in final month with mean value

birth_fit <- glm(births ~ factor(month) + factor(year) + state, family = "poisson", data = births)
birth_pred <- round(predict(birth_fit, newdata = data.frame(month = c(rep(12, 3), rep(11, 3)), year = rep(2024, 6), state = rep(c("ACT", "NSW", "QLD"), 2)), type = "response"), digits = 0)

births <- rbind(births,
                data.frame(month = c(rep(12, 3), rep(11, 3)), year = rep(2024, 6), state = rep(c("ACT", "NSW", "QLD"), 2), births = birth_pred) |>
                      mutate(cohort_birth_months = zoo::as.yearmon(make_date(year, month, 1)),
                             cohort = paste0(state, " ", cohort_birth_months),
                             ind_states = rep(c(1, 2, 3), 2)
                             )
                )

QLD_births <- subset(births, state == "QLD") |> arrange(year, month) |>
  mutate(ind_cohorts = row_number())

saveRDS(QLD_births, file = "QLD_births.rds")

ACT_births <- subset(births, state == "ACT") |> arrange(year, month) |>
  mutate(ind_cohorts = row_number())

saveRDS(ACT_births, file = "ACT_births.rds")

NSW_births <- subset(births, state == "NSW") |> arrange(year, month) |>
  mutate(ind_cohorts = row_number())

saveRDS(NSW_births, file = "NSW_births.rds")

births_all <- births |> arrange(state, year, month) |>
  mutate(ind_cohorts = row_number()) # index to identify the cohorts when running the model

saveRDS(births_all, file = "births_all.rds")

#########################
##### coverage data #####
#########################

dose_df <- read_excel("nirsevimab_doses/F Voight extract 20251114.xlsx", sheet = "clean_data") |>
  mutate(date = as.Date(date, format = "%Y-%m-%d"),
         week = week(date),
         month = month(date),
         year = year(date),
         age_in_months = round(age_in_months, digits = 1),
         cohort_birth_months = zoo::as.yearmon(date %m-% months(age_in_months))
  ) |>
  as.data.frame()

names <- colnames(dose_df)[which(colnames(dose_df) == "totals"):which(colnames(dose_df) == "X")]

for(i in 1:length(names)){
  dose_df[,names[i]] <- suppressWarnings(as.numeric(dose_df[,names[i]]))
  dose_df[is.na(dose_df[,names[i]]), names[i]] <- 0
}

dose_df_sum_age_months <- dose_df |> group_by(month, year, cohort_birth_months) |>
  summarise(ACT = sum(ACT_F) + sum(ACT_M),
            NSW = sum(NSW_F) + sum(NSW_M) + sum(NSW_Unknown),
            QLD = sum(QLD_F) + sum(QLD_M) + sum(QLD_Unknown),
            .groups = "drop") |>
  ungroup() |>
  mutate(start_month = lubridate::make_date(year, month, 1)) |> arrange(year, month, cohort_birth_months)

dose_df_sum_age_months_long <- dose_df_sum_age_months  |>
  pivot_longer(names_to = "state", values_to = "doses", cols = c("ACT", "NSW", "QLD")) |>
  mutate(ind_states = case_when(state == "ACT" ~ 1,
                                 state == "NSW" ~ 2,
                                 state == "QLD" ~ 3),
         cohort = paste0(state, " ", cohort_birth_months))

# checks
if(sum(dose_df_sum_age_months_long$doses) != (sum(dose_df$ACT_F) + sum(dose_df$ACT_M) + sum(dose_df$NSW_F) + sum(dose_df$NSW_M) + sum(dose_df$NSW_Unknown) +
                                              sum(dose_df$QLD_F) + sum(dose_df$QLD_M) + sum(dose_df$QLD_Unknown))){
  stop("incorrect total number of doses in dose_df_sum_age_months_long")
}

saveRDS(dose_df_sum_age_months, "dose_df_sum_age_months.rds")
saveRDS(dose_df_sum_age_months_long, "dose_df_sum_age_months_long.rds")

#############################
##### notification data #####
#############################

# only includes data for those less than 5 years old

# Queensland (QLD)
df_QLD <- read_excel("QLD/QLD data_analysis_original.xlsx")

colnames(df_QLD) <- c("case_number", "date_of_rsv_episode", "week_of_rsv_episode", "year_of_rsv_episode", "week_and_year_of_rsv_episode", "week_and_year_of_rsv_episode_by_age_group",
                      "week_of_birth", "year_of_birth", "sex", "hospitalization",
                      "statistical_area_level_2_of_residence",
                      "age_in_weeks_at_time_of_rsv_episode", "age_in_months_at_time_of_rsv_episode",
                      "age_in_weeks_at_time_of_program_rollout", "age_in_months_at_time_of_program_rollout",
                      "age_6_months_groups_at_time_of_rsv_episode", "binary_age_groups_less_than_6m_vs_geq_6m_at_time_of_rsv_episode", "running_count_of_weekly_notifications", "cumulative_count_of_weekly_notifications", "running_count_of_weekly_notifications_by_age_group",
                      "cumulative_count_of_weekly_notifications_by_age_group")

min_date <- min(df_QLD$date_of_rsv_episode)
max_date <- max(df_QLD$date_of_rsv_episode)

df_QLD <- df_QLD |> mutate(week = week(date_of_rsv_episode),
                           year = year(date_of_rsv_episode),
                           month = month(as.Date(date_of_rsv_episode, format = "%m %Y")),
                           week_year = paste0(week,"-",year),
                           treatment = if_else(week >= 15 & year >= 2024, 1, 0),
                           cohort_birth_start_week = make_date(year_of_birth, 1, 1) + weeks(week_of_birth - 1),
                           cohort_birth_start_month = floor_date(cohort_birth_start_week, unit = "month"),
                           cohort_birth_months = zoo::as.yearmon(make_date(year_of_birth, 1, 1) + weeks(week_of_birth - 1)),
                           state = "QLD",
                           age_rsv_months = interval(cohort_birth_start_month, date_of_rsv_episode) %/% months(1),
                           # model 1
                           age_rsv_months_m1 = interval(cohort_birth_start_week, date_of_rsv_episode) %/% months(1),
                           age_group = if_else(age_rsv_months_m1 <= 8, "age_l_8", "age_g_8")
                           ) |>
  arrange(date_of_rsv_episode) |>
  filter(age_rsv_months >= 0 & age_rsv_months < (5 * 12))

total_QLD <- nrow(df_QLD)

# Australian Capital territories (ACT)

df_ACT <- read_excel("ACT/ACT RSV Data less than 5 years 2023 and 2024.xlsx") |>
  janitor::clean_names()

# notification with ndms_id week of birth is greater than episode date
# similar ndms_ids were entered with similar RSV notification dates
# so, assumed the year was incorrectly entered and should be 2022 not 2023
df_ACT[df_ACT$ndms_id == 1152306, "week_of_birth"] <- "52-2022"

df_ACT <- df_ACT |>
  mutate(week = week(episode_date),
         year = year(episode_date),
         month = month(as.Date(episode_date, format = "%m %Y")),
         rsv_start_month = make_date(year, month, day = 1),
         ndays_per_month = days_in_month(rsv_start_month),
         year_of_birth = as.numeric(str_extract(week_of_birth, "[0-9]{4}$")),
         week_of_birth = as.numeric(str_extract(week_of_birth, "^[0-9]{1,2}")),
         cohort_birth_start_week = make_date(year_of_birth, 1, 1) + weeks(week_of_birth - 1),
         cohort_birth_start_month = floor_date(cohort_birth_start_week, unit = "month"),
         cohort_birth_months = zoo::as.yearmon(make_date(year_of_birth, 1, 1) + weeks(week_of_birth - 1)),
         age_rsv_months = interval(cohort_birth_start_month, rsv_start_month) %/% months(1),
         state = "ACT") |>
  filter(age_rsv_months >= 0 & age_rsv_months < (5 * 12))

total_ACT <- nrow(df_ACT)

# New South Wales (NSW)
# RSV notifications with a calculated onset date between 1 January 2023 to 31 December 2024.
# Source: Notifiable Conditions Records for Epidemiology and Surveillance, NSW Ministry of Health (data extracted 05 March 2026).

df_NSW <- read_excel("NSW/RSV notifications 2023-2024_20260305.xlsx") |>
  janitor::clean_names() |>
  rename("month" = month_of_condition_onset, "year" = year_of_condition_onset) |>
  mutate(rsv_start_month = make_date(year, month, day = 1),
         ndays_per_month = days_in_month(rsv_start_month),
         cohort_birth_start_month = my(year_and_month_of_birth),
         cohort_birth_months = zoo::as.yearmon(cohort_birth_start_month),
         state = "NSW",
         age_rsv_months = interval(cohort_birth_start_month, rsv_start_month) %/% months(1)) |>
  filter(age_rsv_months >= 0 & age_rsv_months < (5 * 12))

total_NSW <- sum(df_NSW$total_notifications)

#############################
##### combined datasets #####
#############################

nirsevimab_duration_months <- 6

# filling in the zero counts
# assumes that the surveillance system is always active post 01/01/2023
# maximum age with surveillance is assumed to be 72 months

all_cohorts <- expand.grid("month" = seq(1, 12), "year" = seq(2018, 2024)) |>
  mutate(cohort_birth_start_month = make_date(year = year, month = month, day = 1))

full_combs <- expand.grid("month" = seq(1, 12),
                          "year" = sort(unique(df_QLD$year)),
                          "cohort_birth_start_month" = unique(all_cohorts$cohort_birth_start_month),
                           "state" = c("ACT", "NSW", "QLD")) |>
  mutate(rsv_start_month = make_date(year = year, month = month, day = 1),
         cohort_birth_months = zoo::as.yearmon(cohort_birth_start_month),
         age_rsv_months = interval(cohort_birth_start_month, rsv_start_month) %/% months(1)
         ) |>
  filter(age_rsv_months >= 0 & age_rsv_months < (5 * 12)) # only includes data for those less than 5 years old

obs_combs <- rbind(df_QLD |> summarise(inc = n(), .by = c(state, cohort_birth_start_month, month, year)),
                   df_NSW |> summarise(inc = sum(total_notifications), .by = c(state, cohort_birth_start_month, month, year)),
                   df_ACT |> summarise(inc = n(), .by = c(state, cohort_birth_start_month, month, year))
                   )

if(sum(obs_combs$inc) != (total_ACT + total_NSW + total_QLD)){
  stop("incorrect total notifications in obs_combs")
}

model_2_data <- left_join(full_combs,
                            obs_combs,
                            by = c("cohort_birth_start_month", "month", "year", "state")) |>
  replace_na(list("inc" = 0)) |>
  mutate(min_date_nirsevimab = lubridate::make_date(year, month, 1) %m-% months(nirsevimab_duration_months),
         offset_months = days_in_month(rsv_start_month),
         cohort = paste0(state, " ", cohort_birth_months),
         cohort_month_year = paste0(cohort,"-", month,"-", year),
         ind_states = case_when(state == "ACT" ~ 1,
                                 state == "NSW" ~ 2,
                                 state == "QLD" ~ 3))

# checks
if(nrow(model_2_data) != nrow(full_combs)){
  stop("incorrect number of data rows in model_2_data")
}

if(sum(model_2_data$inc, na.rm = TRUE) != (total_ACT + total_NSW + total_QLD) |
   sum(subset(model_2_data, state == "ACT")$inc) != total_ACT |
   sum(subset(model_2_data, state == "NSW")$inc) != total_NSW |
   sum(subset(model_2_data, state == "QLD")$inc) != total_QLD){
  stop("incorrect notifications in model_2_data")
}

# checking cohorts are the same in birth and notification datasets
if(!all(unique(model_2_data$cohort) %in% unique(births_all$cohort)) |
   !all(unique(births_all$cohort) %in% unique(model_2_data$cohort))){
  stop("cohorts are not the same in births_all and model_2_data")
}

##### data for all states #####

process_model_data <- function(model_data, births_data){

  return(left_join(model_data,
            data.frame("year" = sort(unique(model_data$year))) |> mutate(ind_years = row_number()),
            by = "year") |>
    left_join(births_data[,c("cohort", "ind_cohorts")], by = "cohort") |>
    left_join(data.frame("cohort_month_year" = unique(model_data$cohort_month_year)) |> mutate(ind_observations = row_number()), by = "cohort_month_year")
    )

  }

model_2_data_all <- process_model_data(model_data = model_2_data, births_data = births_all)

# check
if(length(unique(model_2_data_all$ind_observations)) != nrow(model_2_data_all)){
  stop("each row in model_2_data is not coded as an independent obervation")
}

model_2_data_ACT <- process_model_data(subset(model_2_data, state == "ACT"), ACT_births)
model_2_data_NSW <- process_model_data(subset(model_2_data, state == "NSW"), NSW_births)
model_2_data_QLD <- process_model_data(subset(model_2_data, state == "QLD"), QLD_births)

############################
##### Nirsevimab doses #####
############################

# sums the doses in the previous 6 months for each cohort birth

dose_df <- dose_df_sum_age_months_long |>
  filter(cohort_birth_months %in% unique(model_2_data$cohort_birth_months) & doses > 0)

# assume that all Nirsevimab distributions started on 01-01-2024
# dose_df |> group_by(state) |> summarise(min_start_month = min(start_month))
min_month <- as.Date("01-01-2024", format = "%d-%m-%Y")
max_month <- as.Date("01-12-2024", format = "%d-%m-%Y")

# checks
# missing cohorts
rm_cohorts <- dose_df_sum_age_months_long |> filter((cohort_birth_months %in% unique(model_2_data$cohort_birth_months) == 0)) |> pull(cohort_birth_months) |>
  unique() |> sort() # none of these cohorts are included in our analyses

if(sum(dose_df_sum_age_months_long$doses) !=
   sum(dose_df$doses) + dose_df_sum_age_months_long |> filter(cohort_birth_months %in% rm_cohorts) |> select(doses) |> sum()){
  stop("doses do not sum to the correct number once cohorts have been removed")
}

check_doses <- sum(subset(dose_df, year == 2024)$doses)

# getting the possible combinations once nirsevimab distribution has started
# filling in the 0 counts
# only fill zeroes up to 24-months
dose_df <- left_join(
  unique(model_2_data[,c("cohort_birth_months", "cohort_birth_start_month", "cohort", "state", "ind_states", "rsv_start_month", "month", "year")]) |>
    subset(rsv_start_month >= min_month) |> select(-rsv_start_month),
                         dose_df |> subset(start_month <= max_month),
  by = c("month", "year", "cohort_birth_months", "state", "ind_states", "cohort")
  ) |>
  mutate(start_month = as.Date(ifelse(is.na(start_month), lubridate::make_date(year, month, 1), start_month)),
         age_months = interval(cohort_birth_start_month, start_month) %/% months(1),
         doses = ifelse(is.na(doses), 0, doses)) |>
  filter(!(age_months > 24 & doses == 0)) # remove added zeroes if older than 24 months

# checks
if(sum(dose_df$doses) != check_doses){
  stop("the total doses in doses_df is not correct")
}

saveRDS(dose_df, file = "dose_df.rds")

# fixing zero nirsevimab prior to dose data times

process_doses <- function(model_data, dose_data){

  if(length(unique(dose_data$ind_states)) == 1){
    dose_data$ind_states <- rep(1, nrow(dose_data))
  }

  doses_on_inds <- sort(which(model_data$rsv_start_month >= min_month))

  N_doses_on_inds <- length(doses_on_inds)

  N_doses <- nrow(dose_data)

  # matrix calculation

  doses_mat <- matrix(data = 0, nrow = N_doses_on_inds, ncol = N_doses)

  model_data_doses_ind <- lapply(1:N_doses_on_inds,
                                     function(i){
                                       return(
                                         which(dose_data$cohort == model_data[doses_on_inds[i], "cohort"][[1]] &
                                                 dose_data$start_month > model_data[doses_on_inds[i], "min_date_nirsevimab"][[1]] &
                                                 dose_data$start_month <= model_data[doses_on_inds[i], "rsv_start_month"][[1]] &
                                                 dose_data$state == model_data[doses_on_inds[i], "state"]
                                         )
                                       )
                                     })

  for(i in 1:N_doses_on_inds){
    doses_mat[i, model_data_doses_ind[[i]]] <- 1
  }

  N_doses_mat_one <- sum(doses_mat)

  # checks
  if(!all(unique(colSums(doses_mat)) %in% seq(1, 6, 1))){
    stop("incorrect number of months aggregated for the cumulative doses")
  }
  dose_check <- rep(0, nrow(model_data))
  dose_check[doses_on_inds] <- as.vector(doses_mat %*% dose_data$doses)
  dose_check_sum <- rep(0, nrow(model_data))
  for(i in 1:nrow(model_data[doses_on_inds,])){
    dose_check_sum[doses_on_inds[i]] <- sum(dose_data$doses[model_data_doses_ind[[i]]])
  }
  sparse_doses_mat <- as(Matrix::Matrix(doses_mat, sparse = TRUE), "RsparseMatrix")
  dose_check_sm <- as.numeric(sparse_doses_mat %*% as.matrix(dose_data$doses))

  if(!all(dose_check[doses_on_inds] == dose_check_sm)){
    stop("sparse matrix dose multiplication different to non-sparse matrix dose multiplication")
  }

  if(!all(dose_check_sum[doses_on_inds] == dose_check_sm)){
    stop("manual sum of doses does not equal sparse matrix multiplication")
  }

  return(list("N_doses" = N_doses,
              "dose_data" = dose_data$doses,
              "dose_data_state_index" = dose_data$ind_states,
              "N_doses_on_inds" = N_doses_on_inds,
              "doses_on_inds" = doses_on_inds,
              "doses_mat" = doses_mat,
              "N_doses_mat_one" = N_doses_mat_one))

}

model_2_doses_all <- process_doses(model_2_data_all, dose_df)

model_2_doses_ACT <- process_doses(model_2_data_ACT, subset(dose_df, state == "ACT"))
model_2_doses_NSW <- process_doses(model_2_data_NSW, subset(dose_df, state == "NSW"))
model_2_doses_QLD <- process_doses(model_2_data_QLD, subset(dose_df, state == "QLD"))

### formating the data for Stan model

model_2_data_all$age_rsv_months <- model_2_data_all$age_rsv_months + 0.01

saveRDS(model_2_data_all, "model_2_data_all.rds")

model_2_data_ACT$age_rsv_months <- model_2_data_ACT$age_rsv_months + 0.01
saveRDS(model_2_data_ACT, "model_2_data_ACT.rds")

model_2_data_NSW$age_rsv_months <- model_2_data_NSW$age_rsv_months + 0.01
saveRDS(model_2_data_NSW, "model_2_data_NSW.rds")

model_2_data_QLD$age_rsv_months <- model_2_data_QLD$age_rsv_months + 0.01
saveRDS(model_2_data_QLD, "model_2_data_QLD.rds")

model_3_data_in <- c(list(
  "N" = nrow(model_2_data_all),
  "N_months" = 12,
  "N_years" = length(unique(model_2_data_all$year)),
  "N_states" = length(unique(model_2_data_all$state)),
  "N_cohorts" = length(unique(model_2_data_all$cohort)),
  "N_observations" = length(unique(model_2_data_all$cohort_month_year)),
  "y" = model_2_data_all$inc,
  "sample_cohort_births" = births_all$births,
  "cohort_births_state_index" = births_all$ind_states,
  "months" = model_2_data_all$month,
  "ind_years" = model_2_data_all$ind_years,
  "ind_cohorts" = model_2_data_all$ind_cohorts,
  "ind_observations" = model_2_data_all$ind_observations,
  "ind_states" = model_2_data_all$ind_states,
  "age_months" = model_2_data_all$age_rsv_months,
  "offset_months" = as.vector(model_2_data_all$offset_month),
  "prior_in" = 1.5,
  "sd_age_months" = sd(model_2_data_all$age_rsv_months),
  "prior_in_cohort_births_mean_mean" = c(5000/12, 90000/12, 60000/12), # estimated from https://www.abs.gov.au/statistics/people/population/population-clock-pyramid
  "prior_in_cohort_births_sd_mean" = c(25, 250, 150),
  "prior_in_mu_doses_mean_mean" = c(2.5, 15, 75.0),
  "prior_in_mu_doses_sd_mean" = c(0.5, 2.5, 10.0)),

  model_2_doses_all
)

saveRDS(model_3_data_in, "model_3_data_in.rds")

data_in_individual_state_fun <- function(model_2_data_df, births_data_df, doses_list_in, state){
  model_2_data_ACT_in <- c(list(
  "N" = nrow(model_2_data_df),
  "N_months" = 12,
  "N_years" = length(unique(model_2_data_df$year)),
  "N_cohorts" = length(unique(model_2_data_df$cohort)),
  "N_observations" = length(unique(model_2_data_df$cohort_month_year)),
  "y" = model_2_data_df$inc,
  "sample_cohort_births" = births_data_df$births,
  "months" = model_2_data_df$month,
  "ind_years" = model_2_data_df$ind_years,
  "ind_cohorts" = model_2_data_df$ind_cohorts,
  "ind_observations" = model_2_data_df$ind_observations,
  "age_months" = model_2_data_df$age_rsv_months,
  "offset_months" = as.vector(model_2_data_df$offset_month),
  "prior_in" = 1.5,
  "sd_age_months" = sd(model_2_data_df$age_rsv_months),
  "prior_in_cohort_births_mean_mean" = if(state == "ACT"){5000/12} else if(state == "NSW"){90000/12} else{60000/12}, # estimated from https://www.abs.gov.au/statistics/people/population/population-clock-pyramid
  "prior_in_cohort_births_sd_mean" = if(state == "ACT"){25} else if(state == "NSW"){250} else{150},
  "prior_in_mu_doses_mean_mean" = if(state == "ACT"){2.5} else if(state == "NSW"){15} else{75.0},
  "prior_in_mu_doses_sd_mean" = if(state == "ACT"){0.5} else if(state == "NSW"){2.5} else{10.0}),
  doses_list_in
  )
}

model_2_data_ACT_in <- data_in_individual_state_fun(model_2_data_ACT, ACT_births, model_2_doses_ACT, "ACT")
saveRDS(model_2_data_ACT_in, "model_2_data_ACT_in.rds")

model_2_data_NSW_in <- data_in_individual_state_fun(model_2_data_NSW, NSW_births, model_2_doses_NSW, "NSW")
saveRDS(model_2_data_NSW_in, "model_2_data_NSW_in.rds")

model_2_data_QLD_in <- data_in_individual_state_fun(model_2_data_QLD, QLD_births, model_2_doses_QLD, "QLD")
saveRDS(model_2_data_QLD_in, "model_2_data_QLD_in.rds")

# checks
if(sum(model_3_data_in$y) != sum(c(model_2_data_ACT_in$y, model_2_data_NSW_in$y, model_2_data_QLD_in$y))){stop("data_in notifications not correct")}
if(sum(model_3_data_in$dose_data) != sum(c(model_2_data_ACT_in$dose_data, model_2_data_NSW_in$dose_data, model_2_data_QLD_in$dose_data))){stop("data_in doses not correct")}

########################
##### model 1 data #####
########################

##### Queensland data ######

#QLD_model_1_data <- df_QLD |> group_by(week, year, week_year, age_group, treatment) |> summarise(inc = n()) |> ungroup()

QLD_model_1_data <- df_QLD |> mutate(epi_week = epiweek(date_of_rsv_episode), epi_year = epiyear(date_of_rsv_episode)) |>
  summarise(inc = n(), .by = c(epi_week, epi_year, age_group))

if(sum(QLD_model_1_data$inc) != total_QLD){
  stop("QLD_model_1_data_inc is not correct")
}

epi_weeks_df <- data.frame(
  date = seq(ymd("2023-01-01"), ymd("2024-12-31"), by = "day")
) |> mutate(
  epi_year = epiyear(date),
  epi_week = epiweek(date)
  ) |> summarise(offset_week = n(), .by = c(epi_week, epi_year)) |>
  mutate(week_cont = row_number(),
         week_year = paste0(epi_week,"-",epi_year),
         treatment = if_else(as.Date(paste(epi_year, epi_week, "0", sep = "-"), "%Y-%U-%w") >= as.Date("2024-04-08", format = "%Y-%m-%d"), 1, 0))

# QLD_model_1_data <- left_join(QLD_model_1_data,
#                               data.frame("date" = seq(min_date, max_date, by = "days")) |>
#                                 mutate(week = week(date), year = year(date)) |>
#                                 group_by(week, year) |> summarise(offset_week = n())
#                               )

QLD_model_1_data <- epi_weeks_df |>
  left_join(QLD_model_1_data, by = c("epi_week", "epi_year"))
  #left_join(unique(QLD_model_1_data[,c("week", "year", "week_year")]) |> arrange(year, week) |> ungroup() |> mutate(week_cont = row_number()))

QLD_model_1_data <- QLD_model_1_data |>
  pivot_wider(names_from = age_group, values_from = inc) |>
  replace_na(list("age_g_8" = 0, "age_l_8" = 0)) |>
  rename("inc_less_8m" = age_l_8, "inc_greater_8m" = age_g_8) |>
  mutate(IRR = inc_less_8m/inc_greater_8m)

# checking the total infections are still correct
if(total_QLD != (sum(QLD_model_1_data$inc_greater_8m) + sum(QLD_model_1_data$inc_less_8m))){
  stop("QLD_model_1_data notifications not correct")
}

saveRDS(QLD_model_1_data, file = "QLD_model_1_data.rds")

QLD_model_1_data_in <- list("N" = nrow(QLD_model_1_data),
                            "y_t" = as.integer(QLD_model_1_data$inc_less_8m),
                            "y_c" = as.integer(QLD_model_1_data$inc_greater_8m),
                            "treatment" = as.integer(QLD_model_1_data$treatment),
                            "offset_weeks" = QLD_model_1_data$offset_week,
                            "mu_c_gq" = 1:max(QLD_model_1_data$inc_greater_8m),
                            "offset_weeks_gq" = rep(7, max(QLD_model_1_data$inc_greater_8m)),
                            "N_gq" = max(QLD_model_1_data$inc_greater_8m),
                            "treatment_gq" = rep(1, max(QLD_model_1_data$inc_greater_8m)),
                            "ind_week_year" = as.integer(QLD_model_1_data$week_cont),
                            "N_week_year" = length(sort(unique(QLD_model_1_data$week_cont))),
                            "prior_in" = 1.5,
                            "sd_y_c" = sd(QLD_model_1_data$inc_greater_8m))

saveRDS(QLD_model_1_data_in, file = "QLD_model_1_data_in.rds")

if(sum(model_2_data_QLD$inc) != (sum(QLD_model_1_data$inc_less_8m) + sum(QLD_model_1_data$inc_greater_8m))){
  stop("QLD model 1 and 2 notifications are not equal")
}

##### ACT data #####
ACT_model_1_data <- epi_weeks_df |>
  left_join(df_ACT |> mutate(epi_week = epiweek(episode_date),
                             epi_year = epiyear(episode_date),
    age_group = if_else(age_rsv_months >= 6, "inc_greater_6m", "inc_less_6m")) |>
              summarise(inc = n(), .by = c(epi_week, epi_year, age_group)) |>
              pivot_wider(values_from = inc, names_from = age_group), by = c("epi_week", "epi_year")) |>
  replace_na(list("inc_less_6m" = 0, "inc_greater_6m" = 0))

if((sum(ACT_model_1_data$inc_greater_6m) + sum(ACT_model_1_data$inc_less_6m)) != nrow(df_ACT)){
  stop("notifications in ACT_model_1_data does not sum correctly")
}

if(nrow(unique(ACT_model_1_data[,c("epi_week", "epi_year")])) != nrow(QLD_model_1_data)){
  stop("not all weeks and years have data in ACT_model_1_data")
}

saveRDS(ACT_model_1_data, file = "ACT_model_1_data.rds")

ACT_model_1_data_in <- list("N" = nrow(ACT_model_1_data),
                            "y_t" = as.integer(ACT_model_1_data$inc_less_6m),
                            "y_c" = as.integer(ACT_model_1_data$inc_greater_6m),
                            "treatment" = as.integer(ACT_model_1_data$treatment),
                            "offset_weeks" = ACT_model_1_data$offset_week,
                            "mu_c_gq" = 1:max(ACT_model_1_data$inc_greater_6m),
                            "offset_weeks_gq" = rep(7, max(ACT_model_1_data$inc_greater_6m)),
                            "N_gq" = max(ACT_model_1_data$inc_greater_6m),
                            "treatment_gq" = rep(1, max(ACT_model_1_data$inc_greater_6m)),
                            "ind_week_year" = as.integer(ACT_model_1_data$week_cont),
                            "N_week_year" = length(sort(unique(ACT_model_1_data$week_cont))),
                            "prior_in" = 1.5,
                            "sd_y_c" = sd(ACT_model_1_data$inc_greater_6m))

saveRDS(ACT_model_1_data_in, "ACT_model_1_data_in.rds")

if(sum(model_2_data_ACT$inc) != (sum(ACT_model_1_data$inc_less_6m) + sum(ACT_model_1_data$inc_greater_6m))){
  stop("ACT model 1 and 2 notifications are not equal")
}

##### NSW data #####

NSW_model_1_data <-
  epi_weeks_df |> left_join(
    read_excel(path = "NSW/RSV notifications and hospitalisations 2023-2024_20250923.xlsx", sheet = 1, .name_repair = "universal")[1:210,] |>
      select(Epi.week, Age.group, Total.notifications) |>
      rename("inc" = Total.notifications, "age_group" = Age.group, "week_year" = Epi.week) |>
      mutate(age_group = ifelse(age_group == "0-<6mths", "inc_less_6m", "inc_greater_6m")) |>
      separate_wider_delim(week_year, delim = "-", names = c("epi_week", "epi_year"), cols_remove = FALSE) |>
      mutate(epi_week = as.numeric(epi_week), epi_year = as.numeric(epi_year)) |>
      pivot_wider(values_from = inc, names_from = age_group),
  by = c("epi_week", "epi_year", "week_year")
  )

if(sum(model_2_data_NSW$inc) != (sum(NSW_model_1_data$inc_less_6m) + sum(NSW_model_1_data$inc_greater_6m))){
  print("NSW model 1 and 2 notifications are not equal")
}

saveRDS(NSW_model_1_data, file = "NSW_model_1_data.rds")

NSW_model_1_data_in <- list("N" = nrow(NSW_model_1_data),
                            "y_t" = as.integer(NSW_model_1_data$inc_less_6m),
                            "y_c" = as.integer(NSW_model_1_data$inc_greater_6m),
                            "treatment" = as.integer(NSW_model_1_data$treatment),
                            "offset_weeks" = NSW_model_1_data$offset_week,
                            "mu_c_gq" = 1:max(NSW_model_1_data$inc_greater_6m),
                            "offset_weeks_gq" = rep(7, max(NSW_model_1_data$inc_greater_6m)),
                            "N_gq" = max(NSW_model_1_data$inc_greater_6m),
                            "treatment_gq" = rep(1, max(NSW_model_1_data$inc_greater_6m)),
                            "ind_week_year" = as.integer(NSW_model_1_data$week_cont),
                            "N_week_year" = length(sort(unique(NSW_model_1_data$week_cont))),
                            "prior_in" = 1.5,
                            "sd_y_c" = sd(NSW_model_1_data$inc_greater_6m))

saveRDS(NSW_model_1_data_in, file = "NSW_model_1_data_in.rds")

################################
##### QLD old model 2 data #####
################################
#
# full_combs <- expand.grid("month" = sort(unique(df_QLD$month)),
#                           "year" = sort(unique(df_QLD$year)),
#                           "cohort_birth_start_month" = sort(unique(df_QLD$cohort_birth_start_month)))
#
# obs_combs <- df_QLD |> group_by(cohort_birth_start_month, month, year) |> summarise(inc = n()) |> ungroup()
#
# nrow(obs_combs) == nrow(unique(df_QLD[,c("cohort_birth_start_month", "month", "year")]))
#
# QLD_model_2_data <- left_join(full_combs, obs_combs, by = c("cohort_birth_start_month", "month", "year")) |>
#   replace_na(list("inc" = 0)) |>
#   mutate(min_date_nirsevimab = lubridate::make_date(year, month, 1) %m-% months(nirsevimab_duration_months),
#          rsv_start_month = lubridate::make_date(year, month, 1),
#          age_rsv_months = interval(cohort_birth_start_month, rsv_start_month) %/% months(1),
#          cohort_birth_months = zoo::as.yearmon(cohort_birth_start_month),
#          cohort_month_year = paste0(cohort_birth_months,"-", month,"-", year),
#          offset_months = days_in_month(rsv_start_month)) |>
#   filter(age_rsv_months >= 0 & age_rsv_months <= (6 * 12)) # only includes data for those equal to or less than 5 years old
#
# if(sum(QLD_model_2_data$inc) != nrow(df_QLD)){stop("notifications sum in model 2 not correct")}
#
# max(QLD_model_2_data$age_rsv_months)
#
# total == sum(QLD_model_2_data$inc)
#
# QLD_model_2_data <- left_join(QLD_model_2_data,
#                               data.frame("year" = unique(QLD_model_2_data$year) |> sort()) |> mutate(ind_years = row_number())
#                               )
#
# unique(QLD_model_2_data$cohort_birth_months %in% QLD_births$cohort_birth_months)
#
# unique(QLD_births$cohort_birth_months %in% QLD_model_2_data$cohort_birth_months)
#
# QLD_births <- QLD_births |> mutate(ind_cohort_births = row_number())
#
# QLD_model_2_data <- left_join(QLD_model_2_data,
#                               QLD_births[,c("cohort_birth_months", "ind_cohort_births")]
#                               )
#
# QLD_model_2_data <- left_join(QLD_model_2_data,
#                               data.frame("cohort_month_year" = unique(QLD_model_2_data$cohort_month_year) |> sort()) |> mutate(ind_cohort_month_years = row_number())
# )
#
# # sums the doses in the previous 6 months for each cohort birth
# dose_df_QLD <- dose_df_sum_age_months |> filter(cohort_birth_months %in% unique(QLD_model_2_data$cohort_birth_months) & QLD > 0) |> select(-ACT, -NSW)
#
# sum(dose_df_sum_age_months$QLD)
# sum(dose_df_QLD$QLD) + dose_df_sum_age_months |> filter((cohort_birth_months %in% unique(QLD_model_2_data$cohort_birth_months) == 0)) |> select(QLD) |> sum()
#
# # getting the possible combinations once nirsevimab distribution has started
# # filling in the 0 counts
#
# dose_df_QLD <- full_join(unique(QLD_model_2_data[,c("cohort_birth_months", "month", "year")]) |> subset(year >= 2024 & month >= 1),
#                          dose_df_QLD |> subset(year <= max(QLD_model_2_data$year))) |>
#   mutate(start_month = as.Date(ifelse(is.na(start_month), lubridate::make_date(year, month, 1), start_month)),
#          QLD = ifelse(is.na(QLD), 0, QLD))
#
# sum(dose_df_QLD$QLD)
#
# saveRDS(dose_df_QLD, file = "dose_df_QLD.rds")
#
# # fixing zero nirsevimab prior to dose data times
# doses_on_inds <- sort(which(QLD_model_2_data$rsv_start_month >= min(dose_df_QLD$start_month)))
#
# N_doses_on_inds <- length(doses_on_inds)
# N_doses <- nrow(dose_df_QLD)
#
# # matrix calculation
# doses_mat <- matrix(data = 0, nrow = N_doses_on_inds, ncol = N_doses)
#
# QLD_model_2_data_doses_ind <- lapply(1:nrow(QLD_model_2_data[doses_on_inds,]),
#                                       function(i){
#                                         return(
#                                           which(dose_df_QLD$cohort_birth_months == QLD_model_2_data[doses_on_inds[i], "cohort_birth_months"][[1]] &
#                                                   dose_df_QLD$start_month >= QLD_model_2_data[doses_on_inds[i], "min_date_nirsevimab"][[1]] &
#                                                   dose_df_QLD$start_month <= QLD_model_2_data[doses_on_inds[i], "rsv_start_month"][[1]]
#                                                 )
#                                         )
#                                       })
#
# for(i in 1:N_doses_on_inds){
#   doses_mat[i, QLD_model_2_data_doses_ind[[i]]] <- 1
# }
#
# N_doses_mat_one = sum(doses_mat)
#
# dose_check <- rep(NA, nrow(QLD_model_2_data))
# dose_check[doses_on_inds] <- as.vector(doses_mat %*% dose_df_QLD$QLD)
# unique(dose_check)
#
# QLD_model_2_data$age_rsv_months <- QLD_model_2_data$age_rsv_months + 0.01
#
# saveRDS(QLD_model_2_data, "QLD_model_2_data.rds")
#
# QLD_model_2_data_in <- list(
#   "N" = nrow(QLD_model_2_data),
#   "N_months" = 12,
#   "N_years" = length(unique(QLD_model_2_data$year)),
#   "N_cohort_births" = length(unique(QLD_model_2_data$cohort_birth_months)),
#   "N_cohort_month_years" = length(unique(QLD_model_2_data$cohort_month_year)),
#   "y" = QLD_model_2_data$inc,
#   "sample_cohort_births" = QLD_births$births,
#   "months" = QLD_model_2_data$month,
#   "ind_years" = QLD_model_2_data$ind_years,
#   "ind_cohort_births" = QLD_model_2_data$ind_cohort_births,
#   "ind_cohort_month_years" = QLD_model_2_data$ind_cohort_month_years,
#   "age_months" = QLD_model_2_data$age_rsv_months,
#   "offset_months" = as.vector(QLD_model_2_data$offset_month),
#   "N_doses" = N_doses,
#   "dose_data" = dose_df_QLD$QLD,
#   "N_doses_on_inds" = N_doses_on_inds,
#   "doses_on_inds" = doses_on_inds,
#   "doses_mat" = doses_mat,
#   "N_doses_mat_one" = N_doses_mat_one,
#   "prior_in" = 0.5,
#   "sd_y_doses" = 1,#sd(dose_check)
#   "sd_age_months" = sd(QLD_model_2_data$age_rsv_months),
#   "prior_predictive_check" = 0
#   )
#
# saveRDS(QLD_model_2_data_in, "QLD_model_2_data_in.rds")
