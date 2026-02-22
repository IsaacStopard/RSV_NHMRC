library(tidyverse)
library(readxl)
library(lubridate)
library(stringr)
library(patchwork)
# library(readabs) # Sys.setenv(R_READABS_PATH = "/Users/z5238812/Documents/GitHub/RSV_NHMRC/src/1_data_cleaning")

###########################
##### population data #####
###########################

# births data

# https://www.data.qld.gov.au/dataset/births-by-month
# read_births_csv <- function(year){
#   read.csv(file = paste0("QLD-births-by-month-",year,".csv")) |>
#     rename("month" = Month, "births" = Transactions) |>
#     mutate(year = year, cohort_birth_months = zoo::as.yearmon(make_date(year, month, 1)))
# }
#
# QLD_births <- rbind(read_births_csv(2018), read_births_csv(2019), read_births_csv(2020), read_births_csv(2021), read_births_csv(2022), read_births_csv(2023), read_births_csv(2024))
#
# saveRDS(QLD_births, file = "QLD_births.rds")

births <- read.csv(file = "ABS_BIRTHS_MONTH_OCCURRENCE.csv") |>
  select(Region, Month.of.occurence, TIME_PERIOD, OBS_VALUE) |> rename("state" = Region, "month" = Month.of.occurence, "year" = TIME_PERIOD, "births" = OBS_VALUE) |>
  mutate(state = case_when(state == "New South Wales" ~ "NSW",
                           state == "Australian Capital Territory" ~ "ACT",
                           state == "Queensland" ~ "QLD"),
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
         births = ifelse(month == 12 & year == 2024 | month == 11 & year == 2024, NA, births)
    ) |> na.omit()

# fill in final month with mean value

birth_fit <- glm(births ~ factor(month) + factor(year) + state, family = "poisson", data = births)

birth_pred <- round(predict(birth_fit, newdata = data.frame(month = c(rep(12, 3), rep(11, 3)), year = rep(2024, 6), state = rep(c("ACT", "NSW", "QLD"), 2)), type = "response"), digits = 0)

births <- rbind(births,
                data.frame(month = c(rep(12, 3), rep(11, 3)), year = rep(2024, 6), state = rep(c("ACT", "NSW", "QLD"), 2), births = birth_pred) |>
                      mutate(cohort_birth_months = zoo::as.yearmon(make_date(year, month, 1)),
                             cohort = paste0(state, " ", cohort_birth_months)
                             )
                )

QLD_births <- subset(births, state == "QLD" & year >= 2018) |> arrange(year, month)

saveRDS(QLD_births, file = "QLD_births.rds")

#########################
##### coverage data #####
#########################

dose_df <- read_excel("F Voight extract 20251114.xlsx", sheet = "clean_data") |>
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
  mutate(total = ACT + NSW + QLD,
         start_month = lubridate::make_date(year, month, 1)) |> arrange(year, month, cohort_birth_months)

saveRDS(dose_df, file = "dose_df.rds")
saveRDS(dose_df_sum_age_months, "dose_df_sum_age_months.rds")

#######################
##### Queensland ######
#######################

df_QLD <- read_excel("QLD data_analysis_original.xlsx")

total <- nrow(df_QLD)

colnames(df_QLD) <- c("case_number", "date_of_rsv_episode", "week_of_rsv_episode", "year_of_rsv_episode", "week_and_year_of_rsv_episode", "week_and_year_of_rsv_episode_by_age_group",
                      "week_of_birth", "year_of_birth", "sex", "hospitalization",
                      "statistical_area_level_2_of_residence",
                      "age_in_weeks_at_time_of_rsv_episode", "age_in_months_at_time_of_rsv_episode",
                      "age_in_weeks_at_time_of_program_rollout", "age_in_months_at_time_of_program_rollout",
                      "age_6_months_groups_at_time_of_rsv_episode", "binary_age_groups_less_than_6m_vs_geq_6m_at_time_of_rsv_episode", "running_count_of_weekly_incidence", "cumulative_count_of_weekly_incidence", "running_count_of_weekly_incidence_by_age_group",
                      "cumulative_count_of_weekly_incidence_by_age_group")

min_date <- min(df_QLD$date_of_rsv_episode)
max_date <- max(df_QLD$date_of_rsv_episode)

df_QLD <- df_QLD |> mutate(week = week(date_of_rsv_episode),
                           year = year(date_of_rsv_episode),
                           month = month(as.Date(date_of_rsv_episode, format = "%m %Y")),
                           week_year = paste0(week,"-",year),

                           ndays_per_month = days_in_month(make_date(year, month, day = 1)),

                           treatment = if_else(week >= 15 & year >= 2024, 1, 0),
                           cohort_birth_start_week = make_date(year_of_birth, 1, 1) + weeks(week_of_birth - 1),
                           cohort_birth_start_month = floor_date(cohort_birth_start_week, unit = "month"),
                           cohort_birth_months = zoo::as.yearmon(make_date(year_of_birth, 1, 1) + weeks(week_of_birth - 1)),

                           # model 1
                           age_rsv_months = interval(cohort_birth_start_week, date_of_rsv_episode) %/% months(1),
                           age_group = if_else(age_rsv_months <= 8, "age_l_8", "age_g_8"),
                           ) |>
  arrange(date_of_rsv_episode)

##### model 1 data #####

QLD_model_1_data <- df_QLD |> mutate() |> group_by(week, year, week_year, age_group, treatment) |> summarise(inc = n()) |> ungroup()

sum(QLD_model_1_data$inc) == total

QLD_model_1_data <- left_join(QLD_model_1_data,
                              data.frame("date" = seq(min_date, max_date, by = "days")) |>
                                mutate(week = week(date), year = year(date)) |>
                                group_by(week, year) |> summarise(offset_week = n())
                              )

QLD_model_1_data <- QLD_model_1_data |> left_join(unique(QLD_model_1_data[,c("week", "year", "week_year")]) |>
                                                    arrange(year, week) |> ungroup() |>
                                                    mutate(week_cont = row_number()))

QLD_model_1_data <- QLD_model_1_data |>
  pivot_wider(names_from = age_group, values_from = inc) |>
  replace_na(list("age_g_8" = 0, "age_l_8" = 0)) |>
  rename("inc_less_8m" = age_l_8, "inc_greater_8m" = age_g_8) |>
  mutate(IRR = inc_less_8m/inc_greater_8m)

total == (sum(QLD_model_1_data$inc_greater_8m) + sum(QLD_model_1_data$inc_less_8m)) # checking the total infections are still correct

# checking all weeks and years have data
nrow(unique(df_QLD[,c("week", "year")])) == nrow(expand.grid(week = unique(df_QLD$week), year = unique(df_QLD$year)))

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
                            "prior_in" = 0.75,
                            "sd_y_c" = sd(QLD_model_1_data$inc_greater_8m))

saveRDS(QLD_model_1_data_in, file = "QLD_model_1_data_in.rds")

##### model 2 data #####

nirsevimab_duration_months <- 6

# filling in the zero counts
# assumes that the surveillance system is always active post 01/01/2023
# maximum age with surveillance is assumed to be the maximum age with RSV observed

full_combs <- expand.grid("month" = sort(unique(df_QLD$month)),
                          "year" = sort(unique(df_QLD$year)),
                          "cohort_birth_start_month" = sort(unique(df_QLD$cohort_birth_start_month)))

obs_combs <- df_QLD |> group_by(cohort_birth_start_month, month, year) |> summarise(inc = n()) |> ungroup()

nrow(obs_combs) == nrow(unique(df_QLD[,c("cohort_birth_start_month", "month", "year")]))

QLD_model_2_data <- left_join(full_combs, obs_combs, by = c("cohort_birth_start_month", "month", "year")) |>
  replace_na(list("inc" = 0)) |>
  mutate(min_date_nirsevimab = lubridate::make_date(year, month, 1) %m-% months(nirsevimab_duration_months),
         rsv_start_month = lubridate::make_date(year, month, 1),
         age_rsv_months = interval(cohort_birth_start_month, rsv_start_month) %/% months(1),
         cohort_birth_months = zoo::as.yearmon(cohort_birth_start_month),
         cohort_month_year = paste0(cohort_birth_months,"-", month,"-", year),
         offset_months = days_in_month(rsv_start_month)) |>
  filter(age_rsv_months >= 0 & age_rsv_months <= (6 * 12)) # only includes data for those less than 5 years old

if(sum(QLD_model_2_data$inc) != nrow(df_QLD)){stop("incidence sum in model 2 not correct")}

max(QLD_model_2_data$age_rsv_months)

total == sum(QLD_model_2_data$inc)

QLD_model_2_data <- left_join(QLD_model_2_data,
                              data.frame("year" = unique(QLD_model_2_data$year) |> sort()) |> mutate(ind_years = row_number())
                              )

unique(QLD_model_2_data$cohort_birth_months %in% QLD_births$cohort_birth_months)

unique(QLD_births$cohort_birth_months %in% QLD_model_2_data$cohort_birth_months)

QLD_births <- QLD_births |> mutate(ind_cohort_births = row_number())

QLD_model_2_data <- left_join(QLD_model_2_data,
                              QLD_births[,c("cohort_birth_months", "ind_cohort_births")]
                              )

QLD_model_2_data <- left_join(QLD_model_2_data,
                              data.frame("cohort_month_year" = unique(QLD_model_2_data$cohort_month_year) |> sort()) |> mutate(ind_cohort_month_years = row_number())
)

#QLD_model_2_data <- left_join(QLD_model_2_data,
#                              offset_weeks_df <- data.frame("date" = seq(min(df_QLD$date_of_rsv_episode), max(df_QLD$date_of_rsv_episode), by = "days")) |>
#                                mutate(week = week(date), year = year(date)) |>
#                                group_by(week, year) |> summarise(offset_week = n()/7)
#                              )

# sums the doses in the previous 6 months for each cohort birth
dose_df_QLD <- dose_df_sum_age_months |> filter(cohort_birth_months %in% unique(QLD_model_2_data$cohort_birth_months) & QLD > 0) |> select(-ACT, -NSW, -total)

sum(dose_df_sum_age_months$QLD)
sum(dose_df_QLD$QLD) + dose_df_sum_age_months |> filter((cohort_birth_months %in% unique(QLD_model_2_data$cohort_birth_months) == 0)) |> select(QLD) |> sum()

# getting the possible combinations once nirsevimab distribution has started
# filling in the 0 counts

dose_df_QLD <- full_join(unique(QLD_model_2_data[,c("cohort_birth_months", "month", "year")]) |> subset(year >= 2024 & month >= 1),
                         dose_df_QLD |> subset(year <= max(QLD_model_2_data$year))) |>
  mutate(start_month = as.Date(ifelse(is.na(start_month), lubridate::make_date(year, month, 1), start_month)),
         QLD = ifelse(is.na(QLD), 0, QLD))

sum(dose_df_QLD$QLD)

saveRDS(dose_df_QLD, file = "dose_df_QLD.rds")

# fixing zero nirsevimab prior to dose data times
doses_on_inds <- sort(which(QLD_model_2_data$rsv_start_month >= min(dose_df_QLD$start_month)))

N_doses_on_inds <- length(doses_on_inds)
N_doses <- nrow(dose_df_QLD)

# matrix calculation
doses_mat <- matrix(data = 0, nrow = N_doses_on_inds, ncol = N_doses)

QLD_model_2_data_doses_ind <- lapply(1:nrow(QLD_model_2_data[doses_on_inds,]),
                                      function(i){
                                        return(
                                          which(dose_df_QLD$cohort_birth_months == QLD_model_2_data[doses_on_inds[i], "cohort_birth_months"][[1]] &
                                                  dose_df_QLD$start_month >= QLD_model_2_data[doses_on_inds[i], "min_date_nirsevimab"][[1]] &
                                                  dose_df_QLD$start_month <= QLD_model_2_data[doses_on_inds[i], "rsv_start_month"][[1]]
                                                )
                                        )
                                      })

for(i in 1:N_doses_on_inds){
  doses_mat[i, QLD_model_2_data_doses_ind[[i]]] <- 1
}

N_doses_mat_one = sum(doses_mat)

dose_check <- rep(NA, nrow(QLD_model_2_data))
dose_check[doses_on_inds] <- as.vector(doses_mat %*% dose_df_QLD$QLD)
unique(dose_check)

#
# dose_indices <- unlist(QLD_model_2_data_doses_ind)
# group_lengths <- sapply(QLD_model_2_data_doses_ind, length)
#
# group_starts <- c(1, cumsum(group_lengths) + 1)
# group_starts <- group_starts[1:nrow(QLD_model_2_data[doses_on_inds,])] # remove last element
#
# dose_check <- dose_comp <- rep(NA, nrow(QLD_model_2_data))
#
# dose_check[doses_off_inds] <- dose_comp[doses_off_inds] <- 0
#
# for(i in 1:nrow(QLD_model_2_data[doses_on_inds,])){
#   dose_check[doses_on_inds[i]] <- sum(dose_df_QLD$QLD[QLD_model_2_data_doses_ind[[i]]])
# }
#
# for(i in 1:length(doses_on_inds)){
#   dose_comp[doses_on_inds[i]] = sum(dose_df_QLD$QLD[dose_indices[group_starts[i]:(group_starts[i] + group_lengths[i] - 1)]])
# }

QLD_model_2_data$age_rsv_months <- QLD_model_2_data$age_rsv_months + 0.01

saveRDS(QLD_model_2_data, "QLD_model_2_data.rds")

QLD_model_2_data_in <- list(
  "N" = nrow(QLD_model_2_data),
  "N_months" = 12,
  "N_years" = length(unique(QLD_model_2_data$year)),
  "N_cohort_births" = length(unique(QLD_model_2_data$cohort_birth_months)),
  "N_cohort_month_years" = length(unique(QLD_model_2_data$cohort_month_year)),
  "y" = QLD_model_2_data$inc,
  "sample_cohort_births" = QLD_births$births,
  "months" = QLD_model_2_data$month,
  "ind_years" = QLD_model_2_data$ind_years,
  "ind_cohort_births" = QLD_model_2_data$ind_cohort_births,
  "ind_cohort_month_years" = QLD_model_2_data$ind_cohort_month_years,
  "age_months" = QLD_model_2_data$age_rsv_months,
  "offset_months" = as.vector(QLD_model_2_data$offset_month),
  "N_doses" = N_doses,
  "dose_data" = dose_df_QLD$QLD,
  "N_doses_on_inds" = N_doses_on_inds,
  "doses_on_inds" = doses_on_inds,
  "doses_mat" = doses_mat,
  "N_doses_mat_one" = N_doses_mat_one,
  "prior_in" = 0.5,
  "sd_y_doses" = 1,#sd(dose_check)
  "sd_age_months" = sd(QLD_model_2_data$age_rsv_months),
  "prior_predictive_check" = 0
  )

saveRDS(QLD_model_2_data_in, "QLD_model_2_data_in.rds")

###########################
##### New South Wales #####
###########################

##########################################
##### Australian Capital Territories #####
##########################################

# only includes data for those less than 5 years old
