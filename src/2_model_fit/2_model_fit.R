library(cmdstanr)
library(tidyverse)

orderly::orderly_dependency("1_data_cleaning", "latest()",
                            c("QLD_model_1_data_in.rds",
                              "ACT_model_1_data_in.rds",
                              "NSW_model_1_data_in.rds",
                              "model_3_data_in.rds",
                              "model_2_data_ACT_in.rds",
                              "model_2_data_NSW_in.rds",
                              "model_2_data_QLD_in.rds"))

QLD_model_1_data_in <- readRDS(file = "QLD_model_1_data_in.rds")
ACT_model_1_data_in <- readRDS(file = "ACT_model_1_data_in.rds")
NSW_model_1_data_in <- readRDS(file = "NSW_model_1_data_in.rds")

model_3_data_in <- readRDS(file = "model_3_data_in.rds")
model_2_data_ACT_in <- readRDS(file = "model_2_data_ACT_in.rds")
model_2_data_NSW_in <- readRDS(file = "model_2_data_NSW_in.rds")
model_2_data_QLD_in <- readRDS(file = "model_2_data_QLD_in.rds")

n_iter <- 5000

###################
##### model 1 #####
###################

model_1 <- cmdstan_model(stan_file = "poisson_glm.stan")

fit_1 <- model_1$sample(data = QLD_model_1_data_in,
                        seed = 123,
                        iter_sampling = round(n_iter/2),
                        iter_warmup = round(n_iter/2),
                        chains = 4,
                        init = 2,
                        parallel_chains = 4)

fit_1$save_object(file = "fit_1.rds")

fit_1_ACT <- model_1$sample(data = ACT_model_1_data_in,
                            seed = 123,
                            iter_sampling = round(n_iter/2),
                            iter_warmup = round(n_iter/2),
                            chains = 4,
                            init = 2,
                            parallel_chains = 4)

fit_1_ACT$save_object(file = "fit_1_ACT.rds")

fit_1_NSW <- model_1$sample(data = NSW_model_1_data_in,
                            seed = 123,
                            iter_sampling = round(n_iter/2),
                            iter_warmup = round(n_iter/2),
                            chains = 4,
                            init = 2,
                            parallel_chains = 4)

fit_1_NSW$save_object(file = "fit_1_NSW.rds")

###################
##### model 2 #####
###################

model_2 <- cmdstan_model(stan_file = "poisson_glm_cov_individual_state.stan")

fit_2_ACT <- model_2$sample(data = model_2_data_ACT_in,
                            seed = 123,
                            iter_sampling = round(n_iter/2),
                            iter_warmup = round(n_iter/2),
                            init = 0,
                            chains = 4,
                            parallel_chains = 4,
                            max_treedepth = 12,
                            adapt_delta = 0.85)

fit_2_ACT$save_object(file = "fit_2_ACT.rds")

fit_2_NSW <- model_2$sample(data = model_2_data_NSW_in,
                            seed = 123,
                            iter_sampling = round(n_iter/2),
                            iter_warmup = round(n_iter/2),
                            init = 0,
                            chains = 4,
                            parallel_chains = 4,
                            max_treedepth = 12,
                            adapt_delta = 0.85)

fit_2_NSW$save_object(file = "fit_2_NSW.rds")

fit_2_QLD <- model_2$sample(data = model_2_data_QLD_in,
                            seed = 123,
                            iter_sampling = round(n_iter/2),
                            iter_warmup = round(n_iter/2),
                            init = 0,
                            chains = 4,
                            parallel_chains = 4,
                            max_treedepth = 12,
                            adapt_delta = 0.85)

fit_2_QLD$save_object(file = "fit_2_QLD.rds")

###################
##### model 3 #####
###################

model_3 <- cmdstan_model(stan_file = "poisson_glm_cov_state.stan")

fit_3 <- model_3$sample(data = model_3_data_in,
                            seed = 123,
                            iter_sampling = round(n_iter/2),
                            iter_warmup = round(n_iter/2),
                            init = 0,
                            chains = 4,
                            parallel_chains = 4,
                            max_treedepth = 12,
                            adapt_delta = 0.85)

fit_3$save_object(file = "fit_3.rds")

