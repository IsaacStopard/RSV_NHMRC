library(rstan)
library(tidyverse)

orderly::orderly_dependency("1_data_cleaning", "latest()",
                            c("QLD_model_1_data_in.rds",
                              "QLD_model_2_data_in.rds",
                              "QLD_model_1_data.rds",
                              "QLD_model_2_data.rds"))

QLD_model_1_data_in <- readRDS(file = "QLD_model_1_data_in.rds")

QLD_model_2_data_in <- readRDS(file = "QLD_model_2_data_in.rds")

n_iter <- 5000

###################
##### model 1 #####
###################

model_1 <- stan_model(file = "poisson_glm.stan")

set.seed(123)

fit_1 <- sampling(model_1, data = QLD_model_1_data_in, iter = n_iter, warmup = round(n_iter/2), chains = 4, cores = 4)

saveRDS(fit_1, file = "fit_1.rds")

###################
##### model 2 #####
###################

model_2 <- stan_model(file = "poisson_glm_cov.stan")

set.seed(123)

fit_2 <- sampling(model_2, data = QLD_model_2_data_in, iter = n_iter, warmup = round(n_iter/2), init = 0, chains = 4, cores = 4, control = list(max_treedepth = 15, adapt_delta = 0.85))

saveRDS(fit_2, file = "fit_2.rds")
