//
data {
  int<lower=0> N;
  real<lower=0> N_months;
  int<lower=0> N_years;
  int<lower=0> N_cohorts;
  int<lower=0> N_observations;

  array[N_cohorts] int<lower=0> sample_cohort_births;

  array[N] int<lower=0> y;

  vector<lower=1, upper=N_months>[N] months;
  array[N] int<lower=1, upper=N_years> ind_years;
  array[N] int<lower=1, upper=N_cohorts> ind_cohorts;
  array[N] int<lower=1, upper=N_observations> ind_observations;

  vector<lower=0>[N] age_months;

  vector<lower=0>[N] offset_months;

  // raw counts of the number of doses administered each month
  int<lower=0> N_doses;
  array[N_doses] int<lower=0> dose_data;

  int<lower=0, upper=N> N_doses_on_inds;
  array[N_doses_on_inds] int<lower=1, upper=N> doses_on_inds; // indicator of whether Nirsevimab distribution has started
  matrix<lower=0,upper=1>[N_doses_on_inds, N_doses] doses_mat;
  int<lower=0> N_doses_mat_one;

  real prior_in;
  real<lower=0> sd_age_months;

  real<lower=0> prior_in_cohort_births_mean_mean;
  real<lower=0> prior_in_cohort_births_sd_mean;

  real<lower=0> prior_in_mu_doses_mean_mean;
  real<lower=0> prior_in_mu_doses_sd_mean;

  }

transformed data{
  vector[N] log_offset_months = log(offset_months);
  vector[N] age_months_sq = age_months^2;
  vector[N] age_months_log = log(age_months);
  vector[N] sin_months = sin(months * 2 * pi() / N_months);
  vector[N] cos_months = cos(months * 2 * pi() / N_months);

  vector[N_doses_mat_one] csr_w = csr_extract_w(doses_mat);
  array[N_doses_mat_one] int csr_v = csr_extract_v(doses_mat);
  array[N_doses_on_inds + 1] int csr_u = csr_extract_u(doses_mat);

}

parameters {

  // fixed effects
  real intercept;
  vector[N_years-1] coef_year_raw;
  real coef_age_s;
  real coef_age_l;
  real coef_doses;
  real coef_month_s;
  real coef_month_c;

  real<lower=0> mu_doses_prior_mean;
  real<lower=0> mu_doses_prior_sd;

  vector<lower=0>[N_doses] mu_doses;

  real<lower=0> offset_cohort_births_mean;
  real<lower=0> offset_cohort_births_sd;
  // all cohorts are present so the births must be greater than 1
  vector<lower=0>[N_cohorts] offset_cohort_births;

  // hierarchical parameters
  vector[N_cohorts] coef_cohorts_raw;
  real<lower=0> sigma_cohorts;

  vector[N_observations] coef_observations_raw;
  real<lower=0> sigma_observations;
}

transformed parameters{

  vector[N_years] coef_year = append_row(0, coef_year_raw);

  real<lower=0> mu_doses_prior_shape = square(mu_doses_prior_mean) / square(mu_doses_prior_sd);
  real<lower=0> mu_doses_prior_rate = mu_doses_prior_mean / square(mu_doses_prior_sd);

  real<lower=0> offset_cohort_births_shape = square(offset_cohort_births_mean) / square(offset_cohort_births_sd);
  real<lower=0> offset_cohort_births_rate = offset_cohort_births_mean / square(offset_cohort_births_sd);

  vector[N_cohorts] coef_cohorts = coef_cohorts_raw * sigma_cohorts;
  vector[N_observations] coef_observations = coef_observations_raw * sigma_observations;
}

model {

  // calculations
  vector[N_doses_on_inds] doses = csr_matrix_times_vector(N_doses_on_inds, N_doses, csr_w, csr_v, csr_u, mu_doses);

  vector[N_cohorts] log_offset_cohort_births = log(offset_cohort_births);

  vector[N] lp_mu = log_offset_months +
                    log_offset_cohort_births[ind_cohorts] +
                    intercept +
                    coef_year[ind_years] +
                    coef_month_s * sin_months +
                    coef_month_c * cos_months +
                    (coef_age_s * age_months_log - coef_age_l * age_months) +
                    coef_cohorts[ind_cohorts] +
                    coef_observations[ind_observations];

  lp_mu[doses_on_inds] += (doses ./ (offset_cohort_births[ind_cohorts[doses_on_inds]])) * coef_doses;

  // likelihoods
  dose_data ~ poisson(mu_doses);
  sample_cohort_births ~ poisson(offset_cohort_births);
  y ~ poisson_log(lp_mu);

  // priors
  // measurement error priors

  offset_cohort_births_mean ~ gamma(75.0, 75.0 ./ prior_in_cohort_births_mean_mean);
  offset_cohort_births_sd ~ gamma(20.0, 20.0 ./ prior_in_cohort_births_sd_mean);

  mu_doses_prior_mean ~ gamma(75.0, 75.0 ./ prior_in_mu_doses_mean_mean);
  mu_doses_prior_sd ~ gamma(20.0, 20.0 ./ prior_in_mu_doses_sd_mean);

  mu_doses ~ gamma(mu_doses_prior_shape, mu_doses_prior_rate);

  offset_cohort_births ~ gamma(offset_cohort_births_shape, offset_cohort_births_rate);

  // regression priors
  intercept ~ normal(0, prior_in);
  coef_doses ~ normal(0, prior_in);
  coef_year_raw ~ normal(0, prior_in);
  coef_age_s ~ normal(0, prior_in / sd_age_months);
  coef_age_l ~ normal(0, prior_in / sd_age_months);
  coef_month_s ~ normal(0, prior_in);
  coef_month_c ~ normal(0, prior_in);

  coef_cohorts_raw ~ normal(0, 1);
  sigma_cohorts ~ exponential(10);

  coef_observations_raw ~ normal(0, 1);
  sigma_observations ~ exponential(5);
}

generated quantities{

  array[N] int<lower=0> y_pos_dis;

  real amplitude = sqrt(square(coef_month_c) + square(coef_month_s));
  real phase_shift_cos_months = atan2(coef_month_s, coef_month_c) * N_months / (2*pi());

  vector[N] lp_mu_gq;
  vector[N_cohorts] log_offset_cohort_births_gq = log(offset_cohort_births);
  vector[N_doses_on_inds] doses_gq = csr_matrix_times_vector(N_doses_on_inds, N_doses, csr_w, csr_v, csr_u, mu_doses);

  lp_mu_gq = log_offset_months +
             log_offset_cohort_births_gq[ind_cohorts] +
             intercept +
             coef_year[ind_years] +
             coef_month_s * sin_months +
             coef_month_c * cos_months +
             (coef_age_s * age_months_log - coef_age_l * age_months) +
             coef_cohorts[ind_cohorts] +
             coef_observations[ind_observations];

  lp_mu_gq[doses_on_inds] += (doses_gq ./ (offset_cohort_births[ind_cohorts[doses_on_inds]])) * coef_doses;

  y_pos_dis = poisson_log_rng(lp_mu_gq);

}
