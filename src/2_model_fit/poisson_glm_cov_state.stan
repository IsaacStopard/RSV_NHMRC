//
functions{
 //vector lp_fun(int N, int N_months, vector log_offset_months, vector log_offset_cohort_births, int[] ind_years, vector sin_months, vector cos_months, int[] ind_cohort_births,
 //              int[] ind_cohort_month_years, vector age_months, vector age_months_sq, int[] doses_on_inds, vector doses, real intercept, vector coef_year, real coef_month_s, real coef_month_c,
 //              real coef_age_l, real coef_age_s, real coef_doses, vector offset_cohort_births, vector coef_cohort_births, vector coef_cohort_month_years){


    // outcome monthly notifications

  //  vector[N] lp_mu = log_offset_months +
  //                    log_offset_cohort_births[ind_cohort_births] +
  //                    intercept +
  //                    coef_year[ind_years] +
  //                    coef_month_s * sin_months +
  //                    coef_month_c * cos_months +
  //                    (coef_age_l * age_months) +
  //                    (coef_age_s * age_months_sq) +
  //                    coef_cohort_births[ind_cohort_births] +
  //                    coef_cohort_month_years[ind_cohort_month_years];


  //  lp_mu[doses_on_inds] += (doses ./ (offset_cohort_births[ind_cohort_births[doses_on_inds]]) * coef_doses);

  //  return lp_mu;

  //  }
}

data {
  int<lower=0> N;
  int<lower=0> N_months;
  int<lower=0> N_years;
  int<lower=0> N_cohort_births;
  int<lower=0> N_cohort_month_years;

  int<lower=0> y[N];

  int<lower=0> sample_cohort_births[N_cohort_births];

  vector<lower=1, upper=N_months>[N] months;
  int<lower=1, upper=N_years> ind_years[N];
  int<lower=1, upper=N_cohort_births> ind_cohort_births[N];
  int<lower=1, upper=N_cohort_month_years> ind_cohort_month_years[N];

  vector<lower=0>[N] age_months;

  vector<lower=0>[N] offset_months;

  int<lower=0> N_doses;
  int<lower=0> dose_data[N_doses];

  int<lower=0> N_doses_on_inds;
  int<lower=0> doses_on_inds[N_doses_on_inds];
  matrix<lower=0,upper=1>[N_doses_on_inds, N_doses] doses_mat;
  int<lower=0> N_doses_mat_one;

  real prior_in;
  real<lower=0> sd_y_doses;
  real<lower=0> sd_age_months;

  int prior_predictive_check;
  }

transformed data{
  vector[N] log_offset_months = log(offset_months);
  vector[N] age_months_sq = age_months^2;
  vector[N] age_months_log = log(age_months);
  vector[N] sin_months = sin(months * 2 * pi() / N_months);
  vector[N] cos_months = cos(months * 2 * pi() / N_months);

  vector[N_doses_mat_one] csr_w = csr_extract_w(doses_mat);
  int csr_v[N_doses_mat_one] = csr_extract_v(doses_mat);
  int csr_u[N_doses_on_inds + 1] = csr_extract_u(doses_mat);

}

parameters {

  real intercept;
  vector[N_years-1] coef_year_raw;
  real coef_age_s;
  real coef_age_l;
  real coef_doses;
  real coef_month_s;
  real coef_month_c;

  real<lower=0> mu_doses_prior_shape;
  real<lower=0> mu_doses_prior_mean;
  vector<lower=0>[N_doses] mu_doses;

  real<lower=0> offset_cohort_births_shape;
  real<lower=0> offset_cohort_births_mean;

  // all cohorts are present so the births must be greater than 1
  vector<lower=1.0>[N_cohort_births] offset_cohort_births;

  // hierarchical month parameters
  // vector[N_months] coef_month;
  // real mu_month;
  // real<lower=0> sigma_month;

  vector[N_cohort_births] coef_cohort_births;
  real<lower=0> sigma_cohort_births;

  vector[N_cohort_month_years] coef_cohort_month_years;
  real<lower=0> sigma_cohort_month_years;
}

transformed parameters{

  vector[N_years] coef_year = append_row(0, coef_year_raw);
  real<lower=0> mu_doses_prior_rate = mu_doses_prior_shape / mu_doses_prior_mean;
  real<lower=0> offset_cohort_births_rate = offset_cohort_births_shape / offset_cohort_births_mean;

}

model {

  // calculations
  vector[N_doses_on_inds] doses = csr_matrix_times_vector(N_doses_on_inds, N_doses, csr_w, csr_v, csr_u, mu_doses);

  vector[N_cohort_births] log_offset_cohort_births = log(offset_cohort_births);

  vector[N] lp_mu = log_offset_months +
                    log_offset_cohort_births[ind_cohort_births] +
                    intercept +
                    coef_year[ind_years] +
                    coef_month_s * sin_months +
                    coef_month_c * cos_months +
                    (coef_age_s * age_months_log - coef_age_l * age_months) +
                    coef_cohort_births[ind_cohort_births] +
                    coef_cohort_month_years[ind_cohort_month_years];

  lp_mu[doses_on_inds] += (doses ./ (offset_cohort_births[ind_cohort_births[doses_on_inds]])) * coef_doses;

  // likelihoods
  dose_data ~ poisson(mu_doses);
  sample_cohort_births ~ poisson(offset_cohort_births);
  y ~ poisson_log(lp_mu);

  // priors
  // measurement error prior
  mu_doses_prior_shape ~ gamma(10, 10.0/100.0);
  mu_doses_prior_mean ~ gamma(10, 10.0/5.0);
  mu_doses ~ gamma(mu_doses_prior_shape, mu_doses_prior_rate);

  offset_cohort_births_shape ~ gamma(10, 10.0/100.0);
  offset_cohort_births_mean ~ gamma(100, 100.0/5000.0);
  offset_cohort_births ~ gamma(offset_cohort_births_shape, offset_cohort_births_rate);

  // regression priors
  intercept ~ normal(0, prior_in);
  coef_doses ~ normal(0, prior_in);
  coef_year_raw ~ normal(0, prior_in);
  coef_age_s ~ normal(0, prior_in / sd_age_months);
  coef_age_l ~ normal(0, prior_in / sd_age_months);
  coef_month_s ~ normal(0, prior_in);
  coef_month_c ~ normal(0, prior_in);

  //coef_month ~ normal(mu_month, sigma_month);
  //mu_month ~ normal(0, 0.1);
  //sigma_month ~ exponential(10);

  coef_cohort_births ~ normal(0, sigma_cohort_births);
  sigma_cohort_births ~ exponential(10);

  coef_cohort_month_years ~ normal(0, sigma_cohort_month_years);
  sigma_cohort_month_years ~ exponential(10);
}

generated quantities{

  int<lower=0> y_pos_dis[N];

  real amplitude = sqrt(coef_month_c^2 + coef_month_s^2);
  real phase_shift_cos_months = -atan2(coef_month_s, coef_month_c) * N_months / (2*pi());

  vector[N] lp_mu_gq;
  vector[N_cohort_births] log_offset_cohort_births_gq = log(offset_cohort_births);
  vector[N_doses_on_inds] doses_gq = csr_matrix_times_vector(N_doses_on_inds, N_doses, csr_w, csr_v, csr_u, mu_doses);

  lp_mu_gq = log_offset_months +
             log_offset_cohort_births_gq[ind_cohort_births] +
             intercept +
             coef_year[ind_years] +
             coef_month_s * sin_months +
             coef_month_c * cos_months +
             (coef_age_s * age_months_log - coef_age_l * age_months) +
             coef_cohort_births[ind_cohort_births] +
             coef_cohort_month_years[ind_cohort_month_years];

  lp_mu_gq[doses_on_inds] += (doses_gq ./ (offset_cohort_births[ind_cohort_births[doses_on_inds]])) * coef_doses;

  y_pos_dis = poisson_log_rng(lp_mu_gq);

  if(prior_predictive_check == 1){
    real intercept_prior = normal_rng(0, prior_in);
    real coef_doses_prior = normal_rng(0, prior_in);
    vector[(N_years-1)] coef_year_raw_prior;
    vector[N_years] coef_year_prior;
    real coef_age_s_prior = normal_rng(0, prior_in / sd_age_months);
    real coef_age_l_prior = normal_rng(0, prior_in / sd_age_months);
    vector[N] lp_prior;
    int y_prior_dis[N];
    real coef_month_s_prior = normal_rng(0, prior_in);
    real coef_month_c_prior = normal_rng(0, prior_in);
    vector[N_cohort_births] coef_cohort_births_prior;
    real sigma_cohort_births_prior = exponential_rng(10);
    vector[N_cohort_month_years] coef_cohort_month_years_prior;
    real sigma_cohort_month_years_prior = exponential_rng(10);

    for(i in 1:N_cohort_births){
      coef_cohort_births_prior[i] = normal_rng(0, sigma_cohort_births_prior);
      }

      for(i in 1:N_cohort_month_years){
        coef_cohort_month_years_prior[i] = normal_rng(0, sigma_cohort_month_years_prior);
        }

        for(i in 1:(N_years-1)){
          coef_year_raw_prior[i] = normal_rng(0, prior_in);
          }

    coef_year_prior = append_row(0, coef_year_raw_prior);

    lp_prior = log_offset_months +
               log_offset_cohort_births_gq[ind_cohort_births] +
               intercept_prior +
               coef_year_prior[ind_years] +
               coef_month_s_prior * sin_months +
               coef_month_c_prior * cos_months +
               (coef_age_s_prior * age_months_log - coef_age_l_prior * age_months) +
               coef_cohort_births_prior[ind_cohort_births] +
               coef_cohort_month_years_prior[ind_cohort_month_years];

    lp_prior[doses_on_inds] += (doses_gq ./ offset_cohort_births[ind_cohort_births[doses_on_inds]]) * coef_doses_prior;
    y_prior_dis = poisson_log_rng(lp_prior);
  }

}
