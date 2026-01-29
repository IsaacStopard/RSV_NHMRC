//
functions{

  vector lp_fun(int N, real intercept, real coef_inc_c, real coef_treat, vector mu_c, vector treatment, vector offset_weeks, int[] ind_week_year, vector coef_week_year){

    vector[N] lp_out = log(offset_weeks) + intercept + coef_inc_c * mu_c ./ offset_weeks + coef_treat * treatment + coef_week_year[ind_week_year];

    return lp_out;
  }
}

data {
  int<lower=0> N;
  int<lower=0> y_t[N];
  int<lower=0> y_c[N];
  vector<lower=0>[N] treatment;
  int<lower=1> ind_week_year[N];

  int<lower=0> N_week_year;
  vector<lower=0>[N] offset_weeks;

  int<lower=0> N_gq;
  vector<lower=0>[N_gq] mu_c_gq;
  vector<lower=0>[N_gq] treatment_gq;
  vector<lower=0>[N_gq] offset_weeks_gq;

  real prior_in;
  real<lower=0> sd_y_c;
}


parameters {

  real<lower=0> mu_c_prior_shape;
  real<lower=0> mu_c_prior_rate;

  vector<lower=0>[N] mu_c;
  real intercept;
  real coef_inc_c;
  real coef_treat;

  // hierarchical week parameters
  vector[N_week_year] coef_week_year;
  real<lower=0> sigma_week_year;
}

transformed parameters{

  vector[N] lp_t = lp_fun(N, intercept, coef_inc_c, coef_treat, mu_c, treatment, offset_weeks, ind_week_year, coef_week_year);

}

model{
  y_c ~ poisson(mu_c);
  y_t ~ poisson_log(lp_t);

  // measurement error prior
  mu_c_prior_shape ~ gamma(2, 1);
  mu_c_prior_rate ~ gamma(2, 4);
  mu_c ~ gamma(mu_c_prior_shape, mu_c_prior_rate);

  // regression priors
  intercept ~ normal(0, prior_in);
  coef_inc_c ~ normal(0, prior_in / sd_y_c);
  coef_treat ~ normal(0, prior_in);

  coef_week_year ~ normal(0, sigma_week_year);
  sigma_week_year ~ exponential(10);
}

generated quantities{

  real<lower=0> irr_pos_dis[N];
  int<lower=0> y_t_pos_dis[N];

  real<lower=0> irr_treat;
  vector[N_gq] lp_treat_gq;
  vector[N_gq] lp_untreat_gq;
  vector[N_gq] mu_treat_gq;
  vector[N_gq] mu_untreat_gq;

  real intercept_prior = normal_rng(0, prior_in);
  real coef_inc_c_prior = normal_rng(0, prior_in / sd_y_c);
  real coef_treat_prior = normal_rng(0, prior_in);

  vector[N] lp_t_prior;
  int<lower=0> y_t_prior_dis[N];

  vector[1] week_year_prior;
  real sigma_week_year_prior;
  vector[1] week_pos_r;

  sigma_week_year_prior = exponential_rng(10);
  week_year_prior[1] = normal_rng(0, sigma_week_year_prior);
  week_pos_r[1] = normal_rng(0, sigma_week_year);

  irr_treat = exp(coef_treat);
  lp_treat_gq = lp_fun(N_gq, intercept, coef_inc_c, coef_treat, mu_c_gq, rep_vector(1, N_gq), offset_weeks_gq, rep_array(1, N_gq), week_pos_r);
  lp_untreat_gq = lp_fun(N_gq, intercept, coef_inc_c, coef_treat, mu_c_gq, rep_vector(0, N_gq), offset_weeks_gq, rep_array(1, N_gq), week_pos_r);
  mu_treat_gq = exp(lp_treat_gq);
  mu_untreat_gq = exp(lp_untreat_gq);

  lp_t_prior = lp_fun(N, intercept_prior, coef_inc_c_prior, coef_treat_prior, mu_c, treatment, offset_weeks, rep_array(1, N), week_year_prior);
  y_t_pos_dis = poisson_log_rng(lp_t);
  y_t_prior_dis = poisson_log_rng(lp_t_prior);

  for(i in 1:N){
    if (y_c[i] > 0) {
       irr_pos_dis[i] = (y_t_pos_dis[i] * 1.0) / (y_c[i] * 1.0);
    } else {
       irr_pos_dis[i] = positive_infinity();
    }
  }
}
