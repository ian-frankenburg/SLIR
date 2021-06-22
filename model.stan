functions {
    real[] slir(real t, real[] y, real[] theta,
             real[] x_r, int[] x_i) {
      real dS_dt;
      real dI_dt;
      real dL_dt;
      real dR_dt;
      real S = y[1];
      real L = y[2];
      real I = y[3];
      real R = y[4];
      int N = x_i[1];
      real gamma = theta[2];
      real beta = gamma*theta[1];
      real a = theta[3];
      real b = theta[4];
      dS_dt = -beta  * (S/N) * I - a * S + b * L;
      dL_dt =  a * S - b * L;
      dI_dt = beta * (S/N) * I - gamma * I;
      dR_dt = gamma * I;
      return {dS_dt,dL_dt,dI_dt,dR_dt};
    }
}
data {
  int n_obs;
  int cases[n_obs];
  int N;
  real movement_data[n_obs];
  int t0;
}
transformed data {
  real x_r[1] = {1.0};
  int x_i[1] = {N};
  real y0[4];
  real ts[n_obs];
  real ts_future[n_obs];
  for(i in 1:n_obs){ts[i] = t0+i;}
  for(i in 1:n_obs){ts_future[i] = t0+i;}
  y0[1] = N-cases[1]-0;
  y0[2] =0;
  y0[3] = cases[1];
  y0[4] = 0;
}
parameters {
  real<lower=0> R0;
  real<lower=0> phi_inv;
  real<lower=0> gamma;
  real<lower=0, upper=1> a;
  real<lower=0, upper=1> b;
  real<lower=0> lambda;
}
transformed parameters{
  real ode[n_obs,4] = rep_array(0.0, n_obs,4);
  real phi=1/phi_inv;
  vector[n_obs] R_effective;
  ode = integrate_ode_rk45(slir, y0, t0, ts, {R0, gamma, a, b}, x_r, x_i);
  R_effective = R0*to_vector(ode[,1])/N;
}
model {
  phi_inv ~ inv_gamma(1,1);
  lambda ~ inv_gamma(.1,.1);
  gamma ~ lognormal(0,1);
  a ~ beta(1,5);
  b ~ uniform(0,1);
  R0 ~ lognormal(0,1);
  cases ~ neg_binomial_2(to_vector(ode[,3])+1e-5, phi);
  movement_data ~ beta(lambda*to_vector(ode[,2])/N, lambda*(1-to_vector(ode[,2])/N));
}
generated quantities {
  real pred_cases[n_obs];
  real pred_cases_prior[n_obs];
  real pred_movement[n_obs];
  real y_pred[n_obs,4] = rep_array(0.0, n_obs, 4);
  real y_pred_prior[n_obs,4] = rep_array(0.0, n_obs, 4);;
  y_pred = integrate_ode_rk45(slir, y0, t0, ts_future, {R0, gamma,a, b}, x_r, x_i);
  pred_cases = neg_binomial_2_rng(to_vector(y_pred[,3])+1e-5, phi);
  pred_movement =  beta_rng(lambda*to_vector(y_pred[,2])/N, lambda*(1-to_vector(y_pred[,2])/N));
}
