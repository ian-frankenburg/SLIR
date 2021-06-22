functions {
    real[] sir(real t, real[] y, real[] theta,
             real[] x_r, int[] x_i) {
      real dS_dt;
      real dI_dt;
      real dR_dt;
      real S = y[1];
      real I = y[2];
      real R = y[3];
      int N = x_i[1];
      real beta = theta[1];
      real gamma = theta[2];
      dS_dt = -beta  * (S/N) * I;
      dI_dt = beta * (S/N) * I - gamma * I;
      dR_dt = gamma * I;
      return {dS_dt, dI_dt, dR_dt};
    }
}
data {
  int n_obs;
  int cases[n_obs];
  int N;
  int t0;
}
transformed data {
  real x_r[1] = {1.0};
  int x_i[1] = {N};
  real ts[n_obs];
  real y0[3];
  for(i in 1:n_obs){ts[i] = t0+i;}
  y0[1] = N-cases[1];
  y0[2] = cases[1];
  y0[3] = 0;
}
parameters {
  real<lower=0> R0;
  real<lower=0> phi_inv;
  real<lower=0, upper=1> gamma;
}
transformed parameters{
  real ode[n_obs, 3] = rep_array(0.0, n_obs, 3);
  real phi = 1 / phi_inv; 
  ode = integrate_ode_rk45(sir, y0, t0, ts, {R0 * gamma, gamma}, x_r, x_i);
}
model {
  phi_inv ~ inv_gamma(1,1);
  gamma ~ uniform(0,1);
  R0 ~ lognormal(0,1);
  cases ~ neg_binomial_2(to_vector(ode[,2]) + 1e-5, phi);
}
generated quantities {
  real pred_cases[n_obs];
  real y_pred[n_obs, 3] = rep_array(0.0, n_obs, 3);
  y_pred = integrate_ode_rk45(sir, y0, t0, ts, {R0 * gamma, gamma}, x_r, x_i);
  pred_cases = neg_binomial_2_rng(to_vector(y_pred[,2]) + 1e-5, phi);
}
