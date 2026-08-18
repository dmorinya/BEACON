functions {
  vector sum_to_zero_helmert(vector raw) {
    int K = num_elements(raw) + 1;
    vector[K] z = rep_vector(0.0, K);
    real sum_w = 0;
    for (ii in 1 : (K - 1)) {
      int i = K - ii;
      real n = i;
      real w = raw[i] * inv_sqrt(n * (n + 1));
      sum_w += w;
      z[i] += sum_w;
      z[i + 1] -= w * n;
    }
    return z;
  }
}

data {
  // ---- dimensions and indices ------------------------------------------
  int<lower=1> N;                                   // observations
  int<lower=1> A;                                   // areas (ABS)
  int<lower=2> S;                                   // seasons
  int<lower=1> M;                                   // area-season series
  int<lower=1> W;                                   // weeks per season (52)

  array[N] int<lower=1, upper=A> area;
  array[N] int<lower=1, upper=S> season;
  array[N] int<lower=1, upper=M> series;            // area-season id      
  array[N] int<lower=1, upper=W> week_in_season;    // 1..W, for estimands

  array[M] int<lower=1, upper=N> series_start;
  array[M] int<lower=1> series_len;

  // ---- outcome and covariates ------------------------------------------
  array[N] int<lower=0> y;
  vector[N] log_population;
  vector[N] week;                                   // calendar week index

  // Pivot about which the speed rescaling acts. Setting this to mean(week)
  // rather than leaving it implicitly at zero places the pivot in the middle of
  // the season instead of 26 weeks outside it, which markedly reduces the
  // posterior correlation between shift and log_speed.                 
  real week_pivot;
  vector<lower=0>[N] immunisation_intensity;
  vector<lower=0, upper=1>[N] treated;
  vector[N] ses;
  vector[N] lag_log_y;                              // standardised lag

  // Constants used to rebuild the lag from a simulated count under the
  // dynamic g-formula: L = (log1p(Y) - lag_center) / lag_scale.         
  real lag_center;
  real<lower=0> lag_scale;
  real lag_init;                                    // L in week 1 of a series

  // ---- model configuration ---------------------------------------------
  int<lower=1, upper=6> H;                          // Fourier harmonics   
  real<lower=0> harmonic_decay;                     // prior decay in h   
  int<lower=0, upper=1> intensity_quadratic;        // curvature in dose  

  // ---- prior configuration (defaults in run_beacon_v2.R) ---------------
  real prior_mean_beta_intensity;                   // 0 by default       
  real<lower=0> prior_sd_beta_intensity;
  real<lower=0> prior_sd_shift;                     // weeks
  real<lower=0> prior_sd_speed;                     // log scale
  real<lower=0> prior_sd_ascertainment;             // 0 disables         

  // ---- output switches (control draw file size) ------------------------
  // prior_only = 1 drops the likelihood, so the model samples the prior. Used
  // for the prior predictive check; y_rep then draws from the prior predictive.
  int<lower=0, upper=1> prior_only;
  int<lower=0, upper=1> compute_log_lik;                                // 
  int<lower=0, upper=1> compute_y_rep;
  int<lower=0, upper=1> dynamic_gformula;                               // 
}

transformed data {
  real omega = 2 * pi() / 52.0;
  real stz_shift_scale = inv_sqrt(1 - inv(S));   // unit marginals after
  vector[H] harmonic_scale;
  real eta_cap = 15.0;                           // rng overflow guard   

  // Grid over one period of latent epidemic time, used in generated quantities
  // to locate the maximum of the baseline curve.                        
  int n_grid = 521;
  vector[n_grid] tau_grid;

  for (h in 1 : H) {
    harmonic_scale[h] = pow(h, -harmonic_decay);
  }
  for (g in 1 : n_grid) {
    tau_grid[g] = 52.0 * (g - 1) / n_grid;
  }
}

parameters {
  real alpha;

  // hierarchical intercepts (non-centred)
  vector[A] z_area;
  vector[S] z_season;
  vector[M] z_series;                              // area-season effect 
  real<lower=0> sigma_area;
  real<lower=0> sigma_season;
  real<lower=0> sigma_series;

  // epidemic clock: raw vectors of length S-1 mapped to exact zero-sum  
  vector[S - 1] shift_raw;
  vector[S - 1] log_speed_raw;
  real<lower=0> sigma_shift;
  real<lower=0> sigma_speed;

  // baseline epidemic curve
  vector[H] a_raw;
  vector[H] b_raw;
  real<lower=0> sigma_f;

  // regression coefficients
  real beta_lag;
  real beta_ses;
  real beta_intensity;                             // cycle-average log RR
  real beta_intensity_cos;                         // free-phase modulation
  real beta_intensity_sin;
  real beta_intensity_sq;                          // used iff flag = 1    

  // ascertainment sensitivity parameter, not identified by the data      
  real lambda_raw;

  real<lower=1e-5> reciprocal_phi;                    // phi = 1/reciprocal  
}

transformed parameters {
  vector[A] area_re = sigma_area * z_area;
  vector[S] season_re = sigma_season * z_season;
  vector[M] series_re = sigma_series * z_series;

  vector[S] shift = sigma_shift * sum_to_zero_helmert(shift_raw);
  vector[S] log_speed = sigma_speed * sum_to_zero_helmert(log_speed_raw);

  vector[H] a_f = sigma_f * (a_raw .* harmonic_scale);
  vector[H] b_f = sigma_f * (b_raw .* harmonic_scale);

  real lambda = prior_sd_ascertainment * lambda_raw;
  real phi = inv(reciprocal_phi + 1e-7);

  // Latent epidemic time tau is deliberately NOT declared here. It is an
  // N-vector, and transformed parameters are written to the draws file on every
  // iteration; at N in the tens of thousands that dominates output size. It is
  // recomputed locally in the model and generated quantities blocks instead.
}

model {
  // Latent epidemic time. Strictly increasing in t for every season.
  vector[N] tau = (week - week_pivot - shift[season]) .* exp(log_speed[season]);

  // ---- priors -----------------------------------------------------------
  alpha ~ normal(-6, 2);

  z_area ~ std_normal();
  z_season ~ std_normal();
  z_series ~ std_normal();
  sigma_area ~ normal(0, 0.5);
  sigma_season ~ normal(0, 0.5);
  sigma_series ~ normal(0, 0.5);

  shift_raw ~ normal(0, stz_shift_scale);
  log_speed_raw ~ normal(0, stz_shift_scale);
  sigma_shift ~ normal(0, prior_sd_shift);
  sigma_speed ~ normal(0, prior_sd_speed);

  a_raw ~ std_normal();
  b_raw ~ std_normal();
  sigma_f ~ normal(0, 1);

  beta_lag ~ normal(0, 0.5);
  beta_ses ~ normal(0, 0.3);
  beta_intensity ~ normal(prior_mean_beta_intensity, prior_sd_beta_intensity);
  beta_intensity_cos ~ normal(0, 0.25);
  beta_intensity_sin ~ normal(0, 0.25);
  beta_intensity_sq ~ normal(0, 0.25);

  lambda_raw ~ std_normal();     // prior on the bias parameter

  reciprocal_phi ~ exponential(1);

  // ---- linear predictor --------------------------------------------------
  vector[N] epidemic = rep_vector(0.0, N);
  for (h in 1 : H) {
    epidemic += a_f[h] * sin(h * omega * tau) + b_f[h] * cos(h * omega * tau);
  }

  vector[N] dose = treated .* immunisation_intensity;

  vector[N] intervention = dose
                           .* (beta_intensity
                               + beta_intensity_cos * cos(omega * tau)
                               + beta_intensity_sin * sin(omega * tau));

  if (intensity_quadratic == 1) {
    intervention += beta_intensity_sq * square(dose);
  }

  vector[N] eta = log_population + alpha + area_re[area] + season_re[season]
                  + series_re[series] + epidemic + beta_lag * lag_log_y
                  + beta_ses * ses + intervention;

  // Differential ascertainment enters the likelihood but NOT the estimands,
  // so that its uncertainty propagates into the causal contrast.      
  if (!prior_only) {
    y ~ neg_binomial_2_log(eta - lambda * dose, phi);
  }
}

generated quantities {
  // ---- scalar and low-dimensional estimands only, by default ---------
  real rr_reference = exp(beta_intensity);       // cycle-average RR at I = 1
  real prevented_total;
  real prevented_fraction;
  vector[S] prevented_by_season;

  vector[S] peak_attenuation;      // 1 - peak(observed) / peak(counterfactual)
  vector[S] peak_shift_weeks;      // argmax week, observed minus counterfactual
  vector[S] epidemic_peak_week;    // implied calendar week of the peak     
  real tau_peak;                   // latent time of the baseline maximum

  real prevented_total_dynamic = not_a_number();
  real prevented_fraction_dynamic = not_a_number();

  // Season x week aggregates over ALL observations, used for the national
  // weekly trajectory. Aggregating inside Stan preserves posterior dependence
  // across areas, which summing marginal interval endpoints afterwards would
  // destroy.                                                            
  matrix[S, W] mu_fit_sw;
  matrix[S, W] mu_cf_sw;

  vector[compute_log_lik ? N : 0] log_lik;
  array[compute_y_rep ? N : 0] int y_rep;

  {
    // ---------- controlled contrast (observed lagged history held fixed) ---
    vector[N] tau = (week - week_pivot - shift[season]) .* exp(log_speed[season]);
    vector[N] epidemic_gq = rep_vector(0.0, N);
    for (h in 1 : H) {
      epidemic_gq += a_f[h] * sin(h * omega * tau)
                     + b_f[h] * cos(h * omega * tau);
    }

    vector[N] dose = treated .* immunisation_intensity;

    vector[N] base_no_lag = log_population + alpha + area_re[area]
                            + season_re[season] + series_re[series]
                            + epidemic_gq + beta_ses * ses;

    vector[N] m_dose = dose
                       .* (beta_intensity
                           + beta_intensity_cos * cos(omega * tau)
                           + beta_intensity_sin * sin(omega * tau));

    if (intensity_quadratic == 1) {
      m_dose += beta_intensity_sq * square(dose);
    }

    vector[N] eta_burden = base_no_lag + beta_lag * lag_log_y;
    vector[N] mu_obs = exp(eta_burden + m_dose);
    vector[N] mu_cf = exp(eta_burden);

    // Estimands are accumulated over programme observations only; the weekly
    // aggregates are accumulated over every observation.
    vector[S] num = rep_vector(0.0, S);
    vector[S] den = rep_vector(0.0, S);
    matrix[S, W] agg_obs = rep_matrix(0.0, S, W);
    matrix[S, W] agg_cf = rep_matrix(0.0, S, W);

    for (n in 1 : N) {
      int s = season[n];
      int wk = week_in_season[n];
      agg_obs[s, wk] += mu_obs[n];
      agg_cf[s, wk] += mu_cf[n];
      if (treated[n] > 0) {
        num[s] += mu_obs[n];
        den[s] += mu_cf[n];
      }
    }

    mu_fit_sw = agg_obs;
    mu_cf_sw = agg_cf;

    prevented_by_season = den - num;
    prevented_total = sum(prevented_by_season);
    prevented_fraction = sum(den) > 0 ? 1 - sum(num) / sum(den) : not_a_number();

    for (s in 1 : S) {
      real max_obs = max(agg_obs[s]);
      real max_cf = max(agg_cf[s]);
      if (max_cf > 0) {
        int arg_obs = 1;
        int arg_cf = 1;
        for (w in 1 : W) {
          if (agg_obs[s, w] >= agg_obs[s, arg_obs]) arg_obs = w;
          if (agg_cf[s, w] >= agg_cf[s, arg_cf]) arg_cf = w;
        }
        peak_attenuation[s] = 1 - max_obs / max_cf;
        peak_shift_weeks[s] = arg_obs - arg_cf;
      } else {
        peak_attenuation[s] = not_a_number();
        peak_shift_weeks[s] = not_a_number();
      }
    }

    // Latent time at which the baseline epidemic curve is maximised, located by
    // grid search over one period, and the calendar week to which it maps in
    // each season: t*_s = week_pivot + delta_s + tau* exp(-kappa_s).
    {
      vector[n_grid] f_grid = rep_vector(0.0, n_grid);
      int arg_max = 1;
      for (h in 1 : H) {
        f_grid += a_f[h] * sin(h * omega * tau_grid)
                  + b_f[h] * cos(h * omega * tau_grid);
      }
      for (g in 1 : n_grid) {
        if (f_grid[g] > f_grid[arg_max]) {
          arg_max = g;
        }
      }
      tau_peak = tau_grid[arg_max];
      for (s in 1 : S) {
        epidemic_peak_week[s] = week_pivot + shift[s]
                                + tau_peak * exp(-log_speed[s]);
      }
    }

    if (compute_log_lik == 1) {
      for (n in 1 : N) {
        log_lik[n] = neg_binomial_2_log_lpmf(y[n] | eta_burden[n] + m_dose[n]
                                                    - lambda * dose[n], phi);
      }
    }
    if (compute_y_rep == 1) {
      for (n in 1 : N) {
        y_rep[n] = neg_binomial_2_log_rng(fmin(eta_burden[n] + m_dose[n]
                                               - lambda * dose[n], eta_cap),
                                          phi);
      }
    }

    // ---------- dynamic contrast (history propagated under each arm) -- 
    if (dynamic_gformula == 1) {
      real num_d = 0;
      real den_d = 0;

      for (m in 1 : M) {
        real y_prev_obs = 0;
        real y_prev_cf = 0;

        for (k in 1 : series_len[m]) {
          int n = series_start[m] + k - 1;
          real l_obs = k == 1 ? lag_init
                              : (log1p(y_prev_obs) - lag_center) / lag_scale;
          real l_cf = k == 1 ? lag_init
                             : (log1p(y_prev_cf) - lag_center) / lag_scale;

          real e_obs = base_no_lag[n] + beta_lag * l_obs + m_dose[n];
          real e_cf = base_no_lag[n] + beta_lag * l_cf;

          if (treated[n] > 0) {
            num_d += exp(e_obs);
            den_d += exp(e_cf);
          }

          // Propagate a predictive draw so that the recursion carries the
          // full observation-level variability, not just the conditional mean.
          y_prev_obs = neg_binomial_2_log_rng(fmin(e_obs, eta_cap), phi);
          y_prev_cf = neg_binomial_2_log_rng(fmin(e_cf, eta_cap), phi);
        }
      }

      prevented_total_dynamic = den_d - num_d;
      prevented_fraction_dynamic = den_d > 0 ? 1 - num_d / den_d
                                             : not_a_number();
    }
  }
}
