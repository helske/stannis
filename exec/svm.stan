data {
  int<lower=0> n;
  // # time points (equally spaced)
    vector[n] y;
  // mean corrected return at time t
}
parameters {
  real mu;
  // mean log volatility
  real<lower=-1,upper=1> phi; // persistence of volatility
  real<lower=0> sigma;// white noise shock scale
  // log volatility at time t
  vector[n] h_std;
}

transformed parameters {
  vector[n] h;
  // log volatility at time t
  h = h_std * sigma;
  // now h ~ normal(0, sigma)
  h[1] = h[1] / sqrt(1 - phi * phi); // rescale h[1]
  h = h + mu;
  for (t in 2:n)
    h[t] = h[t] + phi * (h[t-1] - mu);
}

model {
  phi ~ uniform(-1, 1);
  sigma ~ normal(0, 5);
  mu ~ normal(0, 10);
  h_std ~ normal(0, 1);
  y ~ normal(0, exp(h / 2));
}
