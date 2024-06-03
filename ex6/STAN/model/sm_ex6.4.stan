data {
  int<lower=0> n;
  int<lower=0> y;
}

parameters {
  real<lower=0, upper=1> p; 
}

model {
    // a prior for p
    p ~ beta(1, 10);
    
    // Likelihood function
    y ~ binomial(n, p);
}

