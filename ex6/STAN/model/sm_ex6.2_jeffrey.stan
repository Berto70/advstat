//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  int y[N]; // y should be integers for Poisson distribution
}

// The parameters accepted by the model. Our model
// accepts one parameter 'lambda'.
parameters {
  real<lower=0> lambda; // lambda must be positive for Poisson distribution
}

// The model to be estimated. We model the output
// 'y' to be Poisson distributed with rate 'lambda'.
model {
    // a prior for lambda
    lambda ~ gamma(0.00001, 0.00001);
    
    // Likelihood function
    y ~ poisson(lambda);
}


