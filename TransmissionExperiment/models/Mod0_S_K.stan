// Model 1
// Morph same, tree sp same, strains come from grand mean and sd

data {
  int<lower=1> N;  		  	// number of observations
  int virus[N];     			// vector of infected individuals
  int alive[N];   			// vector of total individuals
  vector[N] pathogen;         		// amount of pathogen
}
parameters {
  real<lower=0.00000001,upper=0.1> nu_bar;  		// transmission rate
  real<lower=0.00000001,upper=5> C;			// Heterogeneity
}
transformed parameters {
  vector<lower=0, upper=1>[N] prob;

  for(n in 1:N){
    prob[n] = 1 - (1 + nu_bar* C^2 * pathogen[n] * 7.0*0.032)^(-1/C^2); 	// p_infection
  }

}
model {
  // priors
  nu_bar ~ normal(0.001,0.01);
  C ~ normal(3,1);

  // model
  for(n in 1:N){
    virus[n] ~  binomial(alive[n], prob[n]);	// log-likelihood
  }
}
generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_lpmf(virus[n] | alive[n], prob[n]);
  }
}
