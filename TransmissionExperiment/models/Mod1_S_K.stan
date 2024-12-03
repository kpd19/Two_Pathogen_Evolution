// Model 1
// Morph same, tree sp same, strains come from grand mean and sd

data {
  int<lower=1> N;  		  	// number of observations
  int<lower=1> I;			// number of strains
  int virus[N];     			// vector of infected individuals
  int alive[N];   			// vector of total individuals
  int isolate[N];			// vector of isolate id
  vector[N] pathogen;         		// amount of pathogen
}
parameters {
  vector<lower=0.00000001,upper=0.1>[I] nu_bar;  		// transmission rate
  vector<lower=0.00000001,upper=5>[I] C;			// Heterogeneity
  real<lower = 0, upper = 0.1> mu;
  real<lower = 0, upper = 5> muC;
  real<lower = 0, upper = 4> sigma;
  real<lower = 0, upper = 5> sigmaC;
}
transformed parameters {
  vector<lower=0, upper=1>[N] prob;

  for(n in 1:N){
    prob[n] = 1 - (1 + nu_bar[isolate[n]]* C[isolate[n]]^2 * pathogen[n] * 7.0*0.032)^(-1/C[isolate[n]]^2); 	// p_infection
  }

}
model {
  // priors
  mu ~ normal(0.001,0.01);
  muC ~ normal(3,1);
  sigma ~ lognormal(-3,1);
  sigmaC ~ lognormal(0,0.5);

  for(i in 1:I){
    nu_bar[i] ~ normal(mu,sigma);
    C[i] ~ normal(muC,sigmaC);
  }
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
