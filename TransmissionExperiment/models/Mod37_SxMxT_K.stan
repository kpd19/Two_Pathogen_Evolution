// Model 2
// Morph same, tree sp diff, strains come from grand mean and sd

data {
  int<lower=1> N;  		  	// number of observations
  int<lower=1> I;			// number of isolates
  int<lower=1> T;			// number of tree species
  int virus[N];     			// vector of infected individuals
  int alive[N];   			// vector of total individuals
  int isolate[N];			// vector of isolate id
  int trees[N];				// vector of tree id
  int morph[N];				// vector of morph id
  vector[N] pathogen;         		// amount of pathogen
}
parameters {
  matrix<lower=0.00000001,upper=0.1>[I,T] nu_bar_S;  	// tree
  vector<lower=0.00000001,upper=5>[I] C_S;  		// no
  vector<lower = 0, upper = 0.1>[T] mu_S;		// tree
  real<lower = 0, upper = 5> muC_S;			// no
  vector<lower = 0, upper = 4>[T] sigma_S;		// tree
  real<lower = 0, upper = 5> sigmaC_S;			// no

  vector<lower=0.00000001,upper=0.1>[I] nu_bar_M;  	// no
  vector<lower=0.00000001,upper=5>[I] C_M;  		// no
  real<lower = 0, upper = 0.1> mu_M;			// no
  real<lower = 0, upper = 5> muC_M;			// no
  real<lower = 0, upper = 4> sigma_M;			// no
  real<lower = 0, upper = 5> sigmaC_M;			// no

}
transformed parameters {
  vector<lower=0, upper=1>[N] prob;

  for(n in 1:N){
    if (morph[n] == 1) // if multi-capsid, C varies by tree, nu_bar does 
      prob[n] = 1 - (1 + nu_bar_M[isolate[n]]* C_M[isolate[n]]^2 * pathogen[n] * 7.0*0.032)^(-1/C_M[isolate[n]]^2); 	// p_infection
    else // if single-capsid, C varies by tree, nu_bar does 
      prob[n] = 1 - (1 + nu_bar_S[isolate[n],trees[n]]* C_S[isolate[n]]^2 * pathogen[n] * 7.0*0.032)^(-1/C_S[isolate[n]]^2); 	// p_infection
  }

}
model {
  // priors

  mu_M ~ normal(0.001,0.01);
  sigma_M ~ lognormal(-3,1);

  muC_M ~ normal(3,1);
  sigmaC_M ~ lognormal(0,0.5);

  muC_S ~ normal(3,1);
  sigmaC_S ~ lognormal(0,0.5);

  for(t in 1:T){
    mu_S[t] ~ normal(0.001,0.01);
    sigma_S[t] ~ lognormal(-3,1); 
  }
  
  for(i in 1:I){
    nu_bar_M[i] ~ normal(mu_M,sigma_M);
    C_M[i] ~ normal(muC_M,sigmaC_M);
    C_S[i] ~ normal(muC_S,sigmaC_S);
  }

  for(t in 1:T) {  
    for(i in 1:I){
      nu_bar_S[i,t] ~ normal(mu_S[t],sigma_S[t]);
    }
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
