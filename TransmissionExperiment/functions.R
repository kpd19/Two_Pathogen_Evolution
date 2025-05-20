# functions

# input_data for the functions


get_logC <- function(x,nu_bar,C,time){
  
  log_val = (1/C^2)*log(1 + nu_bar*C^2*time*x*rho)
  
  return(log_val)
}


plot_summary <- function(mod,nu_bars,Cs){
  p_mod <- as.array(mod)
  np_mod <- nuts_params(mod)
  lp_mod <- log_posterior(mod)
  
  if(is.na(Cs)[1] == FALSE){
    poi <- c(nu_bars,Cs)
  } else{
    poi <- c(nu_bars)
  }
  
  color_scheme_set("brewer-Paired")
  
  plt1 <- mcmc_trace(mod, pars = c('lp__',poi),
                     np=np_mod, n_warmup = 5000, window = c(5000,10000)) + 
    theme_classic()

  plt2 <- mcmc_nuts_divergence(np_mod, lp_mod) + 
    theme_classic()
  
  rhats <- rhat(mod)
  plt3 <- mcmc_rhat(rhats, alpha = 0.5)
  
  plt4 <- stan_hist(mod, pars = nu_bars, inc_warmup = FALSE, fill = '#a6cee3', alpha = 1)
  
  print(plt1)
  print(plt2)
  print(plt3)
  print(plt4)
  
  if(is.na(Cs)[1] == FALSE){
    plt5 <- stan_hist(mod, pars = Cs, inc_warmup = FALSE, fill = '#a6cee3', alpha = 1)
    print(plt5)
  } else{
    print("no het")
  }
}

loo_extract <- function(stan_obj) {
  log_lik <- extract_log_lik(stan_obj,
                             merge_chains = F)
  r_eff <- relative_eff(exp(log_lik))
  loo_obj <- loo(log_lik, r_eff = r_eff)
  
  loo_est <- loo_obj$estimates
  return(loo_est)
}


