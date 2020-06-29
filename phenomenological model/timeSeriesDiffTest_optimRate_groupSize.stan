data {
  int All;
	int T[All];
	int maxT;
	int N_indv;
	int Indv[All];
	int Y[All];
	int indiv_group[All];
	real GroupSize_st[All];
	real smallGroupSize; // n = 4
	real largeGroupSize; // n = 12
	real veryLargeGroupSize; // n = 20
}

parameters {
  real mu0;
  real di0;
  vector[2] sizeEffect;
  vector<lower=-pi()/2, upper=pi()/2>[maxT-1] mu_unif;
  vector<lower=-pi()/2, upper=pi()/2>[maxT-1] di_unif;
  vector[N_indv] ri;
  real<lower=0> s_mu;
  real<lower=0> s_di;
  real<lower=0> s_i;
}

transformed parameters {
  vector[maxT] mu;
  vector[maxT] di;
  vector[All] q;
  vector[N_indv] x_ri;
  mu[1] = mu0;
  di[1] = di0;
  x_ri = ri;
  
  // mu walks randomly with a Caychy probability distribution (reparametarized)
  for(t in 2:maxT)
    mu[t] = mu[t-1] + s_mu*tan(mu_unif[t-1]); 
  // di walks randomly with a Caychy probability distribution (reparametarized)
  for (t in 2:maxT)
    di[t] = di[t-1] + s_di*tan(di_unif[t-1]); 
  // Probability to choose the best option q[i] depends on mu and di
  // indiv_group[i] is a binary indicator (indiv=0/group=1)
  // Also, q[i]'s intercept depends on individual and group rondom efffects
  for(idx in 1:All) {
    if(T[idx]<=40) {
      q[idx] = inv_logit( mu[T[idx]] + indiv_group[idx]*(di[T[idx]] + sizeEffect[1]*GroupSize_st[idx] ) + x_ri[Indv[idx]]  );
    } else {
      q[idx] = inv_logit( mu[T[idx]] + indiv_group[idx]*(di[T[idx]] + sizeEffect[2]*GroupSize_st[idx] ) + x_ri[Indv[idx]]  );
    }
  }
}

model {
  // Random effects. Standard deviation is a hiper parameter
  target += normal_lpdf(ri | 0, s_i);
	// Y is a binary data about choice performance (optimal choice = 1/otherwise 0)
	Y ~ bernoulli( q );
}

generated quantities {
  vector[N_indv] log_lik;
  vector[70] p_indiv;
  vector[70] p_group_small;
  vector[70] p_group_large;
  vector[70] p_group_veryLarge;
  vector[70] p_group;
  log_lik = rep_vector(0, N_indv); // initial values for log_lik
  p_indiv = inv_logit(mu);
  p_group = inv_logit(mu + di);
  for(t in 1:70) {
    if(t <= 40) {
      p_group_small[t] = inv_logit(mu[t] + di[t] + sizeEffect[1]*smallGroupSize);
      p_group_large[t] = inv_logit(mu[t] + di[t] + sizeEffect[1]*largeGroupSize);
      p_group_veryLarge[t] = inv_logit(mu[t] + di[t] + sizeEffect[1]*veryLargeGroupSize);
    } else {
      p_group_small[t] = inv_logit(mu[t] + di[t] + sizeEffect[2]*smallGroupSize);
      p_group_large[t] = inv_logit(mu[t] + di[t] + sizeEffect[2]*largeGroupSize);
     p_group_veryLarge[t] = inv_logit(mu[t] + di[t] + sizeEffect[2]*veryLargeGroupSize);
    }
  }
  for(idx in 1:All) {
    // Individual log-likelihood
    if(Y[idx] > 0)
      log_lik[Indv[idx]] = log_lik[Indv[idx]] + q[Indv[idx]];
  }
}
