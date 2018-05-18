data {
  int All;
  int test_size;
  int N_group;
  vector[All] taskDifficulty;
	vector[All] groupSize_standerized;
	vector[All] age_standardized;
	vector[All] female;
	int groupID[All];
	vector[All] conformity_exponent;
	vector[All] soc_change;
	vector[test_size] test_group_size;
}

parameters {
    vector[4] beta;
    real gamma;
  	vector[N_group] r_group;
  	real<lower=0> sigma_e;
  	real<lower=0, upper=5> sigma_r;
}

model {
  real mu;
  real sigma;
  // priors
  beta[1] ~ cauchy(0, 4); //
  for(i in 2:4) {
    // prior for the fixed effects following Gelman 2008
    beta[i] ~ cauchy(0, 2.5); 
  }
  r_group ~ normal(0, sigma_r);
  // likelihood
  for(i in 1:All) {
    mu =  beta[1]+
          beta[2]*conformity_exponent[i] + 
          beta[3]*taskDifficulty[i] + 
          //beta[4]*age_standardized[i] +
          //beta[5]*female[i] +
          beta[4]*conformity_exponent[i]*taskDifficulty[i] +
          //beta[7]*conformity_exponent[i]*age_standardized[i] +
          //beta[8]*conformity_exponent[i]*female[i] +
          //beta[9]*taskDifficulty[i]*age_standardized[i] +
          //beta[10]*taskDifficulty[i]*female[i] +
          //beta[11]*age_standardized[i]*female[i] +
          r_group[groupID[i]];
    sigma = sigma_e + gamma*conformity_exponent[i];
    soc_change[i] ~ normal(mu, sigma);
  }
}

generated quantities {
  vector[test_size] soc_change_pred_25;
  vector[test_size] soc_change_pred_50;
  vector[test_size] soc_change_pred_78;
  for(gs in 1:test_size){
    soc_change_pred_25[gs] = 
    normal_rng(
        beta[1]+
        beta[2]*conformity_exponent[gs] + 
        beta[3]*0 +//difficulty
        beta[4]*conformity_exponent[gs]*0
      , sigma_e + gamma*conformity_exponent[gs]);
    soc_change_pred_50[gs] = 
    normal_rng(
        beta[1]+
        beta[2]*conformity_exponent[gs] + 
        beta[3]*0.5 +//difficulty
        beta[4]*conformity_exponent[gs]*0.5
      , sigma_e + gamma*conformity_exponent[gs]);
    soc_change_pred_78[gs] = 
    normal_rng(
        beta[1]+
        beta[2]*conformity_exponent[gs] + 
        beta[3]*1 +//difficulty
        beta[4]*conformity_exponent[gs]*1
      , sigma_e + gamma*conformity_exponent[gs]);
  }
}
