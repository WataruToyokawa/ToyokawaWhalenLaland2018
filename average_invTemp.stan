data {
  int All;
  int test_size;
  int N_group;
  vector[All] taskDifficulty;
	vector[All] groupSize_standerized;
	vector[All] age_standardized;
	vector[All] female;
	int groupID[All];
	vector[All] average_invTemp;
	vector[test_size] test_group_size;
}

parameters {
    vector[11] beta;
  	vector[N_group] r_group;
  	real<lower=0> sigma_e;
  	real<lower=0, upper=5> sigma_r;
}

model {
  real mu;
  // priors
  beta[1] ~ cauchy(0, 4); //
  for(i in 2:11) {
    // prior for the fixed effects following Gelman 2008
    beta[i] ~ cauchy(0, 2.5); 
  }
  r_group ~ normal(0, sigma_r);
  // likelihood
  for(i in 1:All) {
    mu =  beta[1]+
          beta[2]*groupSize_standerized[i] + 
          beta[3]*taskDifficulty[i] + 
          beta[4]*age_standardized[i] +
          beta[5]*female[i] +
          beta[6]*groupSize_standerized[i]*taskDifficulty[i] +
          beta[7]*groupSize_standerized[i]*age_standardized[i] +
          beta[8]*groupSize_standerized[i]*female[i] +
          beta[9]*taskDifficulty[i]*age_standardized[i] +
          beta[10]*taskDifficulty[i]*female[i] +
          beta[11]*age_standardized[i]*female[i] +
          r_group[groupID[i]];
    average_invTemp[i] ~ normal(mu, sigma_e);
  }
}

generated quantities {
  vector[test_size] average_invTemp_pred_25;
  vector[test_size] average_invTemp_pred_50;
  vector[test_size] average_invTemp_pred_78;
  for(gs in 1:test_size){
    average_invTemp_pred_25[gs] = 
    normal_rng(
       beta[1]+
        beta[2]*test_group_size[gs] + 
        beta[3]*0 +//difficulty
        beta[4]*0 +//mode age
        beta[5]*0.5 +//mode female
        beta[6]*test_group_size[gs]*0 +
        beta[7]*test_group_size[gs]*0 +
        beta[8]*test_group_size[gs]*0.5 +
        beta[9]*0*0 +
        beta[10]*0*0.5 +
        beta[11]*0*0.5
      , sigma_e);
    average_invTemp_pred_50[gs] = 
    normal_rng(
        beta[1]+
        beta[2]*test_group_size[gs] + 
        beta[3]*0.5 +//difficulty
        beta[4]*0 +//mode age
        beta[5]*0.5 +//mode female
        beta[6]*test_group_size[gs]*0.5 +
        beta[7]*test_group_size[gs]*0 +
        beta[8]*test_group_size[gs]*0.5 +
        beta[9]*0.5*0 +
        beta[10]*0.5*0.5 +
        beta[11]*0*0.5
      , sigma_e);
    average_invTemp_pred_78[gs] = 
    normal_rng(
        beta[1]+
        beta[2]*test_group_size[gs] + 
        beta[3]*1 +//difficulty
        beta[4]*0 +//mode age
        beta[5]*0.5 +//mode female
        beta[6]*test_group_size[gs]*1 +
        beta[7]*test_group_size[gs]*0 +
        beta[8]*test_group_size[gs]*0.5 +
        beta[9]*1*0 +
        beta[10]*1*0.5 +
        beta[11]*0*0.5
      , sigma_e);
  }
}
