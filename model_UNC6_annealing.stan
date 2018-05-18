// Model UNC6_annealing
// AL: Softmax choice rule
// Copy when uncertainty is high (positive frequency dependent)
functions {
  real UNC4_lpmf(int Y, real soc, int trial, vector q, vector q_f, vector f) {
    if(trial == 1 || sum(f)==3e-10) {
      return categorical_lpmf(Y | exp(q));
    }
    else {
      return log_sum_exp( log(1- soc ) + categorical_lpmf(Y | exp(q)),
                          log( soc ) + categorical_lpmf(Y | exp(q_f))
                        );
    }
  }
  /*real function_uncertainty(vector Q, int Size) {
    vector[Size] Q2;
    real Q2sum;
    real entropy;
    real uncertainty;
    Q2 = Q - min(Q);
    Q2 = Q2 + 1e-10; // any elements in Q2 should be > 0
    Q2sum = 0;
    entropy = 0;
    for(c in 1:Size)
        Q2sum = Q2sum + Q2[c];
    for(c in 1:Size)
      entropy = entropy -(Q2[c]/Q2sum)*log(Q2[c]/Q2sum);
    return entropy/(-log(1.0/Size));
  }*/
}

data {
    int<lower = 0> All;  // number of observations
    int<lower = 0> Nsub;  // number of subjects
    int<lower = 0> Ncue;  // number of cues
    int<lower = 1> Ncon;  // number of experimental groups 
    int<lower = 0> Ntrial;  // number of trials per subject
    int<lower = 0> sub[All];  // subject index 
    int<lower = 0> Y[All];  // index of chosen option: 0 => missing
    int<lower = 0> trial[All];  // trial number
    real outcome[All];  // 
    real F[Nsub,Ncue,Ntrial];
    int<lower = 1> condition[Nsub];  // group assignment for each subject
    int<lower = 1> MaxTrial[Nsub];
}

transformed data {
  //vector[Nsub] theta;
  matrix[Ncue, Ntrial] f[Nsub]; // frequency information (transfered)
  //theta = rep_vector(theta_raw, Nsub); // ** THETA is FIXED **
  vector<lower = 0>[All] trial_real;
  for(idx in 1:All)
    trial_real[idx] = trial[idx];
  for(idx in 1:All) {
    for(c in 1:Ncue) {
        if(F[sub[idx]][c,trial[idx]] < 0)
          f[sub[idx]][c,trial[idx]] = 1e-01;
        else
          f[sub[idx]][c,trial[idx]] = F[sub[idx]][c,trial[idx]] + 1e-01;
    }
  }
}

parameters {
    vector[Nsub] alpha_raw_raw;  // learning rate -- raw valuable
    vector[Nsub] beta_raw;  // softmax parameter
    vector[Nsub] soc_raw_raw; // copying rate -- raw
    vector[Nsub] theta_raw;  // conformity
    //vector[Nsub] soc_reduction_raw;  // 
    vector[Nsub] annealing_raw;  // softmax parameter
    
    vector[Ncon] mu_alpha;
    vector<lower=0>[Ncon] mu_beta;
    vector[Ncon] mu_soc;
    vector[Ncon] mu_theta;
    //vector[Ncon] mu_soc_reduction;
    vector[Ncon] mu_annealing;
    
    vector<lower=0>[Ncon] s_alpha;
    vector<lower=0>[Ncon] s_beta;
    vector<lower=0>[Ncon] s_soc;
    vector<lower=0>[Ncon] s_theta;
    //vector<lower=0>[Ncon] s_soc_reduction;
    vector<lower=0>[Ncon] s_annealing;
}

transformed parameters {
  matrix[Ncue, Ntrial] Q[Nsub];  // value function for each target
  matrix[Ncue, Ntrial] q[Nsub]; // softmax choice (log) probability
  matrix[Ncue, Ntrial] q_f[Nsub]; // frequency dependent copying (log) probability
  vector[Ntrial] Fsum[Nsub];
  simplex[Ncue] q_f_simplex[Nsub];
  //vector<lower=0, upper=1>[Ntrial] uncertainty[Nsub];
  
  vector[Nsub] alpha_raw;  // learning rate -- raw valuable
  vector[Nsub] beta;  // softmax parameter
  vector[Nsub] soc_raw; // copying rate -- raw
  vector[Nsub] theta;  // softmax parameter
  //vector[Nsub] soc_reduction;  // softmax parameter
  vector[Nsub] annealing;  // softmax parameter
  vector<lower = 0, upper = 1>[Nsub] alpha;  // asocial learning rate
  vector<lower = 0, upper = 1>[Nsub] soc;  // social learning rate
  vector[Nsub] counter;
  
  counter = rep_vector(0, Nsub); // tracking each subject's first experience timing
  
  for(i in 1:Nsub) {
    alpha_raw[i] = mu_alpha[condition[i]] + s_alpha[condition[i]] * alpha_raw_raw[i];
    beta[i] = mu_beta[condition[i]] + s_beta[condition[i]] * beta_raw[i];
    soc_raw[i] = mu_soc[condition[i]] + s_soc[condition[i]] * soc_raw_raw[i];
    theta[i] = mu_theta[condition[i]] + s_theta[condition[i]] * theta_raw[i];
    annealing[i] = mu_annealing[condition[i]] + s_annealing[condition[i]] * annealing_raw[i];
    //soc_reduction[i] = mu_soc_reduction[condition[i]] + s_soc_reduction[condition[i]] * soc_reduction_raw[i];
  }
  
  for(i in 1:Nsub) {
    alpha[i] = 1/(1+exp(-alpha_raw[i]));
    soc[i] = 1/(1+exp(-soc_raw[i]));
  }
  
  for(idx in 1:All) {
    // SETTING INITIAL Q-VALUE
    if(trial[idx] == 1) {
      for(c in 1:Ncue) {
        Q[sub[idx]][c, trial[idx]] = 1e-10; // Q[t==1] = 1e-10
      }
    }
    // SETTING INITIAL Q-VALUE -- END
    
    // CALCULATION UNCERTAINTY   
    /*uncertainty[sub[idx]][trial[idx]] = 
      function_uncertainty(Q[sub[idx]][,trial[idx]], Ncue);
    if(uncertainty[sub[idx]][trial[idx]] > 1)
      uncertainty[sub[idx]][trial[idx]] = 1;*/
    // CALCULATION UNCERTAINTY -- END
    
    // DENOMINATOR OF SOCIAL FREQUENCY INFORMATION
    Fsum[sub[idx]][trial[idx]] = 0;
    for(c in 1:Ncue) {
      Fsum[sub[idx]][trial[idx]] = Fsum[sub[idx]][trial[idx]] + f[sub[idx]][c,trial[idx]]^theta[sub[idx]];
    }
    // DENOMINATOR OF SOCIAL FREQUENCY INFORMATION -- END
    
    // ASOCIAL AND SOCIAL CHOICE PROBABILITY
    for(c in 1:Ncue) {
        // Choice probability by Asocial sofimax rule
      q[sub[idx]][c, trial[idx]] = log_softmax(Q[sub[idx]][, trial[idx]]*( (beta[sub[idx]] + annealing[sub[idx]]*(trial_real[idx]/70)) ))[c];
        // social frequency influence
      q_f_simplex[sub[idx]][c] = exp(theta[sub[idx]]*log(f[sub[idx]][c,trial[idx]])-log(Fsum[sub[idx]][trial[idx]]));
    }
    for(c in 1:Ncue)
      q_f[sub[idx]][c, trial[idx]] = log(q_f_simplex[sub[idx]][c]); // log transform
    // ASOCIAL AND SOCIAL CHOICE PROBABILITY -- END
    
    // Q-VALUE UPDATE
    if(trial[idx] < Ntrial) {
      for(c in 1:Ncue)
        Q[sub[idx]][c, (trial[idx]+1)] = Q[sub[idx]][c, trial[idx]];
      // update chosen option
      if(Y[idx]>0) {
        if(counter[sub[idx]]==0) 
          { // Updating all Q-values at the 1st experience
            for(c in 1:Ncue)
              Q[sub[idx]][c, (trial[idx]+1)] = (1-alpha[sub[idx]])*Q[sub[idx]][c, trial[idx]] + alpha[sub[idx]]*outcome[idx];
            counter[sub[idx]] = 1; 
          }
        else
          { // Updating chosen option's Q-value
            Q[sub[idx]][Y[idx], (trial[idx]+1)] = 
              (1-alpha[sub[idx]])*Q[sub[idx]][Y[idx], trial[idx]] + alpha[sub[idx]]*outcome[idx];
          }
      }
    }
    // Q-VALUE UPDATE -- END
  }
}

model {
  mu_alpha ~ normal(0, 5); // prior -3 ~ 3
  mu_beta ~ normal(0, 5); // prior 1 ~ 5
  mu_soc ~ normal(0, 5); // prior -3 ~ 3
  mu_theta ~ normal(0, 5); // prior -3 ~ 7
  //mu_soc_reduction ~ normal(0, 5); //normal(-2, 3); // prior -4 ~ 2
  mu_annealing ~ normal(0, 5);
  
  s_alpha ~ normal(0, 3); // prior
  s_beta ~ normal(0, 3); // prior
  s_soc ~ normal(0, 3); // prior
  s_theta ~ normal(0, 3); // prior
  //s_soc_reduction ~ normal(0, 3); // prior 
  s_annealing ~ normal(0, 3); // prior
  
  alpha_raw_raw ~ student_t(4, 0, 1); // prior
  beta_raw ~ student_t(4, 0, 1); // prior
  soc_raw_raw ~ student_t(4, 0, 1); // prior
  theta_raw ~ student_t(4, 0, 1); // prior
  //soc_reduction_raw ~ student_t(4, 0, 1); // prior
  annealing_raw ~ student_t(4, 0, 1); // prior
  
  for(idx in 1:All) {
    if(Y[idx] > 0) {
      target += 
        //1/log(All) *
        UNC4_lpmf(Y[idx] | soc[sub[idx]], trial[idx], q[sub[idx]][,trial[idx]], q_f[sub[idx]][,trial[idx]], f[sub[idx]][,trial[idx]] );
    }
  }
}

generated quantities {
  vector[Nsub] log_lik;
  vector[Ntrial] netBeta[Nsub];
  log_lik = rep_vector(0, Nsub); // initial values for log_lik
  for (i in 1:Nsub) {
    for (t in 1:Ntrial) {
      //soc[i][t] = 1/(1+exp(-(soc_raw[i]+(soc_reduction[i]*t/70))));
      netBeta[i][t] = (beta[i] + annealing[i]*t/70);
    }
  }
  for(idx in 1:All) {
    if(Y[idx] > 0) {
      if(trial[idx] == 1 || sum(f[sub[idx]][,trial[idx]])==3e-10)
        log_lik[sub[idx]] = log_lik[sub[idx]] + q[sub[idx]][Y[idx],trial[idx]];
      else
        log_lik[sub[idx]] = log_lik[sub[idx]] + 
          log( (1-soc[sub[idx]])*exp(q[sub[idx]][Y[idx],trial[idx]])+
                soc[sub[idx]]*exp(q_f[sub[idx]][Y[idx],trial[idx]]) 
          );
          //log( (1-1/(1+exp(-(soc_raw[sub[idx]]+(soc_reduction[sub[idx]]*trial[idx]/70)))) )*exp(q[sub[idx]][Y[idx],trial[idx]])+
                //1/(1+exp(-(soc_raw[sub[idx]]+(soc_reduction[sub[idx]]*trial[idx]/70))))*exp(q_f[sub[idx]][Y[idx],trial[idx]]) );
    }
  }
}

