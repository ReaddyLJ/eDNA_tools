// This is a model file for calculating the qPCR results for Skagit eDNA sampling in 2017.

data { /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Standards data
    // Number of observations in various categories
    int N_stand ;   // Number of observations for binomial part of the standards model
    int N_stand_pos ; // Number of observations for count part of the standards model
    
    // Observations
    array[N_stand] int bin_obs_stand    ;
    vector[N_stand_pos] ct_obs_stand ;
    
    // Concentrations
    vector[N_stand] K;
    vector[N_stand] log_K;
    
    // Indexes
    int N_plate;
    array[N_stand] int plate_stand_idx;
    array[N_stand_pos] int pos_stand_idx;
    
    // Field Samples 
    // Counters
    int N_obs_samp;
    int N_obs_samp_pos;
    int N_samp;
               
    // Observations             
    array [N_obs_samp] int bin_obs_samp ;
    vector[N_obs_samp_pos] ct_obs_samp ; 
    
    //Indexes
    array[N_obs_samp_pos] int pos_idx ;
    array[N_obs_samp] int plate_samp_idx;
    array[N_obs_samp] int samp_idx;
    
    // Prior Distributions and constraints
    vector[2] beta0_prior;
    vector[2] beta1_prior ;
    real beta1_max ; // hard constrainton maximum PCR efficiency. Cannot exceed 3.
    vector[2] gamma0_prior ;
    vector[2] gamma1_prior ;
    vector[2] logit_phi_prior ;
    vector[2] log_D_prior;
    
    // Predictions
    int N_pred; // Number of predictions
    vector[N_pred] log_K_pred; // 
    vector[N_pred] K_pred; // 
}
transformed data{
}
parameters { /////////////////////////////////////////////////////////////////////////////////////////////

       vector[N_plate] beta_0;
       vector<upper = beta1_max>[N_plate] beta_1 ;
      
        real<lower=0> sigma_0 ;
      
       real gamma_0 ;
       real gamma_1 ;
       
       real logit_phi ; 
       
       //Unknown field samples.
       vector[N_samp] log_D;
       
}
transformed parameters { ////////////////////////////////////////////////////////////////////////////////
      // Predicted prob of detection 
      vector[N_stand] theta_stand;
      
      // Predicted Mean and SD 
      vector[N_stand] mu_C ; 
      vector[N_stand] sigma_C ; 
      
      mu_C = beta_0[plate_stand_idx] + beta_1[plate_stand_idx] .* (log(2)+log_K);
      sigma_C = pow(pow(sigma_0,2) + pow(exp(gamma_0 + gamma_1 * (log(2)+log_K)),2), 0.5) ;
      theta_stand = 1 - exp(-2*K*inv_logit(logit_phi)) ;
      
      // predictions for each field sample
      vector[N_obs_samp] log_mu_samp;
      vector[N_obs_samp] sigma_samp;
      vector[N_obs_samp] theta_samp;
      
      log_mu_samp =  beta_0[plate_samp_idx] + beta_1[plate_samp_idx] .* (log(2) + log_D[samp_idx]) ;
      sigma_samp = pow(pow(sigma_0,2) + pow(exp(gamma_0 + gamma_1 * (log(2) + log_D[samp_idx])),2), 0.5) ;
      theta_samp = 1 - exp(-2*exp(log_D[samp_idx])*inv_logit(logit_phi)) ;
           
}
model {////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Likelihood components
    //Standards
      // Detection 
      bin_obs_stand   ~ bernoulli(theta_stand) ;
      // Positive observations
      ct_obs_stand   ~ normal(mu_C[pos_stand_idx], sigma_C[pos_stand_idx]) ;
    
    // Field Samples
      // Detection 
      bin_obs_samp   ~ bernoulli(theta_samp) ;
      // Positive observations
      ct_obs_samp ~ normal(log_mu_samp[pos_idx],sigma_samp[pos_idx]);
      
      // Priors
      beta_0  ~ normal(beta0_prior[1],beta0_prior[2]) ;
      beta_1  ~ normal(beta1_prior[1],beta1_prior[2]) ;

      sigma_0 ~ normal(0,1) ;
      gamma_0 ~ normal(gamma0_prior[1],gamma0_prior[2]) ;
      gamma_1 ~ normal(gamma1_prior[1],gamma1_prior[2]) ;

      logit_phi ~ normal(logit_phi_prior[1],logit_phi_prior[2]) ;
      
      //prior for DNA concentration.
      log_D ~ normal(log_D_prior[1],log_D_prior[2]);
}
generated quantities{
  vector[N_pred] pred_SD  ;
  vector[N_pred] pred_theta  ;
  vector[N_samp] D;
  
  real phi ;
  for(i in 1:N_pred){
     pred_SD[i] = normal_rng( 0 , 
                             pow(pow(sigma_0,2) + pow(exp(gamma_0 + gamma_1 * (log(2)+log_K_pred[i])),2), 0.5));
    }
  pred_theta = 1 - exp(-2*K_pred*inv_logit(logit_phi)) ;
  phi = inv_logit(logit_phi) ;
  
  D = exp(log_D);
}
