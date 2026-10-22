// This is a model file for calculating the qPCR results for Skagit eDNA sampling in 2017.

data { /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Number of observations in various categories
    int N_stand_bin ;   // Number of observations for binomial part of the standards model
    int N_stand_pos ; // Number of observations for count part of the standards model
    
    // Observations
    array[N_stand_bin] int bin_obs     ;
    vector[N_stand_pos] ct_obs ;
    
    // Concentrations
    vector[N_stand_bin] K_bin;
    vector[N_stand_pos] log_K_pos;
    
    // Prior Distributions and constraints
    vector[2] beta0_prior;
    vector[2] beta1_prior ;
    real beta1_max ; // hard constrainton maximum PCR efficiency.
    vector[2] gamma0_prior ;
    vector[2] gamma1_prior ;
    vector[2] logit_phi_prior ;
    
    // Predictions
    int N_pred; // Number of predictions
    vector[N_pred] log_K_pred; // 
    vector[N_pred] K_pred; // 
}
transformed data{
}
parameters { /////////////////////////////////////////////////////////////////////////////////////////////

       real beta_0;
       real<upper = beta1_max> beta_1 ;
      
        //real<lower=0> sigma_0 ;
      
       real gamma_0 ;
       real gamma_1 ;
       
       real logit_phi ; 
}
transformed parameters { ////////////////////////////////////////////////////////////////////////////////
      // Predicted prob of detection 
      vector[N_stand_bin] theta_stand;
      
      // Predicted Mean and SD 
      vector[N_stand_pos] mu_C ; 
      vector[N_stand_pos] sigma_C ; 
      
      mu_C = beta_0 + beta_1 * (log(2)+log_K_pos);
      sigma_C = exp(gamma_0 + gamma_1 * (log(2)+log_K_pos));
      
      theta_stand = 1 - exp(-2*K_bin*inv_logit(logit_phi)) ;
}
model {////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Likelihood components
    // Detection
    //bin_obs   ~ bernoulli(theta_stand) ;
    // Positive observations
    ct_obs   ~ normal(mu_C, sigma_C) ;
    
      // Priors
      beta_0  ~ normal(beta0_prior[1],beta0_prior[2]) ;
      beta_1  ~ normal(beta1_prior[1],beta1_prior[2]) ;

      //sigma_0 ~ normal(0,1) ;
      gamma_0 ~ normal(gamma0_prior[1],gamma0_prior[2]) ;
      gamma_1 ~ normal(gamma1_prior[1],gamma1_prior[2]) ;

      logit_phi ~ normal(logit_phi_prior[1],logit_phi_prior[2]) ;
      
}
generated quantities{

  vector[N_pred] pred_theta;  
  vector[N_pred] pred_C;
  for(i in 1:N_pred){
    pred_C[i] = normal_rng(beta_0 + beta_1 * (log(2)+log_K_pred[i]), 
                            exp(gamma_0 + gamma_1 * (log(2)+log_K_pred[i])));
    pred_theta[i] = 1 - exp(-2*K_pred[i]*inv_logit(logit_phi)) ;
  }
}
