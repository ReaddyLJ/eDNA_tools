// This is a model file for calculating the qPCR results for Skagit eDNA sampling in 2017.

data { /////////////////////////////////////////////////////////////////////////////////////////////////////
    
    // Number of observations in various categories
    int N_stand_bin ;   // Number of observations for binomial part of the standards model
    int N_stand_pos ; // Number of observations for count part of the standards model
    
    // Observations
    array[N_stand_bin] int bin_obs   ;
    vector[N_stand_pos] pos_stand ;
    
    // Concentrations
    vector[N_stand_bin] K_stand_bin;
    vector[N_stand_pos] K_stand_pos;
    
    // Prior Distributions
    vector[2] beta0_prior;
}
transformed data{
}
parameters { /////////////////////////////////////////////////////////////////////////////////////////////

       real beta_0;
       real beta_1
      
       real gamma_0 ;
       real gamma_1 ;
       
       real phi ;

}
transformed parameters { ////////////////////////////////////////////////////////////////////////////////
      // Predicted Mean and SD 
      vector[N_stand_pos] mu_C ; 
      vector[N_stand_pos] sigma_C ; 
      
}
model {////////////////////////////////////////////////////////////////////////////////////////////////////
    { //LOCAL VARIABLES DECLARATION START
      vector[N_stand_bin] theta_stand ;
      vector[N_obs_bin] theta_obs ;
      
      vector[N_stand_pos] kappa_stand ;
      vector[N_obs_pos] kappa_obs ;
      
    // Presence-Absence component of model.
    for(i in 1:N_stand_bin){
       theta_stand[i] = phi_0[pcr_stand_bin_idx[i]] + phi_1[pcr_stand_bin_idx[i]] *  (D_bin_stand[i] - OFFSET) ;
    }
    
    for(i in 1:N_obs_bin){
       theta_obs[i]  = phi_0[pcr_obs_bin_idx[i]] + 
                      phi_1[pcr_obs_bin_idx[i]] * (D_contam[sample_bin_idx[i]] + 
                                                      bin_log_dilution_obs[i] - 
                                                      OFFSET);
    }

    // Positive Comonent of the model
    for(i in 1:N_stand_pos){
       kappa_stand[i] = beta_0[pcr_stand_pos_idx[i]] + beta_1[pcr_stand_pos_idx[i]] * (D_pos_stand[i] - OFFSET) ;
    }
    
    for(i in 1:N_obs_pos){
      kappa_obs[i]  = beta_0[pcr_obs_pos_idx[i]] + 
                        beta_1[pcr_obs_pos_idx[i]] * (D_contam[sample_pos_idx[i]] + 
                                                      pos_log_dilution_obs[i] - 
                                                      OFFSET) ;
    }
    
    // Likelihood components
    bin_stand     ~ bernoulli( inv_logit(theta_stand) ) ;
    bin_obs       ~ bernoulli( inv_logit(theta_obs) ) ;

    pos_stand   ~ normal(kappa_stand, sigma_all_stand) ;
    pos_obs     ~ normal(kappa_obs, sigma_all_samp) ;

    } //LOCAL VARIABLES DECLARATION END 


    // Parameters for the gBLOCKS:
     // Random effects
      D ~ normal(0,3);
      D_error   ~ normal(mu_contam_fix,sigma_contam_fix); // Unobserved contamination... D_error is an estimated latent variable.
      
      // Priors
      sigma_stand_int ~ gamma(1,1) ;
      sigma_pcr ~ gamma(1,1) ;

      beta_0 ~ normal(40,5) ;
      beta_1 ~ normal(-3.32,0.1) ;

      phi_0 ~ normal(2, 2) ;
      phi_1 ~ normal(4, 2) ;

      delta ~ normal(0,tau_sample) ;
      tau_sample ~ gamma(2,2) ;
}
