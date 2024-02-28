#This is a negative log likelihood function for the 1D Ising model
#with memory. 

negative_log_likelihood_transition = function(params, S_t, S_t1) {
  J = params[1]
  K = params[2]
  
  log_likelihood = 0
  S_t = t(S_t)
  S_t1 = t(S_t1)

  n_cols = ncol(S_t)

  i = 1
  for (j in 1:n_cols) {

    # Get the indices of the nearest neighbors (assuming periodic boundary conditions)
    neighbors_i = c( i, i)
    neighbors_j = c(j - 1, j + 1)
    neighbors_j = ifelse(neighbors_j < 1, n_cols, neighbors_j)
    neighbors_j = ifelse(neighbors_j > n_cols, 1, neighbors_j)
    
    # Extract the spins of the neighbors
    neighbors = S_t[ cbind(neighbors_i, neighbors_j) ]
    
    # Compute the effective field
    h_eff = K * S_t[i,j] + J * sum( neighbors )
  
    #Probability of flipping the spin
    pft = exp(h_eff*(-S_t[i,j]) ) / (2*cosh(h_eff))
    
    # Ensure that probabilities are finite
    pft  = ifelse(is.finite(pft), pft, 1e-10)
    #pft  = ifelse(pft  > 0, pft, 1e-10)
    
    #print(c(j,pft))
    #Log probability of whether it flipped or not: 
    log_pft = log(ifelse( S_t1[i,j] == S_t[i,j] , 1-pft, pft))
    log_likelihood = log_likelihood + log_pft
 
  }
  
  return(-log_likelihood)  # return the negative log-likelihood for minimization
}


#This function iterates the parameter inference over every subsequent pair
#of time steps. It is a wrapper for nlm.
fit_params_nlm = function( initial_params, final_grid ){

  # Iterate over adjacent time steps and perform optimization for each pair
  num_pairs =nrow(final_grid) - 1
  parameter_estimates = matrix(NA, nrow = num_pairs, ncol = 2)
  parameter_covariances = matrix(NA, nrow = num_pairs, ncol = 2)

  for (t in 1:num_pairs) {
    # Extract spin configurations for the current pair of time steps
    S_t = final_grid[t, ]
    S_t1 = final_grid[(t+1), ]

    #Run nlm to get the result
    result = nlm(f = negative_log_likelihood_transition, p = initial_params, 
                  S_t =S_t, S_t1 = S_t1, hessian=TRUE)

    # Store parameter estimates
    parameter_estimates[t, ] = result$estimate
   
   # Compute covariance matrix from Hessian matrix
    if(sum(result$hessian)>0){
      # Compute standard errors from covariance matrices
      parameter_covariances[t, ] = sqrt(diag(solve(result$hessian)))  
    }
  }

  parm_output = NULL
  parm_output$estimate = parameter_estimates
  parm_output$cov = parameter_covariances

  return(parm_output)
}