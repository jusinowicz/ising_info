#==============================================================================
#These are functions used throughout the code to:
#   Model Ricker metapopulation dynamics on a 1D lattice
#   Model Ising dynamics with memory on a 1D lattice
#   Infer the parameters of an Ising model from a given lattice.
#==============================================================================

#==============================================================================
#Ricker models
#==============================================================================

# Function to simulate one time step of the Ricker metapopulation model for a 1D lattice
ricker_step_1d = function(X, L, kappa, r, sigma) {
  # Define the Ricker operator
  R = function(X) X * exp(r * (1 - X))
  
  # Define the dispersal operator
  D = function(X, kappa, neighbors) {
    (1 - kappa) * X + kappa * sum(neighbors) / length(neighbors)
  }
  
  # Define the environmental noise operator
  N = function(X, sigma) X * exp(sigma * rnorm(1))
  
  # Initialize the updated lattice
  X_new = numeric(L)
  
  # Loop through each patch on the lattice
  for (i in 1:L) {
    # Get the neighboring patches
    neighbors = c(X[i], X[ifelse(i - 1 > 0, i - 1, L)], X[ifelse(i + 1 <= L, i + 1, 1)])
    
    # Update the patch based on the Ricker model with dispersal and noise
    X_new[i] = N(D(R(X[i]), kappa, neighbors), sigma)
  }
  
  return(X_new)
}

# Function to simulate the Ricker metapopulation model over multiple time steps for a 1D lattice
simulate_ricker_1d = function(L, kappa, r, sigma, time_steps) {
  # Initialize the lattice with random initial population densities
  X = runif(L)
  
  # Initialize a matrix to store the lattice at each time step
  lattice_matrix = matrix(0, nrow = L, ncol = time_steps + 1)
  lattice_matrix[, 1] = X
  
  # Simulate the model over the specified number of time steps
  for (t in 1:time_steps) {
    X = ricker_step_1d(X, L, kappa, r, sigma)
    lattice_matrix[, t + 1] = X
  }
  
  return(lattice_matrix)
}

#==============================================================================
#Ising model
#==============================================================================
# Function to update a single spin
update_spin <- function(spin, neighbors, J, K) {
  # Compute the effective field
  h_eff = K * spin + J * sum(neighbors)
  
  # Probability of flipping the spin
    prob_flip = exp(h_eff*(-spin) ) / (2*cosh(h_eff))

  # Update the spin with probability prob_flip
  new_spin = ifelse(runif(1) < prob_flip, -spin, spin)
  
  return(new_spin)
}

# Function to perform a single time step evolution of the Ising model
evolve_ising_model = function(grid, J, K) {
  # Get the dimensions of the grid
  grid=t(grid)
  n_rows = 1
  n_cols = ncol(grid)
  
  # Create a copy of the grid to store the updated spins
  new_grid = matrix(0, n_rows, n_cols)
  
  # Loop over each spin in the grid
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      # Get the indices of the nearest neighbors (assuming periodic boundary conditions)
      neighbors_i = c( i, i)
      neighbors_j = c(j - 1, j + 1)
      neighbors_j = ifelse(neighbors_j < 1, n_cols, neighbors_j)
      neighbors_j = ifelse(neighbors_j > n_cols, 1, neighbors_j)
      
      # Extract the spins of the neighbors
      neighbors = grid[ cbind(neighbors_i, neighbors_j) ]
      
      # Update the spin using a probability function
      new_grid[i, j] = update_spin(grid[i, j], neighbors, J, K)
    }
  }
  
  return(new_grid)
}

#==============================================================================
#Inference
#==============================================================================
# Function to simulate the dynamical Ising model over a specified number of time steps
simulate_ising_model = function(initial_grid, J, K, num_steps) {
  current_grid = initial_grid
  
  for (step in 1:(num_steps-1) ) {
    current_grid[step+1,] = evolve_ising_model(current_grid[step,], J, K)
  }
  
  return(current_grid)
}

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