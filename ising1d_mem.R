
#==============================================================================
#Libraries
#==============================================================================
source("./info_theory_functions/inference.R")
#==============================================================================

#==============================================================================
#Define some functions for the simulation
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

# Function to simulate the dynamical Ising model over a specified number of time steps
simulate_ising_model = function(initial_grid, J, K, num_steps) {
  current_grid = initial_grid
  
  for (step in 1:(num_steps-1) ) {
    current_grid[step+1,] = evolve_ising_model(current_grid[step,], J, K)
  }
  
  return(current_grid)
}

#==============================================================================
#Main body
#==============================================================================
# Set the parameters (these are dimensionless, scaled to temp?)
J = 0.2    # Nearest neighbor coupling
K = 1    # Self-interaction strength

# Set the size of the grid
#n_rows = 25
n_cols = 50

# Number of time steps to simulate
num_steps = 400

# Initialize the grid with random spins (-1 or 1)
initial_grid = matrix(sample(c(-1, 1), n_cols, replace = TRUE), 1, n_cols)
initial_grid = rbind(initial_grid, matrix(0,num_steps,n_cols) )

# Simulate the dynamical Ising model
final_grid = simulate_ising_model(initial_grid, J, K, num_steps)

# Display the initial and final grids
par(mfrow=c(2,1))
plot3D::image2D(initial_grid, main = "Initial Grid")
plot3D::image2D(final_grid, main = "Final Grid")

#==============================================================================
#Inference
#==============================================================================
#Try inferring the parameters from the final spatial grid
initial_params = c(0.5, 1.1)
parameter_estimates = fit_params_nlm(initial_params, final_grid)

#Remove highly uncertain samples with large SE on the estimates: 
pe_keep = parameter_estimates[parameter_covariances[,1] <
            matrix(1,nrow(parameter_covariances), 1 ),]
se_keep = parameter_covariances[parameter_covariances[,1] <
            matrix(1,nrow(parameter_covariances), 1 ),]

# Compute average parameter estimates and standard errors over all pairs
ave_pe = apply(pe_keep, 2, mean,na.rm=T)
ave_se = apply(se_keep, 2, mean,na.rm=T)

# Print results
cat("Average Parameter Estimates:\n")
cat("J:", ave_pe[1], "\n")
cat("K:", ave_pe[2], "\n\n")

cat("Average Standard Errors:\n")
cat("Standard Error of J:", ave_se[1], "\n")
cat("Standard Error of K:", ave_se[2], "\n")
