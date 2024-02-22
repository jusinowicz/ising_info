# Function to update a single spin using Glauber dynamics
update_spin <- function(spin, neighbors, J, K, temperature) {
  # Compute the effective field
  h_eff = K * spin + J * sum(neighbors)
  
  # Probability of flipping the spin
  # prob_flip = 1 / (1 + exp(-2 * h_eff / temperature))
    prob_flip = exp(h_eff*(-spin) ) / 2*cosh(1 + exp(-2 * h_eff / temperature))

  # Update the spin with probability prob_flip
  new_spin = ifelse(runif(1) < prob_flip, -spin, spin)
  
  return(new_spin)
}

# Function to perform a single time step evolution of the Ising model
evolve_ising_model = function(grid, J, K, temperature) {
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
      new_grid[i, j] = update_spin(grid[i, j], neighbors, J, K, temperature)
    }
  }
  
  return(new_grid)
}

# Function to simulate the dynamical Ising model over a specified number of time steps
simulate_ising_model = function(initial_grid, J, K, temperature, num_steps) {
  current_grid = initial_grid
  
  for (step in 1:(num_steps-1) ) {
    current_grid[step+1,] = evolve_ising_model(current_grid[step,], J, K, temperature)
  }
  
  return(current_grid)
}

# Example usage
# Set the parameters
J = 1      # Nearest neighbor coupling
K = 0.5    # Self-interaction strength
temperature = 1.0  # Temperature

# Set the size of the grid
#n_rows = 25
n_cols = 25

# Number of time steps to simulate
num_steps = 100

# Initialize the grid with random spins (-1 or 1)
initial_grid = matrix(sample(c(-1, 1), n_cols, replace = TRUE), 1, n_cols)
initial_grid = rbind(initial_grid, matrix(0,num_steps,n_cols) )

# Simulate the dynamical Ising model
final_grid = simulate_ising_model(initial_grid, J, K, temperature, num_steps)

# Display the initial and final grids
par(mfrow=c(2,1))
plot3D::image2D(initial_grid, main = "Initial Grid")
plot3D::image2D(final_grid, main = "Final Grid")