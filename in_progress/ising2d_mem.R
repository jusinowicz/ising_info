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
evolve_ising_model <- function(grid, J, K, temperature) {
  # Get the dimensions of the grid
  n_rows = nrow(grid)
  n_cols = ncol(grid)
  
  # Create a copy of the grid to store the updated spins
  new_grid = matrix(0, n_rows, n_cols)
  
  # Loop over each spin in the grid
  for (i in 1:n_rows) {
    for (j in 1:n_cols) {
      # Get the indices of the nearest neighbors (assuming periodic boundary conditions)
      neighbors_i = c(i - 1, i + 1, i, i)
      neighbors_j = c(j, j, j - 1, j + 1)
      neighbors_i = ifelse(neighbors_i < 1, n_rows, neighbors_i)
      neighbors_i = ifelse(neighbors_i > n_rows, 1, neighbors_i)
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
  
  for (step in 1:num_steps) {
    current_grid = evolve_ising_model(current_grid, J, K, temperature)
  }
  
  return(current_grid)
}

# Example usage
# Set the parameters
J = 1      # Nearest neighbor coupling
K = 0.5    # Self-interaction strength
temperature = 1.0  # Temperature

# Set the size of the grid
n_rows = 10
n_cols = 10

# Initialize the grid with random spins (-1 or 1)
initial_grid = matrix(sample(c(-1, 1), n_rows * n_cols, replace = TRUE), n_rows, n_cols)

# Number of time steps to simulate
num_steps = 100

# Simulate the dynamical Ising model
final_grid = simulate_ising_model(initial_grid, J, K, temperature, num_steps)

# Display the initial and final grids
print("Initial Grid:")
print(initial_grid)

print("Final Grid:")
print(final_grid)
