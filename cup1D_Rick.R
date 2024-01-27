# Function to simulate one time step of the Ricker metapopulation model for a 1D lattice
ricker_step_1d <- function(X, L, kappa, r, sigma) {
  # Define the Ricker operator
  R <- function(X) X * exp(r * (1 - X))
  
  # Define the dispersal operator
  D <- function(X, kappa, neighbors) {
    (1 - kappa) * X + kappa * sum(neighbors) / length(neighbors)
  }
  
  # Define the environmental noise operator
  N <- function(X, sigma) X * exp(sigma * rnorm(1))
  
  # Initialize the updated lattice
  X_new <- numeric(L)
  
  # Loop through each patch on the lattice
  for (i in 1:L) {
    # Get the neighboring patches
    neighbors <- c(X[i], X[ifelse(i - 1 > 0, i - 1, L)], X[ifelse(i + 1 <= L, i + 1, 1)])
    
    # Update the patch based on the Ricker model with dispersal and noise
    X_new[i] <- N(D(R(X[i]), kappa, neighbors), sigma)
  }
  
  return(X_new)
}

# Function to simulate the Ricker metapopulation model over multiple time steps for a 1D lattice
simulate_ricker_1d <- function(L, kappa, r, sigma, time_steps) {
  # Initialize the lattice with random initial population densities
  X <- runif(L)
  
  # Initialize a matrix to store the lattice at each time step
  lattice_matrix <- matrix(0, nrow = L, ncol = time_steps + 1)
  lattice_matrix[, 1] <- X
  
  # Simulate the model over the specified number of time steps
  for (t in 1:time_steps) {
    X <- ricker_step_1d(X, L, kappa, r, sigma)
    lattice_matrix[, t + 1] <- X
  }
  
  return(lattice_matrix)
}

# Example usage for a 1D lattice
# Parameters
L = 32
time_steps = 100
kappa = 0.15
r = 7.2
sigma = 0.26

final_states = simulate_ricker_1d(L, kappa, r, sigma, time_steps)

# Reset the plotting layout
par(mfrow = c(1, 1))

# Plot the lattice as an image
plot3D::image2D(final_states) 

# Plot the final population distribution
hist(final_states, main = paste("Time Step:", t), xlab = "Population Density", col = "lightblue", border = "black")

#Plot population trajectory of a single site: 
xi = 16
plot(final_states[xi,],t="l")
