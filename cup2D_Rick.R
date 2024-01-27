#==============================================================================
# Example R code for a discrete Ricker 
# metapopulation model implemented on a 2D coupled map lattice. 
#==============================================================================
# Function to simulate one time step of the Ricker metapopulation model
ricker_step = function(X, L, kappa, r, sigma) {
  # Define the Ricker operator
  R = function(X) X * exp(r * (1 - X))
  
  # Define the dispersal operator
  D = function(X, kappa, neighbors) {
    (1 - kappa) * X + kappa * sum(neighbors) / length(neighbors)
  }
  
  # Define the environmental noise operator
  N = function(X, sigma) X * exp(sigma * rnorm(1))
  
  # Initialize the updated lattice
  X_new = matrix(0, nrow = L, ncol = L)
  
  # Loop through each patch on the lattice
  for (i in 1:L) {
    for (j in 1:L) {
      # Get the neighboring patches
      neighbors = c(X[i, j], X[i, ifelse(j - 1 > 0, j - 1, L)], 
                     X[i, ifelse(j + 1 <= L, j + 1, 1)],
                     X[ifelse(i - 1 > 0, i - 1, L), j],
                     X[ifelse(i + 1 <= L, i + 1, 1), j])
      
      # Update the patch based on the Ricker model with dispersal and noise
      X_new[i, j] = N(D(R(X[i, j]), kappa, neighbors), sigma)
    }
  }
  
  return(X_new)
}

# Function to simulate the Ricker metapopulation model over multiple time steps
simulate_ricker = function(L, kappa, r, sigma, time_steps) {
  # Initialize the lattice with random initial population densities
  X = matrix(runif(L * L), nrow = L, ncol = L)
  
  # Initialize a 3D array to store the lattice at each time step
  lattice_array = array(0, dim = c(L, L, time_steps + 1))
  lattice_array[, , 1] = X
  
  # Simulate the model over the specified number of time steps
  for (t in 1:time_steps) {
    X = ricker_step(X, L, kappa, r, sigma)
    lattice_array[, , t + 1] = X
  }
  
  return(lattice_array)
}

# Parameters
L = 32
time_steps = 100
kappa = 0.15
r = 7.2
sigma = 0.26

final_states = simulate_ricker(L, kappa, r, sigma, time_steps)

# Reset the plotting layout
par(mfrow = c(1, 1))

# Plot the lattice as an image
plot3D::image2D(final_states[, , time_steps]) 

# Plot the final population distribution
hist(final_states[, , time_steps], main = paste("Time Step:", t), xlab = "Population Density", col = "lightblue", border = "black")

#Plot population trajectory of a single site: 
xi = 16
yi = 16
plot(final_states[xi,yi,],t="l")