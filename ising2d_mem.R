# Function to initialize a random spin configuration on an L x L lattice
initialize_spins= function(L) {
  return(matrix(sample(c(-1, 1), L * L, replace = TRUE), nrow = L, ncol = L))
}

# Function to update a single spin using Glauber dynamics
glauber_update= function(current_spin, neighbors, J, K) {
  total_field= K * current_spin + J * sum(neighbors)
  probability_flip= 1 / (1 + exp(-2 * total_field))
  
  # Flip the spin with probability determined by Glauber dynamics
  new_spin= ifelse(runif(1) < probability_flip, -current_spin, current_spin)
  return(new_spin)
}

# Function to perform a single Monte Carlo step for the dynamical Ising model
monte_carlo_step= function(lattice, J, K) {
  L= nrow(lattice)
  new_lattice= lattice
  
  for (i in 1:L) {
    for (j in 1:L) {
      # Get the indices of neighboring spins (periodic boundary conditions)
      neighbors= c(lattice[i %% L + 1, j],
                     lattice[i, j %% L + 1],
                     lattice[(i - 2) %% L + 1, j],
                     lattice[i, (j - 2) %% L + 1])
      
      # Update the spin using Glauber dynamics
      new_lattice[i, j]= glauber_update(lattice[i, j], neighbors, J
