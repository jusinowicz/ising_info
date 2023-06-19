#==============================================================================
#Classic 2D Ising model. 
# 
# I have assumed a square lattice with periodic boundary conditions.
# You can adjust the parameters n, temperature, J, and num_steps as needed. 
# The final spin configuration is stored in the final_spins variable, and 
# the energies at each step are stored in the final_energies variable.
#
#
#==============================================================================


# Function to initialize the spin configuration
initialize_spins = function(n) {
  spins = matrix(1, nrow = n, ncol = n)
  spins[sample(1:(n*n), n*n/2)] = -1
  return(spins)
}

# Function to calculate the energy of the spin configuration
calculate_energy = function(spins, J) {
  n = nrow(spins)
  energy = 0
  for (i in 1:n) {
    for (j in 1:n) {
      energy = energy - J * spins[i, j] * (
        spins[i %% n + 1, j] + spins[i, j %% n + 1]
      )
    }
  }
  return(energy)
}

# Function to perform a Monte Carlo update using the Metropolis algorithm
monte_carlo_update = function(spins, temperature, J) {
  n = nrow(spins)
  i = sample(1:n, 1)
  j = sample(1:n, 1)
  delta_energy = 2 * J * spins[i, j] * (
    spins[i %% n + 1, j] + spins[(i + n - 2) %% n + 1, j] +
    spins[i, j %% n + 1] + spins[i, (j + n - 2) %% n + 1]
  )
  
  if (delta_energy <= 0 || runif(1) < exp(-delta_energy / temperature)) {
    spins[i, j] = -spins[i, j]
  }
  
  return(spins)
}

# Function to simulate the Ising model
ising_model_simulation = function(n, temperature, J, num_steps) {
  spins = initialize_spins(n)
  energy = calculate_energy(spins, J)
  
  energies = numeric(num_steps)
  
  for (step in 1:num_steps) {
    spins = monte_carlo_update(spins, temperature, J)
    energy = calculate_energy(spins, J)
    energies[step] = energy
  }
  
  return(list(spins = spins, energies = energies))
}

# Set the simulation parameters
n = 20  # Size of the spin lattice
temperature = 2.0  # Temperature of the system
J = 1.0  # Interaction strength
num_steps = 1000  # Number of Monte Carlo steps

# Run the simulation
simulation_result = ising_model_simulation(n, temperature, J, num_steps)

# Access the final spin configuration and energies
final_spins = simulation_result$spins
final_energies = simulation_result$energies
