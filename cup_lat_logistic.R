# Function to generate a noisy coupled map
generate_noisy_couple_map= function(L, iterations, alpha, beta, noise_strength, epsilon = 1e-10) {
  # Initialize the lattice
  lattice1 = matrix(runif(L * L), nrow = L, ncol = L)
  lattice2 = matrix(runif(L * L), nrow = L, ncol = L)
  
  # Iterate over time steps
  for (t in 1:iterations) {
    # Update lattice values with a small constant added to prevent -Inf
    lattice1 = alpha * (lattice1 + epsilon) * (1 - (lattice1 + epsilon)) + beta * lattice2 + noise_strength * rnorm(L * L)
    lattice2 = alpha * (lattice2 + epsilon) * (1 - (lattice2 + epsilon)) + beta * lattice1 + noise_strength * rnorm(L * L)
  }
  
  
  return(list(lattice1 = lattice1, lattice2 = lattice2))
}

# Set parameters
L= 10       # Lattice size
iterations= 100  # Number of iterations
alpha= 1.2   # Parameter alpha
beta= 0.8    # Parameter beta
noise_strength= 0.1  # Strength of noise

# Generate noisy coupled maps
result= generate_noisy_couple_map(L, iterations, alpha, beta, noise_strength)

# Access the generated lattices
lattice1= result$lattice1
lattice2= result$lattice2

# Print the generated lattices
print("Lattice 1:")
print(lattice1)
print("Lattice 2:")
print(lattice2)
