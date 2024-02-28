
#==============================================================================
#Libraries
#==============================================================================
source("./info_theory_functions/inference.R")
#==============================================================================

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
pe = fit_params_nlm(initial_params, final_grid)

#Remove highly uncertain samples with large SE on the estimates: 
pe_keep = pe$estimate[pe$cov[,1] <
            matrix(1,nrow(pe$cov), 1 ),]
se_keep = pe$cov[pe$cov[,1] <
            matrix(1,nrow(pe$cov), 1 ),]

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
