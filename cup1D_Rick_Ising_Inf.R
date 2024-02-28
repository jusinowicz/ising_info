#==============================================================================
# This code is to run a 1D Ricker model, infer the matching
# 1D Ising parameters J (interaction strength) and K (self-interaction), 
# simulate the resulting Ising model,
# then run the Dynamic Information Theoretic metrics on both sets 
# and compare their behavior.
#==============================================================================
source("./info_theory_functions/info_theory_functions.R")
source("./info_theory_functions/inference.R")
library(dplyr)
library(tidyverse)
library(viridis)
#==============================================================================
#Run the Ricker model
#==============================================================================
# simple 2-cycle oscillations with no dispersal or noise: 
# kappa =0 
# r = 2.1 
# sigma = 0 
#==============================================================================
# Parameters
#Width
L = 50
time_steps = 500
#Dispersal
kappa = 0# 0.26
#Intrinsic growth
r = 2.1 #7.1
#Noise
sigma = 0.16

#Run the model
final_states = simulate_ricker_1d(L, kappa, r, sigma, time_steps)

# Plot the final population distribution
par(mfrow =c(3,1))
hist(final_states, xlab = "Population Density", col = "lightblue", border = "black")

#Plot population trajectory of a single site: 
xi = 16
plot(final_states[xi,],t="l")

#First difference of each oscillator
fd_states = sign(apply(final_states, 1, diff))

#==============================================================================
#Ising Inference
#==============================================================================
#Try inferring the parameters from the final spatial grid of differences, which
#represent the oscillator state.

#Initial parameter guesses for c(J,K)
initial_params = c(1, 0.1)
pe = fit_params_nlm(initial_params, fd_states)

#Remove highly uncertain samples with large SE on the estimates: 
pe_keep = pe$estimate[pe$cov[,1] <
            matrix(1,nrow(pe$cov), 1 ),]
se_keep = pe$cov[pe$cov[,1] <
            matrix(1,nrow(pe$cov), 1 ),]

# Compute average parameter estimates and standard errors over all pairs
ave_pe = apply(pe_keep, 2, mean,na.rm=T)
ave_se = apply(se_keep, 2, mean,na.rm=T)

#==============================================================================
#Run the Ising model with parameters
#==============================================================================
# Set the parameters (these are dimensionless, scaled to temp?)
J = ave_pe[1]    # Nearest neighbor coupling
K = ave_pe[2]    # Self-interaction strength

# Set the size of the grid
n_cols = L

# Number of time steps to simulate
num_steps = time_steps

# Initialize the grid with random spins (-1 or 1)
initial_grid = matrix(sample(c(-1, 1), n_cols, replace = TRUE), 1, n_cols)
initial_grid = rbind(initial_grid, matrix(0,num_steps,n_cols) )

# Simulate the dynamical Ising model
final_grid = simulate_ising_model(initial_grid, J, K, num_steps)

#Plot the Ising model, the oscillator state of the Ricker model,
#and the model itself to compare.

par(mfrow=c(3,1))
plot3D::image2D(final_grid, main = "Ising")
plot3D::image2D(fd_states, main = "Ricker Oscillator")
plot3D::image2D(final_states, main = "Ricker Populations")

#==============================================================================
#Dynamic information theory (DIT) to get computational
#dynamics of all three lattices for comparison. 
#==============================================================================
k =5 #Temporal history length
nlevel = 64

#Ricker model
#Information storage (excess entropy):
#ee_lat = get_ee1D(t(ceiling(final_states)),k=k)
ais_Rick = get_ais1D(t(ceiling(final_states)),k=k)
te_left_Rick = get_te1D_left(t(ceiling(final_states)),k=k)
te_right_Rick = get_te1D_right(t(ceiling(final_states)),k=k)
si_left_Rick = get_si1D_left(t(ceiling(final_states)),k=k)
si_right_Rick = get_si1D_right(t(ceiling(final_states)),k=k)

#Ricker Oscillator
ais_osc = get_ais1D(t(ceiling(fd_states)),k=k)
te_left_osc = get_te1D_left(t(ceiling(fd_states)),k=k)
te_right_osc = get_te1D_right(t(ceiling(fd_states)),k=k)
si_left_osc = get_si1D_left(t(ceiling(fd_states)),k=k)
si_right_osc = get_si1D_right(t(ceiling(fd_states)),k=k)

#Ising 
ais_ising = get_ais1D(t(ceiling(final_grid)),k=k)
te_left_ising = get_te1D_left(t(ceiling(final_grid)),k=k)
te_right_ising = get_te1D_right(t(ceiling(final_grid)),k=k)
si_left_ising = get_si1D_left(t(ceiling(final_grid)),k=k)
si_right_ising = get_si1D_right(t(ceiling(final_grid)),k=k)


#Plot the information theory
# z1 = 25:30
# z2 = 750:760

# z1 = 1:L
# z2 = 1:time_steps

s_upper = L
t_upper = L

z1 = 1:t_upper
z2 = 1:s_upper

tick = list(x = seq(0, t_upper, by = t_upper/10), 
			y = seq(0, s_upper, by = s_upper/5))

#plot3D::image2D((ee_lat)) 
#plot3D::image2D(t(ee_lat))
#fig.name = paste("cup1d_rick_Donly",".pdf",sep="")
#fig.name = paste("cup1d_rick_noDnoN",".pdf",sep="")
#fig.name = paste("cup1d_rick_Nonly",".pdf",sep="")

fig.name = paste("rick_ising_Nonly_ais",".pdf",sep="")
#fig.name = paste("rick_ising_both_ais",".pdf",sep="")

pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow = c(1,3))
global_range = range(ais_osc , ais_ising)
plot3D::image2D(rotateMatrix(t(ais_osc [z2,z1]),180), zlim = global_range, main = "AIS Oscillator", xlab = "Time", 
  ylab = "Space",col=viridis(nlevel), cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(ais_ising[z2,z1]),180), zlim = global_range,main = "AIS Ising", xlab = "Time", 
  ylab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(ais_Rick[z1,z2]),180), main = "AIS Ricker", xlab = "Time", 
  ylab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 


dev.off()


fig.name = paste("rick_ising_Nonly_teleft",".pdf",sep="")
#fig.name = paste("rick_ising_both_teleft",".pdf",sep="")

pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow = c(1,3))
global_range = range(ais_osc , ais_ising)
plot3D::image2D(rotateMatrix(t(te_left_osc [z2,z1]),180), zlim = global_range, main = "TE Oscillator", xlab = "Time", 
  ylab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(te_left_ising[z2,z1]),180), zlim = global_range,main = "TE Ising", xlab = "Time", 
  ylab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(te_left_Rick[z1,z2]),180), main = "TE Ricker", xlab = "Time", 
  ylab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 

dev.off()


fig.name = paste("rick_ising_Nonly",".pdf",sep="")
#fig.name = paste("rick_ising_both",".pdf",sep="")

pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow = c(1,3))
global_range = range(ais_osc , ais_ising)
plot3D::image2D(rotateMatrix(t(fd_states [z2,z1]),180), zlim = global_range, main = "Oscillator", xlab = "Time", 
  ylab = "Space",col=viridis(nlevel), cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(final_grid[z2,z1]),180), zlim = global_range,main = "Ising", xlab = "Time", 
  ylab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(final_states[z1,z2]),180), main = "Ricker", xlab = "Time", 
  ylab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 

dev.off()