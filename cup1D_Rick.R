#==============================================================================
# Example R code for a discrete Ricker 
# metapopulation model implemented on a 1D coupled map lattice. 
#==============================================================================
source("./info_theory_functions/info_theory_functions.R")
source("./info_theory_functions/inference.R")
#==============================================================================
# Example usage for a 1D lattice
#
# simple 2-cycle oscillations with no dispersal or noise: 
# kappa =0 
# r = 2.1 
# sigma = 0 
#==============================================================================
# Parameters
L = 32
time_steps = 1000
kappa =  0.26
r = 2.1 #7.1
sigma = 0#0.16

#Run the model
final_states = simulate_ricker_1d(L, kappa, r, sigma, time_steps)

# Plot the final population distribution
par(mfrow =c(3,1))
hist(final_states, xlab = "Population Density", col = "lightblue", border = "black")

#Plot population trajectory of a single site: 
xi = 16
plot(final_states[xi,],t="l")

#First difference of each oscillator
fd_states = apply(final_states, 1, diff)

#Plot first difference trajectory of a single site: 
xi = 16
plot(fd_states[xi,],t="l")

#==============================================================================
# Classic statistical quantities.
#==============================================================================
#Oscillator state
m_it = matrix( (-1)^(1:time_steps),time_steps, L)*fd_states / 2
m_it_sign = sign(m_it)

#Instantaneous order parameter
iop = rowSums( matrix( 1/(L),time_steps, L) * m_it )

#Synchronization order parameter
sop = 1/(time_steps - time_steps/2 -1 ) * 
          sum(iop[ (time_steps/2):time_steps] )

#==============================================================================
#Dynamic information theory metrics
#==============================================================================
k =5 #Temporal history length
nlevel = 64

#Information storage (excess entropy):
#ee_lat = get_ee1D(t(ceiling(final_states)),k=k)
ais_lat = get_ais1D(t(ceiling(final_states)),k=k)
te_left_lat = get_te1D_left(t(ceiling(final_states)),k=k)
te_right_lat = get_te1D_right(t(ceiling(final_states)),k=k)
si_left_lat = get_si1D_left(t(ceiling(final_states)),k=k)
si_right_lat = get_si1D_right(t(ceiling(final_states)),k=k)

#Plot the information theory
z1 = 25:30
z2 = 750:760
z1 = 1:L
z2 = 1:time_steps
z1 = 1:L
z2 = 1:30
#plot3D::image2D((ee_lat)) 
#plot3D::image2D(t(ee_lat))
fig.name = paste("cup1d_rick_Donly",".pdf",sep="")
#fig.name = paste("cup1d_rick_noDnoN",".pdf",sep="")
#fig.name = paste("cup1d_rick_Nonly",".pdf",sep="")
#fig.name = paste("cup1d_rick_both",".pdf",sep="")


pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow = c(2,3))

plot3D::image2D((final_states[z1,z2]), main = "Population", ylab = "Time", 
  xlab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(ais_lat[z2,z1]),180), main = "Information Storage", ylab = "Time", 
  xlab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(te_left_lat[z2,z1]),180),main = "Transfer L <----", ylab = "Time", 
  xlab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(te_right_lat[z2,z1]),180),main = "Transfer R ---->", ylab = "Time", 
  xlab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(si_left_lat[z2,z1]),180), main = "Modification L <----", ylab = "Time", 
  xlab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5) 
plot3D::image2D(rotateMatrix(t(si_right_lat[z2,z1]),180),main = "Modification R ---->", ylab = "Time", 
  xlab = "Space",col=viridis(nlevel),cex.main=1.3,cex.lab=1.5)

dev.off()



#Dynamic information theory metrics on m_it
#Information storage (excess entropy):
#ee_lat = get_ee1D(t(ceiling(final_states)),k=5)
ais_fd = get_ais1D((ceiling(fd_states)),k=k)
te_left_fd = get_te1D_left((ceiling(fd_states)),k=k)
te_right_fd = get_te1D_right((ceiling(fd_states)),k=k)
si_left_fd = get_si1D_left((ceiling(fd_states)),k=k)
si_right_fd = get_si1D_right((ceiling(fd_states)),k=k)

#Plot the information theory
par(mfrow = c(2,3))
z1 = 25:30
z2 = 750:760
z1 = 1:L
z2 = 1:time_steps
#plot3D::image2D((ee_lat)) 
#plot3D::image2D(t(ee_lat))
plot3D::image2D(t(fd_states[z2,z1])) 
plot3D::image2D(t(ais_fd[z2,z1])) 
plot3D::image2D(t(te_left_fd[z2,z1])) 
plot3D::image2D(t(te_right_fd[z2,z1])) 
plot3D::image2D(t(si_left_fd[z2,z1])) 
plot3D::image2D(t(si_right_fd[z2,z1]))

#Dynamic information theory metrics on the sign of m_it
#Information storage (excess entropy):
#ee_lat = get_ee1D(t(ceiling(final_states)),k=k)
ais_mits = get_ais1D((ceiling(m_it_sign)),k=k)
te_left_mits = get_te1D_left((ceiling(m_it_sign)),k=k)
te_right_mits = get_te1D_right((ceiling(m_it_sign)),k=k)
si_left_mits = get_si1D_left((ceiling(m_it_sign)),k=k)
si_right_mits = get_si1D_right((ceiling(m_it_sign)),k=k)

#Plot the information theory
par(mfrow = c(2,3))
z1 = 25:30
z2 = 750:760
z1 = 1:L
z2 = 1:time_steps
#plot3D::image2D((ee_lat)) 
#plot3D::image2D(t(ee_lat))
plot3D::image2D(t(m_it_sign[z2,z1])) 
plot3D::image2D(t(ais_mits[z2,z1])) 
plot3D::image2D(t(te_left_mits[z2,z1])) 
plot3D::image2D(t(te_right_mits[z2,z1])) 
plot3D::image2D(t(si_left_mits[z2,z1])) 
plot3D::image2D(t(si_right_mits[z2,z1]))


