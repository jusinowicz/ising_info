#Simple function to get h_it for each site at a single time step
sum_hit1D = function (S_it){
  N = max( dim(S_it)[1], dim(S_it)[2] ) #Size of lattice
  kernel1D = c(1,1) #Nearest-neighbor kernel
  
  #Use a convolution to get the sum efficiently
  conv_result = convolve(S_it, kernel1D, 
                  type = "filter", boundary = "wrap")


  return(conv_result)

}

# #Not tested yet! 
# sum_hit2D = function (S_it){

#   N1 = dim(S_it)[1] 
#   N2 = dim(S_it)[2] #Size of lattice
#   kernel2D = matrix(c(1,1,1,1), nrow=2, byrow =T) #Nearest-neighbor kernel
  
#   #Use a convolution to get the sum efficiently
#   conv_result = conv2(S_it, kernel2D, 
#                   boundary = "wrap", conj = FALSE)


#   return(conv_result)
# }

# Define the log-likelihood function based on Eq. (A9)
log_likelihood = function(params, nk, nkf) {
  # Extract parameters J and K
  J = params[1]
  K = params[2]
  
  # Calculate Pf(k) for each bin
  Pf_k = ising1(J, K, nk)
  
  # Calculate the log-likelihood based on Eq. (A9)
  log_likelihood_value = sum(nkf * log(Pf_k) + (nk - nkf) * log(1 - Pf_k))
  
  return(log_likelihood_value)
}

# Define a function to calculate Pf(k) based on parameters J and K
ising1 = function(J, K, nk) {
  # Calculate Pf(k) based on parameters J and K
  Pf_k = exp(J * h_it + K * S_it) * S_it / (2 * cosh(J * h_it + K * S_it))
  
  return(Pf_k)
}


# Perform optimization to maximize the log-likelihood
# Example optimization function: optim()
# Example usage:
# optim(par = c(initial_guess_J, initial_guess_K), fn = log_likelihood, nk = nk_values, nk_f = nk_f_values)

# After optimization, you'll get the inferred parameters J and K
# Additionally, you can calculate the error bars using the bootstrap method
# Example code for bootstrap method:
# bootstrapped_params = function(data, num_iterations) {
#   sampled_params = replicate(num_iterations, {
#     sampled_data = sample(data, replace = TRUE)
#     # Perform optimization on sampled data
#     optim(par = c(initial_guess_J, initial_guess_K), fn = log_likelihood, nk = nk_values_sampled, nk_f = nk_f_values_sampled)$par
#   })
#   return(sampled_params)
# }
# bootstrapped_params = bootstrapped_params(data = nk_data, num_iterations = 100)
# Calculate error bars based on bootstrapped parameters

# Finally, plot the inferred results, possibly including error bars