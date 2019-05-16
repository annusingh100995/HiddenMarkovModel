


# RIBOFLAVIN_MODEL




L_riboflavin = rbind(c(1,1,0,0,0,0,0,0,0,0), c(0,0,1,0,0,0,0,0,0,0), c(0,0,1,0,0,0,0,0,0,0),
                     c(0,0,0,1,1,0,0,0,0,0), c(0,0,0,0,0,1,0,0,0,0), c(0,0,0,0,0,1,0,0,0,0),
                     c(0,0,0,0,0,0,1,1,0,0), c(0,0,0,0,0,0,0,0,1,0), c(0,0,0,0,0,0,0,0,1,0))


R_riboflavin = rbind(c(0,0,1,0,0,0,0,0,0,0), c(0,1,0,0,0,0,0,0,0,0), c(0,1,0,1,0,0,0,0,0,0),
                     c(0,0,0,0,0,1,0,0,0,0), c(0,0,0,0,1,0,0,0,0,0), c(0,0,0,0,1,0,1,0,0,0),
                     c(0,0,0,0,0,0,0,0,1,0), c(0,0,0,0,0,0,0,1,0,0), c(0,0,0,0,0,0,0,1,0,1))

v_riboflavin = R_riboflavin - L_riboflavin

riboflavin_run_1 = c(200,100,0,0,100,0,0,100,0,0)
riboflavin_run_2 = c(2000,1000,0,0,1000,0,0,1000,0,0)
riboflavin_run_3 = c(2000,1000,0,500,1000,0,500,1000,0,0)


stochastic_rate_constant_riboflavin = c( 0.001660302,0.000100000,0.100000000,
                                         0.001660302,0.000100000,0.100000000,
                                         0.001660302,0.000100000,0.100000000)

Reaction_list_riboflavin = list(L_riboflavin, R_riboflavin)

#### ssa 
current_sys_state = riboflavin_run_3

updated_system = list() 
updated_system[[1]] = SystemState
updated_system[[1]] = current_sys_state

#updated_system[[1]] = rbind(c(25000,25000,0,0),c(0,0,30000,0),c(0,0,1,0))


#the reactions
L = Reaction_list_riboflavin[[1]]
R = Reaction_list_riboflavin[[2]]

# the state change vector
V = R-L

#number of reactions
M = nrow(L)
number_of_reaction = nrow(L)
number_of_species = ncol(L)

#vector to store propensity
propensity_vector = c("numeric", M)

#vector to store the probability for each reaction 
reaction_prob = c()

# vector to store the reaction firing
firing_reaction = c()

#list to store the updated system at each simulation time
#updated_system = array(0, c(number_of_reaction, N, NULL))



# time vector to store the time of simulation
time_vector = c()
change = c()

time= 0
# tau in SSA, delta t
tau = 0 

#the initial stage of the system
x = L
#updated_system[[1]] = L
#updated_system[[1]] = matrix(c(250:500), nrow = 3, ncol = 4)

# updated_system[[1]] = the number of molecules for each species


#system_state_at_timestep = list()
# variable to store number of iterations performed
i = 1
Simulation_time = 10000
while(tau < Simulation_time){
  
  propensity_vector = Calculate_Propensity(updated_system[[i]],L_riboflavin,stochastic_rate_constant_riboflavin)
  print(propensity_vector)
  #propensity_vector =  runif(3,0,3)
  propensity_sum = sum(propensity_vector)
  reaction_prob = propensity_vector/propensity_sum
  
  #sampling delta t
  #time = rexp(1,propensity_sum)
  time = exp(propensity_sum^ -1)
  if((tau+time) > Simulation_time){
    break
  }
  
  # random sampling of the reaction firing
  firing_reaction[i] = sample(c(1:number_of_reaction),1,FALSE,reaction_prob)
  
  tau = time +tau 
  time_vector[i] = tau
  #print(tau)
  
  change = updated_system[[i]]+V[firing_reaction[i],]
  
  updated_system[[i+1]] = change
  #updated_system[[i+1]][firing_reaction[i],] = change
  
  
  #change =  as.vector(updated_system[[i]][firing_reaction[i],]+V[firing_reaction[i],])
  #change = as.vector(updated_system[firing_reaction[i],,i]+V[firing_reaction[i],])
  #updated_system[firing_reaction[i],,i+1] =     
  i = i+1 
}
###############################################################################################
#Plotting the graph
##
Riboflavin_model_molecules = matrix(0,ncol(L_riboflavin),length(time_vector-1))
for(t in 1:length(time_vector-1)){
  for(i in 1:ncol(L_riboflavin)){
    # takig the number of moecules from the first change 
    Riboflavin_model_molecules[i,t] = updated_system[[t]][i]
  }
}

for(i in 1:nrow(Riboflavin_model_molecules)){
  #par(mfrow=c(5,2))
  plot(Riboflavin_model_molecules[i,] ~ time_vector, pch=16, cex=.65, col=i,
       main=paste('Species', i, ' vs simulation time'))
}
##################################################################################################


for(i in 1:nrow(Riboflavin_model_molecules)){
  #par(mfrow=c(5,2))
  plot(Riboflavin_model_molecules[i,] ~ time_vector, type="l", col=i,
       main=paste('Species', i, ' vs simulation time'))
}
###





