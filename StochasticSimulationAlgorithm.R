

# Defining Parameters

#Species in order = E,S,ES,P

volume = 1e-15
na = 6.023e23


reaction_rate_constant = c()
#reaction_rate[1] = 1e6/na*v
#reaction_rate[2] = 1e-4
#reaction_rate[3] = 0.1/na*v
reaction_rate_constant[1] = 1e6
reaction_rate_constant[2] = 1e-4
reaction_rate_constant[3] = 0.1


stochastic_rate_constant = c()

stochastic_rate_constant[1] = reaction_rate_constant[1]/(volume*na)
stochastic_rate_constant[2] = reaction_rate_constant[2]
stochastic_rate_constant[3] = reaction_rate_constant[3]

intial_conc = c()
#E 
intial_conc[1] = 5e-7
#S
intial_conc[2] = 2e-7
#ES
intial_conc[3] = 0
#P
intial_conc[4] = 0

number_of_molecules = c()
#E 
number_of_molecules[1]= intial_conc[1]*volume*na
#S
number_of_molecules[2]= intial_conc[2]*volume*na
#ES
number_of_molecules[3]= intial_conc[3]*volume*na
#P
number_of_molecules[4]= intial_conc[4]*volume*na

number_of_molecules = round(number_of_molecules)

run_1 = number_of_molecules
run_2 = number_of_molecules*100
run_3 = number_of_molecules*200
run_4 = c(2000,1000,0,0)
run_5 = c(20000,10000,0,0)
# N = number of species in the model
# M = number of reactions in the model

# L = Intial state before a reaction fires
L = matrix(0, nrow = M , ncol = N)

# R = Final state after a reaction fires
R = matrix(0, nrow = M, ncol = N)

# V = state change reaction 
V = matrix(0, nrow = M, ncol = N)

k = vector("numeric", length = M)



L1 = rbind(c(1,1,0,0),c(0,0,1,0),c(0,0,1,0))
R1 = rbind(c(0,0,1,0), c(1,0,0,0), c(1,0,0,1))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
V1 = R1 - L1
#X = c(x1,x2,x3,x4)
rate_kinetic = c()

Reaction_list_1 = list()

Reaction_list = list(L1, R1)
#####################################################################################################################
######################################################################################################################
#input = current_state , a vetor containg thr number of all elements
# L : the left hand sid fo the reactions
# ReactionRate = vector containing the the reaction rates/stochastic rate constants

Calculate_Propensity = function(CurrentState, L, ReactionRate){
  
  
  # to convert the number of moelcues given in current states as a vector into a matrix
  # By converting it into a matrix, I can calculate prpensity for each reaction.
  
  # Here i'm converting the current state of the system in a matrix c(a,b,c,d)
  # will be converted into [a,b,0,0]
  #                        [0,0,c,0]            
  #                        [0,0,0,d]
  # For reaction 1: the propensity is c1*a*b
  # For reaction 2: the propensity is c2*c
  # For reaction 3: the propensity is c3*d
  state_count_matrix = matrix(0,nrow(L), ncol(L))
  for(i in 1:nrow(L)){
    for(j in 1:ncol(L)){
      state_count_matrix[i,] = CurrentState*L[i,]
    }
  }
  
  temp_vector = c()
  temp_vect_for_prod = c()
  propensity_test = c()  
  
  for(i in 1:nrow(L)){
    temp_vector = state_count_matrix[i,]
    
    temp_vect_for_prod = temp_vector[temp_vector != 0]
    
    propensity_test[i] = prod(temp_vect_for_prod)*ReactionRate[i] 
    
    
  }
  return(propensity_test)  
}
###################################################################################################################
test_prop = Calculate_Propensity(run_1,L1,stochastic_rate_constant)
###################################################################################################################### 
#updated_system = list()  
##########################################################################################
# giev the initila number of molecues for each species

#updated_system[[1]] = rbind(c(25000,25000,0,0),c(0,0,30000,0),c(0,0,1,0))

#Input = SystemState = a vector containgn the initial valus of each species
#Reaction list:  the L and R reaction matrics
# simulaltion time ; the time for simultion

#SSA = function(SystemState,Reaction_list,Simulation_time, stochastic_rate_constant){
  
  #current_sys_state = SystemState
  current_sys_state = run_5
  
  updated_system = list() 
  updated_system[[1]] = SystemState
  updated_system[[1]] = current_sys_state
  
  #updated_system[[1]] = rbind(c(25000,25000,0,0),c(0,0,30000,0),c(0,0,1,0))
  
  
  #the reactions
  L = Reaction_list[[1]]
  R = Reaction_list[[2]]
  
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
  Simulation_time = 500
  while(tau < Simulation_time){
    
    propensity_vector = Calculate_Propensity(updated_system[[i]],L,stochastic_rate_constant)
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
#  return(list(updated_system, time_vector, firing_reaction))
#}
#######################################################################################################################


PLOT_NUMBER_OF_MOLECULES = function(UpdatedStateForEachTime, SimulationTimeVector, NumberOfSpecies){
  
  
  # vector to store the number of molecues for each species
  Number_of_molecules = matrix(0,nrow = NumberOfSpecies , ncol = length(SimulationTimeVector))
  
  # for loop to coun the number of molecule for each species
  
  for(t in 1:length(SimulationTimeVector-1)){
    for(i in 1: NumberOfSpecies){
      # takig the number of moecules from the first change 
      Number_of_molecules[i,t] = UpdatedStateForEachTime[[t]][i]
    }
  }
  
  
  # plotting each species
  # for all in one plot
  #par(mfrow=c(NumberOfSpecies,1))
  
  # for single plots
  for(i in 1:nrow(Number_of_molecules)){
    plot(Number_of_molecules[i,] ~ SimulationTimeVector, pch=16, cex=.65, col=i,
         main=paste('Species', i, ' vs simulation time'))
  }
  
  
}
PLOT_NUMBER_OF_MOLECULES(updated_system, time_vector,4)
#######################################################################################################################
# Note: the SSA function gave error while running, I ran all simulation by feeding in parametrs manually.

#Main Code


PLOT_NUMBER_OF_MOLECULES(run_1_state,run_1_time,4)


PLOT_NUMBER_OF_MOLECULES(updated_system,time_vector,4)





