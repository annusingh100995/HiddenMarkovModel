library("abind", lib.loc="H:/R/win-library/3.4")
library("matlab", lib.loc="H:/R/win-library/3.4")
#############################################


# equilibrium probabilities
equil=function(P){
  e=eigen(t(P))$vectors[,1]
  Re(e/sum(e))
}
# computes the distance between two vectors
# returns -1 if the two vectors have different length
# the number of different elements otherwise
dist=function(r,s){
  if (length(r) != length(s)) { return(-1)}
  d = which (r!=s)
  return (length(d))
}


######################################################################################################
# Transition matrix for 3 states

P1 = rbind(c(0.6,0.1,0.1,0.2),c(0.25,0.3,0.2,0.25),c(0.1,0.6,0.2,0.1),c(0.3,0.2,0.1,0.4))

P2 = rbind(c(0.3,0.1,0.4,0.2),c(0.2,0.4,0.2,0.2),c(0.3,0.3,0.2,0.2),c(0.2,0.5,0.2,0.1))

P3 = rbind(c(0.4,0.3,0.2,0.1),c(0.3,0.4,0.1,0.2),c(0.25,0.2,0.25,0.3),c(0.1,0.4,0.2,0.3))

Transition_matrix = abind(P1,P2,P3,along = 3)

Transition_matrix_list = list(P1,P2,P3)

Segmentation_matrix = rbind(c(0.999,0.0005,0.0005),c(0.001,0.9980,0.001),c(0.005,0.005,0.99))
#########################################################################################################


######################################################################################################
# Transition matrix for 4 states

P1_4_states = rbind(c(0.6,0.1,0.1,0.2),c(0.25,0.3,0.2,0.25),c(0.1,0.6,0.2,0.1),c(0.3,0.2,0.1,0.4))

P2_4_states = rbind(c(0.3,0.1,0.4,0.2),c(0.2,0.4,0.2,0.2),c(0.3,0.3,0.2,0.2),c(0.2,0.5,0.2,0.1))

P3_4_states = rbind(c(0.4,0.3,0.2,0.1),c(0.3,0.4,0.1,0.2),c(0.25,0.2,0.25,0.3),c(0.1,0.4,0.2,0.3))

P4_4_states = rbind(c(0.3,0.3,0.2,0.2),c(0.2,0.4,0.2,0.2),c(0.2,0.3,0.2,0.3),c(0.3,0.2,0.1,0.4))

Transition_matrix_4_states = abind(P1_4_states,P2_4_states,P3_4_states,P4_4_states,along = 3)

#Transition_matrix_list = list(P1,P2,P3)

Segmentation_matrix_4_states = rbind(c(0.999,0.0004,0.0004,0.0002),c(0.0001,0.999,0.0001,0.0008),c(0.005,0.003,0.99,0.002),c(0.002,0.006,0.002,0.99))
#########################################################################################################
######################################################################################################
# Transition matrix for 5 states

P1_5_states = rbind(c(0.6,0.1,0.1,0.2),c(0.25,0.3,0.2,0.25),c(0.1,0.6,0.2,0.1),c(0.3,0.2,0.1,0.4))

P2_5_states = rbind(c(0.3,0.1,0.4,0.2),c(0.2,0.4,0.2,0.2),c(0.3,0.3,0.2,0.2),c(0.2,0.5,0.2,0.1))

P3_5_states = rbind(c(0.4,0.3,0.2,0.1),c(0.3,0.4,0.1,0.2),c(0.25,0.2,0.25,0.3),c(0.1,0.4,0.2,0.3))

P4_5_states = rbind(c(0.3,0.3,0.2,0.2),c(0.2,0.4,0.2,0.2),c(0.2,0.3,0.2,0.3),c(0.3,0.2,0.1,0.4))

P5_5_states = rbind(c(0.3,0.4,0.2,0.1),c(0.1,0.4,0.4,0.1),c(0.3,0.1,0.2,0.4),c(0.5,0.1,0.1,0.3))

Transition_matrix_5_states = abind(P1_5_states,P2_5_states,P3_5_states,P4_5_states,P5_5_states,along = 3)

#Transition_matrix_list = list(P1,P2,P3)

Segmentation_matrix_5_states = rbind(c(0.999,0.0004,0.0002,0.0002,0.0002),c(0.0001,0.999,0.0001,0.0004,0.0004),c(0.002,0.003,0.003,0.99,0.002),c(0.002,0.001,0.005,0.002,0.99),c(0.001,0.003,0.004,0.002,0.99))
#########################################################################################################

# simulate hidden sequence
# n = length of the requested sequence
# lambda = segmentation matrix
hssim_FMC=function(n,Segmentation_Matrix)
{
  r=dim(Segmentation_Matrix)[1]
  q = dim(Segmentation_Matrix)[1]
  states=1:q
  s=vector("numeric",n)
  
  equil_seg_mat = c()
  
  #Equlibrium probabities for segmentation matirx
  equil_seg_mat = equil(Segmentation_Matrix)
  
  s[1]=sample(states,1,FALSE,equil_seg_mat)
  for (t in 2:n) {
    s[t]=sample(states,1,FALSE,prob=Segmentation_Matrix[s[t-1],])
  }
  s
}

Hidden_seq = hssim_FMC(40000,Segmentation_matrix)
Hidden_seq_4_states = hssim_FMC(40000,Segmentation_matrix_4_states)
#########################################################################
#Generation of observed sequence from the hidden state


simulated_obs_state_FMC = function(Hidden_sequence, Transition_matrices){
  
  n = length(Hidden_sequence)
  r = dim(Transition_matrices)[1]
  q = dim(Transition_matrices)[1]
  obs_states = 1:q
  
  obs_seq = vector("numeric",n)
  
  emission_prob = equil(Transition_matrices[,,Hidden_sequence[1]])
  obs_seq[1] = sample(obs_states,1,FALSE,emission_prob)
  #print(emission_prob)
  
  for(t in 2:n){
    #print(Transition_matrices[obs_seq[t-1],,Hidden_sequence[t]])
    obs_seq[t] = sample(obs_states,1,FALSE,prob = Transition_matrices[obs_seq[t-1],,Hidden_sequence[t]])
  }
  obs_seq
}


obs_seq = simulated_obs_state_FMC(Hidden_seq,Transition_matrix)
obs_seq_4_states = simulated_obs_state_FMC(Hidden_seq,Transition_matrix)
##################################################################################
######################################################################################
FORWARD_BACKWARD_ALGORITHM = function(TranstionMatrix, SegmantationMatrix, ObsSeq){
  
  r = dim(SegmantationMatrix)[1]
  n = length(ObsSeq)
  
  s = vector("numeric",n)
  filtered_prob = matrix(0,nrow=r,ncol=n)
  backward_prob = matrix(0,nrow=r,ncol=n)
  pred_prob = vector("numeric",r)
  b = matrix(0,nrow=r,ncol=n)
  m = matrix(0,nrow=r,ncol=n)
  
  
  m_new = matrix(0,nrow=r,ncol=n)
  s_new = vector("numeric",n)
  
  #Equlibrium probabilites for transition matrices
  b_new = array(0,c(r,r,n))
  b_new_2 = array(0,c(r,r,n))
  temp_equil_matrix = matrix(0,nrow = r, ncol = 4)
  for(i in 1:r){
    temp_equil_matrix[i,] = equil(TranstionMatrix[,,i])
  }
  
  
  
  #First filtered probability  
  filtered_prob[,1] = temp_equil_matrix[,ObsSeq[1]]*equil(SegmantationMatrix)
  filtered_prob[,1]=filtered_prob[,1]/sum(filtered_prob[,1])
  
  #calculation of filtered probabilites for n = 2: n
  for(i in 2:n){
    for (k in 1:r){
      pred_prob[k] = sum(SegmantationMatrix[,k] * filtered_prob[,i-1])
    }
    filtered_prob[,i] = TranstionMatrix[ObsSeq[i-1],ObsSeq[i],]*pred_prob
    filtered_prob[,i] = filtered_prob[,i]/sum(filtered_prob[,i])
    
  }
  ###############################################################################################
  ##backward probabilities
  #####################################
  # probabity that s(t+1) = k given that S(t) = j at y = t , where t belongs to 1:n
  
  for(t in 1:n){
    for(k in 1:r){
      for(j in 1:r){
        b_new_2[j,k,t] = SegmantationMatrix[j,k]*filtered_prob[j,t]
      }
    }
  }
  
  
  
  #########################################################################################################
  #Second normalisation method
  b_new_3 = array(0,c(r,r,n))
  for(t in 1:n){
    for(k in 1:r){
      temp_row_sum = 0
      # for(j in 1:r){
      #temp_row_sum = sum(b_new_2[,k,t])
      
      # }
      b_new_3[,k,t] = b_new_2[,k,t]/sum(b_new_2[,k,t])
    }
  }
  ###################################################################################################################
  
  #######################################################################################
  ##### #####             estimated hidden seq    ############
  #############################################################################################
  m_new = matrix(0,r,n)
  m_new[,n] = filtered_prob[,n]
  for(t in (n-1):1){
    for(j in 1:r){
      temp_b = 0
      temp_b = sum((b_new_3[j,,t])*m_new[,(t+1)])
      m_new[j,t] = temp_b
    }
  }
  estmated_hidden_seq = c()
  for(i in 1:n){
    estmated_hidden_seq[i] = which.max(m_new[,i])
  }
  plot(estmated_hidden_seq)
  
  return(estmated_hidden_seq)
}

########################################################################################
########################################################################################

#############################################################################################
#############################################################################





#MAIN CODE
# 3 states
Hidden_seq = hssim_FMC(40000,Segmentation_matrix)
obs_seq = simulated_obs_state_FMC(Hidden_seq,Transition_matrix)

Estimated_hidden_seq_3_states = FORWARD_BACKWARD_ALGORITHM_TEST(Transition_matrix, Segmentation_matrix, obs_seq)
Sequence_plot(Hidden_seq, Estimated_hidden_seq_3_states)

#LENGTH = 60k


Hidden_seq_60k = hssim_FMC(60000,Segmentation_matrix)
obs_seq_60k = simulated_obs_state_FMC(Hidden_seq_60k,Transition_matrix)

Estimated_hidden_seq_3_states_60k = FORWARD_BACKWARD_ALGORITHM(Transition_matrix, Segmentation_matrix, obs_seq_60k)
Sequence_plot(Hidden_seq_60k, Estimated_hidden_seq_3_states_60k)


##4 states
Hidden_seq_4_states = hssim_FMC(40000,Segmentation_matrix_4_states)

obs_seq_4_states = simulated_obs_state_FMC(Hidden_seq_4_states,Transition_matrix_4_states)

HS_4_state = FORWARD_BACKWARD_ALGORITHM(Transition_matrix_4_states,Segmentation_matrix_4_states, obs_seq_4_states)

Sequence_plot(Hidden_seq_4_states, HS_4_state)

#lenght = 60k

Hidden_seq_4_states_60k = hssim_FMC(60000,Segmentation_matrix_4_states)

obs_seq_4_states_60k = simulated_obs_state_FMC(Hidden_seq_4_states_60k,Transition_matrix_4_states)

HS_4_state_60k = FORWARD_BACKWARD_ALGORITHM(Transition_matrix_4_states,Segmentation_matrix_4_states, obs_seq_4_states_60k)

Sequence_plot(Hidden_seq_4_states_60k, HS_4_state_60k)


##5 states
Hidden_seq_5_states = hssim_FMC(40000,Segmentation_matrix_5_states)

obs_seq_5_states = simulated_obs_state_FMC(Hidden_seq_5_states,Transition_matrix_5_states)

Estimated_HS_5_state = FORWARD_BACKWARD_ALGORITHM(Transition_matrix_5_states,Segmentation_matrix_5_states, obs_seq_5_states)

Sequence_plot(Hidden_seq_5_states, Estimated_HS_5_state)

#LENGTH = 60K 
Hidden_seq_5_states_60k = hssim_FMC(60000,Segmentation_matrix_5_states)

obs_seq_5_states_60k = simulated_obs_state_FMC(Hidden_seq_5_states_60k,Transition_matrix_5_states)

Estimated_HS_5_state_60k = FORWARD_BACKWARD_ALGORITHM(Transition_matrix_5_states,Segmentation_matrix_5_states, obs_seq_5_states_60k)

Sequence_plot(Hidden_seq_5_states_60k, Estimated_HS_5_state_60k)
