#!/usr/bin/python3

# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 20:38:09 2018

@author: rafael
"""

import numpy as np
import functielijst as f
import sys
#import matplotlib.pyplot as plt

np.seterr(all='raise')

N = 100 #20
T = 2 #24
v_i_num = 10
v_i = np.ones((N), dtype=int) * v_i_num
v_i_j = np.ones((N, T)) * v_i_num / T         #The matrix holding all our cells, 
                                              #with the genes in the X direction and the different cells in the Y                                                                  
v_max = 20 #24
decay = 0.001
period_output=10*N

replication_count=0
#build a_j nested array                                                        
a_j = np.zeros(shape=(N, T))                  #We need this array to determine the probability a different gene is chosen, and to change replicase affinities to for example 0.9.
for n in range(N):                                                              
    for j in range(T):
        if j == 0:                                               
            a_j[n][j] = 0.8#1.0
        else:
            a_j[n][j] = 1.0
            
R_icell = np.ones(N)                                                                   
d = 0 
count = v_i_num  

for n in range(N):
        R_icell[n]  = f.R_i(v_max, T, v_i_j[n], v_i[n])
                                                           
while d <= 500000000:                #S: not enough time
    
    #constantly calculate the RI   S: <- no need, just recalculate those that you modified
    for n in range(N):
        R_icell[n]  = f.R_i(v_max, T, v_i_j[n], v_i[n])
    #First, calculate R_i values for N cells
    '''                              
    if sum of prob == 0 : 
      print simulation exitnt 
      bye
      exit
    '''
    
    sum_R_i = np.sum(R_icell)
    if np.sum(sum_R_i)<0.000001: 
        sys.stderr.write("Extinction, simulation terminated\n")
        sys.exit(1)
    prob1 = R_icell / sum_R_i   #list with probability of being chosen for each cell
    #Then, randomly pick a cell and return R_i value and cell index
    ind_Cell = np.random.choice(range(N), size=None, replace=True, p=prob1) #S: what happens if all zero?
    
    joint_prob2 = np.multiply(a_j[ind_Cell], v_i_j[ind_Cell])  #multiplies the a_j value for each type of replicator with the number of replicators
    # print (joint_prob2)
    # sys.exit(1)
    prob2 = joint_prob2 / np.sum(joint_prob2)
    ind_Chosen_Gene = np.random.choice(range(T), size=None, replace=True, p=prob2) #chose a gene
    #Increase chosen cell ribozome count by 1
    v_i[ind_Cell] = v_i[ind_Cell] + 1                                                
    
    v_i_j[ind_Cell][ind_Chosen_Gene] = v_i_j[ind_Cell][ind_Chosen_Gene] + 1
    
    #Now, make if statement to determine if cells die/duplicate
    if v_i[ind_Cell] >= v_max:
        #select random cell to destroy
        rand_Cell = np.random.choice(range(N))
        if rand_Cell == ind_Cell:
            while rand_Cell == ind_Cell:
                rand_Cell = np.random.choice(range(N))
        else:
            rand_Cell = rand_Cell # S: no need
        
        
        new_Cell1 = np.ones(T)
        
        
        for t in range(T):
            new_Cell1[t] = np.random.binomial(v_i_j[ind_Cell][t], 0.5)
            
        new_Cell2 = v_i_j[ind_Cell] - new_Cell1
        
        v_i_j[rand_Cell] = np.copy(new_Cell1)  # S: is assigning array like this just makes reference: VERY RISKY !!!
        v_i_j[ind_Cell] = new_Cell2[:]
        
        v_i[rand_Cell] = np.sum(new_Cell1)
        v_i[ind_Cell] = np.sum(new_Cell2)
        
        replication_count +=1
        
    for ind_Cell in range(N):
      for ind_Chosen_Gene in range(T):
        howmanydecay=np.random.binomial(v_i_j[ind_Cell][ind_Chosen_Gene], decay)
        v_i[ind_Cell] = v_i[ind_Cell] - howmanydecay
        v_i_j[ind_Cell][ind_Chosen_Gene] = v_i_j[ind_Cell][ind_Chosen_Gene] - howmanydecay
        # print("decay in cell,gene, new number: ",ind_Cell,ind_Chosen_Gene, v_i_j[ind_Cell][ind_Chosen_Gene] )
    
    if d % period_output == 0: 
          viable_count = list(R_icell).count(0)
          viable_frac = (N - np.array(viable_count)) / N
          #message = 'Viable ribocell fraction '
          #print(message),
          try:
            mean_health = np.mean([ x[0]/float(sum(x)) if sum(x)!=0 else 0. for x in v_i_j  ])
            # std_health = np.std([ x[0]/float(sum(x)) if sum(x)!=0 else 0. for x in v_i_j ])
          except:
            # print( [ x[0]/float(sum(x)) for x in v_i_j ])
            print( [ float(sum(x)) for x in v_i_j ])
            sys.exit(1)
            
          print(d/N, viable_frac , N*replication_count/period_output , [(str(x[0])+','+str(x[1])) for x in v_i_j])
          sys.stdout.flush()
          # print( v_i_j)
          # print("repl per time step: ", replication_count/period_output)
          replication_count=0
    d = d + 1 #S: What is this doing here? Might be actually correct????
          


