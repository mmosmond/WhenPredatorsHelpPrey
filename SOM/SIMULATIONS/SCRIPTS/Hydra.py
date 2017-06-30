# -*- coding: utf-8 -*-
"""
Created on Thu Aug 27 10:08:21 2015

@author: mmosmond
"""

######################################################################
##MODULES##
######################################################################s

import numpy as np
import random as rng
import math
import time

######################################################################
##HELPER FUNCTIONS##
######################################################################

def open_output_files(n,s):
    """
    This function opens the output files and returns file
    handles to each.
    """

    outfile_A = open("P_k%d_r%d.dat" %(n+1,s+1),"w")

    return [outfile_A]

def write_data_to_output(fileHandles, time, data):
    """
    This function writes a (time, data) pair to the
    corresponding output file. We write densities
    not abundances.
    """

    for i in range(0,len(fileHandles)):
        fileHandles[i].write("%5.4e  %8.7f\n" %(time, data[i]))
    
def close_output_files(fileHandles):
    """
    This function closes all output files.
    """

    for i in range(0,len(fileHandles)):
        fileHandles[i].close()

######################################################################
##PARAMETERS##
######################################################################

Rtot = 1000 #total resource abundance (R_tot)
mu = 0.01 #max prey birth rate (b_max)
mp = 0.1 #prey mortality rate (m)
c = 1 #per capita consumption rate by generalist predator (k)

gamma = 1.0 #selection on prey (gamma)
sdmut = 0.05 #standard deviation in prey segregation (alpha)

kappas = [0.000]  #rates of environmental change (kappa)
NumSims = 1 #number of replicates

maxTime = 10 #total time
outputFreq = 1000 #how many iterations (birth and death events) between recording
recsteps = 10 #number of timesteps used to calculate simulation means 

######################################################################
##INITIAL VALUES##
######################################################################

P0 = int(Rtot - (mp + c) / mu) #initial prey abundance (N)
R0 = Rtot - P0 #initial resource abundance

zoptP0 = 0 #intial abiotic optimum (theta)
zP0 = zoptP0 #intial mean prey trait value (mean z)

expsd = ( 2 * sdmut**2 )**0.5 #initial standard deviation in prey traits (V_z**0.5)

rng.seed(1) #random number seed

######################################################################
##SIMULATION##
######################################################################

def main():
    
    n = 0
    while n < len(kappas):
        
        kappa = kappas[n] #set rate of environmental change

        s = 0    
        while s < NumSims:    
        
            #initialize    
            time = 0
            iteration = 0
            abundances = [R0, P0] #abundances (R,P)
            zP = np.random.normal(zP0, expsd, abundances[1])   #phenotype distribution
            
            # open output files
            fileHandles = open_output_files(n,s)
            zPdat = open("zP_k%d_r%d.dat" %(n+1,s+1),"w")     
            
            # simulate until maximum simulation time
            while time < maxTime:
                
                zo = zoptP0 + kappa * time #optimal prey trait        
                
                #calculate propensities for reactions     
                abp = abundances[0] * mu * np.exp( - (gamma / 2) * (zP - zo)**2 ) #birth propensity for each prey (Gaussian)
                adp = (mp + c) * abundances[1] #total prey death propensity     
                a_0 = sum(abp) + adp #total propensity
               
                # pick a random number, compute the time increment
                # and update t (equation 21a in Gillespie's paper)
                rand_1 = rng.random()
                tau    = 1.0 / a_0 * math.log(1 / rand_1)
                time  += tau
                
                # find the reaction to execute (equation 21b in
                # Gillespie's paper). We need a second random number
                # here.  
                rand_2    = rng.random()
                threshold = a_0 * rand_2
                
                summation = 0
                count     = 0
                while threshold > summation:
                    if count < len(abp):            
                        summation += abp[count]
                    else:
                        summation += adp
                    count += 1
                
                #birth, death, and mutation
                if count - 1 < len(abp): #if birth event
                    new_zP = ( zP[count - 1] + zP[rng.randrange(0, len(zP))] ) / 2.0 #random mating with focal prey
                    new_zP = np.random.normal(new_zP, sdmut) #add random segregation effect to offspring trait value
                    zP = np.append(zP, new_zP) #add offpsring to list of trait values
                else: #if death event
                    zP = np.delete(zP, rng.randrange(0, len(zP))) #remove (random) dead individual's trait value from list
                
                #update abundances
                abundances[0] = Rtot - len(zP)
                abundances[1] = len(zP)
                
                #end simulation if prey go extinct        
                if abundances[1] == 0: 
                    print("Prey extinct")
                    write_data_to_output(fileHandles, time, abundances[1:])
                    xP = (np.mean(zP), np.std(zP)**2)
                    zPdat.write(",".join(map(str, xP)) + "\n")                    
                    break 
                
                #otherwise continue
                # dump data every outputFreq iteration
                # also print a short progess message 
                if (iteration % outputFreq) == 0:
                    write_data_to_output(fileHandles, time, abundances[1:])
                    xP = (np.mean(zP), np.std(zP)**2)
                    zPdat.write(",".join(map(str, xP)) + "\n")
                    print("kappa %.4f   simulation %d   iteration %d   time %5.4g" % (kappa, s+1, iteration, time))   
            
                iteration += 1    
                
            # cleanup
            close_output_files(fileHandles)
            zPdat.close()            
            
            #next replicate
            s += 1
            
        #next rate of change    
        n += 1
        
######################################################################
##RUNNING##
######################################################################    
    
#run (with timer)
start = time.time()
main()
end = time.time()
print(end-start)  
