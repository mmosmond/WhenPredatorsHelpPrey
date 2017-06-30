# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 16:10:44 2015

@author: mmosmond
"""

######################################################################
##MODULES##
######################################################################

import random as rng
import math
import numpy as np
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
    outfile_B = open("Z_k%d_r%d.dat" %(n+1,s+1),"w")

    return [outfile_A, outfile_B]

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
mp = 0.1 #min prey death rate (m_min)
e = 1 #predator conversion efficiency
c = 0.01 #max predator capture rate (k)
mz = 4 #predator death rate (m_P)

gammaP = 0.03 #selection on prey from death (gamma)
gammaZ = 0.5 #selection on prey and predator from predation (gamma_k)
sdmutP = 0.05 #standard deviation in prey segregation (alpha)
sdmutZ = 0.05 #standard deviation in predator segregation (alpha_P)

kappas = [0.0]  #rates of environmental change (delta)
NumSims = 1 #number of replicates

maxTime = 10 #total time of simulation
outputFreq = 1000 #how many iterations (birth and death events) between recording
recsteps = 10 #number of timesteps used to calculate simulation means 

######################################################################
##INITIAL VALUES##
######################################################################

P0 = int(mz / (c * e)) #initial prey abundance (N)
Z0 = int((mu * (Rtot - mz / (c * e)) - mp) / (mu + c)) #initial predator abundance (P)
R0 = Rtot - P0 - Z0 #intial resource abundance

zoptP0 = 0.1 #initial prey optimal (theta)
zP0 = zoptP0 #initial mean prey trait value (mean z)
zZ0 = 0 #initial mean predator trait (mean z_P)

expsdP = (2 * sdmutP ** 2)**0.5 #initial standard dev. in prey traits (V_z**0.5)
expsdZ = (2 * sdmutZ ** 2)**0.5 #initial standard dev. in predator traits (V_z,P**0.5)

rng.seed(1) #random number seed (for repeatability)

######################################################################
##SIMULATION##
######################################################################

def main():
    
    n = 0
    while n < len(kappas):    
        
        kappa = kappas[n]
        s = 0
        while s < NumSims:  
            
            time = 0
            iteration = 0
            abundances = [R0, P0, Z0]
            zP = np.random.normal(zP0, expsdP, abundances[1])   #initialize prey traits           
            zZ = np.random.normal(zZ0, expsdZ, abundances[2])   #initialize predator traits
        
            # open output files
            fileHandles = open_output_files(n,s)
            zPdat = open("zP_k%d_r%d.dat" %(n+1,s+1),"w")      
            zZdat = open("zZ_k%d_r%d.dat" %(n+1,s+1),"w")      
            
            # we simulate until we hit our maximum simulation time
            while time < maxTime:
                
                zoptP = zoptP0 + kappa * time #optimal prey trait                 
                
                if abundances[2] > 0: #if predators still alive
                    zZmean = np.mean(zZ) #record mean trait value
                else:
                    zZmean = 0

                #calculate propensities for reactions        
                a_bz = e * c * np.clip(1 - (gammaZ / 2) * (np.mean(zP) - zZ)**2, 0, 1) * abundances[1] #birth propensity for each predator                 
                a_dp = mp + (gammaP / 2) * (zP - zoptP)**2 + abundances[2] * c * np.clip(1 - (gammaZ / 2) * (zP - zZmean)**2, 0, 1) #background death propensity for each prey                
                a_dz = mz * abundances[2] #total predator death propensity
                a_bp = mu * abundances[0] * abundances[1] #total prey birth propensity
                a_i = [a_dz, a_bp]        
                a_0 = sum(a_bz) + sum(a_dp) + sum(a_i) #total propensity  
                
                # pick a random number, compute the time increment
                # and update t (equation 21a in Gillespie's paper)
                rand_1 = rng.random()
                tau    = 1.0/a_0 * math.log(1/rand_1)
                time  += tau
        
                # find the reaction to execute (equation 21b in
                # Gillespie's paper). We need a second random number
                # here.  
                rand_2    = rng.random()
                threshold = a_0 * rand_2
                
                summation = 0
                count     = 0
                while threshold > summation:
                    if count < len(a_i):            
                        summation += a_i[count]
                    elif count < len(a_i) + len(a_dp):
                        summation += a_dp[count - len(a_i)]
                    else:
#                        summation += a_pred[int((count - (len(a_i) + len(a_dp))) / len(a_pred[0]))][(count - (len(a_i) + len(a_dp))) % len(a_pred[0])]
                        summation += a_bz[count - (len(a_i) + len(a_dp))]
                    count += 1
        
                #predator death 
                if count - 1 == 0:
                    zZ = np.delete(zZ, rng.randrange(0, len(zZ))) #remove random predator from list
                    abundances[2] -= 1 #remove a predator  
                    abundances[0] += 1 #add available resource
                #prey birth
                elif count - 1 == 1:
                    new_zP = ( zP[rng.randrange(0, len(zP))] + zP[rng.randrange(0, len(zP))] ) / 2.0 #random mating (both parents random)
                    new_zP = np.random.normal(new_zP, sdmutP) #segregation and recombination
                    zP = np.append(zP, new_zP) #add offpsring to list of trait values
                    abundances[0] -= 1 #remove a resource
                    abundances[1] += 1 #add a prey
                #prey death
                elif count - 1 < len(a_i) + len(a_dp):
                    zP = np.delete(zP, count - 1 - len(a_i)) #remove dead prey from list
                    abundances[1] -= 1 #remove a prey
                    abundances[0] += 1 #add a resource
                #predator birth
                else:
                    new_zZ = ( zZ[count - 1 - (len(a_i) + len(a_dp))] + zZ[rng.randrange(0, len(zZ))] ) / 2.0 #random mating with focal
                    new_zZ = np.random.normal(new_zZ, sdmutZ) #segregation and recombination
                    zZ = np.append(zZ, new_zZ) #add offpsring to list of trait values
                    abundances[2] += 1 #add a predator 
                    abundances[0] -= 1 #remove a resource

                #end simulation (and record) if prey go extinct        
                if abundances[1] == 0: 
                    print("Prey extinct")
                    write_data_to_output(fileHandles, time, abundances[1:])
                    xP = (np.mean(zP), np.std(zP)**2)
                    zPdat.write(",".join(map(str, xP)) + "\n")                    
                    break

                #checks
                if sum(abundances) != Rtot:
                    print("Error: System not closed!")
                    return
                for i in range(len(abundances)):
                    if abundances[i] < 0:
                        print("Error: NEGATIVE ABUNDANCE!")                    
                        return
                
                # dump data every outputFreq iteration
                # we also print a short progess message 
                if (iteration % outputFreq) == 0:
                    write_data_to_output(fileHandles, time, abundances[1:])
                    xP = (np.mean(zP), np.std(zP)**2)
                    zPdat.write(",".join(map(str, xP)) + "\n")  
                    if abundances[2] > 0:
                        xZ = (np.mean(zZ), np.std(zZ)**2)                      
                        zZdat.write(",".join(map(str, xZ)) + "\n")
                    print("kappa %.4f simulation %d iteration %d time %5.4g" %(kappa, s+1, iteration, time))   
                
                iteration += 1    
    
            # cleanup
            close_output_files(fileHandles)
            zPdat.close()            
            zZdat.close()
            
            s += 1
             
        n += 1
    
######################################################################
##RUNNING##    
######################################################################
    
#run (with timer)
start = time.time()
main()
end = time.time()
print(end-start)   
