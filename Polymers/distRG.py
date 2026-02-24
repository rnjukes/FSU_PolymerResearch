# Ryan Jukes
# Distribution Func for Radius of Gyration
# 12-3-25

import gsd.hoomd
import numpy as np
import matplotlib.pyplot as plt 

inputFile = gsd.hoomd.open('Ltrajectoryr0-2.0.gsd', 'rb')
confirmlf = inputFile[-1]
frames = len(inputFile)
positions = confirmlf.particles.position      #Get (x,y,z) position array from gsd 
degreePol = 192   #Aka (N)
nChains = 800
vals = []


for chainNumb in range(nChains):
    # First calculate the chain center of mass via R_COM=1/N \Sum r_i 
    # where N is the degree of polymerization (192) 
    # and r_i is the position of the ith bead in the chain

    #Rg_z=1/N \sum (r_i,z-R_COM,z)^2

    start = chainNumb * degreePol
    end = start + degreePol 
    r = positions[start:end]

    #Loop and add each bead vector value to the running vector, then divide by degree of poly to average position vector (com)
    cmvector = np.array([0.0, 0.0, 0.0])   
    for newBead in r: 
        cmvector += newBead
    centOfMass = cmvector/degreePol   

    zcm = centOfMass[2]  #gets only the z val from the array
    count = 0.0 
    for newBead in r: 
        z_i = newBead[2]
        count += (z_i - zcm)**2      #(bead postion from COM of z)^2
    gyration = count/degreePol       #Accounts for the 1/N before the sigma notation 

    vals.append(gyration)

vals = np.array(vals)

#Plot em
plt.hist(vals, bins=100)
plt.xlabel("Radial Gyration Z")
plt.ylabel("Count")
plt.show()
