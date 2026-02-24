#Ryan Jukes
#4C for week 6 LJ sim 
#11-10-25

# The mean value of <r^2> as a function of beads along the chain, ie 2 beads, 3 beads, 4 beads, ... , 31 beads, 32 beads.
# If done correctly, the data should reach a plateau according to <r^2>=Nb^2 where b is called the statistical segment length.

import gsd.hoomd
import numpy as np
import matplotlib.pyplot as plt 

nChains = 32**2
nBeadsPerChain = 32

inputFile = gsd.hoomd.open('unwrappedTrajectory.gsd', 'r')
frames = len(inputFile)

r2vals = []
r2tracker = []
pltvals = []
count = 1

for i in range(nBeadsPerChain-1):
    r2vals.append(0.0)
    r2tracker.append(0)


for frame in inputFile:
    print (f"Frame #{count}")                     

    for chainNumb in range(nChains):       #Loop though specific chains in the frame           
        start = chainNumb * nBeadsPerChain        
        end = (chainNumb + 1) * nBeadsPerChain
        chainPosition = frame.particles.position[start:end] 

        for x in range(nBeadsPerChain):        #Loops over reference bead position

            for y in range(1, nBeadsPerChain-x):              #y is the distance in beads between beads
                r = chainPosition[x+y] - chainPosition[x]     
                r2 = np.dot(r,r)                              
                r2vals[y-1] += r2         #Adds calculated r^2 value to list with proper perameters 
                r2tracker[y-1] += 1       #Tracks amt of calcs done for y value with +1

    count += 1 

xvals = []
yvals = []

for y in range(1, nBeadsPerChain):
    xvals.append(y)
    f = (r2vals[y-1]/r2tracker[y-1])
    yvals.append(f)

xvals = np.array(xvals)
yvals = np.array(yvals)

plt.figure(figsize=(10, 6))
plt.plot(xvals,yvals)
plt.xlabel('Bead Distance')
plt.ylabel('r^2')
#plt.show()
