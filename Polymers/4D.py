#Ryan Jukes
#4D for week 5 LJ sim 
#11-10-25

# Calculate the end to end time correlation function <r(0)r(t)>. 
# At time equal to 0, <r(0) dot r(t)> becomes <r^2(0)> and should be the same as the end to end distance. 
# As time passes, the vectors become less and less correlated and the function should decay to zero.

import gsd.hoomd
import numpy as np
import matplotlib.pyplot as plt

nChains = 32**2
nBeadsPerChain = 32

inputFile = gsd.hoomd.open('unwrappedTrajectory.gsd', 'r')
frames = len(inputFile)

pltvals = []          #2D array for my plot that stores the time value and the final time correlation function value (x,y)
r0vectors = []
r2vals = []
r0rtdotprods = []       #r(0)r(t) dot products array (done for each vector in the loop
count = 1

for frame in inputFile:
    print (f"Frame #{count}")

    for chainNumb in range(nChains):
        start = chainNumb * nBeadsPerChain
        end = (chainNumb + 1) * nBeadsPerChain
        chainPosition = frame.particles.position[start:end]        # chainPosition - numpy array storing all particle vals in the vector

        if count == 1:
            r0vectors.append(chainPosition[-1] - chainPosition[0])

        r0vector = r0vectors[chainNumb]
        rt = chainPosition[-1] - chainPosition[0]
        r0rtdotprod = np.dot(r0vector,rt)
        r0rtdotprods.append(r0rtdotprod)

        r2dotprod = np.dot(r0vector,r0vector)
        r2vals.append(r2dotprod)


    avgr2 = np.mean(r2vals)
    avgr0rtdotprod = np.mean(r0rtdotprods)

    Qrt = (avgr0rtdotprod/avgr2)
    pltvals.append((count-1,Qrt))
    #pltvals.append((count-1,avgr0rtdotprod))    #Def from textbook or README file
    count += 1

pltvals = np.array(pltvals)          #If its not a numpy array it throws "list indices must be integers or slices, not tuple"
plt.figure(figsize=(10, 6))
plt.plot(pltvals[:,0],pltvals[:,1])
plt.xlabel('Time')
plt.ylabel('Qrt')
plt.show()
