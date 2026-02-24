import gsd.hoomd
import numpy as np

traj = gsd.hoomd.open("Ltrajectoryr0-4.0.gsd", "r")
frame = traj[-1] #Last frame of the pulled back sim
position = frame.particles.position
bonds = frame.bonds.group
box = frame.configuration.box

x,y,z = box[:3] #Getting box side lengths
lengths = [] #Bond lengths

for i,j in bonds: #Loops over each pair found (i to j)
    dx, dy, dz = position[j] - position[i] #Find and log the displacement between the bonded beads in the chain
    dist = np.sqrt(dx**2 + dy**2 + dz**2)
    lengths.append(dist)

lengths = np.array(lengths)

maxVal=1.55
minVal=1.45
above = np.sum(lengths > maxVal)
below = np.sum(lengths < minVal)
total = len(lengths)

pctAbove = 100*above/total
pctBelow = 100*below/total

print("Mean bond length: ",np.mean(lengths))
print("Std dev:          ",np.std(lengths))
print("Min distance:     ",np.min(lengths))
print("Max distance:     ",np.max(lengths))
print("Percent above",maxVal,"= ",pctAbove)
print("Percent below",minVal,"= ",pctBelow)