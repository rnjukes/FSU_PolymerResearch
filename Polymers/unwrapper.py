import gsd.hoomd
import numpy as np

#Open input and output trajectories
traj=gsd.hoomd.open('trajectory.gsd','r')
outTraj=gsd.hoomd.open('unwrappedTrajectory.gsd','w')

#Alter these for your system
nChains=800
nBeadsPerChain=192
nFrames=len(traj)

i=0
for frame in traj:

    i+=1
    #Make a new array of zeros based on the current positions
    newXYZ=np.zeros((nChains*nBeadsPerChain,3))
    beadIndex=0
    boxLengths=frame.configuration.box[0:3]

    for chainIndex in range(nChains):
        #The first position for each chain is always that of the first bead
        PS = chainIndex * nBeadsPerChain
        PF = (chainIndex + 1) * nBeadsPerChain
        positions = frame.particles.position[PS:PF]
        newXYZ[PS, :] = positions[0, :]
        r0 = positions[0,:]

        for beadIndex in range(1,nBeadsPerChain):
            #Calculate the nearest neighbor displacement

            currentBeadPosition=positions[beadIndex,:]
            deltaR=currentBeadPosition-r0

            #Calculate the nearest neighbor distance
            deltaR=deltaR-boxLengths*np.round(deltaR/boxLengths)

            #Calculate the true bead position
            r=deltaR+r0
            newXYZ[PS + beadIndex,:]=r

            #Set up the next r
            r0=r
            beadIndex+=1

    frame.particles.position=newXYZ
    outTraj.append(frame)