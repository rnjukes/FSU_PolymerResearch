import hoomd
import itertools
import numpy as np
import math
import gsd.hoomd
import sys
import time

#Set the degree of polymerization

DoP=192

degreeOfPolymerization=DoP
Nchains=800

#Choose either GPU or CPU. If no GPU is present will produce an error

deviceChoice=hoomd.device.GPU()
#deviceChoice=hoomd.device.CPU()

simulation=hoomd.Simulation(device=deviceChoice, seed=1)
simulation.timestep = 0
simulation.create_state_from_gsd(filename="Ltrajectoryr0-4.0.gsd",frame=-1)

integrator=hoomd.md.Integrator(dt=0.005)
cell=hoomd.md.nlist.Cell(buffer=0.4)

#Define forces
#Note! The cutoff here is different than before!
lj=hoomd.md.pair.LJ(nlist=cell,mode="shift")
lj.params[("A","A")]=dict(epsilon=3.000,sigma=1.0)
lj.r_cut[("A","A")]=2.5
lj.params[("B","B")]=dict(epsilon=1.000,sigma=1.0)
lj.r_cut[("B","B")]=2**(1/6)
lj.params[("A","B")]=dict(epsilon=4.500,sigma=1.0)
lj.r_cut[("A","B")]=2**(1/6)


#FENE bonds use a special potential which diverges at r0 and prevents the bond from stretching
#However FENE bonds can be unstable

#Adjust r0 here
r0Base=1.5
r0New=4.0

#bonded=hoomd.md.bond.FENEWCA()
#bonded.params["polymer"] = dict(k=30.0/(r0New/r0Base)**2, r0=r0New, epsilon=1.0, sigma=1.0, delta=0)

#integrator.forces.append(lj)
#integrator.forces.append(bonded)
langevin=hoomd.md.methods.Langevin(filter=hoomd.filter.All(),kT=1.0,default_gamma=1.0)
#npt=hoomd.md.methods.ConstantPressure(filter=hoomd.filter.All(), S=5.0,
                                   #   tauS=1.0,
                                   #   couple="xyz",thermostat=hoomd.md.methods.thermostats.Bussi(kT=1.0))
integrator.methods.append(langevin)

simulation.operations.integrator=integrator

#Custom update to scramble some bonds

#simulation.operations.updaters.append(box_resize)

#Write trajectory file
configSteps=10000
config=hoomd.write.GSD(hoomd.trigger.Periodic(configSteps),"trajectory.gsd",mode='wb',dynamic=['property','topology'])
simulation.operations.writers.append(config)

#Thermodynamic data
loggerSim=hoomd.logging.Logger(categories=['scalar','sequence'])
thermodynamic_properties=hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
logSteps=1000
loggerSim.add(thermodynamic_properties,quantities=['kinetic_temperature','pressure','kinetic_energy','potential_energy','pressure_tensor'])
performanceFile=open('thermodynamicData','w')
table=hoomd.write.HDF5Log(trigger=hoomd.trigger.Periodic(logSteps),filename='thermoData.h5',mode='w',logger=loggerSim)
simulation.operations.writers.append(table)

#Simulation progress
statusSteps=10000
simStatus=hoomd.logging.Logger(categories=['scalar'])
thermodynamic_properties=hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
simulation.operations.computes.append(thermodynamic_properties)
simStatus.add(simulation,
              quantities=['timestep','tps','walltime'])
simStatus.add(thermodynamic_properties,
              quantities=['kinetic_temperature','kinetic_energy','potential_energy'])
statusOut=hoomd.write.Table(trigger=hoomd.trigger.Periodic(statusSteps),logger=simStatus)
simulation.operations.writers.append(statusOut)

#Run simulation

simulation.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=1.0)
print('Precheck finished commencing simulation')

#Starting with an already pulled file so dont use or start the FENE bonds?
#simulation.run(3e6)
#simulation.run(5e6)

bondedHarmonic=hoomd.md.bond.Harmonic()
bondedHarmonic.params["polymer"] = dict(k=30.0, r0=0)

steps = int(1e6)
updateFreq = 10000
totalUpdates = steps//updateFreq  #int division

kBase = 1000.0
kStart = 1.0 * kBase
kFinal = 5.0 * kBase
#Start at 20% of the max?

rStart=4.0
rFinal=1.5

#New integrator vals (slow down pulling back together)
dtStart = 1e-4
dtEnd = 0.005

integrator.forces.append(bondedHarmonic)

for i in range(0,totalUpdates):
    scale = (i+1)/totalUpdates

    #Fraction of the loop completed

    #integrator.dt = dtStart + scale * (dtEnd-dtStart)
    #Main increasing pull function that will scale via how much of the process 
    #has been completed already

    kStep = kStart + scale * (kFinal-kStart)
    rStep = rStart + scale * (rFinal-rStart)
    #k is bond pulling strength in a way - attractive force strength  
    #use same formula as integrator but for k? 

    bondedHarmonic.params["polymer"] = dict(
    #Make this edit harmonic bonds and not the FENE from the main sim
        k=kStep,
        r0=rStep,
        )
        #Harmonics only take k and r0 values as inputs here

    simulation.run(updateFreq)

    if i % 20 == 0:
        print("Ramp", i+1, "/", totalUpdates, "dt=", integrator.dt, "k=", kStep)

hoomd.write.GSD.write(state=simulation.state,filename="finalConfig.gsd",mode='wb')

