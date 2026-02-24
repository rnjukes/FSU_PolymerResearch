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
simulation.create_state_from_gsd(filename="startConfigLow.gsd",frame=-1)

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
r0New=1.5

bonded=hoomd.md.bond.FENEWCA()
bonded.params["polymer"] = dict(k=30.0/(r0New/r0Base)**2, r0=r0New, epsilon=1.0, sigma=1.0, delta=0)

integrator.forces.append(lj)
integrator.forces.append(bonded)
langevin=hoomd.md.methods.Langevin(filter=hoomd.filter.All(),kT=1.0,default_gamma=0.02)
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

simulation.run(3e6)
#simulation.run(5e6)

hoomd.write.GSD.write(state=simulation.state,filename="finalConfig.gsd",mode='wb')