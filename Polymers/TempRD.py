# Ryan Jukes
# 10/20/25
# Week 1 - Temperature based radial distribution

import gsd.hoomd
import freud
import numpy as np
import matplotlib.pyplot as plt

temp = 0.0
plt.figure(figsize=(10,6))                         #make graph and define on screen size

while temp <= 1.0 + 0.0001:                        #1.0 is max value. 0.0001 to counter possible floating point rounding issues
    inFile = f"trajectoryT{temp:.1f}.gsd"          #gain input. {temp:.1f} to get 1 decimal point when looking for correct file
    print("Reading",inFile)

    f = gsd.hoomd.open(name = inFile, mode = 'r')  #Makes f a hoomd object containing the file data. r=read mode
    frame = f[0]                                   #Takes first frame from the data to analyze and use graphically
    fp = frame.particles.position                  #fp object stores positional data in NumPy array (# particles, 3). x, y, z rows (3) 
    sideLength = float(frame.configuration.box[0]) #x dimension
    box = freud.box.Box.cube(sideLength)           #makes the freud box object
    print(sideLength)

    rdf = freud.density.RDF(bins = 75, r_max = sideLength/2-1e-6)    #Makes RDF obj. Small number ensures it is smaller than half the box
    rdf.compute(system=(box, fp))                                    #Computes values of RDF 

    plt.plot(rdf.bin_centers, rdf.rdf, label=f"T={temp:.1f}")        #Graphs computed values (x, y, label for legend)

    temp += 0.1
    f.close()


#Graphical setup with matplotlib

plt.title("RDF for temp")
plt.xlabel("r")
plt.ylabel("g(r)")
plt.legend()
plt.grid(True)
plt.show()