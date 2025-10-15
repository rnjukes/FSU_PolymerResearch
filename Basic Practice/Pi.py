#Ryan Jukes

valpi=0
runpi=0

for i in range(10000):
 valpi += ((-1)**i) / (2*i + 1)
 runpi+=valpi

print(valpi*4)
