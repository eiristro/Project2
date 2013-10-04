import numpy as np
from matplotlib.pylab import *

data = []
datar = []

filename = "zdata500.dat"
filenamer = "zdatar100.dat"

ofile = open(filename, "r")
for line in ofile:
    linedata = line.split("/n")
    for elem in linedata:
        data.append(float(elem))
ofile.close()

ofile = open(filenamer, "r")
for line in ofile:
    linedata = line.split("/n")
    for elem in linedata:
        datar.append(float(elem))
ofile.close()

rho = np.linspace(0, 5, len(data))

plot(rho, [data[i]**2 for i in range(len(data))])
title("With repulsion")
#figure()
#plot(rho, [datar[i]**2 for i in range(len(data))])
#title("With repulsion")
show()

#raw_input("enter to close")
