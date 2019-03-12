import numpy as np
import matplotlib.pyplot as plt
import sys
import os

dir_path = os.path.dirname(os.path.realpath(__file__))
filename = sys.argv[1]

data = np.genfromtxt(dir_path + '/data/' + filename+'.dat', delimiter=' ')

Y = data[:,1:]
X = data[:,0]

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title(filename)
plt.ylabel(r"$s []$")
plt.xlabel(r"T $[K]$")

#plt.ylabel(r"$E_{pot} [\frac{eV}{Volume}]$")
#plt.xlabel(r"Unit cell volume $[\AA ^3]$")

ax1.plot(X,Y)

plt.savefig('plots/'+filename+'.png')

if "show" in sys.argv:
    plt.show()
