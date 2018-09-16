import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


datosRK1 = np.genfromtxt("datosRK1.dat")
datosRK2 = np.genfromtxt("datosRK2.dat")
datosRK4 = np.genfromtxt("datosRK4.dat")
datosRK3 = np.genfromtxt("datosRK3.dat")
datosRK5 = np.genfromtxt("datosRK5.dat")
datosRK6 = np.genfromtxt("datosRK6.dat")

plt.figure(figsize=(10,7))
plt.plot(datosRK1[:,1],datosRK1[:,2],label="Starts at (0.8,0.0)")
plt.plot(datosRK2[:,1],datosRK2[:,2],label="Starts at (2.0,0.0)")
#plt.plot(datosRK3[:,1],datosRK3[:,2],label="Starts at (-2.0,0.0)")
plt.plot(datosRK4[:,1],datosRK4[:,2],label="Starts at (-1.0,0.6)")
plt.plot(datosRK5[:,1],datosRK5[:,2],label="Starts at (-2.0,-2.0)",color='k')
plt.plot(datosRK6[:,1],datosRK6[:,2],label="Starts at (0.05,0.0)",color='r')
plt.legend()
plt.savefig("graph1.pdf")
plt.clf()
plt.close()
