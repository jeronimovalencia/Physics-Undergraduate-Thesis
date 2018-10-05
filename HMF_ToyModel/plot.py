import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


datosRK = np.genfromtxt("datosRK.dat")

x_0 = np.cos(datosRK[0,:])
y_0 = np.sin(datosRK[0,:])

x_f = np.cos(datosRK[1,:])
y_f = np.sin(datosRK[1,:])

theta = np.linspace(0,2*np.pi,100)

plt.figure(figsize=(16,7))
plt.subplot(121)
plt.title("Initial distribution",size=20)
plt.scatter(x_0,y_0, label="Initial", color="red",s=10)
plt.plot(np.cos(theta), np.sin(theta),color="gray",linewidth=1)
plt.xlabel('x',size=22)
plt.ylabel('y',size=22)
plt.subplot(122)
plt.title("Final distribution",size=20)
plt.scatter(x_f,y_f, label="Final", color="black",s=10)
plt.plot(np.cos(theta), np.sin(theta),color="gray",linewidth=1)
plt.xlabel('x',size=22)
plt.ylabel('y',size=22)
plt.savefig("graficas.pdf")
plt.clf()
plt.close()

