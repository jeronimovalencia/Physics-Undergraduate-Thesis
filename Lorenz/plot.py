import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D


datosRK = np.genfromtxt("datosRK.dat")

plt.figure(figsize=(10,10))
plt.title("Lorentz System")
plt.subplot(221)
plt.title("X vs. t")
plt.plot(datosRK[:,0],datosRK[:,1])
plt.subplot(222)
plt.title("Y vs. t")
plt.plot(datosRK[:,0],datosRK[:,2])
plt.subplot(223)
plt.title("Z vs. t")
plt.plot(datosRK[:,0],datosRK[:,3])
plt.subplot(224)
plt.title("Y vs. X")
plt.plot(datosRK[:,1],datosRK[:,2])
plt.savefig("graficas.pdf")
plt.clf()
plt.close()

fig2 = plt.figure()
ax = fig2.gca(projection = '3d')
ax.plot(datosRK[:,1],datosRK[:,2],datosRK[:,3],linewidth=0.5, markersize=12)
plt.savefig("graficas_3D.pdf")
plt.clf()
plt.close()
