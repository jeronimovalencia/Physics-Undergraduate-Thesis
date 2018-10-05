import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import scipy.special as scisp
import scipy.optimize as sciop


datosLyap = np.genfromtxt("datosLyap.dat")

datosLyap2 = np.genfromtxt("datosLyap2.dat")

x = datosLyap[:,0]
y = datosLyap[:,1]

x1 = datosLyap2[:,0]
y1 = datosLyap2[:,1]


#Autoconsistencia
beta = np.linspace(1,300,3000)
roots = np.array([])

def f(x,lel):
	return (x/lel)-(scisp.iv(1,x)/scisp.iv(0,x))


for b in beta:
	root = sciop.newton(f,b,args=(b,))
	roots = np.append(roots,root)

U = 1/(2*beta) + 0.5*(1-(roots/(beta))*(roots/(beta)))

x2 = np.linspace(0.75,3.0,3000)

plt.figure(figsize=(7,7))
plt.scatter(x,y, color="black",label='Simulation')
plt.xlabel('U', size=22)
plt.ylabel('|M|', size=22)
plt.plot(U,roots/(beta),color='silver',label='Theoretical')
plt.plot(x2,0*x2,color='silver')
plt.legend(fontsize = 'large')
plt.savefig("graficaMagnetizacion.pdf")
plt.clf()
plt.close()
