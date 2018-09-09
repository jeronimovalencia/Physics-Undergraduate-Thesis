import numpy as np
import matplotlib.pyplot as plt

data1 = np.genfromtxt("datos1.dat")
data2 = np.genfromtxt("datos2.dat")


plt.figure(figsize=(10,10))
plt.plot(data1[:,0],data1[:,1])
plt.plot(data2[:,0],data2[:,1])
plt.title("Lyapunov Exponent vs. Number of terms in the sum")
plt.savefig("convergence.pdf")
plt.ylim(0,2)
plt.clf()
plt.close()
