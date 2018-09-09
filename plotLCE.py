import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("datos3.dat")

dataX = data[:,0]
dataY = data[:,1]

dataNeg = dataX[dataY <= 0]
dataYNeg = dataY[dataY <= 0]
dataPos = dataX[dataY > 0]
dataYPos = dataY[dataY > 0]

plt.figure(figsize=(8,5))
plt.title("Maximum Lyapunov Exponent (Lorenz System)")
plt.plot(dataX,dataY)
#plt.plot(dataNeg, dataYNeg)
#plt.plot(dataPos, dataYPos)
plt.xlabel(r"$\gamma$")
plt.ylabel(r"$\lambda_1$")
plt.savefig("LorenzLCE.pdf")
