import matplotlib.pyplot as plt
import numpy as np

datos1 = np.genfromtxt("datos1.dat") 
S=130
M = 10
plt.figure()
for i in range(M):
	plt.scatter(datos1[i,1], datos1[i,2],s=S)
	plt.arrow(datos1[i,1], datos1[i,2],datos1[i+1,1]-datos1[i,1], datos1[i+1,2]-datos1[i,2],head_width=0.03,ec='k',fc='k', linestyle='dotted')

i=M
plt.scatter(datos1[i,1], datos1[i,2],s=S)
#plt.legend()
plt.savefig("henon-graph.pdf")

