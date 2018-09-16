all : graph1.pdf

graph1.pdf : datosRK1.dat datosRK2.dat datosRK3.dat datosRK4.dat plot.py
	python3 plot.py

datosRK1.dat datosRK2.dat datosRK3.dat datosRK4.dat : datosRK.exe
	./datosRK.exe 

datosRK.exe : trajectory.cpp
	g++ -O3 trajectory.cpp -o datosRK.exe
