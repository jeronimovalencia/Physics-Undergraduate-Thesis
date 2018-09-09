all : graficas.pdf graficas_3D.pdf

graficas.pdf graficas_3D.pdf : datosRK.dat plot.py
	python3 plot.py
	

datosRK.dat : datosRK.exe
	./datosRK.exe 

datosRK.exe : trajectory.cpp
	g++ -O3 trajectory.cpp -o datosRK.exe
