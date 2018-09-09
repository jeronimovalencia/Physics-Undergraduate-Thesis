all : graficas.pdf graficas_3D.pdf

graficas.pdf graficas_3D.pdf : datosRK.dat plot.py
	python3 plot.py
	

datosRK.dat : datosRK.exe
	./datosRK.exe 
	rm datosRK.exe

datosRK.exe : trajectory.cpp
	g++ trajectory.cpp -o datosRK.exe
