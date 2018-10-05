all : graficas.pdf

graficas.pdf : datosRK.dat plot.py
	python3 plot.py
	
datosRK.dat : datosRK.exe
	./datosRK.exe 

datosRK.exe : trajectory2.cpp
	g++ -O3 trajectory2.cpp -o datosRK.exe
