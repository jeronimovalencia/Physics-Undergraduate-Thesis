all : convergence.pdf

convergence.pdf : datos1.dat datos2.dat plot_convergence.py
	python3 plot_convergence.py
	
datos1.dat datos2.dat : convergence.exe
	./convergence.exe
	convergence.exe

convergence.exe : series_convergence.cpp
	g++ -O3 series_convergence.cpp -o convergence.exe
