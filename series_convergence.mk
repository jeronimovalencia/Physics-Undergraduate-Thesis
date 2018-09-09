all : convergence.pdf

convergence.pdf : datos1.dat datos2.dat plot_convergence.py
	python3 plot_convergence.py
	
datos1.dat datos2.dat : convergence.exe
	./convergence.exe
	rm convergence.exe

convergence.exe : series_convergence.cpp
	g++ series_convergence.cpp -o convergence.exe
