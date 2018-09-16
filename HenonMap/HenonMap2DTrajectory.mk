all : henon-graph.pdf

henon-graph.pdf : plot.py datos1.dat datos2.dat
	python3 plot.py

datos1.dat datos2.dat : henon.exe
	./henon.exe
	
henon.exe : HenonMap2DTrajectory.cpp
	g++ -O3 HenonMap2DTrajectory.cpp -o henon.exe
