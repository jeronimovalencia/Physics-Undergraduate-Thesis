all : LorenzLCE.pdf

LorenzLCE.pdf : datos3.dat plotLCE.py
	python3 plotLCE.py

datos3.dat : LCE.exe
	./LCE.exe

LCE.exe : LCE.cpp
	g++ -O3 LCE.cpp -o LCE.exe
