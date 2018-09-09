all : LorenzLCE.pdf

LorenzLCE.pdf : datos3.dat plotLCE.py
	python3 plotLCE.py

datos3.dat : LCE.exe
	./LCE.exe
	rm LCE.exe

LCE.exe : LCE.cpp
	g++ LCE.cpp -o LCE.exe
