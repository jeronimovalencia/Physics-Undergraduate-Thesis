all : graficaLyap.pdf

graficaLyap.pdf : datosLyap.dat plotLyapMagnetization.py plotLyapExponents.py
	python3 plotLyapMagnetization.py
	python3 plotLyapExponents.py

datosLyap.dat : Lyap.exe
	./Lyap.exe

Lyap.exe : Lyapunov2.cpp
	g++ -O3 Lyapunov2.cpp -o Lyap.exe
