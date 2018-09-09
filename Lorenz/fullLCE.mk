all : fullLCE.exe
	./fullLCE.exe

fullLCE.exe : fullLCE.cpp
	g++ -O3 fullLCE.cpp -o fullLCE.exe
