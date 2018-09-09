all : fullLCE.exe
	./fullLCE.exe
	rm fullLCE.exe

fullLCE.exe : fullLCE.cpp
	g++ fullLCE.cpp -o fullLCE.exe
