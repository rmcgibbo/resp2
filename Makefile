all:
	g++ -Wall -c -fPIC -I. dotsphere.cc
	g++ -Wall -c -fPIC -I. vdwsurface.cc
	g++ -shared -o dotsphere.so dotsphere.o vdwsurface.o
