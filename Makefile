all:
	g++ -c -fPIC -I. dotsphere.cc
	g++ -c -fPIC -I. vdwsurface.cc
	g++ -shared -o dotsphere.so dotsphere.o vdwsurface.o
