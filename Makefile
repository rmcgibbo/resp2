all:
	g++ -c -fPIC -I. dotsphere.cc
	g++ -shared -o dotsphere.so dotsphere.o
