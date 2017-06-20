all : main.cpp mtwist.c
	g++ main.cpp mtwist.c -o test -O3 -g
