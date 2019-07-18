CC = g++
CFLAGS = -Wall -O3 --std=c++11

Ex02.1.exe : Ex02.1.o random.o
	$(CC) random.o Ex02.1.o -o Ex02.1.exe
Ex02.1.o : Ex02.1.cpp
	$(CC) -c Ex02.1.cpp -o Ex02.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Ex02.1.exe seed.out
