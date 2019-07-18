CC = g++
CFLAGS = -Wall -O3 --std=c++11

Metropolis.exe : Metropolis.o random.o
	$(CC) random.o Metropolis.o -o Metropolis.exe
Metropolis.o : Metropolis.cpp
	$(CC) -c Metropolis.cpp -o Metropolis.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Metropolis.exe seed.out
