CC = g++
CFLAGS = -Wall -O3 --std=c++11

Ex03.1.exe : Ex03.1.o random.o
	$(CC) random.o Ex03.1.o -o Ex03.1.exe
Ex03.1.o : Ex03.1.cpp
	$(CC) -c Ex03.1.cpp -o Ex03.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Ex03.1.exe seed.out
