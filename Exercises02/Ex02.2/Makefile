CC = g++
CFLAGS = -Wall -O3 --std=c++11

Ex02.2.exe : Ex02.2.o random.o 
	$(CC) random.o Ex02.2.o -o Ex02.2.exe
Ex02.2.o : Ex02.2.cpp
	$(CC) -c Ex02.2.cpp -o Ex02.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Ex02.2.exe seed.out
