CC = g++
CFLAGS = -Wall -O3 --std=c++11

Ex01.1.exe : Ex01.1.o random.o
	$(CC) random.o Ex01.1.o -o Ex01.1.exe
Ex01.1.o : Ex01.1.cpp
	$(CC) -c Ex01.1.cpp -o Ex01.1.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Ex01.1.exe seed.out
