CC = g++
CFLAGS = -Wall -O3 --std=c++11

Ex01.3.exe : Ex01.3.o random.o
	$(CC) random.o Ex01.3.o -o Ex01.3.exe
Ex01.3.o : Ex01.3.cpp
	$(CC) -c Ex01.3.cpp -o Ex01.3.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o Ex01.3.exe seed.out
