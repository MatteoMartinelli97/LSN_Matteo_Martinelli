CC = g++
CFLAGS = -Wall -O3 --std=c++11

Ex01.2.exe : Ex01.2.o random.o histo.o
	$(CC) random.o Ex01.2.o histo.o -o Ex01.2.exe
Ex01.1.o : Ex01.1.cpp
	$(CC) -c Ex01.2.cpp -o Ex01.2.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
histo.o : histo.cpp histo.h
	$(CC) -c histo.cpp -o histo.o $(CFLAGS)
clean :
	rm *.o Ex01.2.exe seed.out
