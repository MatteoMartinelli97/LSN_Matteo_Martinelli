CC = g++
CFLAGS = -Wall -O3 --std=c++11

Genetic_Opt.exe : Genetic_Opt.o random.o path.o
	$(CC) random.o path.o Genetic_Opt.o -o Genetic_Opt.exe
Genetic_Opt.o : Genetic_Opt.cpp
	$(CC) -c Genetic_Opt.cpp -o Genetic_Opt.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
path.o : path.cpp path.h
	$(CC) -c path.cpp -o path.o $(CFLAGS)
clean :
	rm *.o Genetic_Opt.exe seed.out
