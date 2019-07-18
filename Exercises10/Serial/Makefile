CC = g++
CFLAGS = -Wall -O3 --std=c++11

Simulated_Annealing.exe : Simulated_Annealing.o random.o path.o
	$(CC) random.o path.o Simulated_Annealing.o -o Simulated_Annealing.exe
Simulated_Annealing.o : Simulated_Annealing.cpp
	$(CC) -c Simulated_Annealing.cpp -o Simulated_Annealing.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
path.o : path.cpp path.h
	$(CC) -c path.cpp -o path.o $(CFLAGS)
clean :
	rm *.o Simulated_Annealing.exe seed.out
