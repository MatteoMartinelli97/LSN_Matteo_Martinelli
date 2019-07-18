CC = mpicxx

Parallel_Sim_Ann.exe : Parallel_Sim_Ann.o random.o path.o
	$(CC) random.o path.o Parallel_Sim_Ann.o -o Parallel_Sim_Ann.exe
Parallel_Sim_Ann.o : Parallel_Sim_Ann.cpp
	$(CC) -c Parallel_Sim_Ann.cpp -o Parallel_Sim_Ann.o
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o
path.o : path.cpp path.h
	$(CC) -c path.cpp -o path.o
clean :
	rm *.o Parallel_Sim_Ann.exe seed.out
