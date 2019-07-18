===============================================
****************Exercises_10*******************
********Numerical Simulation Laboratory********
**************Matteo Martinelli****************
===============================================

Simulated Anealing

Legend of files:

Serial
	directory with main, header and results for serial simulation

	Simulated Annealing.cpp
		main

	Simulated Annealing.h
		header with all global variables and functions declaration

	Simulated Annealing.exe 
		executable

	High_T
		result from a simulation starting from higher T (see jupyter for more information)

Parallel
	directory with main, header and results for parallel simulation

	Parallel_Sim_Ann.cpp
		main

	Parallel_Sim_Ann.h
		header with all global variables and functions declaration

	Parallel_Sim_Ann.exe 
		executable



In both drectories are

path.h
	header file for class path used in main

path.cpp
	class path used in main

bin.out
	extraction per index in path (just to check whether the algorithm works and samples preferentiall first paths)

cts_pos_out
	map with coordinates of citites ( x \t y )

'Geometry'
	directory with result for simulation in 'geometry'

input.dat
	file from which input values are read
	
		ReadInput >> ncts;					'int' number of cities
		ReadInput >> temp;					'double' starting temperature for simulation
		ReadInput >> nstep;					'int' starting number of steps per block
		ReadInput >> nblk;					'int' number of blocks
		ReadInput >> geometry; 				'bool' select displayement of citites: 1 = square, 0 = circle

seed.in
    random seed for generator

Primes
    other seed for random generator, must have for parallel simulation

Exercises_10.ipynb 
	Analysis of the data of simulation already performed


How to produce data:

Serial
	choose input data
	compile and run the simulation, the output file will appear in the respective directory

Parallel
	choose input data
	compile with makefile
	execution with command mpiexec -np 'n' ./Parallel_Sim_Ann.exe
		where 'n' is the number of cores for parallel run

-make clean
	removes. o
	removes executable 
	removes saved seed

