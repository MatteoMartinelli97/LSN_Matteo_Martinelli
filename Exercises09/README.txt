===============================================
****************Exercises_09*******************
********Numerical Simulation Laboratory********
**************Matteo Martinelli****************
===============================================

Genetic Algorithm

Legend of files:

Genetic_Opt.cpp
	main and functions definition

Genetic_Opt.h
	header with all global variables and functions declaration	

Genetic_Opt.exe 
	executable

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
		ReadInput >> npop;					'int' number of population (suggested ncts*ncts)
		ReadInput >> ngen;					'int' number of generation
		ReadInput >> geometry; 				'bool' select displayement of citites: 1 = square, 0 = circle

seed.in
    random seed for generator

Primes
    other seed for random generator

Exercises_09.ipynb 
	Analysis of the data of simulation already performed


How to produce data:

./
	choose input data
	compile and run the simulation, the output file will appear in the respective directory


-make clean
	removes. o
	removes executable 
	removes saved seed

