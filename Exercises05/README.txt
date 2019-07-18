===============================================
****************Exercises_05*******************
********Numerical Simulation Laboratory********
**************Matteo Martinelli****************
===============================================

Legend of files:
Metropolis
	directory with main and input files

Metropolis.cpp
	main and functions definition

Metropolis.h
	header with all global variables and functions declaration	

Metropolis.exe 
	executable

.out
	files with average outputs

input.dat
	file from which input values are read
	
		ReadInput >> state;				'string' name of state 3 options are possible '1s', '2p', '3d'
		ReadInput >> pdf;				'bool' 0 means uniform, 1 means gaussian probability for Metropolis sampling
		ReadInput >> delta;				'double' parameter of pdf (suggested 0.1 as input, the code will equilibrate)
		ReadInput >> x;					'double' starting position of electron along x axis
		ReadInput >> y;					'double' starting position of electron along y axis
		ReadInput >> z;					'double' starting position of electron along z axis
		ReadInput >> nsteps;			'int' number of total steps 
		ReadInput >> nblock;			'int' number of blocks in which divide 'nstep' (must be 'nstep' = k * 'nblock' or trunkating)

Gauss_T 
	directory with results for gaussian sampling


Exercises_05.ipynb 
	Analysis of the data of simulation already performed


How to produce data:

Metropolis
	choose input data
	run the simulation, the output file will appear in the same directory if uniform, and in Gauss_T otherwise

-clean.sh
	removes all data files and frames
	not necessary since programs do not append, if running you will lose the simulations till now performed

-make clean
	removes. o
	removes executable 
	removes saved seed

