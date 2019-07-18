===============================================
****************Exercises_07*******************
********Numerical Simulation Laboratory********
**************Matteo Martinelli****************
===============================================

Monte Carlo NVT

Legend of files:
MonteCarlo_NVT
	directory with main and input files

Monte_Carlo_NVT.cpp
	main and functions definition

MonteCarlo_NVT.h
	header with all global variables and functions declaration	

MonteCarlo_NVT.exe 
	executable


output_*.dat
	files with average values outputs

input.dat
	file from which input values are read
	
		ReadInput >> Restart; 				'bool' activates restarting option 
		ReadInput >> Realsimulation;		'bool' activates SI measures 
		ReadInput >> temp;					'double' input temperature
		ReadInput >> npart;					'int' number of particles (if changing remember to change also MolDyn.h since there's the declarion of 'const int m_part')
		ReadInput >> rho;					'double' input density
		ReadInput >> rcut;					'double' input cut-off radius
		ReadInput >> delta;					'double' parameter for metropoli probability
		ReadInput >> nblk;					'int' number of blocks
		ReadInput >> nstep;					'int' number of steps per block

		default values are in files input.* where '*' is the phase to simulate

Real_Model.dat
	file from which properties for simulation of real element are read
	
		ReadProperty >> name;				'string' name of element (select in which directory write the outputs)
		ReadProperty >> sigma;				'double' sigma value in Lennard-Jones parametrization (nm)
		ReadProperty >> epsilon				'double' epsilon/k_b value in Lennard-Jones parametrization (K)
		ReadProperty >> mass;				'double' molecular mass for element (uma)

		default values are in files Real_Model.* where '*' is the element to simulate

config.0
	input configuration
	default FCC configuration is in config.fcc

config.final
	output configuration after simulation (needed to restart)

seed.in
    random seed for generator

Primes
    other seed for random generator
    
'Element'/*
	directory with measures in SI units for 'Element' 
	'Element'/'phase'
		Saved values of already produced simulation for 'Element' in phase 'phase'			

frames
	directory with printed position step by step, useful for animations

Autocorrelation
	directory with long run and instantaneous values for computing autocorrelation (see jupyter)


Exercises_07.ipynb 
	Analysis of the data of simulation already performed


How to produce data:

Monte_Carlo_NVT
	choose input data
	run the simulation, the output file will appear in the same directory in reduce units, and in the directory with the name of the element for SI units

-clean.sh
	removes all data files and frames
	not necessary since program calls if the case

-make clean
	removes. o
	removes executable 
	removes saved seed

