/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_NVT.h"

using namespace std;

int main() { 
	Input(); //Inizialization
	Equilibrate();
	int nconf = 1;
	for(int iblk=1; iblk <= nblk; ++iblk) {
	//Simulation 
	    Reset(iblk);   //Reset block averages
		for(int istep=1; istep <= nstep; ++istep) {
			Move();
			Measure();
			Accumulate(); //Update block averages
//			PrintMeasure (nstep*(iblk-1)+istep);
			if(istep%10 == 0){
//			ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
				nconf += 1;
			}
		}
	Averages(iblk);   //Print results for current block
	}
	ConfFinal(); //Write final configuration

	return 0;
}


void Input(void) {
	ifstream ReadInput,ReadConf;

	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;


  
//Read input informations
	ReadInput.open("input.dat");
	
	ReadInput >> Restart;
	ReadInput >> Realsimulation;

	if (Restart == 0) rnd.Initialize("seed.in");
	if (Restart == 1) rnd.Initialize("seed.out");
	else
	if (Realsimulation == 1) {
		ifstream ReadProperty;
		ReadProperty.open("Real_Model.dat");
		ReadProperty >> name;
		cout << endl << "Simulating " << name << endl << endl;

		ReadProperty >> sigma;
		ReadProperty >> epsilon_kb;
		ReadProperty >> mass;
		ReadProperty.close();

	//If first run of simulation clean directory of same element
	//only if realsimulation
		string command = "rm " + name +"/*";
		const char* tcommand = command.c_str();
		system(tcommand);	//Command to delete de file if existing
	
	}

	ReadInput >> temp;
	beta = 1.0/temp;
	cout << "Temperature = " << temp << endl;
	
	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;

	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	box = pow(vol,1.0/3.0);
	cout << "Volume of the simulation box = " << vol << endl;
	cout << "Edge of the simulation box = " << box << endl;

	ReadInput >> rcut;
	cout << "Cutoff of the interatomic potential = " << rcut << endl << endl;
   
//Tail corrections for potential energy and pressure
	vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
	ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
	cout << "Tail correction for the potential energy = " << vtail << endl;
	cout << "Tail correction for the virial           = " << ptail << endl; 

	ReadInput >> delta;

	ReadInput >> nblk;

	ReadInput >> nstep;

	cout << "The program perform Metropolis moves with uniform translations" << endl;
	cout << "Moves parameter = " << delta << endl;
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;
	ReadInput.close();

	

//Prepare arrays for measurements
	iv = 0; //Potential energy
	iw = 1; //Virial
 
	n_props = 2; //Number of observables

//measurement of g(r)
	igofr = 2;
	nbins = 100;
	n_props = n_props + nbins;
	bin_size = (box*0.5)/(double)(nbins);		//just up to L/2 for problems of simmetry in a square

//Select acceptance limits (for gas 50%rule does not work)
	min_acc = 0.45;
	max_acc = 0.55;

	if (temp == 1.2) {	//Gas temperature
		min_acc = 0.53;
		max_acc = 0.63;
	}



//Read initial configuration
	if (Restart == 1) {
		cout << "Read initial configuration from file config.final " << endl << endl;
		ReadConf.open("config.final");
		for (int i=0; i<npart; ++i) {
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = Pbc( x[i] * box );
			y[i] = Pbc( y[i] * box );
			z[i] = Pbc( z[i] * box );
		}
		ReadConf.close();
	} else {

		system ("./clean.sh");
		cout << "Read initial configuration from file config.0 " << endl << endl;
		ReadConf.open("config.0");
		for (int i=0; i<npart; ++i) {
			ReadConf >> x[i] >> y[i] >> z[i];
			x[i] = Pbc( x[i] * box );
			y[i] = Pbc( y[i] * box );
			z[i] = Pbc( z[i] * box );
		}
		ReadConf.close();
	}
//Evaluate potential energy and virial of the initial configuration
	Measure();

//Print initial values for the potential energy and virial
	cout << "Initial potential energy (with tail corrections) = " << walker[iv]/(double)npart + vtail << endl;
	cout << "Virial                   (with tail corrections) = " << walker[iw]/(double)npart + ptail << endl;
	cout << "Pressure                 (with tail corrections) = " << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl << endl;
}


void Move(void) {
	int o;
	double p, energy_old, energy_new;
	double xold, yold, zold, xnew, ynew, znew;


	for(int i=0; i<npart; ++i) {
//Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
		o = (int)(rnd.Rannyu()*npart);

		//Old
		xold = x[o];
		yold = y[o];
		zold = z[o];

		energy_old = Boltzmann(xold,yold,zold,o);

		//New
		xnew = Pbc( x[o] + delta*rnd.Rannyu(-0.5, 0.5) );
		ynew = Pbc( y[o] + delta*rnd.Rannyu(-0.5, 0.5) );
		znew = Pbc( z[o] + delta*rnd.Rannyu(-0.5, 0.5) );

		energy_new = Boltzmann(xnew,ynew,znew,o);

		//Metropolis test
		p = exp(beta*(energy_old-energy_new));
		if(p >= rnd.Rannyu()) {
			//Update
			x[o] = xnew;
			y[o] = ynew;
			z[o] = znew;
    
			accepted = accepted + 1.0;
		}
		attempted = attempted + 1.0;
	}
}

double Boltzmann(double xx, double yy, double zz, int ip) {
	double ene=0.0;
	double dx, dy, dz, dr;

	for (int i=0; i<npart; ++i)  {
		if(i != ip) {
			// distance ip-i in pbc
			dx = Pbc(xx - x[i]);
			dy = Pbc(yy - y[i]);
			dz = Pbc(zz - z[i]);

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			if(dr < rcut) ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
		}
	}

	return 4.0*ene;
}

void Measure() {
	int bin;
	double v = 0.0, w = 0.0;
	double vij, wij;
	double dx, dy, dz, dr;

//reset the hystogram of g(r)
	for (int k=igofr; k<n_props; ++k) walker[k]=0.0;

//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i) {
		for (int j=i+1; j<npart; ++j) {

			// distance i-j in pbc
			dx = Pbc(x[i] - x[j]);
			dy = Pbc(y[i] - y[j]);
			dz = Pbc(z[i] - z[j]); 

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);
//cout << dr << endl;
	//update of the histogram of g(r)
			bin = dr/bin_size;			
//cout << bin << endl;
			walker[igofr+bin] += 2.0;
			
			if(dr < rcut) {
				vij = 1.0/pow(dr,12) - 1.0/pow(dr,6);
				wij = 1.0/pow(dr,12) - 0.5/pow(dr,6);

	// contribution to energy and virial
				v += vij;
				w += wij;
			}
		}
	}
	walker[iv] = 4.0 * v;
	walker[iw] = 48.0 * w / 3.0;
	v = (4.0 * v + vtail)/(double) npart;
	w = rho * temp + (48.0 * w / 3.0 + ptail * (double)npart) / vol;
}

void Reset(int iblk) {
	//Reset block averages
   
	if(iblk == 1) {
		for(int i=0; i<n_props; ++i) {
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i) blk_av[i] = 0;
//	for(int i=0; i<nbins; ++i) bin_counter[i] = 0;
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}


void Accumulate(void) {
	//Update block averages

	for(int i=0; i<n_props; ++i) blk_av[i] = blk_av[i] + walker[i];
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) { //Print results for current block
    
	double r;
	double gdir = 0;
	ofstream Gofr, Gave, Epot, Pres;
	const int wd=12;
   
	cout << "Block number " << iblk << endl;
	cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
	Epot.open("output.epot.dat",ios::app);
	Pres.open("output.pres.dat",ios::app);
	Gofr.open("output.gofr.dat",ios::app);
	Gave.open("output.gave.dat",ios::app);
    
	stima_pot = blk_av[iv]/blk_norm/(double)npart + vtail; //Potential energy
	glob_av[iv] += stima_pot;
	glob_av2[iv] += stima_pot*stima_pot;
	err_pot=Block_Error(iblk, glob_av[iv],glob_av2[iv]);
    
	stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
	glob_av[iw] += stima_pres;
	glob_av2[iw] += stima_pres*stima_pres;
	err_press=Block_Error(iblk, glob_av[iw],glob_av2[iw]);

	//Potential energy per particle
	Epot << iblk << scientific <<  "	" << stima_pot << "	" << glob_av[iv]/(double)iblk << "	" << err_pot << endl;
	//Pressure
	Pres << iblk << scientific <<  "	" << stima_pres << "	" << glob_av[iw]/(double)iblk << "	" << err_press << endl;

//g(r)
	ofstream si_Gdir;
	if (Realsimulation == 1) si_Gdir.open(name + "/" + name + "_output.gofr.dat",ios::app);		//File for realsimulation

	for (int i=igofr; i<n_props; ++i) {
		r = (double)(i-igofr)*bin_size;
		double vr = (4.0/3.0)*pi*(pow(r+bin_size, 3) - pow(r, 3));
		double norm = rho*npart*vr;
		gdir = blk_av[i]/(blk_norm*norm);
		glob_av[i] += gdir ;
		glob_av2[i] += gdir*gdir;
		err_gdir=Block_Error(iblk, glob_av[i],glob_av2[i]);
		r = (double)(i-igofr)*bin_size + bin_size*0.5;
		Gofr << iblk << "	" << i-igofr << "	" << r << "	" << gdir << "	" << glob_av[i]/(double)(iblk) << "	" << err_gdir  << "\t\t\t" << blk_av[i]/blk_norm << endl;
		if (iblk == nblk) Gave << i-igofr << "	" << r << "	" << gdir << "	" << glob_av[i]/(double)(iblk) << "	" << err_gdir << endl;
		if (Realsimulation == 1 and iblk == nblk) si_Gdir << i-igofr << "\t" << setprecision (7) << r*sigma << "	" << glob_av[i]/(double)(iblk) << "\t" << err_gdir << endl; //coorection for IS units
	}

	if (Realsimulation == 1) si_Gdir.close();

	cout << "----------------------------" << endl << endl;

	Epot.close();
	Pres.close();
	Gofr.close();
	Gave.close();

	if (Realsimulation == 1) {
		ofstream si_Epot, si_Press; //International system average measure units files
		double pot, press, si_err_pot, si_err_press;
	
		si_Epot.open(name + "/" + name + "_output.epot.dat",ios::app);
		si_Press.open(name + "/" + name + "_output.press.dat",ios::app);

	//Corrections for IS units
		pot = glob_av[iv]/(double)(iblk) * epsilon_kb*kb; 						//Potential energy
		press = glob_av[iw]/(double)(iblk) * 1E27*(epsilon_kb*kb/pow(sigma,3));	//Pressure (1E27 is a factor to adjust the dimension of sigma that is in nm)

		si_err_pot = Block_Error ( iblk, glob_av[iv]*epsilon_kb*kb, glob_av2[iv]*(epsilon_kb*kb)*(epsilon_kb*kb));	
		si_err_press = Block_Error ( iblk, glob_av[iw]*1E27*(epsilon_kb*kb/pow(sigma,3)), glob_av2[iw]*(1E27*(epsilon_kb*kb/pow(sigma,3)))*(1E27*(epsilon_kb*kb/pow(sigma,3))));

		si_Epot << iblk << "\t" << setprecision (7) << scientific << pot << "\t" << si_err_pot << endl;
		si_Press << iblk << "\t" << setprecision (7) << scientific << press << "\t" << si_err_press << endl;
		
		si_Epot.close();
		si_Press.close();
	}	
}


void ConfFinal(void) {
	ofstream WriteConf, WriteSeed;

	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");
	for (int i=0; i<npart; ++i) WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;

	WriteConf.close();
	rnd.SaveSeed();
}

void ConfXYZ(int nconf) { //Write configuration in .xyz format
	ofstream WriteXYZ;

	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for (int i=0; i<npart; ++i) WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	
	WriteXYZ.close();
}

double Pbc(double r) { //Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r/box);
}

double Block_Error (int block, double sum, double sum2) {
	//block is number of block	
	//sum is cumulative sum
	//sum2 is squred cumulative sum
	if (block == 1) return 0;
	else return sqrt ((sum2/(double)(block) - (sum*sum)/(block*block))/(block-1.));
}

void Equilibrate () {
	Reset(0);
	for (int i=0; i<500; ++i) Move();
//Reset initial conditions for next equilibration or for starting simulation
	ifstream ReadConf;
	ReadConf.open("config.0");
	for (int i=0; i<npart; ++i) {
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = Pbc( x[i] * box );
		y[i] = Pbc( y[i] * box );
		z[i] = Pbc( z[i] * box );
	}
	ReadConf.close();
	if (accepted/attempted < min_acc or accepted/attempted > max_acc) {
		cout << "Equlibrating..." << endl;
		delta += 0.05;
		Equilibrate ();
	}
	else {
		cout << "Equilibrated!" << endl;
		cout << "Using delta = " << delta << endl << endl;
	}

	return;
}

void PrintMeasure(int istep) {
	if (istep == 1) system("rm Autocorrelation/*");
	ofstream WriteEpot, WritePress;
	WriteEpot.open("Autocorrelation/ist_epot.dat", ios::app);
	WritePress.open("Autocorrelation/ist_press.dat", ios::app);
	WriteEpot << istep << "\t" << walker[iv]/(double)npart + vtail << endl;
	WritePress << istep << "\t" << rho * temp + (walker[iw] + (double)npart * ptail) / vol << endl;
	
	WriteEpot.close();
	WritePress.close();
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
