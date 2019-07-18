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
#include <string>
#include <cmath>
#include <iomanip>
#include "Metropolis.h"

using namespace std;
 
int main (){
	Input();
	for (int iblock=1; iblock<=nblock; iblock++) {
		for (int istep=1; istep<=bsteps; istep++) {
			Move();
			PrintPosition();	
		}
		BlockMeasure(iblock);
		cout << "-------------------------------" << endl;
		cout << "Block " << iblock << endl;
		cout << "Acceptance rate = " << acount/(double) bsteps << endl;
		acount = 0;
	}
	
rnd.SaveSeed();
return 0;
}	

void Input() {
	ifstream ReadInput;
	ReadInput.open("input.dat"); //Read input

	cout << "Hydrogen Atom probability distributions" << endl;
	cout << "The program uses Bohr radius as lenght unit" << endl << endl;

	//Initialize random generator
	rnd.Initialize();
	
	ReadInput >> state;

	cout << "Simulating " + state + " state" << endl << endl;

	ReadInput >> pdf;
	ReadInput >> delta;
	ReadInput >> x;
	ReadInput >> y;
	ReadInput >> z;

	r = sqrt(x*x + y*y + z*z);

//Delete old files in order to substitute them with new simulation
	string command;
	if (pdf == 0) command = "rm " + state +"_position.out " + state+"_r.out";
	if (pdf == 1) command = "rm Gauss_T/" + state +"_position.out " + "Gauss_T/"+state+"_r.out";
	const char* tcommand = command.c_str();
	system(tcommand);	//Command to delete de file if existing

//Equilibration
	x_start = x;
	y_start = y;
	z_start = z;
	Equilibrate ();
	
	ReadInput >> nsteps;
	ReadInput >> nblock;

	ReadInput.close();

	bsteps = nsteps/nblock;
	r_sum = 0;
	glob_av = 0;
	glob_av = 0;

}

void Move () {
	double A, Tr, Tx, Ty, Tz, q;
	double px, py;

	if (pdf == 0) {
		Tx = x + delta*rnd.Rannyu(-0.5,0.5);
		Ty = y + delta*rnd.Rannyu(-0.5,0.5);
		Tz = z + delta*rnd.Rannyu(-0.5,0.5);
	}

	if (pdf == 1) {
		Tx = rnd.Gauss(x, delta);
		Ty = rnd.Gauss(y, delta);
		Tz = rnd.Gauss(z, delta);
	}

	Tr = sqrt(Tx*Tx + Ty*Ty + Tz*Tz);

	if (state == "1s") {
		px = exp(-2.0*r);
		py = exp(-2.0*Tr);
	}
	if (state == "2p") { 
		px = r*r*exp(-r)*(z/r)*(z/r);
		py = Tr*Tr*exp(-Tr)*(Tz/Tr)*(Tz/Tr);
	}
	if (state == "3d") {
		px = r*r*exp(-2./3.*r)*pow(3.0 *(z/r)*(z/r) - 1.0,2);
		py = Tr*Tr*exp(-2./3.*Tr)*pow(3.0 *(Tz/Tr)*(Tz/Tr) - 1.0,2);
	}

	q = py/px;
	r_sum += r;
	A = min (1., q);
	if (Accept(A)) {
		x = Tx;
		y = Ty;
		z = Tz;
		r = Tr;
		acount++;
	}

	return;
}


bool Accept (double A) {
	if (A == 1.) return true;
	else {
		double k = rnd.Rannyu();
		if (k <= A) return true;
		else return false;
	}
}

void PrintPosition (void){ //Write final configuration
	ofstream WritePosition;
	if (pdf == 0) WritePosition.open(state+"_position.out",ios::app);
	if (pdf == 1) WritePosition.open("Gauss_T/"+state+"_position.out",ios::app);

	WritePosition << x << "\t" <<  y << "\t" << z << endl;
	WritePosition.close();
	return;
}

void BlockMeasure (int block) {
	ofstream WriteMeasure;
	if (pdf == 0) WriteMeasure.open(state+"_r.out",ios::app);
	if (pdf == 1) WriteMeasure.open("Gauss_T/"+state+"_r.out",ios::app);
	r_mean = r_sum/(double) (bsteps);
	glob_av += r_mean;
	glob_av2+= r_mean*r_mean;
	if (block == 1) r_std = 0;
	else r_std = sqrt ((glob_av2/(double)(block) - (glob_av*glob_av)/(double)(block*block))/(block-1.));
	r_sum = 0;
	WriteMeasure << block << "\t" <<  glob_av/(double)(block) << "\t" << r_std << endl;
	WriteMeasure.close();
	return;
}

void Equilibrate () {
	acount = 0;
	for (int i=0; i<eqsteps; i++) Move();
	x = x_start;
	y = y_start;
	z = z_start;
	r = sqrt(x*x + y*y + z*z);
	if (acount/(double)(eqsteps) < 0.47 or acount/(double)(eqsteps) > 0.53) {
		cout << "Equlibrating..." << endl;
		delta += 0.3;
		Equilibrate ();
	}
	else {
		cout << "Equilibrated!" << endl;
		cout << "Using delta = " << delta << endl;
	}

	return;
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
