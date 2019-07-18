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
#include <numeric>
#include "Simulated_Annealing.h"

using namespace std;

int main (int argc, char *argv[]){

	Input();
	if (Check_val) return -1;

	for (int iblk=1; iblk<=nblk; ++iblk) {		//cycle over blocks
		Reset();
		cout << "Block " << iblk << endl;
		cout << "Block temperature = " << temp << endl;
		cout << "Block steps = " << nstep << endl;
		for (int istep=0; istep<nstep; istep++) {	//cycle over steps in block
			Move(0);
		//Check if everything is ok with new path
			if (Check_val) return -1;
		}

		for (int i=0; i<4; ++i)	cout << "Acceptance Rate of mutation " << i << " = " << accepted[i]/attempted[i] << endl;
		cout << endl << "--------------------------------" << endl;
		double l = Travel_Path.GetLength(dist);
		if ( l < best_l) {
			Best_Path = Travel_Path;
			best_l = l;
		}			
		Print(iblk);
		New_Params();
	}
	cout << "Best path length = " << Best_Path.GetLength(dist) << endl;
	Print_Best_Path_Coord ();

return 0;
}


void Input() {
	ifstream ReadInput;
	
	rnd.Initialize();
	
	cout << "Simulated Annealing		" << endl;
	cout << "Travelling Salesman Problem		" << endl << endl;
	
	ReadInput.open("input.dat");

	ReadInput >> ncts;
	cout << "There are " << ncts << " cities to visit	"  << endl;

	ReadInput >> temp;
	cout << "Starting simulation at temperature " << temp << endl;

	ReadInput >> nstep;
	ReadInput >> nblk;
	cout << "Simulation of " << nblk << " differernt temperatures each with " << nstep << "steps	" << endl;
	cout << "Distances are calculated with L^2 norm" << endl;

	ReadInput >> geometry;
	ReadInput.close();

//Initialize value for checking moves	
	Check_val = false;

//Starting random path
	Travel_Path.Initialize(ncts);
	Check_val = Travel_Path.Check();
	if (Check_val) return;	

//Vectors of cities coordinates
	x = new double [ncts];
	y = new double [ncts];
	ofstream WriteCts;
	if (geometry == 0) {
		cout << "Cities are displayed on a circle" << endl;
		geo_name = "Circle";
		WriteCts.open(geo_name +"/cts_pos.out");
		//Circle dispositions
		for (int i=0; i<ncts; ++i) {	
			double theta = 2.0*rnd.Angle();
			x[i] = cos(theta);
			y[i] = sin(theta);
			WriteCts << i << "\t" << x[i] << "\t" << y[i] << endl;
		}
	} else {
		cout << "Cities are displayed inside a square" << endl;
		geo_name = "Square";
		WriteCts.open(geo_name +"/cts_pos.out");
		//Inside square disposition
		for (int i=0; i<ncts; ++i) {	
			x[i] = rnd.Rannyu(0.0, 1.0);
			y[i] = rnd.Rannyu(0.0, 1.0);
			WriteCts << i << "\t" << x[i] << "\t" << y[i] << endl;
		}
	}
	WriteCts.close();

//If existing file of mean_l and best_l, remove
	string command = "rm " + geo_name + "/*length.out";
	const char* tcommand = command.c_str();
	system(tcommand);

	dist = new double* [ncts];
	for (int i=0; i<ncts; ++i) dist[i] = new double [ncts];
//Fill distances matrix (Dij = distance between city i-j = j-i)
//save all distances to optimize calulation
	for (int i=0; i<ncts; ++i) {
		for (int j=i; j<ncts; ++j) {
			double x2 = (x[i]-x[j])*(x[i]-x[j]);
			double y2 = (y[i]-y[j])*(y[i]-y[j]);
			dist[i][j] = x2 + y2;		//Symmetric matrix
			dist[j][i] = x2 + y2;
		}
	}

	Best_Path = Travel_Path;
	best_l = Best_Path.GetLength(dist);
}


void Move (int imut) {
	Path Trial_Path (ncts);
	Trial_Path = Travel_Path;

	if (imut == 0) {			//Swap 2 random cities
		int cty_1 = rnd.Rannyu(0, ncts);
		int cty_2 = rnd.Rannyu(0, ncts);
		Trial_Path.Swap(cty_1, cty_2);
	}
	if (imut == 1) {			//Swap m random cities 
		int cty_1 = rnd.Rannyu(0, ncts);
		int m = rnd.Rannyu(0, ncts/2);
		int cty_2 = rnd.Rannyu(0, ncts-2*m+1);
		Trial_Path.Swap(cty_1, cty_2, m);
	}
	if (imut == 2) {			//Shift +n, m cities from cty_1
		int cty_1 = rnd.Rannyu(0, ncts);
		int m = rnd.Rannyu(1, ncts+1);
		int plus_shift = rnd.Rannyu(m, ncts);
		Trial_Path.Shift(plus_shift, m, cty_1);
	}
	if (imut == 3) {			//Reverse m random cities from cty1
		int cty_1 = rnd.Rannyu(0, ncts);
		int m = rnd.Rannyu(0, ncts);
		Trial_Path.Reverse(cty_1, m);
	}

	Check_val = Trial_Path.Check();
	if (Check_val) return;

	double length_diff = Trial_Path.GetLength(dist) - Travel_Path.GetLength(dist);
	double q = exp (-beta*length_diff);
	double A = min (1., q);
	if (Accept(A)) {
		Travel_Path = Trial_Path;
		accepted[imut] += 1.0;
	}
	attempted[imut] += 1.0;
	if (imut!=3) Move(imut+1);
	return;
}

void Reset() {//Reset block averages
	for (int i=0; i<4; ++i) {
		attempted[i] = 0.0;
		accepted[i] = 0.0;
	}
	beta = 1.0/temp;
}


bool Accept (double A) {
	if (A == 1.) return true;
	else {
		double k = rnd.Rannyu();
		if (k <= A) return true;
		else return false;
	}
}

void Print_Best_Path_Coord () {
	ofstream WriteBPath;
	WriteBPath.open(geo_name + "/best_path.out");
	for (int i=0; i<ncts; ++i) WriteBPath << Best_Path.GetPath()[i] << "\t" << x[Best_Path.GetPath()[i]] << "\t" << y[Best_Path.GetPath()[i]] << endl;
	WriteBPath << Best_Path.GetPath()[0] << "\t" << x[Best_Path.GetPath()[0]] << "\t" << y[Best_Path.GetPath()[0]] << endl;
	return;
}

void Print (int iblk) {
	ofstream BestLength;
	BestLength.open(geo_name + "/best_length.out", ios::app);	
	BestLength << iblk << setprecision(15) << scientific << "\t" << temp << "\t" << best_l << endl;
	return;
}

void New_Params () {
	temp = temp*t_incr;
	nstep = nstep + s_incr;
	if (nstep > max_steps) {
		nstep = max_steps;
		s_incr = 0;
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
