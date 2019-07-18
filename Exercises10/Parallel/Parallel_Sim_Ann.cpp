/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "mpi.h"
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <numeric>
#include "Parallel_Sim_Ann.h"

using namespace std;

int main (int argc, char *argv[]){

	MPI::Init(argc, argv);
	int size = MPI::COMM_WORLD.Get_size();
	int rank = MPI::COMM_WORLD.Get_rank();

		Input(rank);
		if (Check_val) return -1;
		if (rank == 0) cout << "Running on " << size << " cores" << endl;
		
		dbl_int best_received_l;
		best_received_l.val  = 0;
		best_received_l.rank = -1;
	
		for (int iblk=1; iblk<=nblk; ++iblk) {		//cycle over blocks
			Reset();
			for (int istep=0; istep<nstep; istep++) {	//cycle over steps in block
				Move(0);
			//Check if everything is ok with new path
				if (Check_val) return -1;
			}

		//For all process compute distance and check if better than his best
			//if it is save value of length, rank and best path
			double l = Travel_Path.GetLength(dist);
			if ( l < best_l.val) {
				Best_Path = Travel_Path;
				best_l.val = l;
				best_l.rank = rank;
			}
			
		//Process rank == 0 receives minimum path and which process has it
			best_received_l.val = 0;
			MPI_Reduce(&best_l, &best_received_l, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI::COMM_WORLD);
		
		//Process rank == 0 writes on a file the best length	
			if (rank == 0)	Print(iblk, best_received_l.val);
		//Change parameters and simulate again
			New_Params();
		}
		cout << rank <<  " Best path length = " << Best_Path.GetLength(dist) << endl;

	//All processes receive from 0 the rank of final best path
		int best_rank = -1;
		if (rank == 0)  best_rank = best_received_l.rank;
		MPI_Bcast (&best_rank, 1, MPI_INTEGER, 0, MPI::COMM_WORLD );

	//The process with best path writes on a file the coordinates
		if (rank == best_rank) Print_Best_Path_Coord ();

	MPI::Finalize();

return 0;
}


void Input(int rank) {
	ifstream ReadInput;
	
	rnd.Initialize();
	
ReadInput.open("input.dat");

	ReadInput >> ncts;
	ReadInput >> temp;
	ReadInput >> nstep;
	ReadInput >> nblk;

	if (rank == 0 ) {
		cout << "Simulated Annealing		" << endl;
		cout << "Travelling Salesman Problem		" << endl << endl;

		cout << "There are " << ncts << " cities to visit	"  << endl;

		cout << "Starting simulation at temperature " << temp << endl;

		cout << "Simulation of " << nblk << " different temperatures" << endl;
		cout << "Scaling T by factor " << t_incr << endl;
		cout << "Increasing step +" << s_incr << "per block, with max = " << max_steps << endl;
		cout << "Distances are calculated with L^2 norm" << endl;
	}

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
		if (rank == 0) cout << "Cities are displayed on a circle" << endl;
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
		if (rank == 0)cout << "Cities are displayed inside a square" << endl;
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
	if(rank == 0) {
		string command = "rm " + geo_name + "/*length.out";
		const char* tcommand = command.c_str();
		system(tcommand);
	}

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
	best_l.val = Best_Path.GetLength(dist);
	rnd.Par_Initialize(rank+1);
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

void Print (int iblk, double bestlength) {
	ofstream BestLength;
	BestLength.open(geo_name + "/best_length.out", ios::app);	
	BestLength << iblk << setprecision(15) << scientific  << "\t" << temp << "\t" << bestlength << endl;
}

void New_Params () {
	temp = temp*t_incr;
	nstep = nstep + s_incr;
	if (nstep > max_steps) {
		nstep = max_steps;
		s_incr = 0;
	}
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
