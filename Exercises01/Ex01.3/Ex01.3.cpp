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
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

	Random rnd;
	rnd.Initialize();

/****************************************
Buffon Experiment
****************************************/

	ofstream output;
	output.open("Pi.out");
	double L = 0.75;
	double d = 1.;
	int M = 1000000;	//throws
	int N = 100;		//blocks
	int S = M/N;		//steps per block
	double sum_pi = 0;
	double sum_pi2 = 0;
	double mean_pi = 0;
	double mean_stddev = 0;

	for (int i=0; i<N; i++) {			//blocks cycle 
		int n_hit = 0;	
		for (int j=0; j<S; j++) {		//experiment in each block
			double theta = rnd.Angle();
			double x1 = rnd.Rannyu(1., 10.);
			double x2 = x1 + L*cos(theta);
			int r1 = x1/d;
			int r2 = x2/d;
			if (r1 != r2) n_hit++;
		}
		double pi = (2.*L*(double) (S))/((double) (n_hit)*d);
		sum_pi += pi;
		sum_pi2 += pi*pi;
		mean_pi = sum_pi/(i+1);
		output << setprecision(7) << fixed << (i+1)*S << "\t" << mean_pi << "\t";
		if (i==0) mean_stddev = 0;
		else mean_stddev = sqrt (((sum_pi2/(i+1.)) - (mean_pi*mean_pi))/i);
		output << mean_stddev << endl;
	}

	output.close();
	rnd.SaveSeed();

return 0;
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
