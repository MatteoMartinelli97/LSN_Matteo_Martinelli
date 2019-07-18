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
double Block_StdDev (int block, double cumsum, double cumsum2);

int main (int argc, char *argv[]){

	Random rnd;
	rnd.Initialize();
/*****************************************
1) Uniform sampling for Integral evaluation

*****************************************/
	ofstream output;
	output.open("Integral.out");
	
	int M = 1000000;
	int N = 1000;
	int S = M/N; 	//Steps per block
	double A = 0;		//block result
	double sum_A = 0;	//cumulative average
	double sum_A2 = 0; 	//cumulative square average
	double mean = 0;
	double mean_stddev = 0;

	for (int i=1; i<=N; i++) {
		double sum = 0;
		for (int k=0; k<S; k++)	sum += M_PI*0.5*cos(M_PI*0.5*rnd.Rannyu(0., 1.));	//integrand evaluation
		A = sum / (double) (S);
		sum_A += A;
		sum_A2 += A*A;
		mean = sum_A/(double)(i);
		mean_stddev = Block_StdDev(i, sum_A, sum_A2);
		output << setprecision(7) << i << "\t" << mean << "\t" << mean_stddev << endl;
	}
	output.close();

/*****************************************
2) Importance sampling for Integral evaluation
	
*****************************************/

	A = 0;	//block result
	sum_A = 0;	//cumulative average
	sum_A2 = 0; 	//cumulative square average
	mean = 0;
	mean_stddev = 0;
	output.open("IntegralIS.out");

	for (int i=1; i<=N; i++) {
		double sum = 0;
		for (int k=0; k<S; k++) {						//integrand evaluation with linear distribution of probability
			double j = rnd.RannIS();
			sum += M_PI*0.25/(1.-j)*cos(M_PI*0.5*j);
		}
		A = sum / (double) (S);
		sum_A += A;
		sum_A2 += A*A;
		mean = sum_A/(double)(i);
		mean_stddev = Block_StdDev(i, sum_A, sum_A2);
		output << setprecision(7) << i << "\t" << mean << "\t" << mean_stddev << endl;
	}

	output.close();
	rnd.SaveSeed();

return 0;
}



double Block_StdDev (int block, double cumsum, double cumsum2) {
	//block is number of block	
	//cumsum is cumulative sum
	//cumsum2 is squred cumulative sum
	if (block == 1) return 0;
	else return sqrt ((cumsum2/(double)(block) - (cumsum*cumsum)/(block*block))/(block-1.));
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
