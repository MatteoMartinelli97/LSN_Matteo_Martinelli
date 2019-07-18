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
 
int main (int argc, char *argv[]) {

	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
    	Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

/******************************************
1) mean value of a uniform distribution in [0,1]

+*****************************************/
	ofstream output;
	output.open("Average.out");
	int M = 1000000;	//steps
	int N = 1000;	//blocks
	int S = M/N;	//steps per block
	double A = 0;	//block result
	double sum_A = 0;	//cumulative average
	double sum_A2 = 0; 	//cumulative square average
	double mean = 0;
	double mean_stddev = 0;


	for (int i=0; i<N; i++) {								//cycle for Blocks
		double sum = 0;										//sum for mean in each block
		for (int j=0; j<S; j++) sum += rnd.Rannyu();		//cycle for each block
		A = sum / (double) (S);
		sum_A += A;
		sum_A2 += A*A;	
		mean = sum_A/(i+1.);
		output << fixed << setprecision(7) << (i+1)*S << "	" << mean -0.5<< "	";
		if (i==0) mean_stddev = 0;
		else mean_stddev = sqrt (((sum_A2/(i+1.)) - (mean*mean))/i);
		output << mean_stddev << endl;
	}

	output.close();

/*****************************************
2) fluctuatios of mean value of a uniform distribution in [0,1]

*****************************************/
	output.open("Sigma2.out");
	A = 0;	//block result
	sum_A = 0;	//cumulative average
	sum_A2 = 0; 	//cumulative square average
	mean = 0;
	mean_stddev = 0;


	for (int i=0; i<N; i++) {								//cycle for Blocks
		double sum = 0;										//sum for mean in each block
		for (int j=0; j<S; j++) {							//cycle for each block
			double k = rnd.Rannyu();		
			sum += (k-0.5)*(k-0.5);
		}
		A = sum / (double) (S);
		sum_A += A;
		sum_A2 += A*A;	
		mean = sum_A/(i+1.);
		output << fixed << setprecision(7) << (i+1)*S << "	" << mean -1./12.<< "	";
		if (i==0) mean_stddev = 0;
		else mean_stddev = sqrt (((sum_A2/(i+1.)) - (mean*mean))/i);
		output << mean_stddev << endl;
	}

	output.close();

/*****************************************/
//3) ChiSquared calculation

/*****************************************/

	M = 100;
	int counts [100] = {0};
	int n = 10000;
	output.open("ChiSquared.out");

	for (int l=0; l<100; l++) {
		for (int i=0; i<n; i++) {
			double k = rnd.Rannyu(0., 1.);
			int j = k/(1./double(M)); 						//interval in which is k
			counts[j]++;
		}
		double chi2 = 0;
		for (int i=0; i<M; i++) {
		chi2 += (double (M))*((double (counts[i]) - (double(n)/double(M))) * (double (counts[i]) - (double(n)/double(M))))/(double(n)); 
		counts[i] = 0;
		}
		output << l+1 << "	" << chi2 << endl;
	}
	output.close();
	output.open("ChiSquaredImp.out");
	M = 10000;
	n = 1000000;
	int counts2 [10000] = {0};
	for (int l=0; l<100; l++) {
		for (int i=0; i<n; i++) {
			double k = rnd.Rannyu(0., 1.);
			int j = k/(1./double(M)); 						//interval in which is k
			counts2[j]++;
		}
		double chi2 = 0;
		for (int i=0; i<M; i++) {
		chi2 += (double (M))*((double (counts2[i]) - (double(n)/double(M))) * (double (counts2[i]) - (double(n)/double(M))))/(double(n)); 
		counts2[i] = 0;
		}
		output << l+1 << "	" << chi2 << endl;
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
