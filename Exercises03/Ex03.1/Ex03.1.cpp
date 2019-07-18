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

/****************************************
1) Direct Sampling of final asset price for C = Call & P = Put options

****************************************/
	int M = 100000;
	int N = 100;
	int L = M/N;
	ofstream output;
	output.open("Direct_sampling.out");
	double S0 = 100.;		//asset price at t=0
	double T = 1.;			//delivery time
	double K = 100.;		//strike price
	double r = 0.1 ;		//risk-free interest rate
	double sigma = 0.25;	//volatility
	double C_final_price_sum = 0;
	double C_final_error_sum = 0;
	double P_final_price_sum = 0;
	double P_final_error_sum = 0;

	for (int i=1; i<=N; i++) {
		double C_price_sum = 0;
		double P_price_sum = 0;
		for (int j=1; j<=L; j++) {
			//one simulation
			double W = rnd.Gauss(0., sqrt(T));
			double ST = S0 * exp ((r-0.5*sigma*sigma)*T + sigma * W);
			double C_price = 0;
			double P_price = 0;			
			C_price = exp(-r*T)* max (0., (ST-K));
			P_price = exp(-r*T)* max (0., (K-ST));
			C_price_sum += C_price;
			P_price_sum += P_price;	
		}
		C_final_price_sum += C_price_sum/(double) (L);
		C_final_error_sum += (C_price_sum/(double) (L))*(C_price_sum/(double) (L));
		P_final_price_sum += P_price_sum/(double) (L);
		P_final_error_sum += (P_price_sum/(double) (L))*(P_price_sum/(double) (L));
		double C_final_price = C_final_price_sum/(double)(i);
		double C_final_error = Block_StdDev (i, C_final_price_sum, C_final_error_sum);
		double P_final_price = P_final_price_sum/(double)(i);
		double P_final_error = Block_StdDev (i, P_final_price_sum, P_final_error_sum);
		output << i+1 << setprecision(7) << fixed << "\t" << C_final_price << "\t" << C_final_error;
		output << "\t" << P_final_price << "\t" << P_final_error<< endl;
	}
	output.close();	

/****************************************
2)Sampling discretized final asset price for C = Call & P = Put options

****************************************/
	C_final_price_sum = 0;
	C_final_error_sum = 0;
	P_final_price_sum = 0;
	P_final_error_sum = 0;
	int S = 100;
	double dt = T/(double) (S);
	output.open ("Discretized_sampling.out");

	for (int i=0; i<N; i++) {
		double C_price_sum = 0;
		double P_price_sum = 0;
		for (int j=0; j<L; j++) {
			//one simulation
			double St = S0;
			for (int k=0; k<S; k++) {
				double W = rnd.Gauss(0., 1.);
				St = St*exp ((r-0.5*sigma*sigma)*dt + sigma*W*sqrt(dt));
			}
			double C_price = 0;
			double P_price = 0;			
			C_price = exp(-r*T)* max (0., (St-K));
			P_price = exp(-r*T)* max (0., (K-St));
			C_price_sum += C_price;
			P_price_sum += P_price;	
		}
		C_final_price_sum += C_price_sum/(double) (L);
		C_final_error_sum += (C_price_sum/(double) (L))*(C_price_sum/(double) (L));
		P_final_price_sum += P_price_sum/(double) (L);
		P_final_error_sum += (P_price_sum/(double) (L))*(P_price_sum/(double) (L));
		double C_final_price = C_final_price_sum/(double)(i+1.);
		double C_final_error = Block_StdDev (i+1., C_final_price_sum, C_final_error_sum);
		double P_final_price = P_final_price_sum/(double)(i+1.);
		double P_final_error = Block_StdDev (i+1., P_final_price_sum, P_final_error_sum);
		output << i+1 << setprecision(7) << fixed << "\t" << C_final_price << "\t" << C_final_error;
		output << "\t" << P_final_price << "\t" << P_final_error<< endl;
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
