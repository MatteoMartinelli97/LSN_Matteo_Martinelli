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

/*****************************************
a) Random Walk on a 3D lattice

*****************************************/
	ofstream output;
	output.open("3D_latticeRW.out");
	
	int M = 10000; 	//Simulaiion
	int N = 100;	//Blocks
	int S = 100;	//Random Walk steps
	int L = M/N;	//Steps per block
	const double a = 1.; //step lenght
	int dir = 0; //direction
	double sum_err [S] = {0};		//cumulative for blocks mean
	double sum_mean_r2 [S] = {0};	//cumulative for blocks error
	
		for (int i=1; i<=N; i++) {			//block cycle
			double sum_r2 [S] = {0};
			for (int j=1; j<=L; j++) {			//steps per block
				int position [3] = {0};
				for (int s=0; s<S; s++) {			//one cycle of random walk (100 steps)
				//step on the lattice
					double k = rnd.Rannyu(0., 3.);
					if (k < 1.) dir = 0;		 //x
					else if (k >= 2.) dir = 2;	 //z
					else dir = 1;				 //y
					k = rnd.Rannyu(0., 2.);
					if (k < 1.) position [dir] -= a;	//backward
					else position [dir] += a;			//forward
				//calculation of r^2 
					double r2 = 0;					//for every step "s" need to put r2=0 and calculate again
					for (int p=0; p<3; p++) r2 += position[p]*position[p];
					sum_r2 [s] += r2;				//sum r2 for the L steps in a block (vector in function of step "s" of RW)

					if (j==L) {						//in the last step of a block calculate the mean value (=result of this block)
						sum_mean_r2 [s] += sum_r2 [s]/(double) (L);
						sum_err [s] += (sum_r2 [s] * sum_r2 [s])/(double)(L*L);	//count for the error (= sum_mean_r2^2)
						if (i==N) {											//in the last block(=finished simulation) calculate
							double r_mean = sqrt(sum_mean_r2[s] /(double) (N));	//the mean value of blocks and error for the last block
							double r2_stddev = sqrt (((sum_err [s]/(double)(N)) - (r_mean *r_mean )*(r_mean *r_mean))/(N-1.));	
							double r_mean_stddev = 0.5*r2_stddev/(r_mean);
							output << s+1 << "	" << r_mean << "	" << r_mean_stddev << endl;
						//reset of vectors for continuum random walk							
							sum_err [s] = 0;
							sum_mean_r2 [s] = 0;
						}
					}
				}
			}
		}

	output.close();

/*****************************************
a) Random Walk on the continuum

*****************************************/
	output.open("3D_continuumRW.out");
	
	M = 10000; 	//Simulaiion
	N = 100;	//Blocks
	S = 100;	//Random Walk steps
	L = M/N;	//Steps per block
	
	for (int i=1; i<=N; i++) {				//block cycle
			double sum_r2 [S] = {0};
			for (int j=1; j<=L; j++) {			//steps per block
				double position [3] = {0};
				for (int s=0; s<S; s++) {			//one cycle of random walk (100 steps)
				//step on the lattice
					double theta = rnd.Solid_Angle_Theta();
					double phi = rnd.Solid_Angle_Phi();
					position [0] += a*sin(theta)*cos(phi);	 //x
					position [1] += a*sin(theta)*sin(phi);	 //y
					position [2] += a*cos(theta);			 //z
				//calculation of r^2 
					double r2 = 0;					//for every step "s" need to put r2=0 and calculate again
					for (int p=0; p<3; p++) r2 += position[p]*position[p];
					sum_r2 [s] += r2;				//sum r2 for the L steps in a block (vector in function of step "s" of RW)

					if (j==L) {					//in the last step of a block calculate the mean value (=result of this block)
						sum_mean_r2 [s] += sum_r2 [s]/(double) (L);
						sum_err [s] += (sum_r2 [s] * sum_r2 [s])/(L*L);	//count for the error
						if (i==N) {											//in the last block(=finished simulation) calculate
							double r_mean = sqrt(sum_mean_r2[s] /(double) (N));	//the mean value of blocks and error for the last block
							double r2_stddev = sqrt (((sum_err [s]/(double)(N)) - (r_mean *r_mean )*(r_mean *r_mean))/(N-1.));	
							double r_mean_stddev = 0.5*r2_stddev/(r_mean);
							output << s+1 << "	" << r_mean << "	" << r_mean_stddev << endl;
						}
					}
				}
			}
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
