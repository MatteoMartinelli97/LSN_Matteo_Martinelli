/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __Genetic_Alg__
#define __Genetic_Alg__

#include <vector>
#include <algorithm>
#include <cstring> 
#include <string>
#include "random.h"
#include "path.h"
#include <chrono>

using namespace std;
using namespace std::chrono;

//Random generator
Random rnd;


int ncts;
Path Travel_Path; 	//path
Path Best_Path;
double best_l;

//simulation parameters
const double t_incr = 0.8;
const int max_steps = 10000;
int s_incr = 100;
int nstep, nblk;
double temp;
double beta;
double accepted[4];
double attempted [4];

//geometry parameters
int geometry;
string geo_name;
double** dist;
double* x;
double* y;

//Check parameters
bool Check_val;


void Input();
void Reset();
void Move(int imut);
bool Accept (double A);
void Print_Best_Path_Coord ();
void Print (int iblk);
void New_Params();
#endif
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
