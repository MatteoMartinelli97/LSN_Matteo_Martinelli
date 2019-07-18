#include "random.h"
#include <cstring>

using namespace std;

//random generator
Random rnd;

//M(RT)2 properties
string state;
double delta;
int pdf;
double x, y, z;
double r;

//cumulative and mean for data blocking
double r_sum;
double r_mean, r_std;
double glob_av, glob_av2;

//simulation
int nsteps;
int nblock;
int bsteps;

//equilibration

const int eqsteps = 1000;
int acount;
double x_start, y_start, z_start;

//functions
void Input(void);
void Move(void);
bool Accept (double A);
void PrintPosition (void);
void BlockMeasure (int block);
void Equilibrate();

