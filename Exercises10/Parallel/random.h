/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __Random__
#define __Random__

class Random {

private:
	int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
	// constructors
	Random();
	// destructor
	~Random();
	// methods
	void Initialize();
	void Par_Initialize (int rank);
	void SetRandom(int * , int, int);
	void SaveSeed();
	double Rannyu(void);
	double Rannyu(double min, double max);
	double Gauss(double mean, double sigma);
	double Exp (double lambda);
	double Lorentz (double mean, double gamma);
	double Angle ();		//Uniform angle between [0, pi]
	double Solid_Angle_Theta ();
	double Solid_Angle_Phi();
};

#endif // __Random__

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
