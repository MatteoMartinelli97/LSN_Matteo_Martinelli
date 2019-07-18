/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#ifndef __PATH__
#define __PATH__

#include <vector>
#include <iterator> 
using namespace std;

class Path {

private:
	vector <int> _path;
	int _size;
	double _length;
	int Pbc (int idx);



public:
	Path ();
	Path(int size);		//Constructs a random path of "size" cities
	~Path();
	void Initialize(int size);
	double GetLength (double** distances);	//Returns length of path

	bool Check ();	//Check that all values in the path are different
					//this guarantees no error in genetic moves
					//if two equals city are found returns true
	void Print();	//Print the path
	void WriteGene (int gene, int idx);
//Genetic mutation functions
	void Swap (int idx_1, int idx_2);
	void Swap (int idx_1, int idx_2, int n_cts);
	void Shift (int nshift);
	void Shift (int nshift, int n_cts, int idx_1);
	void Reverse (int idx_1, int n_cts);
	vector <int> GetPath() {return _path;}
	vector <int>::iterator GetStart() {return _path.begin();}
	vector <int>::iterator GetEnd () {return _path.end();}
};

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
