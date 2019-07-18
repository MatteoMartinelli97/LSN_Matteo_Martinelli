/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <vector>
#include <algorithm>
#include <iostream>
#include "path.h"
using namespace std;

//Default constructor
Path::Path() { }

//Constructor of random path of size = size
Path::Path (int size) {
//Vector of size different cities (in order from 0 to size-1)
	_path.resize(size);
	for (int i=0; i<size; ++i) _path[i] = i;
//Shuffling path
	random_shuffle(_path.begin(), _path.end());
	_size = size;
	_length = 0;
}

Path::~Path() { }


double Path::GetLength (double** distances) {
	double length = 0;
	for (int i=0; i<_size; ++i) length += distances[_path[Pbc(i)]][_path[Pbc(i+1)]];
	_length = length;
	return length;
}

bool Path::Check () {

	vector <int> sorted_path (_path);					//Copy the vector and sort 
	sort(sorted_path.begin(), sorted_path.end());		//two equal elements will be adjacent now
	for (int i=0; i<_size-1; ++i) {
		if (sorted_path[i] == sorted_path[i+1]) {
			cerr << "Error! Problem with genetic move... arrested program" << endl;
			Print();
			return true;
		}
	}
	return false;
}

void Path::Print () {
	for (int i=0; i<_size; ++i) cout << _path[i] << "\t";
	cout << endl;
}


int Path::Pbc (int idx) {
	if (idx < 0) Pbc(idx + _size);
	return idx%_size;
}

void Path::Swap (int idx_1, int idx_2) {
	swap(_path[Pbc(idx_1)], _path[Pbc(idx_2)]);
	return;
}

void Path::Swap (int idx_1, int idx_2, int n_cts) {	//swap of n_cts from idx_1 with n_cts from idx_2
	if (n_cts > _size*0.5) n_cts = _size*0.5;
	idx_2 = Pbc(idx_1+n_cts+idx_2);	//Adjust 2nd index in pbc, because it's exctracted in a way that can permit overlaps
	for (int i=0; i<n_cts; ++i) swap(_path[Pbc(idx_1+i)], _path[Pbc(idx_2+i)]);	
	return;
}

void Path::Shift (int nshift) {	//shift of nshift position
	if (nshift > _size) nshift = nshift%_size;
	rotate (_path.rbegin(), _path.rbegin() +nshift, _path.rend());
	return;
}

void Path::Shift (int nshift, int n_cts, int idx_1) {
	idx_1 = Pbc(idx_1);
//	if (nshift < n_cts) nshift = n_cts+1;
	int Start_cts [n_cts];
	for (int i=0; i<n_cts; ++i) Start_cts[i] = _path[Pbc(idx_1 + i)];
	for (int i=idx_1; i<nshift+idx_1; ++i) _path[Pbc(i)] = _path[Pbc(i + n_cts)];
	for (int i=0; i<n_cts; ++i) _path[Pbc(idx_1 + nshift +i)] = Start_cts[i];
	return;
}
		
void Path::Reverse (int idx_1, int n_cts) {	
	if (n_cts > _size) n_cts = _size;
	int swaps = n_cts*0.5; //number of swap to do (reverse = first to last, last to first ecc)
	for (int i=0; i<swaps; ++i) Swap(idx_1+i, idx_1+n_cts-1-i);
	return;
}

void Path::WriteGene (int gene, int idx) {
	_path[idx] = gene;
	return;
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


