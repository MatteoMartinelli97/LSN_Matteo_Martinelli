#include "histo.h"
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

Histo::Histo (int Nbin, const char* filename) {
	N_bin = Nbin;
	bin = new double [N_bin];
	counts = new int [N_bin];
	for (int i=0; i<N_bin; i++) counts[i]=0;
	_max = 0;
	_min = 0;
	N_values = 0;
	ifstream input (filename);
	if (input.is_open()){
		while ( !input.eof() ){
			double value = 0;
			input >> value;
			if (_max < value) _max = value;
			if (_min > value or N_values==0) _min = value;
			N_values++;
		}
	} else cerr << "PROBLEM: Unable to open" << filename << endl;
	input.close();
	_step = (_max-_min)/ (double) (N_bin);
	for (int i=0; i<N_bin; i++) bin[i] = _min+i*_step;
	input.open (filename);
	if (input.is_open()){
		while ( !input.eof() ){
			double value = 0;
			input >> value;
			int j = (value - _min)/_step;
			counts[j]++;
		}
	} else cerr << "PROBLEM: Unable to open" << filename << endl;
}

Histo::Histo (int Nbin, int skip, int size, double* v) {
	N_bin = Nbin;
	bin = new double [N_bin];
	counts = new int [N_bin];
	for (int i=0; i<N_bin; i++) counts[i]=0;
	N_values = size;
	_max = 0;
	_min = 0;
	for (int i=skip; i<size; i++) {
		if (_max < v[i]) _max = v[i];
		if (_min > v[i] or i==0) _min = v[i];
	}	
	_step = (_max-_min)/ (double) (N_bin);
	for (int i=0; i<N_bin; i++) bin[i] = _min+i*_step;
	for (int i=skip; i<size; i++) {
		int j = (v[i] - _min)/_step;
		counts[j]++;
	}		
}


Histo::Histo (int Nbin, int skip, int size, const char* filename) {
	N_bin = Nbin;
	bin = new double [N_bin];
	counts = new int [N_bin];
	for (int i=0; i<N_bin; i++) counts[i]=0;
	_max = 0;
	_min = 0;
	N_values = 0;
	ifstream input (filename);
	if (input.is_open()){
		for (int i=0; i<skip; i++) {
			double value = 0;
			input >> value;
		}
		while (N_values<size && !input.eof()) {
			double value = 0;
			input >> value;
			if (_max < value) _max = value;
			if (_min > value or N_values==0) _min = value;
			N_values++;
		}
	} else cerr << "PROBLEM: Unable to open" << filename << endl;
	input.close();

	_step = (_max-_min)/ (double) (N_bin);
	for (int i=0; i<N_bin; i++) bin[i] = _min+i*_step;
	input.open (filename);
	if (input.is_open()){
		for (int i=0; i<N_values+skip; i++) {
			double value = 0;			
			input >> value;
			if (i<skip) continue;
			int j = 0;
			if (value == _max) j = N_bin-1;
			else j = (value - _min)/_step;
			counts[j]++;
		}
	} else cerr << "PROBLEM: Unable to open" << filename << endl;
}

Histo::Histo (int Nbin, int skip, int size, double range_min, double range_max, double* v) {
	N_bin = Nbin;
	bin = new double [N_bin];
	counts = new int [N_bin];
	for (int i=0; i<N_bin; i++) counts[i]=0;
	N_values = size;
	_max = 0;
	_min = 0;
	for (int i=skip; i<size; i++) {
		if (v[i] < range_min or v[i] > range_max) continue;
		if (_max < v[i]) _max = v[i];
		if (_min > v[i] or i==0) _min = v[i];
	}	
	_step = (_max-_min)/ (double) (N_bin);
	for (int i=0; i<N_bin; i++) bin[i] = _min+i*_step;
	for (int i=skip; i<size; i++) {
		if (v[i] < range_min or v[i] > range_max) continue;
		int j = (v[i] - _min)/_step;
		counts[j]++;
	}		
}

Histo::Histo (int Nbin, int skip, int size, double range_min, double range_max, const char* filename) {
	N_bin = Nbin;
	bin = new double [N_bin];
	counts = new int [N_bin];
	for (int i=0; i<N_bin; i++) counts[i]=0;
	_max = 0;
	_min = 0;
	N_values = 0;
	ifstream input (filename);
	if (input.is_open()){
		for (int i=0; i<skip; i++) {
			double value = 0;
			input >> value;
		}
		while (N_values<size && !input.eof()) {
			double value = 0;
			input >> value;
			if (value < range_min or value > range_max) continue;
			if (_max < value) _max = value;
			if (_min > value or N_values==0) _min = value;
			N_values++;
		}
	} else cerr << "PROBLEM: Unable to open" << filename << endl;
	input.close();

	_step = (_max-_min)/ (double) (N_bin);
	for (int i=0; i<N_bin; i++) bin[i] = _min+i*_step;
	input.open (filename);
	if (input.is_open()){
		for (int i=0; i<size+skip; i++) {
			double value = 0;
			input >> value;
			if (i<skip) continue;
			if (value < range_min or value > range_max) continue;
			int j = 0;
			if (value == _max)j = N_bin-1;
			else j = (value - _min)/_step;
			counts[j]++;
		}
	} else cerr << "PROBLEM: Unable to open" << filename << endl;
}

Histo::Histo (int Nbin, int skip, double range_min, double range_max, const char* filename) {
	N_bin = Nbin;
	bin = new double [N_bin];
	counts = new int [N_bin];
	for (int i=0; i<N_bin; i++) counts[i]=0;
	_max = 0;
	_min = 0;
	N_values = 0;
	ifstream input (filename);
	if (input.is_open()){
		for (int i=0; i<skip; i++) {
			double value = 0;
			input >> value;
		}
		while ( !input.eof() ){
			double value = 0;
			input >> value;
			if (value < range_min or value > range_max) continue;
			if (_max < value) _max = value;
			if (_min > value or N_values==0) _min = value;
			N_values++;
		}
	} else cerr << "PROBLEM: Unable to open" << filename << endl;
	input.close();
	_step = (_max-_min)/ (double) (N_bin);
	for (int i=0; i<N_bin; i++) bin[i] = _min+i*_step;
	input.open (filename);
	if (input.is_open()){
		for (int i=0; i<skip; i++) {
			double value = 0;
			input >> value;
		}
		while ( !input.eof() ){
			double value = 0;
			input >> value;
			if (value < range_min or value > range_max) continue;
			int j = (value - _min)/_step;
			counts[j]++;
		}
	} else cerr << "PROBLEM: Unable to open" << filename << endl;
}
	

Histo::~Histo() {
	delete counts;
	delete bin;
}

void Histo::SetPoint(double value) {
	if (value > _max or value < _min) cerr << "PROBLEM: invalid value" << endl;
	else {
		int j = (value - _min)/_step;
		counts[j]++;	
		N_values++;
	}
return;
}

void Histo::Print (const char* filename) {
	ofstream output (filename, ios_base::out | ios_base::app);
	for (int i=0; i<N_bin; i++) output << fixed << setprecision(8) << bin[i] << "\t" << counts[i] << endl;
	output.close();
}

void Histo::Print_norm (const char* filename) {
	ofstream output (filename, ios_base::out | ios_base::app);
	for (int i=0; i<N_bin; i++) output << fixed << setprecision(8) << bin[i] << "\t" << _norm_counts[i] << endl;
	output.close();
}


void Histo::Normalize () {
	_norm_counts = new double [N_bin];
	for (int i=0; i<N_bin; ++i) _norm_counts[i] = (double)(counts[i]) / (_step*(double)(N_values));
}
