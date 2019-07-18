#ifndef _HISTOGRAM_H
#define _HISTOGRAM_H

class Histo {

private:
	int N_bin;
	int N_values;
	int* counts;
	double* _norm_counts;
	double* bin;
	double _max, _min;
	double _step;
protected:

public:
	//constructors
	Histo (int Nbin, const char* filename);
	Histo (int Nbin, int skip, int size, double* v);
	Histo (int Nbin, int skip, int size, const char* filename);
	Histo (int Nbin, int skip, double range_min, double range_max, const char* filename);
	Histo (int Nbin, int skip, int size, double range_min, double range_max, double* v);
	Histo (int Nbin, int skip, int size, double range_min, double range_max, const char* filename);

	//destructor
	~Histo();
	//methods
	void SetPoint(double value);
	double Max_value () {return _max;};
	double Min_value () {return _min;};
	int* Bin_count () {return counts;};
	double* Bin_start () {return bin;};
	double Step() {return _step;};
	void Print (const char* filename);
	void Print_norm (const char* filename);
	void Normalize();
};

#endif
