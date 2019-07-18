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
#include "histo.h"

using namespace std;
 
int main (int argc, char *argv[]){

	system("rm *histo.out");	//Command to delete de file if existing

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

/*****************************************/


/*****************************************/
	ofstream output;
	ofstream output2;
	ofstream output3;
	output.open("Standard_dice.out");
	output2.open("Exponential_dice.out");
	output3.open("Lorentzian_dice.out");
	int M = 10000;
	int N = 1;
	double std_mean = 0;
	double exp_mean = 0;
	double lor_mean = 0;
	
	for (N=1; N<=100; N++) {
		if (N!=1 && N!=2 && N!=10 && N!=100) continue;		//just some values for N
		for (int i=0; i<M; i++) {
			double std_sum = 0;		//cumulative sums
			double exp_sum = 0;
			double lor_sum = 0;
			for (int j=0; j<N; j++) {
				double std = rnd.Rannyu ();			//standard dice
				double exp = rnd.Exp (1.);			//exponential dice
				double lor = rnd.Lorentz (0., 1.); //lorentzian dice		
				std_sum += std;
				exp_sum += exp;
				lor_sum += lor;
			}
			std_mean = std_sum/(double) (N);
			exp_mean = exp_sum/(double) (N);
			lor_mean = lor_sum/(double) (N);
			output << std_mean << endl;
			output2 << exp_mean << endl;
			output3 << lor_mean << endl;
		}
	}
	
	output.close();
	output2.close();
	output3.close();

	output.open("Steps.out");
//N = 1 Histogram
	Histo* histo = new Histo (50, 0, 10000, "Standard_dice.out");
	histo -> Normalize();
	histo-> Print_norm ("N1_histo.out");
	output << "N=1\t" << histo->Step() << "\t";
	delete histo;
	histo = new Histo (50, 0, 10000, "Exponential_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N1_histo.out");
	output << histo->Step() << "\t";
	delete histo;
	histo = new Histo (500, 0, 10000,-20., 20., "Lorentzian_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N1_histo.out");
	output << histo->Step() << endl;
	delete histo;

//N = 2 Histogram
	histo = new Histo (50, 10000, 10000, "Standard_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N2_histo.out");
	output << "N=2\t" << histo->Step() << "\t";
	delete histo;
	histo = new Histo (50, 10000, 10000, "Exponential_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N2_histo.out");
	output << histo->Step() << "\t";
	delete histo;
	histo = new Histo (500, 10000, 10000,-20., 20., "Lorentzian_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N2_histo.out");
	output << histo->Step() << endl;
	delete histo;

//N = 10 Histogram
	histo = new Histo (50, 20000, 10000, "Standard_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N10_histo.out");
	output << "N=10\t" << histo->Step() << "\t";
	delete histo;
	histo = new Histo (50, 20000, 10000, "Exponential_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N10_histo.out");
	output << histo->Step() << "\t";
	delete histo;
	histo = new Histo (500, 20000, 10000,-20., 20., "Lorentzian_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N10_histo.out");
	output << histo->Step() << endl;
	delete histo;

//N = 100 Histogram
	histo = new Histo (50, 30000, 10000, "Standard_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N100_histo.out");
	output << "N=100\t" << histo->Step() << "\t";
	delete histo;
	histo = new Histo (50, 30000, 10000, "Exponential_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N100_histo.out");
	output << histo->Step() << "\t";
	delete histo;
	histo = new Histo (500, 30000, 10000,-20., 20., "Lorentzian_dice.out");
	histo -> Normalize();
	histo->Print_norm ("N100_histo.out");
	output << histo->Step() << endl;
	delete histo;

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
