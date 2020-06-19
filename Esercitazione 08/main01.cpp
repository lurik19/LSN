#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;

// Variabili
double x;
double step = 2.8; // con questo valore l'accettazione è tra il 40% e il 60%
double mu, sigma;
const int m_props = 1000;
int n_props = 2;
double varia_param;
int bool_varia;
double vec_varia[2];

int ie = 0; // indice energia (valore medio Hamiltoniana)

Random rnd;

// Data blocking
double walker[m_props];
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_ene, err_ene;


// funzioni
void Reset(int);
void Accumulate(void);
void Averages(int);

void Init();
void Metropolis();
double p(double);
double psi(double);
double gauss(double, double, double);
void Measure();
double Error(double,double,int);


int nblk; // Numero di blocchi
int nstep; // Numero di o tiri per blocco

double temp, beta, delta;
int mu_vs_sigma;

int main (int argc, char *argv[]){
	
	Init();

	double E;

	ofstream Energia;
	Energia.open("output.ene_mu_sigma.0", ios::app);

	for(int i=0; i < 1 + 79*bool_varia; i++){		

		sigma += vec_varia[0];
		mu += vec_varia[1];

		for(int iblk=1; iblk <= nblk; ++iblk){

			Reset(iblk);
			for(int istep=1; istep <= nstep; ++istep){
				Metropolis();
				Measure();
				Accumulate();
			}

			Averages(iblk);

		}

		E = glob_av[ie]/(double)nblk;
		
		cout << "Passo " << i+1 << endl;
		cout << "mu = " << mu << "  " << "sigma = " << sigma << "  " << "E = " << E << endl;

		Energia << mu << "  " << sigma << "  " << E << endl;

	}

	Energia.close();

	cout << "Accettazione media = " << accepted/attempted << endl << endl;


	return 0;
}


void Init(){

	// Inizializziamo una variabile di tipo Random
	rnd.SetRandom();

	x = rnd.Rannyu(-step, step);

	ifstream ReadInput;
	ReadInput.open("input.dat");

	ReadInput >> varia_param;
	ReadInput >> mu;
	ReadInput >> sigma;
	ReadInput >> delta;
	
	// Con varia_param decidiamo se variare mu o sigma o nessuno dei due
	if(varia_param == 0){
		cout << "Il programma funziona con mu e sigma fissati\n";
		cout << "sigma = " << sigma << endl;
		cout << "mu    = " << mu << endl;

		bool_varia = 0;
		vec_varia[0] = 0;
		vec_varia[1] = 0;
	} else {
		if(varia_param == 1) {
			cout << "Il programma varia sigma\n";
			cout << "sigma iniziale  = " << sigma << endl;
			cout << "mu	         = " << mu << endl;
			cout << "Lunghezza passi = " << delta << endl;
			bool_varia = 1;
			vec_varia[0] = delta;
			vec_varia[1] = 0;
		} else {
			cout << "Il programma varia mu\n";
			cout << "mu iniziale  = " << mu << endl;
			cout << "sigma        = " << sigma << endl;
			cout << "Lunghezza passi = " << delta << endl;

			bool_varia = 1;
			vec_varia[0] = 0;
			vec_varia[1] = delta;
		}
	}

	ReadInput >> nblk;
	ReadInput >> nstep;

	cout << "Numero di blocchi         = " << nblk << endl;
	cout << "Numero di tiri per blocco = " << nstep << endl;

	cout << endl;

	ReadInput.close();
}


void Reset(int iblk) { // Reset block averages
   
	if(iblk == 1) {
		for(int i=0; i<n_props; ++i) {
			glob_av[i] = 0;
           		glob_av2[i] = 0;
       		}
   	}

	for(int i=0; i<n_props; ++i) {
		blk_av[i] = 0;
	}

	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}

void Accumulate(void) { // Update block averages

	for(int i=0; i<n_props; ++i){
		blk_av[i] = blk_av[i] + walker[i];
	}

	blk_norm = blk_norm + 1.0;

}

void Averages(int iblk) { // Print results for current block

	ofstream Ene, Ene_;
	const int wd=24;

	//cout << "Block number " << iblk << endl;
	//cout << "Acceptance rate " << accepted/attempted << endl << endl;

	Ene.open("output.ene.0", ios::app);

	// Valore medio dell'Hamiltoniana
	stima_ene = blk_av[ie]/double(blk_norm);
	glob_av[ie] += stima_ene;
	glob_av2[ie] += stima_ene * stima_ene;

	err_ene = Error(glob_av[ie], glob_av2[ie], iblk);

	Ene << setw(wd) << iblk <<  setw(wd) << stima_ene << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_ene << endl;
	
	Ene.close();

	//cout << "----------------------------" << endl << endl;


}


void Metropolis(){
	double X = rnd.Rannyu(x - step, x + step);
	double alpha = min(1., p(X)/p(x));
	
	if (rnd.Rannyu() < alpha) {x = X; accepted += 1;}

	attempted += 1;

	// Scrivo i punti in un file
	ofstream WritePoints;
	
	WritePoints.open("points.0", ios::app);
	WritePoints << x << endl;
	WritePoints.close();
}

double p(double x){
	return pow(psi(x), 2);
}

double psi(double x){
	return gauss(x, mu, sigma) + gauss(x, -mu, sigma);
}

double gauss(double x, double mu, double sigma){
	return exp(-(x-mu)*(x-mu)/(2*sigma*sigma));
}

void Measure(){

	double kin, pot;
	ofstream Ene;
	Ene.open("output.ene_ist.0", ios::app);
	
	kin = 1./(2*sigma*sigma)*(1.-1./psi(x)*((x-mu)*(x-mu)/(sigma*sigma)*gauss(x,mu,sigma)
					       +(x+mu)*(x+mu)/(sigma*sigma)*gauss(x,-mu,sigma)));

	pot = (x * x - 2.5) * x * x;

	walker[ie] = (kin + pot);

	Ene << walker[ie] << endl;
	Ene.close();

}

double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
}
