#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Funzioni.h"

using namespace std;

int main (int argc, char *argv[]){

	// ************************* PARTE 1 *************************


	// Lo svolgimento è stato fatto nei file random.h e random.cpp


	// ************************* PARTE 2 *************************
	
	int N[4] = {1,2,10,100};
	int M = 1E6; // Numero di realizzazioni dei diversi S_N

	double *S_N_unif = new double[M*4];    // Vettori che raccolgono S1, S2, S10 e S100
	double *S_N_exp = new double[M*4];     // per tutte le M realizzazioni per tutte e
	double *S_N_lorentz = new double[M*4]; // tre le distribuzioni di probabilità

	// Inizializziamo una variabile di tipo Random
	Random rnd;
	rnd.SetRandom();

	// Apertura dei file necessari per scrivere i dati
	ofstream file1("./Data/Es01.2/S_N_unif.txt");
	ofstream file2("./Data/Es01.2/S_N_exp.txt");
	ofstream file3("./Data/Es01.2/S_N_lorentz.txt");
	int wd = 24;

	// *********************** SVOLGIMENTO ***********************

	for(int i=0; i<M; i++){
		for(int j=0; j<4; j++){
			// Calcoliamo gli S_N per i diversi N e per le diverse distribuzioni
			S_N_unif[j + 4*i] = 0;
			S_N_exp[j + 4*i] = 0;
			S_N_lorentz[j + 4*i] = 0;	
			for(int k=0; k<N[j]; k++){
				S_N_unif[j + 4*i] += rnd.Rannyu();
				S_N_exp[j + 4*i] += rnd.Exp(1.);
				S_N_lorentz[j + 4*i] += rnd.Lorentz(0., 1.);
			}
			S_N_unif[j + 4*i]/=N[j];
			S_N_exp[j + 4*i]/=N[j];
			S_N_lorentz[j + 4*i]/=N[j];
			
			// Scriviamo su file i risultati
			file1 << S_N_unif[j + 4*i] << setw(wd);
			file2 << S_N_exp[j + 4*i] << setw(wd);
			file3 << S_N_lorentz[j + 4*i] << setw(wd);
			
			// Sui file verranno scritti gli S_N in questo ordine:
			// S_1^(0), S_2^(0), S_10^(0), S_100^(0), S_1^(1), S_2^(1), ....
			// dove l'indice tra parentesi rappresenta la i-esima realizzazione
		}
		
		file1 << endl;
		file2 << endl;
		file3 << endl;
	}

	// Salvo il file seed.out

	rnd.SaveSeed();

	// Elimino i vettori allocati precendentemente

	delete[] S_N_unif;
	delete[] S_N_exp;
	delete[] S_N_lorentz;

	// Chiudo i file aperti precedentemente

	file1.close();
	file2.close();
	file3.close();

	return 0;
}
