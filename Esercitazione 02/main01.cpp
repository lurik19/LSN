#include <iostream>
#include <fstream>
#include <cmath>
#include "Funzioni.h"

using namespace std;

int main (int argc, char *argv[]){
	
	// ************************* PARTE 1 *************************
	
	int M = 1E6; // Numero di tiri
	int N = 100; // Numero di blocchi
	int L = int(M/N); // Numero di esperimenti (o tiri) per blocco

	double sum; // Variabile d'appoggio per calcolare il risultato dell'integrazione, I

	// Inizializziamo una variabile di tipo Random
	Random rnd;
	rnd.SetRandom();

	// Variabili che servono per applicare il blocking method
	double *ave = new double[N];
	double *av2 = new double[N];
	double *I_prog = new double[N];
	double *I2_prog = new double[N];
	double *err_prog = new double[N];

	// Apertura dei file necessari per scrivere i dati della parte 1
	ofstream file1("./Data/Es02.1/I_es2.1.a.txt");

	// *********************** SVOLGIMENTO ***********************

	// Stimo I all'interno di ogni blocco
	for(int i=0; i<N; i++){		
		sum = 0;
		for(int j=0; j<L; j++){
			sum += M_PI/2*cos(M_PI/2*rnd.Rannyu());
		}
		sum /= L;
		ave[i] = sum;
		av2[i] = pow(sum,2);
	}

	// Ora ne calcoliamo la media (su i blocchi) e l'incertezza statistica
	for(int i=0; i<N; i++){
		I_prog[i] = 0;
		I2_prog[i] = 0;
		for(int j=0; j<i+1; j++){
			I_prog[i] += ave[j];
			I2_prog[i] += av2[j];
		}
		I_prog[i]/=(i+1);
		I2_prog[i]/=(i+1);
		err_prog[i] = error(I_prog, I2_prog, i);

		// Scriviamo su file i risultati
		file1 << I_prog[i] << "   " << err_prog[i] << endl;
	}


	// ************************* PARTE 2 *************************

	double r; // Variabile che utilizzeremo per generare
		  	  // numeri casuali che seguano la PDF desiderata
		  	  // ovvero una PDF ottimale per il calcolo di I

	// Apertura dei file necessari per scrivere i dati della parte 2
	ofstream file2("./Data/Es02.1/I_es2.1.b.txt");
	
	// *********************** SVOLGIMENTO ***********************

	// Stimo I all'interno di ogni blocco
	for(int i=0; i<N; i++){		
		sum = 0;
		for(int j=0; j<L; j++){
			r = 2./M_PI*asin(rnd.Rannyu());
			sum += 1*M_PI/2*cos(M_PI/2*r)/(M_PI/2*cos(M_PI/2*r));
		}
		sum /= L;
		ave[i] = sum;
		av2[i] = pow(sum,2);
	}

	// Ora ne calcoliamo la media (su i blocchi) e l'incertezza statistica
	for(int i=0; i<N; i++){
		I_prog[i] = 0;
		I2_prog[i] = 0;
		for(int j=0; j<i+1; j++){
			I_prog[i] += ave[j];
			I2_prog[i] += av2[j];
		}
		I_prog[i]/=(i+1);
		I2_prog[i]/=(i+1);
		err_prog[i] = error(I_prog, I2_prog, i);

		// Scriviamo su file i risultati
		file2 << I_prog[i] << "   " << err_prog[i] << endl;
	}

	// Salvo il file seed.out

	rnd.SaveSeed();

	// Elimino i vettori allocati precendentemente	

	delete[] ave;
	delete[] av2;
	delete[] I_prog;
	delete[] I2_prog;
	delete[] err_prog;
	
	// Chiudo i file aperti precedentemente
	
	file1.close();
	file2.close();

	return 0;
}
