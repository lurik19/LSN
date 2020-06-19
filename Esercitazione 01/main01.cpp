#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Funzioni.h"

using namespace std;

int main (int argc, char *argv[]){
	
	int M = 1E6; // Numero di tiri
	int N = 100; // Numero di blocchi
	int L = int(M/N); // Numero di esperimenti (o tiri) per blocco

	double sum; // Variabile d'appoggio per calcolare <r>

	// Inizializziamo una variabile di tipo Random
	Random rnd;
	rnd.SetRandom();

	// Variabili che servono per applicare il blocking method
	double *ave = new double[N];
	double *av2 = new double[N];
	double *sum_prog = new double[N];
	double *su2_prog = new double[N];
	double *err_prog = new double[N];


	// ************************* PARTE 1 *************************

	// Apertura dei file necessari per scrivere i dati della parte 1
	ofstream file1("./Data/Es01.1/sum_es1.1.a.txt");
	int wd = 24;

	// *********************** SVOLGIMENTO ***********************

	// Calcoliamo <r> per ogni blocco
	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<L; j++){
			sum += rnd.Rannyu();
		}
		ave[i] = sum/L;
		av2[i] = pow(ave[i],2);
	}

	// Ora ne calcoliamo la media (su i blocchi) e l'incertezza statistica
	for(int i=0; i<N; i++){
		sum_prog[i] = 0;
		su2_prog[i] = 0;
		for(int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}
		sum_prog[i]/=(i+1);
		su2_prog[i]/=(i+1);
		err_prog[i] = error(sum_prog, su2_prog, i);

		// Scriviamo su file i risultati
		file1 << sum_prog[i] << setw(wd) << err_prog[i] << '\n';
	}

	// ************************* PARTE 2 *************************

	// Apertura dei file necessari per scrivere i dati della parte 2
	ofstream file2("./Data/Es01.1/sum_es1.1.b.txt");

	// *********************** SVOLGIMENTO ***********************

	// Calcoliamo <(r-<r>)^2> per ogni blocco
	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<L; j++){
			sum += pow(rnd.Rannyu()-0.5,2);
		}
		ave[i] = sum/L;
		av2[i] = pow(ave[i],2);
	}

	// Ora ne calcoliamo la media (su i blocchi) e l'incertezza statistica
	for(int i=0; i<N; i++){
		sum_prog[i] = 0;
		su2_prog[i] = 0;
		for(int j=0; j<i+1; j++){
			sum_prog[i] += ave[j];
			su2_prog[i] += av2[j];
		}
		sum_prog[i]/=(i+1);
		su2_prog[i]/=(i+1);
		err_prog[i] = error(sum_prog, su2_prog, i);

		// Scriviamo su file i risultati
		file2 << sum_prog[i] << setw(wd) << err_prog[i] << '\n';
	}

	// ************************* PARTE 3 *************************
	
	M = 100; // Numero di sottointervalli in cui dividiamo [0,1]
	int n = 1E6; // Numero di numeri pseudo-random generati
	int n_throw = 1E4; // Numero di numeri pseudo-random usati per calcolare il chi^2
			   // UNA volta. Scegli n_throw affinchè n/n_throw sia intero

	double rand; // Variabile d'appoggio per salvare numeri random

	int j = 0;
	int *count = new int[M]; // Vettore che conta il numero di eventi in ogni sottointervallo
	double *chi_quadro = new double[n/n_throw];

	// Apertura del file necessario per scrivere i dati della parte 3
	ofstream file3("./Data/Es01.1/es1.1.c.txt");

	// *********************** SVOLGIMENTO ***********************	

	for(int counter=0; counter<n/n_throw; counter++){
		for(int i=0; i<M; i++){ // Pongo a zero il conteggio del numero
			count[i] = 0;   // di eventi in ogni sottointervallo
		}

		for(int i=counter*n_throw; i<counter*n_throw + n_throw; i++){
			j = 0;
			rand = rnd.Rannyu();
			while(rand < double(j)/M || rand >= double(j+1)/M){ // in questo modo trovo l'indice
				j++;					    					// del sottointervallo in cui è 
			}						    						// caduto il mio numero random
			count[j] += 1;
		}

		// Calcolo chi^2
		chi_quadro[counter] = 0;
		for(int i=0; i<M; i++){
			chi_quadro[counter] += pow((count[i]-n_throw/M),2)/(n_throw/M);	
		}
		file3 << chi_quadro[counter] << endl;
	}

	// Salvo il file seed.out

	rnd.SaveSeed();

	// Elimino i vettori allocati precendentemente	

	delete[] ave;
	delete[] av2;
	delete[] sum_prog;
	delete[] su2_prog;
	delete[] err_prog;
	delete[] count;
	delete[] chi_quadro;	
	
	// Chiudo i file aperti precedentemente

	file1.close();
	file2.close();
	file3.close();

	return 0;
}
