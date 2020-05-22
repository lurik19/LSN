// DOMANDA: SE GENERO UN VETTORE RANDOM DA 0 A 3 E LO DIVIDO IN 6 INTERVALLI MI VIENE LA GIUSTA LEGGE
//SE INVECE GENERO DUE VETTORI RANDOM: UNO CHE MI FA SCEGLIERE LA DIREZIONE E UNO BERNOULLI CHE MI FA SCEGLIERE SE ANDARE AVANTI O INDIETRO, ALLORA NON MI VIENE LA LEGGE GIUSTA! 
//VUOL DIRE CHE IL GENEREATORE DI NUMERI CREA CORRELAZIONI TRA QUESTI DUE VETTORI RANDOM?

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Funzioni.h"

using namespace std;

int main (int argc, char *argv[]){
	
	// ************************* PARTE 1 *************************
	
	int N = 100; // Numero di passi per singolo RW
	int L = 1E4; // Numero di simulazioni ( = numero di RW simulati)

	// A queste L simulazioni applichiamo il blocking method
	int n = 100; // Numero di blocchi
	int l = int(L/n); // Numero di RW simulati per blocco

	double a = 1.; // Distanza tra i punti del reticolo in tutte e tre le direzioni
	double rsquared;
	double p = 0.5; // Probabilità che il camminatore faccia un passo a dx (1-p a sx)

	int dir_index; // Variabile che useremo per scegliere
		       // in che direzione compiere il passo

	double pos[3] = {0.,0.,0.}; // Posizione del camminatore (inizia nell'origine)

	// Inizializziamo una variabile di tipo Random
	Random rnd;
	rnd.SetRandom();

	// Variabili che servono per applicare il blocking method
	double *ave = new double[N*n];
	double *av2 = new double[N*n];

	double *Rsquared_prog = new double[N];
	double *Rsquared2_prog = new double[N];
	double *err_prog = new double[N];

	// Apertura dei file necessari per scrivere i dati della parte 1
	ofstream file1("./Data/Es02.2/Rsquared_es2.2.a.txt");
	int wd = 24;

	// *********************** SVOLGIMENTO ***********************

	for(int i=0; i<N*n; i++){
		ave[i] = 0.;
		av2[i] = 0.;
	}

	// Facciamo scorrere l'indice k, che indica che siamo nel k-esimo blocco
	for(int k=0; k<n; k++){
	// In un blocco simulo l random walk da N passi l'uno
		for(int i=0; i<l; i++){
		// Fissato i, ho fissato il RW che sto considerando
			pos[0] = 0.;
			pos[1] = 0.;
			pos[2] = 0.;
			for(int j=0; j<N; j++){
				// Fissato j, ho fissato il numero di passi che considero

				// Decidiamo se il camminatore fa il passo lungo x, y o z
				dir_index = int(rnd.Rannyu(0.,3.)); // Vale: 0, 1 o 2
				// Il camminatore fa un passo + a o -a nella direzione selezionata
				pos[dir_index] += (2 * rnd.Bernoulli(p) - 1) * a;

				rsquared = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
				// Calcoliamo <r^2_j> e il suo quadrato per il k-esimo blocco
				ave[j + k*N] += rsquared / l;
				av2[j + k*N] += pow(rsquared, 2) / l;
			}
		}
	}


	// Calcoliamo (considerando n blocchi) la media e la dev std della media per r^2 al passo k
	for(int k=0; k<N; k++){		
		Rsquared_prog[k] = 0.;
		Rsquared2_prog[k] = 0.;
		for(int i=0; i<n; i++){
			Rsquared_prog[k] += ave[i*N + k];
			Rsquared2_prog[k] += av2[i*N + k];
		}
		Rsquared_prog[k] /= n;
		Rsquared2_prog[k] /= n;

		err_prog[k] = error(Rsquared_prog[k], Rsquared2_prog[k], n);

		file1 << Rsquared_prog[k] << setw(wd) << err_prog[k] << endl;
	}

	// ************************* PARTE 2 *************************

	double theta, phi; // Angoli che compongono l'angolo solido
	
	// Apertura dei file necessari per scrivere i dati della parte 1
	ofstream file2("./Data/Es02.2/Rsquared_es2.2.b.txt");

	// *********************** SVOLGIMENTO ***********************

	for(int i=0; i<N*n; i++){
		ave[i] = 0.;
		av2[i] = 0.;
	}

	// Facciamo scorrere l'indice k, che indica che siamo nel k-esimo blocco
	for(int k=0; k<n; k++){
	// In un blocco simulo l random walk da N passi l'uno
		for(int i=0; i<l; i++){
		// Fissato i, ho fissato il RW che sto considerando
			pos[0] = 0.;
			pos[1] = 0.;
			pos[2] = 0.;
			for(int j=0; j<N; j++){
				// Fissato j, ho fissato il numero di passi che considero
				rnd.AngoloSolido(theta,phi);
				
				pos[0] += a * sin(theta) * cos(phi);
				pos[1] += a * sin(theta) * sin(phi);
				pos[2] += a * cos(theta);				

				rsquared = pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2];
				// Calcoliamo <r^2_j> e il suo quadrato per il k-esimo blocco
				ave[j + k*N] += rsquared / l;
				av2[j + k*N] += pow(rsquared, 2) / l;
			}
		}
	}


	// Calcoliamo (considerando n blocchi) la media e la dev std della media per r^2 al passo k
	for(int k=0; k<N; k++){		
		Rsquared_prog[k] = 0.;
		Rsquared2_prog[k] = 0.;
		for(int i=0; i<n; i++){
			Rsquared_prog[k] += ave[i*N + k];
			Rsquared2_prog[k] += av2[i*N + k];
		}
		Rsquared_prog[k] /= n;
		Rsquared2_prog[k] /= n;

		err_prog[k] = error(Rsquared_prog[k], Rsquared2_prog[k], n);

		file2 << Rsquared_prog[k] << setw(wd) << err_prog[k] << endl;
	}

	// Salvo il file seed.out

	rnd.SaveSeed();
	
	// Elimino i vettori allocati precendentemente	

	delete[] ave;
	delete[] av2;
	delete[] Rsquared_prog;
	delete[] Rsquared2_prog;
	delete[] err_prog;

	// Chiudo i file aperti precedentemente
	
	file1.close();
	file2.close();

	return 0;
}
