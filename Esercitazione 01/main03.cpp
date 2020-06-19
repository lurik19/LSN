#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Funzioni.h"

using namespace std;

int main (int argc, char *argv[]){
	
	// ************************* PARTE 1 *************************

	double d = 1.; // Distanza tra le righe
	double L = 0.7*d; // Lunghezza dell'ago

	double side = 1E3*d; // Supponiamo di avere un "pavimento" su cui lanciamo il nostro ago
		            	 // e questo pavimento è quadrato e di lato "side"

	int M = 1E5; // Numero di tiri
	int N = 100; // Numero di blocchi
	int num_exp = M/N; // Numero di esperimenti per blocco

	int N_hit; // Numero di volte che l'ago finisce su una linea

	double y1; // Coordinata di uno degli estremi dell'ago (che genereremo random)
	double y2; // Coordinata dell'altro estremo dell'ago (che calcoleremo)
		  	   // In teoria dovremmo generare anche le x1 e x2, ma al fine della
		   	   // simulazione sono inutili.

	double theta; // Per ottenere y2 (e x2) da y1 (e x1) genereremo un
		      // angolo theta random in modo uniforme theta tra 0 e 2*PI

	// Inizializziamo una variabile di tipo Random
	Random rnd;
	rnd.SetRandom();

	// Variabili che servono per applicare il blocking method
	double *ave = new double[N];
	double *av2 = new double[N];
	double *pi_prog = new double[N];
	double *pi2_prog = new double[N];
	double *err_prog = new double[N];

	// Apertura dei file necessari per scrivere i dati
	ofstream file1("./Data/Es01.3/pi_es1.3.txt");
	int wd = 24;

	// *********************** SVOLGIMENTO ***********************

	for(int k=0; k<N; k++){
		// Ora facciamo i num_exp = M/N esperimenti di UN blocco		
		N_hit = 0;
		for(int i=k*num_exp; i<k*num_exp + num_exp; i++){
			y1 = rnd.Rannyu(0. + L, side - L); // (0 + L) e (side - L): condizioni poste affinchè
							  				   // l'ago sia tutto all'interno del pavimento
			theta = rnd.Rannyu(0., 2.*M_PI);
			// Generiamo y1 (e x1) in modo random uniforme
			y2 = y1 + L*sin(theta);

			// Controlliamo se l'ago interseca una delle linee
			for(int j=1; j<int(side/d); j++){
				if ((y1<=j*d && y2>=j*d) || (y1>=j*d && y2<=j*d)){
					N_hit++;
					j = int(side/d); // Questo comando interrompe il
				}					 // ciclo quando troviamo l'ago	
			}
		}
		// Calcoliamo A_i e (A_i)^2 di UN blocco
		ave[k] = 2*L*num_exp/(double(N_hit)*d);
		av2[k] = pow(2*L*num_exp/(double(N_hit)*d),2);
	}
	
	// Ora facciamo la media di N_hit e (N_hit)^2 tra i diversi blocchi
	for(int k=0; k<N; k++){
		pi_prog[k] = 0;
		pi2_prog[k] = 0;
		
		for(int j=0; j<k+1; j++){
			pi_prog[k] += ave[j];
			pi2_prog[k] += av2[j];
		}
		// Calcoliamo <A>_k e <A^2>_k e l'incertezza statistica k-esima
		pi_prog[k] /= (k+1);
		pi2_prog[k] /= (k+1);
		err_prog[k] = error(pi_prog, pi2_prog, k);
		
		file1 << pi_prog[k] << setw(wd) << err_prog[k] << '\n';
	}

	// Salvo il file seed.out

	rnd.SaveSeed();

	// Elimino i vettori allocati precendentemente	

	delete[] ave;
	delete[] av2;
	delete[] pi_prog;
	delete[] pi2_prog;
	delete[] err_prog;
	
	// Chiudo i file aperti precedentemente

	file1.close();
	
	return 0;
}
