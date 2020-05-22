#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Funzioni.h"

using namespace std;

int main (int argc, char *argv[]){
	
	// ************************* PARTE 1 *************************
	
	int M = 1E6; // Numero di tiri
	int N = 100; // Numero di blocchi
	int L = int(M/N); // Numero di esperimenti (o tiri) per blocco

	// Parametri del nostro sistema
	int S_0 = 100, K = 100;
	double T = 1., r = 0.1, sigma = 0.25;


	double W_T, S_T;
	double C_0, P_0;

	// Inizializziamo una variabile di tipo Random
	Random rnd;
	rnd.SetRandom();

	// Variabili che servono per applicare il blocking method
	double *ave_C0 = new double[N];
	double *av2_C0 = new double[N];
	double *C0_prog = new double[N];
	double *C02_prog = new double[N];
	double *errC0_prog = new double[N];

	double *ave_P0 = new double[N];
	double *av2_P0 = new double[N];
	double *P0_prog = new double[N];
	double *P02_prog = new double[N];
	double *errP0_prog = new double[N];

	// Apertura dei file necessari per scrivere i dati della parte 1
	ofstream file1("./Data/Es03.1/C0_es3.1.a.txt");
	ofstream file2("./Data/Es03.1/P0_es3.1.a.txt");
	int wd = 24;

	// *********************** SVOLGIMENTO ***********************

	for(int i=0; i<N; i++){
		C_0 = 0.;
		P_0 = 0.;
		for(int j=0; j<L; j++){
			// Campioniamo S(T), il che significa campionare W(T) (Gauss(0,T))
			W_T = rnd.Gauss(0, sqrt(T));	
			S_T = S_0 * exp((r - 0.5 * sigma * sigma) * T + sigma * W_T);

			// C[S(0),0]: ev. temp. al contrario su max(0, S(T) - K)
			C_0 += exp(-r * T) * max(0., S_T - K);

			// P[S(0),0]: ev. temp. al contrario su max(0, K - S(T))
			P_0 += exp(-r * T) * max(0., K - S_T);
		}
		C_0 /= L;
		ave_C0[i] = C_0;
		av2_C0[i] = pow(C_0, 2);

		P_0 /= L;
		ave_P0[i] = P_0;
		av2_P0[i] = pow(P_0, 2);
	}

	// Ora ne calcoliamo la media (su k blocchi) e l'incertezza statistica
	for(int k=0; k<N; k++){
		C0_prog[k] = 0.;
		C02_prog[k] = 0.;

		P0_prog[k] = 0.;
		P02_prog[k] = 0.;

		for(int j=0; j<k+1; j++){
			C0_prog[k] += ave_C0[j];
			C02_prog[k] += av2_C0[j];

			P0_prog[k] += ave_P0[j];
			P02_prog[k] += av2_P0[j];
		}

		C0_prog[k] /= (k+1);
		C02_prog[k] /= (k+1);
		errC0_prog[k] = error(C0_prog, C02_prog, k);

		P0_prog[k] /= (k+1);
		P02_prog[k] /= (k+1);
		errP0_prog[k] = error(P0_prog, P02_prog, k);

		// Scriviamo su file i risultati
		file1 << C0_prog[k] << setw(wd) << errC0_prog[k] << '\n';
		file2 << P0_prog[k] << setw(wd) << errP0_prog[k] << '\n';
	}

	// ************************* PARTE 2 *************************

	double Z_k;
	int N_t = 100; // Numero di intervalli temporali in cui dividiamo [0, T]	
	double S_t;
	
	// Apertura dei file necessari per scrivere i dati della parte 2
	ofstream file3("./Data/Es03.1/C0_es3.1.b.txt");
	ofstream file4("./Data/Es03.1/P0_es3.1.b.txt");


	// *********************** SVOLGIMENTO ***********************

	for(int i=0; i<N; i++){
		C_0 = 0.;
		P_0 = 0.;
		for(int j=0; j<L; j++){			
			S_t = S_0;
			for(int k=0; k<N_t; k++){
				Z_k = rnd.Gauss(0., 1.);
				// notiamo che t_k+1 - t_k = (k + 1 - k) * T / N_t = T / N_t
				S_t *= exp((r - 0.5 * sigma * sigma) * T / N_t + sigma * Z_k * sqrt(T / N_t));		
			}
			
			C_0 += exp(-r * T) * max(0., S_t - K);
			P_0 += exp(-r * T) * max(0., K - S_t);
		}
			C_0 /= L;
			ave_C0[i] = C_0;
			av2_C0[i] = pow(C_0, 2);

			P_0 /= L;
			ave_P0[i] = P_0;
			av2_P0[i] = pow(P_0, 2);
	}

	// Ora ne calcoliamo la media (su k blocchi) e l'incertezza statistica
	for(int k=0; k<N; k++){
		C0_prog[k] = 0.;
		C02_prog[k] = 0.;

		P0_prog[k] = 0.;
		P02_prog[k] = 0.;

		for(int j=0; j<k+1; j++){
			C0_prog[k] += ave_C0[j];
			C02_prog[k] += av2_C0[j];

			P0_prog[k] += ave_P0[j];
			P02_prog[k] += av2_P0[j];
		}

		C0_prog[k] /= (k+1);
		C02_prog[k] /= (k+1);
		errC0_prog[k] = error(C0_prog, C02_prog, k);

		P0_prog[k] /= (k+1);
		P02_prog[k] /= (k+1);
		errP0_prog[k] = error(P0_prog, P02_prog, k);

		// Scriviamo su file i risultati
		file3 << C0_prog[k] << setw(wd) << errC0_prog[k] << '\n';

		file4 << P0_prog[k] << setw(wd) << errP0_prog[k] << '\n';
	}




	// Salvo il file seed.out

	rnd.SaveSeed();

	// Elimino i vettori allocati precendentemente	

	delete[] ave_C0;
	delete[] av2_C0;
	delete[] C0_prog;
	delete[] C02_prog;
	delete[] errC0_prog;

	delete[] ave_P0;
	delete[] av2_P0;
	delete[] P0_prog;
	delete[] P02_prog;
	delete[] errP0_prog;

	// Chiudo i file aperti precedentemente

	file1.close();
	file2.close();

	file3.close();
	file4.close();

	return 0;
}
