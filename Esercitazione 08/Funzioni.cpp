#include "Funzioni.h"

using namespace std;

// Funzione con cui calcoliamo l'incertezza statistica

double error (double* av, double* av2, int n){
	if (n==0)
		return 0;
	else
		return sqrt((av2[n]-pow(av[n],2))/n);
}

double error (double av, double av2, int n){
	if (n==0)
		return 0;
	else
		return sqrt((av2-pow(av,2))/n);
}
