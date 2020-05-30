#include <iostream>
#include <fstream>
#include <cmath>
#include "Funzioni.h"

using namespace std;

int main (int argc, char *argv[]){
	
	
	int M = 1E6; // Numero di tiri
	int nblocks = 100; // Numero di blocchi
	int L = int(M/nblocks); // Numero di esperimenti (o tiri) per blocco


	double X, Y, Z;
	double *x = new double[M];
	double *y = new double[M];
	double *z = new double[M];
	double *r = new double[M];

	double rand, alpha;

	// Inizializziamo una variabile di tipo Random
	Random rnd;
	rnd.SetRandom();

	// Variabili che servono per applicare il blocking method
	double *ave = new double[nblocks];
	double *av2 = new double[nblocks];
	double r_prog, r2_prog, err_prog;

	// *********************** SVOLGIMENTO ***********************

	// Implemento Metropolis per T(x|y) Gaussiana


	cout << "\n****** Orbitale 1s ******\n\n";

	// Apertura dei file necessari
	ofstream fileXYZ;
	ofstream fileBlocks;

	fileXYZ.open("XYZ_1s_gauss.out");
	fileBlocks.open("r_1s_gauss.txt");

	double step = 0.76; // scelto calcolando la media di alpha con diversi valori di step
	double ave_alpha = 0, ave_alpha_tot = 0;
	double max_alpha, min_alpha;

	x[0] = rnd.Gauss(0, step);
	y[0] = rnd.Gauss(0, step);
	z[0] = rnd.Gauss(0, step);

	for(int i=0; i<M; i++){
		X = rnd.Gauss(x[i], step);
		Y = rnd.Gauss(y[i], step);
		Z = rnd.Gauss(z[i], step);
	
		alpha = min(1.,(1./M_PI*exp(-2*sqrt(X*X+Y*Y+Z*Z)))/
			       (1./M_PI*exp(-2*sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]))));

		ave_alpha_tot += alpha;
		ave_alpha += alpha;
		if(i%L == 0 && i > 0){
			ave_alpha = ave_alpha*100/L;
			if (ave_alpha < min_alpha || i == L) min_alpha = ave_alpha;
			if (ave_alpha > max_alpha || i == L) max_alpha = ave_alpha;
			ave_alpha = 0;
		}
			
		rand = rnd.Rannyu();
		if (rand <= alpha) {x[i+1] = X; y[i+1] = Y; z[i+1] = Z;}
		else		   {x[i+1] = x[i]; y[i+1] = y[i]; z[i+1] = z[i];}
	}
	cout << "Max [<A(x_i+1|x_i)>]_block = " << max_alpha << "%" << endl;
	cout << "Min [<A(x_i+1|x_i)>]_block = " << min_alpha << "%" << endl << endl;
	
	ave_alpha_tot /= M;
	cout << "<A(x_i+1|x_i)>_tot = " << ave_alpha_tot*100 << "%" << endl;


	for(int i=0; i<M; i++){
		fileXYZ << x[i] << "   " << y[i] << "   " << z[i] << endl;
	}

	for(int i=0; i<M; i++){
		r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
	}

	// Calcoliamo <r> per ogni blocco
	for(int i=0; i<nblocks; i++){
		ave[i] = 0;
		for(int j=0; j<L; j++){
			ave[i] += r[j + i * L]/L;
		}
		av2[i] = pow(ave[i],2);
	}


	// Ora ne calcoliamo la media (su i blocchi) e l'incertezza statistica
	for(int i=0; i<nblocks; i++){
		r_prog = 0;
		r2_prog = 0;
		for(int j=0; j<i+1; j++){
			r_prog += ave[j];
			r2_prog += av2[j];
		}
		r_prog/=(i+1);
		r2_prog/=(i+1);
		err_prog = error(r_prog, r2_prog, i);

		// Scriviamo su file i risultati
		fileBlocks << r_prog << "   " << err_prog << endl;
	}

	fileXYZ.close();
	fileBlocks.close();

	cout << "\n****** Orbitale 2pz ******\n\n";

	// Apertura dei file necessari
	fileXYZ.open("XYZ_2pz_gauss.out");
	fileBlocks.open("r_2pz_gauss.txt");

	step = 1.87; // scelto calcolando la media di alpha con diversi valori di step
	ave_alpha = 0, ave_alpha_tot = 0;

	x[0] = rnd.Gauss(0, step);
	y[0] = rnd.Gauss(0, step);
	z[0] = rnd.Gauss(0, step);

	for(int i=0; i<M; i++){
		X = rnd.Gauss(x[i], step);
		Y = rnd.Gauss(y[i], step);
		Z = rnd.Gauss(z[i], step);
	
		alpha = min(1., (1./(32*M_PI)*Z*Z*exp(-sqrt(X*X+Y*Y+Z*Z)))/
			(1./(32*M_PI)*z[i]*z[i]*exp(-sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]))));

		ave_alpha_tot += alpha;
		ave_alpha += alpha;
		if(i%L == 0 && i > 0){
			ave_alpha = ave_alpha*100/L;
			if (ave_alpha < min_alpha || i == L) min_alpha = ave_alpha;
			if (ave_alpha > max_alpha || i == L) max_alpha = ave_alpha;
			ave_alpha = 0;
		}
				
		rand = rnd.Rannyu();
		if (rand <= alpha) {x[i+1] = X; y[i+1] = Y; z[i+1] = Z;}
		else		   {x[i+1] = x[i]; y[i+1] = y[i]; z[i+1] = z[i];}
	}
	cout << "Max [<A(x_i+1|x_i)>]_block = " << max_alpha << "%" << endl;
	cout << "Min [<A(x_i+1|x_i)>]_block = " << min_alpha << "%" << endl << endl;
	
	ave_alpha_tot /= M;
	cout << "<A(x_i+1|x_i)>_tot = " << ave_alpha_tot*100 << "%" << endl;


	for(int i=0; i<M; i++){
		fileXYZ << x[i] << "   " << y[i] << "   " << z[i] << endl;
	}

	for(int i=0; i<M; i++){
		r[i] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
	}

	// Calcoliamo <r> per ogni blocco
	for(int i=0; i<nblocks; i++){
		ave[i] = 0;
		for(int j=0; j<L; j++){
			ave[i] += r[j + i * L]/L;
		}
		av2[i] = pow(ave[i],2);
	}


	// Ora ne calcoliamo la media (su i blocchi) e l'incertezza statistica
	for(int i=0; i<nblocks; i++){
		r_prog = 0;
		r2_prog = 0;
		for(int j=0; j<i+1; j++){
			r_prog += ave[j];
			r2_prog += av2[j];
		}
		r_prog/=(i+1);
		r2_prog/=(i+1);
		err_prog = error(r_prog, r2_prog, i);

		// Scriviamo su file i risultati
		fileBlocks << r_prog << "   " << err_prog << endl;
	}

	fileXYZ.close();
	fileBlocks.close();

	// Salvo il file seed.out

	rnd.SaveSeed();

	// Elimino i vettori allocati precendentemente	
	
	delete[] x;
	delete[] y;
	delete[] z;
	delete[] r;
	delete[] ave;
	delete[] av2;

	// Chiudo i file aperti precedentemente

	return 0;
}
