#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include "Funzioni.h"

using namespace std;

int main (int argc, char *argv[]){
	
	
	int M = 1E6; // Numero di tiri
	int nblocks = 100; // Numero di blocchi
	int neq = 0; // Numero di step MC necessari per raggiungere il sampling corretto
	int L = int((M-neq)/nblocks); // Numero di esperimenti (o tiri) per blocco

	double X, Y, Z;
	double *x = new double[M];
	double *y = new double[M];
	double *z = new double[M];
	double *r = new double[M-neq];

	double rand, alpha;

	// Inizializziamo una variabile di tipo Random
	Random rnd;
	rnd.SetRandom();

	// Variabili che servono per applicare il blocking method
	double *ave = new double[nblocks];
	double *av2 = new double[nblocks];
	double r_prog, r2_prog, err_prog;

	// *********************** SVOLGIMENTO ***********************

	// Implemento Metropolis per T(x|y) uniforme

	cout << "\n****** Orbitale 1s ******\n\n";

	// Apertura dei file necessari
	ifstream ReadInput("input.dat");
	ofstream fileXYZ;
	ofstream fileBlocks;

	double step = 1.23; // scelto calcolando la media di alpha con diversi valori di step
	double ave_alpha = 0, ave_alpha_tot = 0;
	double max_alpha, min_alpha;
	
	int point;
	ReadInput >> point;
	string name = to_string(point);

	fileXYZ.open("XYZ_1s_"+name+".out");
	fileBlocks.open("r_1s_"+name+".txt");

	x[0] = rnd.Rannyu(point - step, point + step);
	y[0] = rnd.Rannyu(point - step, point + step);
	z[0] = rnd.Rannyu(point - step, point + step);

	for(int i=0; i<M; i++){
		// Compiamo un passo utilizzando la probabilità di transizione T
		X = rnd.Rannyu(x[i] - step, x[i] + step);
		Y = rnd.Rannyu(y[i] - step, y[i] + step);
		Z = rnd.Rannyu(z[i] - step, z[i] + step);
	
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
		
		// Qui avviene l'accettazione (o meno) di una mossa
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

	M -= neq;
	int i;
	// Calcoliamo le posizioni dei punti (dopo l'equilibrazione)
	for(int j=0; j<M; j++){
		i = j + neq;
		r[j] = sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]);
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
	fileXYZ.open("XYZ_2pz_"+name+".out");
	fileBlocks.open("r_2pz_"+name+".txt");

	step = 2.98; // scelto calcolando la media di alpha con diversi valori di step
	ave_alpha = 0, ave_alpha_tot = 0;

	x[0] = rnd.Rannyu(-step, step);
	y[0] = rnd.Rannyu(-step, step);
	z[0] = rnd.Rannyu(-step, step);

	for(int i=0; i<M; i++){
		X = rnd.Rannyu(x[i] - step, x[i] + step);
		Y = rnd.Rannyu(y[i] - step, y[i] + step);
		Z = rnd.Rannyu(z[i] - step, z[i] + step);
	
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
	ReadInput.close();

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
