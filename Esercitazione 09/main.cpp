#include "GeneticAlg.h"

using namespace std;

// Global variables that we need
int soc; // 0 = circumference, !0 = square	
int ncities, npop, ngenerations;
string filename;

void Input();
double L2_average(vector<double>);

int main (int argc, char *argv[]){

	cout << "\n---------- Inizio del programma ----------\n";
	
	Input();
	
	// Open files
	ofstream cities, L2, L2_ave;
	
	cities.open("cities_"+filename+".dat");
	L2.open("cost_"+filename+".dat");
	L2_ave.open("cost_ave_"+filename+".dat");
	
	Population pop(npop, ncities, soc);
	vector<vector<int>> appop(2);
	vector<vector<int>> popul(npop);
	
	// We write the cities' positions on a file
	for(int i=0; i<ncities; i++) {
		cities << pop.get_ind(0).get_cities()[i].x << "  "<< pop.get_ind(0).get_cities()[i].y << endl;
	}
	cities.close();
	
	cities.open("cities_final_"+filename+".dat");
	
	pop.order_pop();
	cout << "Migliore valore di L2 iniziale: " << pop.get_votes()[npop - 1] << endl;
	
	// We write the best cost function value and the average
	// value of the cost function for the starting population
	L2 << 0 << "  " << pop.get_votes()[npop - 1] << endl;
	L2_ave << 0 << "  " << L2_average(pop.get_votes()) << endl;
	
	// *************** The genetic part of the algorithm starts here ***************
	for(int j=0; j<ngenerations; j++){
		for(int i=0; i<npop/2; i++){
			appop = pop.crossover(); // This gives us two offsprings
			popul[2*i] = appop[0];
			popul[2*i + 1] = appop[1];

		}
		// We write the new population in the object pop
		pop.set_pop(popul, npop);
		pop.mutations();
		
		pop.order_pop();
		
		// We write the best cost function value and the average
		// value of the cost function for the i-th population
		L2 << j+1 << "  " << pop.get_votes()[npop - 1] << endl;
		L2_ave << j+1 << "  " << L2_average(pop.get_votes()) << endl;
	}
	// *************** The genetic part of the algorithm ends here ***************

	cout << "Migliore valore di L2 finale: " << pop.get_votes()[npop-1] << endl;
	
	// We write the final order of the cities for the individual with the best L2 value
	for(int i=0; i<ncities; i++) {
		cities << pop.get_ind(npop-1).get_cities()[i].x << "  "
			   << pop.get_ind(npop-1).get_cities()[i].y << endl;
	}
	
	// We write the first city another time, because
	// it is also the ending point of the journey
	cities << pop.get_ind(npop-1).get_cities()[0].x << "  "
		   << pop.get_ind(npop-1).get_cities()[0].y << endl;

	cout << "\n----------- Fine del programma -----------\n\n";

	// Close files
	cities.close();
	L2.close();
	L2_ave.close();
	
	// Save random numbers' seed
	rnd.SaveSeed();

	return 0;
}

void Input(){

	ifstream ReadInput("input.dat");
	
	ReadInput >> soc; // 0 = circumference, !0 = square
	ReadInput >> ncities;
	ReadInput >> npop;
	ReadInput >> ngenerations;
	
	if(soc == 0){
		cout << "\nGeneriamo le città su una circonferenza\n";
		filename = "circ";
	} else {
		cout << "\nGeneriamo le città in un quadrato\n";
		filename = "square";
	}
	
	cout << "Numero di città = " << ncities << endl;
	cout << "Numero di commessi viaggiatori = " << npop << endl;
	cout << "Numero di generazioni = " << ngenerations << endl;
	
	cout << "---------------------------------\n\n";
	
	ReadInput.close();
	
}

double L2_average(vector<double> votes){

	double ave = 0;
	for(int k=npop-1; k>npop/2; k--){
		ave += votes[k];
	}
	ave /= (npop/2 - 1);
	
	return ave;
}

