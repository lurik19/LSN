#include "GeneticAlg.h"

using namespace std;

// Global variables that we need
int soc; // 0 = circumference, !0 = square	
int ncities, nmoves;
string filename;

void Input();

int main (int argc, char *argv[]){

	cout << "\n---------- Inizio del programma ----------\n";
	
	// Initialization of variables
	Input();
	
	// Open files
	ofstream cities, Length;
	
	cities.open("./Data/cities_"+filename+".dat");
	Length.open("./Data/Length_"+filename+".dat");
	
	Individual TSold(ncities, soc);
	Individual TSnew = TSold;

	// We write the cities' positions on a file
	for(int i=0; i<ncities; i++){
		cities << TSold.get_cities()[i].x << "  "<< TSold.get_cities()[i].y << endl;
	}
	cities.close();
	
	cout << "Energia iniziale = " << TSold.L1() << endl;
	
	cities.open("./Data/cities_final_"+filename+".dat");
	
	// *************** Simulated Annealing ***************
	// Variables that we need for Metropolis algorithm
	double Eold, Enew, p;
	double accepted = 0, attempted = 0;
	double beta;
	
	vector<double> acceptance;

	for(int i=1; i <= 1000; i++){
		accepted = 0;
		attempted = 0;
		
		beta = i;
		
		for(int j=0; j<nmoves; j++){
			TSnew = TSold;
			// We apply mutations to the individual
			TSnew.pair_perm();
			TSnew.shift_cities();
			TSnew.permutation();
			TSnew.inversion();
			
			Eold = TSold.L1();
			Enew = TSnew.L1();
			
			p = min(1., exp(-beta * (Enew - Eold)));

			// Metropolis
			if(rnd.Rannyu() < p) {
				TSold = TSnew;
				accepted += 1;
			}
			
			attempted += 1;	
		}
	
		acceptance.push_back(accepted/attempted);
		
		Length << beta << "  " << TSold.L1() << endl;
	}

	cout << "\nTemperatura tra " << 1./1000 << " e " << 1./1 << endl;
	cout << "Accettazione tra " << *min_element(acceptance.begin(), acceptance.end())*100
		 << "% e " << *max_element(acceptance.begin(), acceptance.end())*100 <<"%\n\n";

	cout << "Energia finale (dopo simulated annealing) = " << TSold.L1() << endl;

	// We write the final order of the cities for the individual with the best L1 value
	for(int i=0; i<ncities; i++) {
		cities << TSold.get_cities()[i].x << "  " << TSold.get_cities()[i].y << endl;
	}
	
	// We write the first city another time, because
	// it is also the ending point of the journey
	cities << TSold.get_cities()[0].x << "  " << TSold.get_cities()[0].y << endl;

	
	cout << "\n----------- Fine del programma -----------\n\n";

	// Close files
	cities.close();
	Length.close();
	
	// Save random numbers' seed
	rnd.SaveSeed();

	return 0;
}

void Input(){

	ifstream ReadInput("input.dat");
	
	ReadInput >> soc; // 0 = circumference, !0 = square
	ReadInput >> ncities;
	ReadInput >> nmoves;
	
	if(soc == 0){
		cout << "\nGeneriamo le città su una circonferenza\n";
		filename = "circ";
	} else {
		cout << "\nGeneriamo le città in un quadrato\n";
		filename = "square";
	}
	
	cout << "Numero di città = " << ncities << endl;
	cout << "Numero di iterazioni per ogni temperatura = " << nmoves << endl;
	
	cout << "---------------------------------\n\n";
	
	ReadInput.close();
	
}
