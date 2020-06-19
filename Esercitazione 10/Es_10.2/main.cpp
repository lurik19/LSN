#include "GeneticAlg.h"
#include "mpi.h"

using namespace std;

// Global variables that we need
int soc; // 0 = circumference, !0 = square	
int ncities, npop, ngenerations, nmigr;
string filename;

// For parallel programming
int size, Rank;

void Input();
double L1_average(vector<double>);

int main (int argc, char *argv[]){

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
	MPI_Status stat;
	
	// We want to print this kind of indications only once 
	// (only for one process) just for aesthetic reasons
	if(Rank == 0) cout << "\n---------- Inizio del programma ----------\n";
	
	// Initialization of variables
	Input();
	
	// Open files
	ofstream cities, L1, L1_ave;
	
	cities.open("./Data/cities_"+filename+"_"+to_string(Rank)+".dat");
	L1.open("./Data/cost_"+filename+"_"+to_string(Rank)+".dat");
	L1_ave.open("./Data/cost_ave_"+filename+"_"+to_string(Rank)+".dat");
	
	Population pop(npop, ncities, soc);
	vector<vector<int>> appop(2);
	vector<vector<int>> popul(npop);
	
	vector<int> appo_ind1(ncities);
	vector<int> appo_ind2(ncities);
	
	// I did't find a more compact way to exchange a vector of structs
	vector<City> Cities(ncities); 
	vector<double> cities_x(ncities);
	vector<double> cities_y(ncities);
	
	vector<int> swaps{0,1,2,3};
	int tag1, tag2;
	
	// We want every population to "have" the same cities, so we
	// copy the cities from the first process to the other three
	if(Rank == 0){
		for(int i=0; i<ncities; i++){
			cities_x[i] = pop.get_ind(0).get_cities()[i].x;
			cities_y[i] = pop.get_ind(0).get_cities()[i].y;
		}	
	}
	
	MPI_Bcast(&cities_x.front(), ncities, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cities_y.front(), ncities, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	for(int i=0; i<ncities; i++){
			Cities[i].x = cities_x[i];
			Cities[i].y = cities_y[i];
	}
	
	pop.set_cities(Cities);
	
	// We write the cities' positions on a file
	for(int i=0; i<ncities; i++) {
		cities << pop.get_ind(0).get_cities()[i].x << "  " << pop.get_ind(0).get_cities()[i].y << endl;
	}
	cities.close();
	
	cities.open("./Data/cities_final_" + filename + "_" + to_string(Rank) + ".dat");
	
	// Now we order the population and we write the starting lowest value of L1
	pop.order_pop();
	cout << "Migliore valore di L1 iniziale processo " << Rank << ": " << pop.get_votes()[npop - 1] << endl;
	
	// We write the best cost function value and the average
	// value of the cost function for the starting population
	L1 << 0 << "  " << pop.get_votes()[npop - 1] << endl;
	L1_ave << 0 << "  " << L1_average(pop.get_votes()) << endl;
	
	// *************** The genetic part of the algorithm starts here ***************
	for(int j=1; j<=ngenerations; j++){
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
		L1 << j+1 << "  " << pop.get_votes()[npop - 1] << endl;
		L1_ave << j+1 << "  " << L1_average(pop.get_votes()) << endl;

		if(j % nmigr == 0){
			// We randomly shuffle the vector swaps, which we will use to
			// swap the best individuals between the processes
			if(Rank == 0) random_shuffle(swaps.begin(), swaps.end());
			MPI_Bcast(&swaps.front(), 4, MPI_INT, 0, MPI_COMM_WORLD);
			
			// We swap the best individuals between the processes with rank swaps[i] and swaps[i+1]
			// (swaps[0] with swaps[1] and swaps[2] with swaps[3])
			for(int k=0; k<2; k++){
				tag1 = 2 * k;
				tag2 = 2 * k + 1;
				if(Rank == swaps[2*k]){
					appo_ind1 = pop.get_ind(npop-1).get_path();
					MPI_Send(&appo_ind1.front(), ncities, MPI_INT, swaps[2*k+1], tag1, MPI_COMM_WORLD);
					MPI_Recv(&appo_ind2.front(), ncities, MPI_INT, swaps[2*k+1], tag2, MPI_COMM_WORLD, &stat);
					pop.set_path(appo_ind2, npop-1);
				}
				
				if(Rank == swaps[2*k+1]){
					appo_ind2 = pop.get_ind(npop-1).get_path();
					MPI_Recv(&appo_ind1.front(), ncities, MPI_INT, swaps[2*k], tag1, MPI_COMM_WORLD, &stat);
					MPI_Send(&appo_ind2.front(), ncities, MPI_INT, swaps[2*k], tag2, MPI_COMM_WORLD);
					pop.set_path(appo_ind1, npop-1);
				}
			}
			// Then we order the populations
			pop.order_pop();
		}
		
	}
	// *************** The genetic part of the algorithm ends here ***************

	//  We write the final lowest value of L1
	cout << "Migliore valore di L1 finale processo " << Rank << ": " << pop.get_votes()[npop-1] << endl;
	
	// We write the final order of the cities for the individual with the best L1 value
	for(int i=0; i<ncities; i++) {
		cities << pop.get_ind(npop-1).get_cities()[i].x << "  "
			   << pop.get_ind(npop-1).get_cities()[i].y << endl;
	}
	
	// We write the first city another time, because
	// it is also the ending point of the journey
	cities << pop.get_ind(npop-1).get_cities()[0].x << "  "
		   << pop.get_ind(npop-1).get_cities()[0].y << endl;

	if(Rank==0) cout << "\n----------- Fine del programma -----------\n\n";

	// Close files
	cities.close();
	L1.close();
	L1_ave.close();
	
	// Save random numbers' seed
	rnd.SaveSeed();

	MPI_Finalize();
	return 0;
}

void Input(){

	rnd.SetRandom(Rank);
	srand(Rank + 1); // We need to have different initializations (in every process)
					 // for rand because I use random_shuffle in the different programs
	
	ifstream ReadInput("input.dat");

	ReadInput >> soc; // 0 = circumference, !0 = square
	ReadInput >> ncities;
	ReadInput >> npop;
	ReadInput >> ngenerations;
	ReadInput >> nmigr;
	
	if(soc == 0){
		if(Rank == 0) cout << "\nGeneriamo le città su una circonferenza\n";
		filename = "circ";
	} else {
		if(Rank == 0) cout << "\nGeneriamo le città in un quadrato\n";
		filename = "square";
	}
	
	if(Rank == 0){
		cout << "Numero di città = " << ncities << endl;
		cout << "Numero di commessi viaggiatori = " << npop << endl;
		cout << "Numero di generazioni = " << ngenerations << endl;
		cout << "Generazioni dopo cui si ha una migrazione = " << nmigr << endl;		
		cout << "---------------------------------\n\n";
	}
	ReadInput.close();
	
}

double L1_average(vector<double> votes){

	double ave = 0;
	for(int k=npop-1; k>npop/2; k--){
		ave += votes[k];
	}
	ave /= (npop/2 - 1);
	
	return ave;
}

