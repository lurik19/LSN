#ifndef GA_H
#define GA_H

#include <iostream>		// cout
#include <fstream>		// ofstream
#include <cmath>		// M_PI, cos, sin, pow
#include <string>		// string
#include <vector>		// vector
#include <algorithm>	// sort, random_shuffle, upper_bound, swap

#include "random.h"

Random rnd;

using namespace std;

struct City{

	double x;
	double y;

};

class Individual {

private:
	int _ncities;
	vector<int> _path;
	vector<City> _cities;
	
	double _pm = 0.05; // Mutation probability

protected:

public:

	// Constructors
	Individual();
	Individual(int, int);
	
	// Copy constructor
	Individual(const Individual&);
	// Destructors
	~Individual();

	// Methods
	vector<int> get_path();
	void print_path(int);
	void set_path(vector<int>);
	
	vector<City> get_cities();
	void set_cities(vector<City>);
	
	bool check();
	double L1();
	
	// Pbc
	int Pbc(int);
	int PbcWO1(int); // Pbc without considering 1 (which has to stay in first position)
					 // (we will use this method in the mutations)
	// Genetic mutations
	void pair_perm();
	void shift_cities();
	void permutation();
	void inversion();

};

class Population {

private:
	int _npop;
	int _ncities;
	vector<Individual> _pop;
	vector<double> _votes;
	
	double _pc = 0.7; // Crossover probability

protected:

public:

	// Constructors
	Population(int, int, int);
	// Destructors
	~Population();

	// Methods
	vector<Individual> get_pop();
	void set_pop(vector<vector<int>>, int);
	void print_pop();
	
	Individual get_ind(int);
	vector<double> get_votes();
	void set_cities(vector<City>);
	void set_path(vector<int>, int);
	
	// Crossover
	void order_pop();
	int select();
	vector<vector<int>> crossover();
	
	// Mutations
	void mutations();

};

// ******************* Implementation of the methods ******************* //

// ******************* Individual ******************* //

Individual :: Individual() {}

Individual :: Individual(int ncities, int soc) : _cities(ncities) {

	_ncities = ncities;

	for(int i=1; i<=_ncities; i++){
		_path.push_back(i);
	}

	random_shuffle(_path.begin() + 1, _path.end());
	
	// We randomly initialize the positions of the cities
	if(soc == 0){ // soc means "square or circumference"
		// On a circumference
		double theta;
		for(int i=0; i<_ncities; i++){
			theta = rnd.Rannyu(0, 2*M_PI);
			_cities[i].x = cos(theta);
			_cities[i].y = sin(theta);
		}
	} else {
		// In a square
		for(int i=0; i<_ncities; i++){
			// Generate _x and _y between -1 and 1
			_cities[i].x = 2 * rnd.Rannyu() - 1.;
			_cities[i].y = 2 * rnd.Rannyu() - 1.;
		}		
	}

}

Individual :: Individual(const Individual& ind){

	_ncities = ind._ncities;
	_path = ind._path;
	_cities = ind._cities;

}

Individual :: ~Individual() {}

vector<int> Individual :: get_path(){

	return _path;

}

void Individual :: print_path(int endline){

	cout << "[";
	for(int i=0; i<_ncities; i++){
		cout << _path[i];
		if(i < _ncities - 1) cout << ",";
	}
	cout << "]";
	if(endline == 1) cout << endl;

}

void Individual :: set_path(vector<int> path){

	_path = path;

}

vector<City> Individual :: get_cities(){
	
	// We return the cities in the order
	// in which the TS visit them
	vector<City> cities(_ncities);
	for(int i=0; i<_ncities; i++){
		cities[i] = _cities[_path[i]-1];
	}

	return cities;

}

void Individual :: set_cities(vector<City> cities){

	_cities = cities;

}

// We check if an individual respects
// the bounds of the problem
bool Individual :: check(){

	bool check=true;
	int count;
	if(_path[0] != 1) check = false;
	else{
		for(int i=1; i<=_ncities; i++){
			count = 0;
			for(int j=0; j<_ncities; j++) if(_path[j] == i) count++;
			if(count > 1) {check = false; break;}
		}
	}

	return check;
}

double Individual :: L1(){

	double L1 = 0;
	for(int i=0; i<_ncities; i++){
		L1 += sqrt(pow(_cities[_path[i]-1].x - _cities[_path[Pbc(i+1)]-1].x, 2)
			 + pow(_cities[_path[i]-1].y - _cities[_path[Pbc(i+1)]-1].y, 2));
	}

	return L1;

}

int Individual :: Pbc(int i){
	
	return i % _ncities;

}

int Individual :: PbcWO1(int i){

	if(i < _ncities) return i;
	else			 return i % _ncities + 1;

}

void Individual :: pair_perm(){

	int i1, i2;
	if(rnd.Rannyu() < _pm){
		// Generate two random integers between 2 and 32 (indeces for the swap)
		i1 = rnd.RandInt(1, _ncities - 1);
		do{
			i2 = rnd.RandInt(1, _ncities - 1);
		} while(i2 == i1);
		swap(_path[i1], _path[i2]);
	}

}

void Individual :: shift_cities(){

	int i1, m, n;
	if(rnd.Rannyu() < _pm){
		i1 = rnd.RandInt(1, _ncities - 1);
		m = rnd.RandInt(1, _ncities - 1);
		
		// n goes from 1 to (m-1). At n=m
		// we go back to the initial individual
		n = rnd.RandInt(1, m);
		for(int i=0; i<n; i++){
			// This for represent a single shift
			for(int j=0; j<m-1; j++){
				swap(_path[PbcWO1(i1 + j)], _path[PbcWO1(i1 + m - 1)]);
			}
		}
		
	}
}

void Individual :: permutation(){

	int m, i1, i2;
	if(rnd.Rannyu() < _pm){
		m = rnd.RandInt(1, (_ncities + 1) / 2 - 1); // 1 <= m < _ncities/2
		i1 = rnd.RandInt(1, _ncities - 1);
		i2 = PbcWO1(rnd.RandInt(i1 + m, _ncities + i1 - m - 1));
		for(int i=0; i<m; i++){
			swap(_path[PbcWO1(i1+i)], _path[PbcWO1(i2+i)]);
		}
	}

}

void Individual :: inversion(){

	int i1, m;
	if(rnd.Rannyu() < _pm){
		i1 = rnd.RandInt(1, _ncities - 1);
		m = rnd.RandInt(1, _ncities);
		for(int i=0; i<(int)m/2; i++){ // The (int) is not necessary
			swap(_path[PbcWO1(i1 + i)], _path[PbcWO1(i1 + (m - 1) - i)]);
		}
	}
}

// ******************* Population ******************* //

Population :: Population(int npop, int ncities, int soc) : _votes(npop, 0) {			
									// There is no need to check the starting	
	_npop = npop;					// population because it has been constructed
	_ncities = ncities;				// in such a way as it is correct
	for(int i=0; i<_npop; i++){
		_pop.push_back(Individual(_ncities, soc));
	}

	vector<City> cities(_ncities);
	// We randomly initialize the positions of the cities
	// NOTE: All the individuals of the population share
	// the same positions for the cities!
	if(soc == 0){ // soc means "square or circumference"
		// On a circumference
		double theta;
		for(int i=0; i<_ncities; i++){
			theta = rnd.Rannyu(0, 2*M_PI);
			cities[i].x = cos(theta);
			cities[i].y = sin(theta);
		}
		for(int i=0; i<_npop; i++){
			_pop[i].set_cities(cities);
		}
	} else {
		// In a square
		for(int i=0; i<_ncities; i++){
			// Generate _x and _y between -1 and 1	
			cities[i].x = 2 * rnd.Rannyu() - 1.;
			cities[i].y = 2 * rnd.Rannyu() - 1.;
		}
		for(int i=0; i<_npop; i++){
			_pop[i].set_cities(cities);
		}
	}
}

Population :: ~Population() {}

Individual Population :: get_ind(int i){

	return _pop[i];

}

vector<Individual> Population :: get_pop(){

	return _pop;

}

void Population :: print_pop(){
	
	for(int i=0; i<_npop; i++){
		_pop[i].print_path(0);
		cout << ", L1 = " << _pop[i].L1() << endl;
	}
	
}

void Population :: set_pop(vector<vector<int>> pop, int n){

	for(int i=0; i<n; i++){
		_pop[i].set_path(pop[i]);
	}

}

vector<double> Population :: get_votes(){

	return _votes;

}

void Population :: set_cities(vector<City> cities){

	for(int i=0; i<_npop; i++){
		_pop[i].set_cities(cities);
	}
}

void Population :: set_path(vector<int> path, int ind){

	_pop[ind].set_path(path);

}

struct Votes {
	Individual ind;
	double vote;
};

void Population :: order_pop(){

	vector<Votes> votes(_npop);
	for(int i=0; i<_npop; i++){
		votes[i].ind = _pop[i];
		votes[i].vote = _pop[i].L1();
	}
	
	// We order the population on a fitness basis:
	// the higher L1, the lower the rank
	sort(votes.begin(), votes.end(),
		 [](Votes const &a, Votes const &b) { return a.vote > b.vote; });

	// We write the population (and the respective votes) in a sorted fashion
	 for(int i=0; i<_npop; i++) {
        _pop[i] = votes[i].ind;
        _votes[i] = votes[i].vote;
    }

}

int Population :: select(){

	double p = 1./5;
	double r = rnd.Rannyu();

	// We return the index of the selected individual
	return int(_npop * pow(r, p));

}

vector<vector<int>> Population :: crossover(){

	vector<vector<int>> new_pop(2);
	vector<int> ind1, ind2;
	int i1, i2, i_cut;
	int start1=1, start2 = 1;
	int elem;
	
	if(_npop > 1){
	
		i1 = select();
		do {
			i2 = select();
		} while(i2 == i1);
	
		if(rnd.Rannyu() < _pc){
			i_cut = rnd.RandInt(1, _ncities - 1); // At this position we cut the individuals
			
			ind1 = _pop[i1].get_path();
			ind2 = _pop[i2].get_path();
			
			for(int i=i_cut; i<_ncities; i++){
				// We use this cycle to write ind1 properly
				for(int j=start1; j<_ncities; j++){
					elem = _pop[i2].get_path()[j];
					if(find(ind1.begin(), ind1.begin() + i_cut, elem) == ind1.begin() + i_cut){
						// if ind2[j] is not in ind1 from 0 to i_cut-1
						ind1[i] = elem;
						start1++;
						break;
					} else { start1++;}
				}
				// We use this cycle to write ind2 properly
				for(int j=start2; j<_ncities; j++){
					elem = _pop[i1].get_path()[j];
					if(find(ind2.begin(), ind2.begin() + i_cut, elem) == ind2.begin() + i_cut){
						// if ind1[j] is not in ind2 from 0 to i_cut-1
						ind2[i] = elem;
						start2++;
						break;
					} else {start2++;}
				}
			}
			
			new_pop[0] = ind1;
			new_pop[1] = ind2;

		} else {
			// If there is no crossover we just copy
			// the selected individuals
			new_pop[0] = _pop[i1].get_path();
			new_pop[1] = _pop[i2].get_path();
		
		}
		
	} else {
			cout << "In order to have crossover, the population should be formed by at least 2 individuals\n";
		}
		
	return new_pop;
}

void Population :: mutations(){

	for(int i=0; i<_npop; i++){
		_pop[i].pair_perm();
		_pop[i].shift_cities();
		_pop[i].permutation();
		_pop[i].inversion();
	}

}

#endif
