/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main()
{

  Input(); //Inizialization

  for(int itemp = 0; itemp <= (2. - 0.5)/temp_step * ReadTemp; itemp++){
	
	beta = 1. / temp;
	
	cout << "*****************************" << endl;
	cout << "Temperature = " << temp << endl;
	cout << "*****************************" << endl << endl;

	for(int iblk=1; iblk <= nblk; ++iblk) //Simulation
	  {
	    Reset(iblk); //Reset block averages
	    for(int istep=1; istep <= nstep; ++istep)
	    {
	      Move(metro);
	      Measure();
	      Accumulate(); //Update block averages
	    }
	    Averages(iblk); //Print results for current block
	  }
	
	temp += temp_step;
	
  }
  ConfFinal(); //Write final configuration

  return 0;
}


void Input(void)
{
  ifstream ReadInput, ReadConf;

  cout << endl << "Classic 1D Ising model" << endl;
  cout << "Monte Carlo simulation" << endl << endl;
  cout << "Nearest neighbour interaction" << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T" << endl;
  cout << "The program uses k_B=1 and mu_B=1 units" << endl << endl;

//Read seed for random numbers
   int p1, p2;
   ifstream Primes("Primes");
   Primes >> p1 >> p2 ;
   Primes.close();

   ifstream input("seed.in");
   input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
   rnd.SetRandom(seed,p1,p2);
   input.close();
  
//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> ReadTemp;
  ReadInput >> temp_step;

  ReadInput >> temp;
  if (ReadTemp == 0){
  	beta = 1./temp; // Dubbio: potrei non metterlo
  	cout << "Temperature = " << temp << endl;
  } else {
        if (ReadTemp == 1) {
		cout << "Temperature will vary from 0.5 to 2.0 during the simulation" << endl;
		temp = 0.5;
		beta = 1./temp; // Dubbio: potrei non metterlo
	} else {
		cerr << " *** Wrong ReadTemp value: it should be 0 or 1 ***" << endl;
		exit(1);
	}
  }

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;
    
  ReadInput >> metro; // if metro=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  ReadInput >> restart;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;

  cout << "Restart = ";
  if(restart == 1) cout << "True" << endl;
  else 		   cout << "False" << endl;
  cout << endl;

  ReadInput.close();


// Prepare arrays for measurements
  iu = 0; // Energy
  ic = 1; // Heat capacity
  im = 2; // Magnetization
  ix = 3; // Magnetic susceptibility
 
  n_props = 4; // Number of observables

// initial configuration

  if (restart != 1) {
	for (int i=0; i<nspin; ++i) {
		if(rnd.Rannyu() >= 0.5) s[i] = 1;
		else s[i] = -1;
	}
  } else {
	cout << "Read previous spin configuration from file config.final " << endl << endl;
	ReadConf.open("config.final");
	for(int i=0; i<m_spin ; i++){
		ReadConf >> s[i];
	}
	ReadConf.close();
  }

// Evaluate energy etc. of the initial configuration
  Measure();

// Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl << endl;
}


void Move(int metro)
{
  int o;
  double deltaE, alpha, p;
  double r;

  for(int i=0; i<nspin; ++i)
  {
  // Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1) // Metropolis
    {
	
	deltaE = -2 * Boltzmann(s[o], o);

	alpha = min(1., exp(-beta*deltaE));
	attempted += 1;
	// Decide whether to accept the move or not
	r = rnd.Rannyu();
	if (r <= alpha) {s[o] *= -1;
			 accepted += 1;
			}
    }
    else // Gibbs sampling
    {

	deltaE = -2 * Boltzmann(s[o], o) * s[o]; // s[o] * s[o] = 1

	p = 1./(1 + exp(-beta * deltaE)); // probability to have p( s_o = 1 | {s_j : j=/=o})

	attempted += 1;
	accepted += 1;

	// Decide what value should assume spin s_o
	r = rnd.Rannyu();
	if (r <= p) {s[o] = +1;}
	else 	    {s[o] = -1;}
    }
  }
}


double Boltzmann(int sm, int ip)
{
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure()
{
  int bin;
  double u = 0.0, m = 0.0;

// cycle over spins
  for (int i=0; i<nspin; ++i)
  {
     u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
     m += s[i]; 
  }


  walker[iu] = u;
  walker[ic] = u * u;
  walker[im] = m;
  walker[ix] = beta * m * m;
}


void Reset(int iblk) // Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;
   attempted = 0;
   accepted = 0;
}


void Accumulate(void) // Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) // Print results for current block
{
    
   ofstream Ene, Heat, Mag, Chi;

   const int wd=12;
    
    //cout << "Block number " << iblk << endl;
    //cout << "Acceptance rate " << accepted/attempted << endl << endl;
    
    Ene.open("output.ene.0",ios::app);
    stima_u = blk_av[iu]/blk_norm/(double)nspin; // Energy
    glob_av[iu]  += stima_u;
    glob_av2[iu] += stima_u*stima_u;
    err_u = Error(glob_av[iu], glob_av2[iu],iblk);
    Ene << setw(wd) << iblk << setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
	Ene.close();

    Heat.open("output.heat.0", ios::app);
    stima_c = beta*beta*(blk_av[ic]/blk_norm - pow(blk_av[iu]/blk_norm,2))/(double)nspin; // Heat capacity
    glob_av[ic]  += stima_c;
    glob_av2[ic] += stima_c*stima_c;
    err_c = Error(glob_av[ic],glob_av2[ic],iblk);
    Heat << setw(wd) << iblk << setw(wd) << stima_c << setw(wd) << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
    Heat.close();

    Mag.open("output.mag.0", ios::app);
    stima_m = blk_av[im]/blk_norm/(double)nspin; // Magnetization
    glob_av[im]  += stima_m;
    glob_av2[im] += stima_m*stima_m;
    err_m = Error(glob_av[im],glob_av2[im],iblk);
    Mag << setw(wd) << iblk << setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
    Mag.close();

    Chi.open("output.chi.0", ios::app);
    stima_x = blk_av[ix]/blk_norm/(double)nspin; // Susceptibility
    glob_av[ix]  += stima_x;
    glob_av2[ix] += stima_x*stima_x;
    err_x = Error(glob_av[ix],glob_av2[ix],iblk);
    Chi << setw(wd) << iblk << setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
    Chi.close();

    string ene = "output.ene_temp", heat = "output.heat_temp",
	   mag = "output.mag_temp", chi = "output.chi_temp";

    string Metro = "_Metropolis.0", Gibbs = "_Gibbs.0";

    if (iblk == nblk){
	
	if(metro == 1){
		ene += Metro;
		heat += Metro;
		mag += Metro;
		chi += Metro;
	}
	else{
		ene += Gibbs;
		heat += Gibbs;
		mag += Gibbs;
		chi += Gibbs;
	}	

	Ene.open(ene, ios::app);	
	Ene << temp << setw(wd) << glob_av[iu]/(double)nblk << setw(wd) << err_u << endl;
	Ene.close();

	Heat.open(heat, ios::app);
	Heat << temp << setw(wd) << glob_av[ic]/(double)nblk << setw(wd) << err_c << endl;
	Heat.close();

	Mag.open(mag, ios::app);
	Mag << temp << setw(wd) << glob_av[im]/(double)nblk << setw(wd) << err_m << endl;
	Mag.close();

	Chi.open(chi, ios::app);
	Chi << temp << setw(wd) << glob_av[ix]/(double)nblk << setw(wd) << err_x << endl;
	Chi.close();
    } else {
	//cout << "----------------------------" << endl << endl;
    }
}

void ConfFinal(void)
{
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i)
  {
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
    if(i >= nspin) i = i - nspin;
    else if(i < 0) i = i + nspin;
    return i;
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
