#include <iostream>
#include <stdlib.h>
#include <armadillo>
#include <fstream>
#include <iomanip>
#include <time.h>
//#include <mpi.h>

using namespace std;
using namespace arma;

// Creating the random initial matrix of spin-values:
void random_init(mat &spin, double &E, double &M, int N){
    //srand(time(NULL));
    srand(78);          //make 78 a random number for each run later
    for (int i = 0;i<N+2;i++){
        for (int j=0;j<N+2;j++){
            int random = (rand()%2)+1;
            if (random == 1){
                spin(i,j) = 1.;
            }
            else{
                spin(i,j) = -1.;
            }
        }
    }
    // Boundary conditions
    for (int i = 0.; i< N+2; i++){
        spin(i,0) = spin(i,N);
        spin(i,N+1) = spin(i,1);
        spin(0,i) = spin(N,i);
        spin(N+1,i) = spin(1,i);
    }

    // Updating the elements of the matrices
    for (int i = 1;i<N+1;i++){
        for (int j = 1;j<N+1;j++){
            E -= spin(i,j)*(spin(i,j+1) + spin(i+1,j));
            M += (double) spin(i,j);
        }
    }
}

void change_spin(mat &spin, int i, int j, int N){

    // flip current spin:
    spin(i,j) = -1*spin(i,j);

    //boundary conditions:
    if (i==1){
        spin(N+1,j) = spin(i,j);
    }

    if (i==N){
        spin(0,j) = spin(i,j);
    }

    if (j==1){
        spin(i,N+1) = spin(i,j);
    }

    if (j==N){
        spin(i,0) = spin(i,j);
    }
}

void Metropolis(mat &spin, int N, double &E, double &M, int &MCcount, double *w){
    for(int atoms = 0; atoms <N*N; atoms++){
        int i = (rand()%N)+1;
        int j = (rand()%N)+1;
        int dE = 2.0*spin(i,j)*(spin(i+1,j)+spin(i,j+1) + spin(i-1,j)+spin(i,j-1));

        if (dE>0){
            double ran = ((double)rand()/(RAND_MAX));
            if (ran <= w[dE+8]){
                change_spin(spin,i,j,N);
                E += (double) dE;
                M += (double) 2*spin(i,j);
                MCcount += 1;
            }
        }
        else{
            change_spin(spin,i,j,N);
            E += (double) dE;
            M += (double) 2*spin(i,j);
            MCcount += 1;

        }
    }
}

void write_to_file(int N, int MCcycles, double T, vec average){

    double norm = 1.0/((double)(MCcycles));
    double E_average = average[0]*norm;
    double E2_average = average[1]*norm;
    double M_average = average[2]*norm;
    double M2_average = average[3]*norm;
    double Mabs_average = average[4]*norm;

    // all expectation values are per spin, divide by N*N

    double Evariance = (E2_average - E_average*E_average)/N/N;
    double Mvariance = (M2_average - Mabs_average*Mabs_average)/N/N;

    double Cv = Evariance/T/T;
    double chi = Mvariance/T;

    //ofile << setiosflags(ios::showpoint|ios:uppercase);
    outfile1 << setw(15) << setprecision(8) << T << endl;
    outfile2 << setw(15) << setprecision(8) << E_average/N/N << endl;
    outfile3 << setw(15) << setprecision(8) << Cv << endl;
    outfile4 << setw(15) << setprecision(8) << M_average/N/N << endl;
    outfile5 << setw(15) << setprecision(8) << chi << endl;
    outfile6 << setw(15) << setprecision(8) << Mabs_average/N/N <<endl;
}

int main(){

    outfile1.open('temperature.txt');
    outfile2.open('E_values.txt');
    outfile3.open('Cv_values.txt');
    outfile4.open('M_values.txt');
    outfile5.open('chi_values.txt');
    outfile6.open('Mabs_values.txt');


    int max_MCcycles = 1e5;
    int N = 2;

    // temperatures in units of kT:
    double initial_T = 1.0;
    double final_T = 1.0;
    double step_T = 0.1;

    // Monte Carlo trials:
    for (double T = initial_T; T<= final_T; T+= step_T){

        double beta = 1.0/T;
        double E = 0.;
        double M = 0.;

        // average vector with zeros as entry:
        double average[5];
        for(int i=0; i<5; i++) average[i] = 0.0;

        mat spin = zeros<mat>(N+2,N+2); // N+2 too make BCs easier
        random_init(spin,E, M, N);      // random initial spins tates

        // precalculating w since this saves CPU time:
        double w[17];
        for(int de=-8; de<=8; de++) w[de+8]= 0;
        for(int de=-8; de<=8; de+=4) w[de+8] = exp(-beta*de);

        // total average vector with zeros as entry:
        //double tot_average[4];
        //for(int i=0; i<=4; i++) tot_average[i] = 0.0;

        for (int cycle =1; cycle<= max_MCcycles; cycle++){
            int MCcount = 0;
            Metropolis(spin, N, E, M, MCcount, w);
            average[0] += E;
            average[1] += E*E;
            average[2] += M;
            average[3] += M*M;
            average[4] += fabs(M);
        }

        //write_to_file(N, MCcycles, T, average);
    }



}
