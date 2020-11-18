# include <chrono>
# include <random>
# include <cstdlib>
# include <sstream>
# include <string>
# include <fstream>
# include <iostream>
# include <iomanip>
# include <mpi.h>

// Attention , ne marche qu'en C++ 11 ou supérieur :
double approximate_pi( unsigned long nbSamples ) 
{
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = beginning.time_since_epoch();
    unsigned seed = d.count();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution <double> distribution ( -1.0 ,1.0);
    unsigned long nbDarts = 0;
    // Throw nbSamples darts in the unit square [-1 :1] x [-1 :1]
    for ( unsigned sample = 0 ; sample < nbSamples ; ++ sample ) {
        double x = distribution(generator);
        double y = distribution(generator);
        // Test if the dart is in the unit disk
        if ( x*x+y*y<=1 ) nbDarts ++;
    }
    // Number of nbDarts throwed in the unit disk
    double ratio = double(nbDarts)/double(nbSamples);
    return 4*ratio;
}

int main(int nargs, char* argv[]){
	MPI_Init( &nargs, &argv );
	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
	int nbp;
	MPI_Comm_size(globComm, &nbp);
	int rank;
	MPI_Comm_rank(globComm, &rank);
	int tag = 1;
	MPI_Status Stat;
	
	double val;
	double pi = 0.0;
	unsigned long nbSamples = 1000;
	
	if(rank == 0){
		std::cout << "Tâche numéro (" << rank << ") ";
		val = approximate_pi(nbSamples);
		std::cout << "Valeur obtenue: " << val << std::endl ;
		pi = pi+val;
		MPI_Send(&pi, 1, MPI_DOUBLE, 1, tag, globComm); }
	else if(rank == nbp-1){
		std::cout << "Tâche maitresse: Numéro (" << rank << ") ";
		val = approximate_pi(nbSamples);
		MPI_Recv(&pi, 1, MPI_DOUBLE, nbp-2, tag, globComm, &Stat);
		std::cout << "Valeur obtenue: " << val << std::endl;
		pi = pi+val;
		std::cout << "Approximation finale de pi: " << pi/nbp << std::endl; }
	else{
		std::cout << "Tâche numéro ("<< rank << ") ";
		double val = approximate_pi(nbSamples);
		std::cout << "Valeur obtenue: " << val << std::endl;
		MPI_Recv(&pi, 1, MPI_DOUBLE, rank-1, tag, globComm, &Stat);
		pi = pi+val;
		MPI_Send(&pi, 1, MPI_DOUBLE, rank+1, tag, globComm);
	}
	
	MPI_Finalize();
	return EXIT_SUCCESS;
	
}

