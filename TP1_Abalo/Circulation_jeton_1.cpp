# include <iostream>
# include <cstdlib>
# include <mpi.h>

int main( int nargs, char* argv[] )
{
	
	int next , prev , val, tag=1;
	MPI_Status stats;

	MPI_Init( &nargs, &argv );
	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
	int nbp;
	MPI_Comm_size(globComm, &nbp);
	int rank;
	MPI_Comm_rank(globComm, &rank);

	prev = rank -1;
	next = rank+1;
	
	if (rank == 0) {
		prev = nbp - 1;
		val = 1;
		std::cout << "Numéro (" << rank << "):" << "J'envoie " << val << std::endl; 
		MPI_Send(&val, 1, MPI_INT, next, tag, globComm);
		MPI_Recv(&val, 1, MPI_INT, prev, tag, globComm, &stats);}
	else { 
		MPI_Recv(&val, 1, MPI_INT, prev, tag, globComm, &stats);
		std::cout << "Numéro (" << rank << "):";
		std::cout << "J'ai reçu " <<val << " et j'envoie " <<val+1 << std::endl; 
	    val =val+1;
		MPI_Send(&val, 1, MPI_INT, (next)%nbp, tag, globComm);
	}

	MPI_Finalize();
	return EXIT_SUCCESS;
}
