# include <iostream>
# include <cstdlib>
# include <mpi.h>

int main( int nargs, char* argv[] )
{
	
	int next , prev , val1, val2, tag=1;
	MPI_Status stats;

	MPI_Init( &nargs, &argv );
	MPI_Comm globComm;
	MPI_Comm_dup(MPI_COMM_WORLD, &globComm);
	int nbp;
	MPI_Comm_size(globComm, &nbp);
	int rank;
	MPI_Comm_rank(globComm, &rank);

	prev = (rank -1)%nbp;
	next = (rank+1)%nbp;
	
	val1= rank*rank;
	std::cout << "Numéro (" << rank << "):" << "J'envoie " << val1 << std::endl; 
	MPI_Send(&val1, 1, MPI_INT, next, tag, globComm);
	MPI_Recv(&val2, 1, MPI_INT, prev, tag, globComm, &stats);
    val2 = val2+1;
	std::cout << "Numéro (" << rank << "):";
	std::cout << "J'ai reçu " << val2 << std::endl; 
	

	MPI_Finalize();
	return EXIT_SUCCESS;
}
