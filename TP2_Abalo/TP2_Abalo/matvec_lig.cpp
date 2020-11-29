// Produit matrice-vecteur, partition par lignes
# include <cassert>
# include <vector>
# include <iostream>
#include <mpi.h>
#include <stdio.h>

// ---------------------------------------------------------------------
class Matrix : public std::vector<double>
{
public:
    Matrix (int dim);
    Matrix( int nrows, int ncols );
    Matrix( const Matrix& A ) = delete;
    Matrix( Matrix&& A ) = default;
    ~Matrix() = default;

    Matrix& operator = ( const Matrix& A ) = delete;
    Matrix& operator = ( Matrix&& A ) = default;
    
    double& operator () ( int i, int j ) {
        return m_arr_coefs[i + j*m_nrows];
    }
    double  operator () ( int i, int j ) const {
        return m_arr_coefs[i + j*m_nrows];
    }
    
    std::vector<double> operator * ( const std::vector<double>& u ) const;
    std::vector<double> operator + ( const std::vector<double>& u ) const;

    std::ostream& print( std::ostream& out ) const
    {
        const Matrix& A = *this;
        out << "[\n";
        for ( int i = 0; i < m_nrows; ++i ) {
            out << " [ ";
            for ( int j = 0; j < m_ncols; ++j ) {
                out << A(i,j) << " ";
            }
            out << " ]\n";
        }
        out << "]";
        return out;
    }
private:
    int m_nrows, m_ncols;
    std::vector<double> m_arr_coefs;
};
// ---------------------------------------------------------------------
inline std::ostream& 
operator << ( std::ostream& out, const Matrix& A )
{
    return A.print(out);
}
// ---------------------------------------------------------------------
inline std::ostream&
operator << ( std::ostream& out, const std::vector<double>& u )
{
    out << "[ ";
    for ( const auto& x : u )
        out << x << " ";
    out << " ]";
    return out;
}
// ---------------------------------------------------------------------
std::vector<double> 
Matrix::operator * ( const std::vector<double>& u ) const
{
    const Matrix& A = *this;
    assert( u.size() == unsigned(m_ncols) );
    std::vector<double> v(m_nrows, 0.);
    for ( int i = 0; i < m_nrows; ++i ) {
        for ( int j = 0; j < m_ncols; ++j ) {
            v[i] += A(i,j)*u[j];
        }            
    }
    return v;
}

// =====================================================================
Matrix::Matrix (int dim) : m_nrows(dim), m_ncols(dim),
                           m_arr_coefs(dim*dim)
{
    for ( int i = 0; i < dim; ++ i ) {
        for ( int j = 0; j < dim; ++j ) {
            (*this)(i,j) = (i+j)%dim;
        }
    }
}
// ---------------------------------------------------------------------
Matrix::Matrix( int nrows, int ncols ) : m_nrows(nrows), m_ncols(ncols),
                                         m_arr_coefs(nrows*ncols)
{
    int dim = (nrows > ncols ? nrows : ncols );
    for ( int i = 0; i < nrows; ++ i ) {
        for ( int j = 0; j < ncols; ++j ) {
            (*this)(i,j) = (i+j)%dim;
        }
    }    
}
// =====================================================================
int main(int argc , char *argv [] ) 
{
   
    MPI_Init(&argc ,&argv ) ;
    int nbp;
	MPI_Comm_size(MPI_COMM_WORLD, &nbp);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int N = 4;
    int N_loc = N/nbp;
    int imin_loc = N_loc*rank;
    std::vector<double> recvbuf (N);

    Matrix A(N);
    Matrix A_loc(N_loc,N);
    for ( int i = 0; i < N_loc; ++i ) {
        for (int j = 0; j < N; ++j) {
            A_loc(i,j) = A(i+imin_loc,j);
        }
    }

    std::vector<double> u( N );
    for (int i = 0; i < N; ++i) u[i] = i+1;  

    //Multiplication partielle pour chaque processus
    std::vector<double> v = A_loc*u;
    std::cout << "Résultat partiel : " << v << "\n" << std::endl;
    // on réunit les vecteurs obtenus  
    MPI_Allgather (v.data(), N_loc, MPI_DOUBLE, recvbuf.data(), N_loc, MPI_DOUBLE, MPI_COMM_WORLD );
    printf ("rank= %d  Results : %f %f %f %f \n" , rank , recvbuf [0] ,recvbuf [1] , recvbuf [2] , recvbuf [3] ) ;

    MPI_Finalize () ;
    return EXIT_SUCCESS;
}