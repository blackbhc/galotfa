#include "barinfo.hpp"
#include <cmath>
#include <mpi.h>
using namespace std;

namespace otf {

/**
 * @brief Calculate the A0 Fourier coefficient
 *
 * @param partNum particle number
 * @param mass masses of partciles
 * @return the A0 value
 */
auto bar_info::A0( const unsigned partNum, const double* masses ) -> double
{
    double A0sum = 0;
    for ( auto i = 0U; i < partNum; ++i )
    {
        A0sum += masses[ i ];
    }
    MPI_Allreduce( MPI_IN_PLACE, &A0sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    return A0sum;
}

/**
 * @brief Calculate the A2 Fourier coefficient
 *
 * @param partNum particle number
 * @param mass masses of partciles
 * @param phi azimuthal angle of the particles
 * @return the A2 value
 */
auto bar_info::A2( const unsigned partNum, const double* masses, const double* phis ) -> double
{
    double A2sumRe = 0;  // real part
    double A2sumIm = 0;  // imaginary part
    for ( auto i = 0U; i < partNum; ++i )
    {
        A2sumRe += masses[ i ] * cos( 2 * phis[ i ] );
        A2sumIm += masses[ i ] * sin( 2 * phis[ i ] );
    }
    MPI_Allreduce( MPI_IN_PLACE, &A2sumRe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( MPI_IN_PLACE, &A2sumIm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    return sqrt( A2sumRe * A2sumRe + A2sumIm * A2sumIm );
}

/**
 * @brief Calculate the bar angle as the phase angle of the m=2 Fourier mode
 *
 * @param partNum particle number
 * @param mass masses of partciles
 * @param phi azimuthal angle of the particles
 * @return the bar angle
 */
auto bar_info::bar_angle( const unsigned partNum, const double* masses,
                          const double* phis ) -> double
{
    double A2sumRe = 0;  // real part
    double A2sumIm = 0;  // imaginary part
    for ( auto i = 0U; i < partNum; ++i )
    {
        A2sumRe += masses[ i ] * cos( 2 * phis[ i ] );
        A2sumIm += masses[ i ] * sin( 2 * phis[ i ] );
    }
    // MPI reduce
    MPI_Allreduce( MPI_IN_PLACE, &A2sumRe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( MPI_IN_PLACE, &A2sumIm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    return atan2( A2sumIm, A2sumRe ) / 2;
}

/**
 * @brief Calculate the buckling strength parameter.
 *
 * @param partNum particle number
 * @param mass masses of partciles
 * @param phi azimuthal angle of the particles
 * @param zed z coordinates of the particles
 * @return the value of buckling strength
 */
auto bar_info::Sbuckle( const unsigned partNum, const double* masses, const double* phis,
                        const double* zeds ) -> double
{
    const double A0value     = A0( partNum, masses );
    double       numeratorRe = 0;  // real part
    double       numeratorIm = 0;  // imaginary part
    for ( auto i = 0U; i < partNum; ++i )
    {
        numeratorRe += masses[ i ] * zeds[ i ] * cos( 2 * phis[ i ] );
        numeratorIm += masses[ i ] * zeds[ i ] * sin( 2 * phis[ i ] );
    }
    MPI_Allreduce( MPI_IN_PLACE, &numeratorRe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( MPI_IN_PLACE, &numeratorIm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    return sqrt( numeratorRe * numeratorRe + numeratorIm * numeratorIm ) / A0value;
}

}  // namespace otf
