#include "eigen.hpp"
#include "gsl/gsl_matrix_double.h"
#include "gsl/gsl_vector_double.h"
#include <gsl/gsl_eigen.h>
#define GSL_RANGE_CHECK_OFF  // turn off the range check in gsl
#define HAVE_INLINE          // inline function for the gsl_xxx_get APIs

/**
 * @brief Calculate the eigenvalues and eigenvectors of a 3x3 symmetric marix
 *
 * @paramin matrixData 1D array of the 3x3 matrix data, in row-major order
 * @paramout eigenValues the eigenvalues of the given array, sorted by their absolute value.
 * @param eigenVectors 1D array of the 3x3 eigenmatrix of the corresponding eigenvectors, in
 * row-major order and each column is an eigenvector.
 * @return
 */
void eigen::eigens_sym_33( double matrixData[ matrixDim ], double eigenValues[ vecDim ],
                           double eigenVectors[ matrixDim ] )
{
    gsl_matrix_view matrixView = gsl_matrix_view_array( matrixData, vecDim, vecDim );

    // allocate the eigenvalues and the eigenvectors
    gsl_vector* eval = gsl_vector_alloc( vecDim );
    gsl_matrix* evec = gsl_matrix_alloc( vecDim, vecDim );

    gsl_eigen_symmv_workspace* workSpace = gsl_eigen_symmv_alloc( vecDim );
    gsl_eigen_symmv( &matrixView.matrix, eval, evec, workSpace );
    gsl_eigen_symmv_free( workSpace );

    // sort the eigenvalues in the ascending order of magnitude
    gsl_eigen_symmv_sort( eval, evec, GSL_EIGEN_SORT_ABS_ASC );

    for ( int i = 0; i < 3; ++i )
    {
        eigenValues[ i ]          = gsl_vector_get( eval, i );
        eigenVectors[ i * 3 + 0 ] = gsl_matrix_get( evec, i, 0 );
        eigenVectors[ i * 3 + 1 ] = gsl_matrix_get( evec, i, 1 );
        eigenVectors[ i * 3 + 2 ] = gsl_matrix_get( evec, i, 2 );
    }

    gsl_vector_free( eval );
    gsl_matrix_free( evec );
}
