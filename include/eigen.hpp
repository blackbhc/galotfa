/**
 * @file eigen.hpp
 * @brief A class to organize the eigenvalue and eigenvector calculations for
 * 3x3 matrices.
 */

#pragma once

/**
 * @class eigen
 * @brief Wrapper class of the eigen-system.
 *
 */
class eigen
{
private:
    static constexpr auto vecDim    = 3;
    static constexpr auto matrixDim = 3;

public:
    // calculate the eigenvalues and eigenvectors of a given 3x3 symmetric
    // matrix.
    static void eigens_sym_33(double matrixData[matrixDim],
                              double eigenValues[vecDim],
                              double eigenVectors[matrixDim]);
};
