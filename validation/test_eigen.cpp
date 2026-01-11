/**
 * @file test_eigen.cpp
 * @brief Test the eigenvalue, eigenvector calculations.
 */

#define DEBUG 1
#include "../include/eigen.hpp"
#include <iostream>
using namespace std;

void print_matrix(double (&matrix)[9])
{
    for (int i = 0; i < 3; ++i)
        cout << "[" << matrix[i * 3 + 0] << ", " << matrix[i * 3 + 1] << ", "
             << matrix[i * 3 + 2] << "]," << endl;
}

void calAndPrint(double (&matrix)[9])
{
    eigen  eigenCalculator;
    double values[3], vectors[9];
    eigenCalculator.eigens_sym_33(matrix, values, vectors);
    cout << "Eigenvalues:" << endl;
    for (int i = 0; i < 3; ++i)
        cout << " " << values[i];
    cout << endl;

    cout << "Eigenmatrix:" << endl;
    print_matrix(vectors);
}

int main()
{
    double matrix1[] = {1, 2, 3, 2, 4, 5, 3, 5, 6};
    print_matrix(matrix1);
    calAndPrint(matrix1);
    double matrix2[] = {7, 3, -7, 3, 11, 2, -7, 2, 5};
    print_matrix(matrix2);
    calAndPrint(matrix2);
    return 0;
}
