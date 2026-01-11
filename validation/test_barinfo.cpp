/**
 * @file test_barinfo.cpp
 * @brief Test the bar info APIs
 */

#define DEBUG 1
#include "../include/barinfo.hpp"
#ifdef DEBUG
#include "../include/myprompt.hpp"
#endif
#include <cassert>
#include <cmath>
#include <cstring>
#include <mpi.h>
using namespace std;
using namespace otf;
#define THRESHOLD 1e-6  // the equal threshold of floating numbers

int main(int argc, char* argv[])
{
    auto returnCode = 0;
    auto floatEq    = [](double a, double b) -> bool {
        return abs(a - b) <= THRESHOLD;
    };

    // NOTE: test in rank 4 mpi process program: 10 coordinate in each rank
    int rank = -1, size = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(size == 4);  // check the mpi size

    double massSend[40] = {
        0.58801452, 0.69910875, 0.18815196, 0.04380856, 0.20501895, 0.10606287,
        0.72724014, 0.67940052, 0.4738457,  0.44829582, 0.01910695, 0.75259834,
        0.60244854, 0.96177758, 0.66436865, 0.60662962, 0.44915131, 0.22535416,
        0.6701743,  0.73576659, 0.25799564, 0.09554215, 0.96090974, 0.25176729,
        0.28216512, 0.76825393, 0.7979234,  0.5440372,  0.38270763, 0.38165095,
        0.28582739, 0.74026815, 0.23898683, 0.4377217,  0.8835387,  0.28928114,
        0.78450686, 0.75895366, 0.41778538, 0.22576877};
    double phiSend[40] = {
        0.,         0.16110732, 0.32221463, 0.48332195, 0.64442926, 0.80553658,
        0.96664389, 1.12775121, 1.28885852, 1.44996584, 1.61107316, 1.77218047,
        1.93328779, 2.0943951,  2.25550242, 2.41660973, 2.57771705, 2.73882436,
        2.89993168, 3.061039,   3.22214631, 3.38325363, 3.54436094, 3.70546826,
        3.86657557, 4.02768289, 4.1887902,  4.34989752, 4.51100484, 4.67211215,
        4.83321947, 4.99432678, 5.1554341,  5.31654141, 5.47764873, 5.63875604,
        5.79986336, 5.96097068, 6.12207799, 6.28318531};
    double zedSend[40] = {
        42.00981389, 6.43636909,  59.64326869, 83.73237223, 89.24863863,
        20.05274439, 50.23952344, 89.53818445, 25.59209314, 86.7232343,
        1.64879348,  55.24969544, 52.79053865, 92.33503869, 24.59484435,
        6.40183762,  90.21047046, 87.40398004, 16.36672906, 99.97413067,
        34.68039703, 31.28781593, 84.7104021,  88.02311026, 67.65586515,
        5.36751543,  55.92137735, 69.45129418, 82.41973027, 31.14286589,
        50.52305408, 84.90037879, 29.35156327, 67.71195507, 42.09064022,
        68.17127136, 22.12279895, 54.89976984, 84.88467195, 73.65669013};
    double massRecv[10] = {0};
    double phiRecv[10]  = {0};
    double zedRecv[10]  = {0};

    MPI_Scatter(massSend, 10, MPI_DOUBLE, massRecv, 10, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    MPI_Scatter(phiSend, 10, MPI_DOUBLE, phiRecv, 10, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    MPI_Scatter(zedSend, 10, MPI_DOUBLE, zedRecv, 10, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);

    /* // check the data are correctly scattered
    for ( int j = 0; j < 4; ++j )
    {
        if ( j == rank )
        {
            print( "Rank [%d]:", rank );
            for ( int i = 0; i < 10; i++ )
            {
                print( "%lf", massRecv[ i ] );
            }
        }
        MPI_Barrier( MPI_COMM_WORLD );
    } */

    double A0target      = 19.631915485490538;
    double getA0         = bar_info::A0(10, massRecv);
    double SbarTarget    = 0.08786394735228277;
    auto   A2            = bar_info::A2(10, massRecv, phiRecv);
    double getSbar       = A2 / getA0;
    double SbuckleTarget = 4.279046929538018;
    double getSbuckle    = bar_info::Sbuckle(10, massRecv, phiRecv, zedRecv);

    if (not floatEq(getA0, A0target))
    {
        MPI_ERROR(rank, "A0: Target is [%lf] but get [%lf].", A0target, getA0);
        returnCode += 1;
    }

    if (not floatEq(getSbar, SbarTarget))
    {
        MPI_ERROR(rank, "Sbar: Target is [%lf] but get [%lf].", SbarTarget,
                  getSbar);
        returnCode += 1;
    }

    if (not floatEq(getSbuckle, SbuckleTarget))
    {
        MPI_ERROR(rank, "Sbuckle: Target is [%lf] but get [%lf].",
                  SbuckleTarget, getSbuckle);
        returnCode += 1;
    }

    MPI_Finalize();
    return returnCode;
}
/* python codes:
import numpy as np
np.random.seed(2024)

mass = np.random.rand(40)
zs = np.random.rand(40) * 100
phis = np.linspace(0, np.pi*2, 40)

A0 = mass.sum()
A2 = np.abs(np.sum(mass * np.exp(phis*2j)))
Sbar = A2/A0
Sbuckle = np.abs(np.sum(mass * zs * np.exp(phis*2j))) / A0
*/
