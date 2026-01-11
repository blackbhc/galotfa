/**
 * @file test_bin1d.cpp
 * @brief Compare the results with the output of python
 */

#define DEBUG 1
#define THRESHOLD 1e-6  // the equal threshold of floating numbers
#include "../include/myprompt.hpp"
#include "../include/statistic.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <mpi.h>
using namespace std;

int main(int argc, char* argv[])
{
    // NOTE: test in rank 4 mpi process program
    int rank = -1;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double        xmin = -0.1, xmax = 10.7;
    unsigned long binNum = 13, dataNum = 100;
    // The following datas and results are generated with the python codes at
    // the end of this file
    double        xs[100] = {
        5.88014519, 6.99108748, 1.8815196,  0.43808564, 2.05018952, 1.06062874,
        7.27240144, 6.79400524, 4.73845703, 4.48295824, 0.19106948, 7.52598337,
        6.02448539, 9.61777575, 6.64368647, 6.06629619, 4.49151315, 2.25354163,
        6.70174297, 7.35766592, 2.57995638, 0.95542154, 9.60909742, 2.51767287,
        2.82165119, 7.68253935, 7.97923397, 5.44037198, 3.82707631, 3.8165095,
        2.85827388, 7.40268153, 2.38986832, 4.37721705, 8.83538703, 2.8928114,
        7.84506857, 7.58953657, 4.17785385, 2.25768768, 4.20098139, 0.64363691,
        5.96432687, 8.37323722, 8.92486386, 2.00527444, 5.02395234, 8.95381844,
        2.55920931, 8.67232343, 0.16487935, 5.52496954, 5.27905386, 9.23350387,
        2.45948444, 0.64018376, 9.02104705, 8.740398,   1.63667291, 9.99741307,
        3.4680397,  3.12878159, 8.47104021, 8.80231103, 6.76558652, 0.53675154,
        5.59213774, 6.94512942, 8.24197303, 3.11428659, 5.05230541, 8.49003788,
        2.93515633, 6.77119551, 4.20906402, 6.81712714, 2.21227989, 5.48997698,
        8.48846719, 7.36566901, 4.99622586, 3.79664993, 7.87520811, 1.68869308,
        5.86358614, 4.31210671, 0.61910185, 2.89454775, 7.34145398, 2.88655455,
        3.90398114, 6.35617324, 8.31148857, 3.19421001, 1.59224788, 7.11664216,
        8.72708643, 5.93156365, 6.94712879, 1.7323332};
    double values[100] = {
        53.17325923, 87.09586207, 84.10902681, 97.205554,   78.22572097,
        19.70305064, 61.06260656, 47.88555059, 61.66369961, 13.99332433,
        41.12358162, 77.76303394, 93.97255196, 10.45794082, 93.84821973,
        79.73871717, 33.08027187, 31.17857505, 29.015382,   17.38895907,
        88.1149556,  82.00564872, 73.66588776, 84.69682305, 60.91156238,
        34.42330056, 22.96689914, 94.46257885, 29.17257109, 41.0044586,
        33.76642754, 94.93285956, 23.09289326, 16.94720391, 49.08976324,
        90.47534478, 88.39996873, 19.49115928, 52.02179666, 17.48993239,
        89.95458788, 37.81993342, 70.28812831, 51.17454521, 63.48510613,
        93.39949409, 21.04902224, 33.62040076, 65.94671821, 41.46033628,
        25.53134401, 87.56022493, 40.17913092, 56.13679289, 17.16830974,
        0.87779352,  40.87887158, 10.34688852, 52.5046826,  45.98609702,
        41.15242062, 75.79761599, 1.77834276,  61.28675777, 65.11424264,
        72.20668012, 57.35854382, 99.68571069, 61.00515511, 50.23623555,
        35.38243353, 14.33942325, 76.33674223, 86.65747519, 46.68150964,
        94.20867848, 69.41974645, 72.11863575, 60.34835808, 87.60981203,
        88.14587996, 91.38249169, 4.11769049,  80.07797289, 99.11918671,
        12.64027624, 14.30259137, 6.68313846,  34.05274713, 79.63163875,
        56.64854964, 74.11762909, 25.72342763, 60.95754285, 93.94392841,
        73.9667158,  36.46318003, 34.87172757, 13.35308954, 84.34558265};
    double xsRecv[100] = {0}, valuesRecv[100] = {0};
    for (auto i = 0; i < 100; ++i)
    {
        xsRecv[i]     = nan("");
        valuesRecv[i] = nan("");
    }
    mpi_print(rank, "Data before MPI_Scatter:");
    for (int i = 0; i < 4; ++i)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == i)
        {
            myprint("Values in rank %d:", rank);
            for (int j = 0; j < 100; ++j)
            {
                printf("%lf ", valuesRecv[j]);
            }
            myprint("\nEnd of values.");

            myprint("xsRecv in rank %d:", rank);
            for (int j = 0; j < 100; ++j)
            {
                printf("%lf ", xsRecv[j]);
            }
            myprint("\nEnd of xsRecv.");
        }
    }
    // Scatter the data to 4 ranks
    MPI_Scatter(xs, 25, MPI_DOUBLE, xsRecv + rank * 25, 25, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    MPI_Scatter(values, 25, MPI_DOUBLE, valuesRecv + rank * 25, 25, MPI_DOUBLE,
                0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    mpi_print(rank, "Data after MPI_Scatter:");
    for (int i = 0; i < 4; ++i)
    {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == i)
        {
            myprint("Values in rank %d:", rank);
            for (int j = 0; j < 100; ++j)
            {
                printf("%lf ", valuesRecv[j]);
            }
            myprint("\nEnd of values.");

            myprint("xsRecv in rank %d:", rank);
            for (int j = 0; j < 100; ++j)
            {
                printf("%lf ", xsRecv[j]);
            }
            myprint("\nEnd of xsRecv.");
        }
    }
    // MPI_Finalize();
    // return 0;

    double targetCount[] = {7., 2.,  11., 13., 5., 8., 8.,
                            7., 14., 7.,  14., 3., 1.};
    double targetSum[]   = {
        289.06747805, 101.70869936, 707.78755557, 790.72305512, 259.36049164,
        326.98267013, 496.25644999, 505.28120004, 890.94505152, 342.09491171,
        551.00055635, 140.26062147, 45.98609702};
    double targetMean[] = {41.29535401, 50.85434968, 64.34432323, 60.82485039,
                           51.87209833, 40.87283377, 62.03205625, 72.18302858,
                           63.63893225, 48.87070167, 39.3571826,  46.75354049,
                           45.98609702};

    double targetStd[] = {30.94756174, 31.15129904, 27.15978424, 25.89924962,
                          21.59708911, 25.37731245, 25.92368752, 20.76489556,
                          29.16292821, 34.36105793, 19.49406309, 26.64388871,
                          0.};

    // test the C++ codes
    auto count = statistic::bin1d(rank, xsRecv, xmin, xmax, binNum,
                                  statistic_method::COUNT, dataNum);
    auto sum   = statistic::bin1d(rank, xsRecv, xmin, xmax, binNum,
                                  statistic_method::SUM, dataNum, values);
    auto mean  = statistic::bin1d(rank, xsRecv, xmin, xmax, binNum,
                                  statistic_method::MEAN, dataNum, values);
    auto std   = statistic::bin1d(rank, xsRecv, xmin, xmax, binNum,
                                  statistic_method::STD, dataNum, values);
    if (rank == 0)  // check the results in the root process
    {
        cout << "results of count:" << endl;
        for (auto i = 0UL; i < binNum; ++i)
        {
            cout << count[i] << " ";
        }
        cout << endl;
        for (auto i = 0UL; i < binNum; ++i)
        {
            if (abs(count[i] - targetCount[i]) >= THRESHOLD)
            {
                cout << "Target value=" << targetCount[i] << endl;
                cout << "Get value=" << count[i] << endl;
                MPI_Finalize();
                return -1;
            }
        }

        cout << "results of sum:" << endl;
        for (auto i = 0UL; i < binNum; ++i)
        {
            cout << sum[i] << " ";
        }
        cout << endl;
        for (auto i = 0UL; i < binNum; ++i)
        {
            if (abs(sum[i] - targetSum[i]) >= THRESHOLD)
            {
                cout << "Target value=" << targetSum[i] << endl;
                cout << "Get value=" << sum[i] << endl;
                MPI_Finalize();
                return -1;
            }
        }

        cout << "results of mean:" << endl;
        for (auto i = 0UL; i < binNum; ++i)
        {
            cout << mean[i] << " ";
        }
        cout << endl;
        for (auto i = 0UL; i < binNum; ++i)
        {
            if (abs(mean[i] - targetMean[i]) >= THRESHOLD)
            {
                cout << "Target value=" << targetMean[i] << endl;
                cout << "Get value=" << mean[i] << endl;
                MPI_Finalize();
                return -1;
            }
        }

        cout << "results of std:" << endl;
        for (auto i = 0UL; i < binNum; ++i)
        {
            cout << std[i] << " ";
        }
        cout << endl;
        for (auto i = 0UL; i < binNum; ++i)
        {
            if (abs(std[i] - targetStd[i]) >= THRESHOLD)
            {
                cout << "Target value=" << targetStd[i] << endl;
                cout << "Get value=" << std[i] << endl;
                MPI_Finalize();
                return -1;
            }
        }
    }
    MPI_Finalize();
    return 0;
}
/* python:
import numpy as np
from scipy.stats import binned_statistic as bin1d
np.random.seed(2024)
xs = np.random.rand(100)*10
values= np.random.rand(100)*100
bin1d(x=xs, values=values, range=[-0.1, 10.7], bins=13,
statistic="count")[0]

bin1d(x=xs, values=values, range=[-0.1, 10.7], bins=13,
statistic="sum")[0]

bin1d(x=xs, values=values, range=[-0.1, 10.7], bins=13,
statistic="mean")[0]

bin1d(x=xs, values=values, range=[-0.1, 10.7], bins=13,
statistic="std")[0]
*/
