/**
 * @file test_bin2d.cpp
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

    double        xmin = -0.1, xmax = 10.7, ymin = -5.5, ymax = 5.9;
    unsigned long xBinNum = 7, yBinNum = 5, dataNum = 100;
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
    double ys[100] = {
        0.31732592,  3.70958621,  3.41090268,  4.7205554,   2.8225721,
        -3.02969494, 1.10626066,  -0.21144494, 1.16636996,  -3.60066757,
        -0.88764184, 2.77630339,  4.3972552,   -3.95420592, 4.38482197,
        2.97387172,  -1.69197281, -1.88214249, -2.0984618,  -3.26110409,
        3.81149556,  3.20056487,  2.36658878,  3.46968231,  1.09115624,
        -1.55766994, -2.70331009, 4.44625788,  -2.08274289, -0.89955414,
        -1.62335725, 4.49328596,  -2.69071067, -3.30527961, -0.09102368,
        4.04753448,  3.83999687,  -3.05088407, 0.20217967,  -3.25100676,
        3.99545879,  -1.21800666, 2.02881283,  0.11745452,  1.34851061,
        4.33994941,  -2.89509778, -1.63795992, 1.59467182,  -0.85396637,
        -2.4468656,  3.75602249,  -0.98208691, 0.61367929,  -3.28316903,
        -4.91222065, -0.91211284, -3.96531115, 0.25046826,  -0.4013903,
        -0.88475794, 2.5797616,   -4.82216572, 1.12867578,  1.51142426,
        2.22066801,  0.73585438,  4.96857107,  1.10051551,  0.02362355,
        -1.46175665, -3.56605767, 2.63367422,  3.66574752,  -0.33184904,
        4.42086785,  1.94197465,  2.21186357,  1.03483581,  3.7609812,
        3.814588,    4.13824917,  -4.58823095, 3.00779729,  4.91191867,
        -3.73597238, -3.56974086, -4.33168615, -1.59472529, 2.96316388,
        0.66485496,  2.41176291,  -2.42765724, 1.09575428,  4.39439284,
        2.39667158,  -1.353682,   -1.51282724, -3.66469105, 3.43455826};
    double values[100] = {
        42.7620314,  78.14583606, 18.9724537,  0.3767982,   88.02606778,
        88.04640679, 58.21273656, 10.01290842, 51.22482089, 48.86088709,
        49.05773617, 75.49584767, 92.30034393, 33.30440777, 61.64161585,
        99.64066784, 78.44932456, 86.67834718, 94.65055783, 43.59819161,
        9.58963269,  75.03290737, 80.3815177,  33.76788689, 15.99791267,
        86.48828747, 29.21358913, 15.22476098, 11.11208716, 20.39822981,
        71.13443183, 93.38274409, 4.02540923,  32.87191043, 87.2596159,
        43.44727151, 71.27224976, 48.982064,   30.79856241, 70.94975755,
        60.27478453, 45.17272994, 40.3284794,  45.67048296, 61.05176795,
        55.62012677, 94.22231418, 74.81408527, 86.44826938, 19.00776983,
        2.78977403,  9.90092896,  84.27080459, 35.34296751, 47.79382155,
        50.72157251, 17.50349893, 66.92980461, 31.55633508, 66.52874536,
        24.68204224, 91.34630229, 25.72381727, 37.50175673, 40.17362384,
        22.04980775, 67.72909211, 43.83255821, 93.36336195, 38.80041004,
        29.04259212, 28.01182603, 66.02144988, 11.19386939, 16.36478308,
        30.15467584, 29.14997735, 2.35583521,  11.90600575, 70.64064622,
        23.87837643, 37.86730729, 17.97593933, 26.40806644, 61.39795082,
        13.56619533, 77.24610208, 11.0063208,  66.29359616, 70.02013443,
        51.37024849, 28.10358977, 41.5205741,  8.15149347,  61.76688924,
        86.40467441, 3.86381641,  75.33962543, 96.851591,   79.21758306};
    double xsRecv[100] = {0}, ysRecv[100] = {0}, valuesRecv[100] = {0};
    for (auto i = 0; i < 100; ++i)
    {
        xsRecv[i]     = nan("");
        ysRecv[i]     = nan("");
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
    MPI_Scatter(ys, 25, MPI_DOUBLE, ysRecv + rank * 25, 25, MPI_DOUBLE, 0,
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

    double targetCount[] = {2., 3., 1., 2., 1., 3., 3., 2., 9., 4., 3., 2.,
                            7., 1., 2., 0., 4., 3., 3., 5., 2., 3., 2., 4.,
                            7., 4., 5., 7., 1., 1., 1., 0., 2., 1., 0.};
    double targetSum[]   = {
        1.27967675e+02, 1.36008911e+02, 4.90577362e+01, 9.70827151e+01,
        3.76798201e-01, 1.29749900e+02, 1.61838188e+02, 4.75542478e+01,
        4.98031889e+02, 1.70423920e+02, 9.52989929e+01, 8.95614117e+01,
        1.90565770e+02, 9.13463023e+01, 9.81420918e+01, 0.00000000e+00,
        2.82875336e+02, 1.61715944e+02, 1.42324982e+02, 2.02702361e+02,
        1.40449783e+02, 2.09926218e+02, 6.82256450e+01, 2.30177736e+02,
        3.88991946e+02, 1.38641387e+02, 2.35900352e+02, 3.12212492e+02,
        6.10517680e+01, 7.12722498e+01, 3.33044078e+01, 0.00000000e+00,
        1.01871713e+02, 8.03815177e+01, 0.00000000e+00};
    double targetMean[] = {
        63.9838373,  45.33630359, 49.05773617, 48.54135756, 0.3767982,
        43.24996664, 53.94606275, 23.77712388, 55.33687655, 42.60598006,
        31.76633095, 44.78070586, 27.22368136, 91.34630229, 49.07104591,
        nan(""),     70.71883408, 53.9053148,  47.44166082, 40.54047222,
        70.2248913,  69.975406,   34.11282249, 57.54443392, 55.57027795,
        34.66034681, 47.18007048, 44.60178458, 61.05176795, 71.27224976,
        33.30440777, nan(""),     50.93585644, 80.3815177,  nan("")};

    double targetStd[] = {
        13.26226478, 34.80606676, 0.,          26.49154981, 0.,
        24.68182686, 35.86509038, 7.77921121,  26.3455048,  20.16987494,
        14.43018933, 33.6686187,  13.44070506, 0.,          11.20373862,
        nan(""),     24.9716356,  10.36749137, 40.03359154, 31.53224069,
        26.62669969, 18.82497676, 24.09991407, 24.10225655, 26.60936488,
        18.9983667,  30.13269624, 30.95461918, 0.,          0.,
        0.,          nan(""),     15.59288893, 0.,          nan("")};

    // test the C++ codes
    auto count =
        statistic::bin2d(rank, xsRecv, xmin, xmax, xBinNum, ysRecv, ymin, ymax,
                         yBinNum, statistic_method::COUNT, dataNum);
    auto sum =
        statistic::bin2d(rank, xsRecv, xmin, xmax, xBinNum, ysRecv, ymin, ymax,
                         yBinNum, statistic_method::SUM, dataNum, values);
    auto mean =
        statistic::bin2d(rank, xsRecv, xmin, xmax, xBinNum, ysRecv, ymin, ymax,
                         yBinNum, statistic_method::MEAN, dataNum, values);
    auto std =
        statistic::bin2d(rank, xsRecv, xmin, xmax, xBinNum, ysRecv, ymin, ymax,
                         yBinNum, statistic_method::STD, dataNum, values);
    if (rank == 0)  // check the results in the root process
    {
        cout << "results of count:" << endl;
        for (auto i = 0UL; i < xBinNum; ++i)
        {
            for (unsigned long j = 0; j < yBinNum; ++j)
                cout << count[i * yBinNum + j] << " ";
            cout << endl;
        }
        for (auto i = 0UL; i < xBinNum * yBinNum; ++i)
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
        for (auto i = 0UL; i < xBinNum; ++i)
        {
            for (unsigned long j = 0; j < yBinNum; ++j)
                cout << sum[i * yBinNum + j] << " ";
            cout << endl;
        }
        for (auto i = 0UL; i < xBinNum * yBinNum; ++i)
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
        for (auto i = 0UL; i < xBinNum; ++i)
        {
            for (auto j = 0UL; j < yBinNum; ++j)
                cout << mean[i * yBinNum + j] << " ";
            cout << endl;
        }
        for (auto i = 0UL; i < xBinNum * yBinNum; ++i)
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
        for (auto i = 0UL; i < xBinNum; ++i)
        {
            for (auto j = 0UL; j < yBinNum; ++j)
                cout << std[i * yBinNum + j] << " ";
            cout << endl;
        }
        for (auto i = 0UL; i < xBinNum * yBinNum; ++i)
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
from scipy.stats import binned_statistic_2d as bin2d
np.random.seed(2024)
xs = np.random.rand(100)*10
ys = np.random.rand(100)*10 - 5
values= np.random.rand(100)*100
bin2d(x=xs, y=ys, values=values, range=[[-0.1, 10.7], [-5.5, 5.9]], bins=[7, 5],
statistic="count")[0].reshape(7*5)

bin2d(x=xs, y=ys, values=values, range=[[-0.1, 10.7], [-5.5, 5.9]], bins=[7, 5],
statistic="sum")[0].reshape(7*5)

bin2d(x=xs, y=ys, values=values, range=[[-0.1, 10.7], [-5.5, 5.9]], bins=[7, 5],
statistic="mean")[0].reshape(7*5)

bin2d(x=xs, y=ys, values=values, range=[[-0.1, 10.7], [-5.5, 5.9]], bins=[7, 5],
statistic="std")[0].reshape(7*5)
*/
