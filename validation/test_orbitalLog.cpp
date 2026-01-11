/**
 * @file test_orbitalLog.cpp
 * @brief Test the orbital log.
 */

#define DEBUG 1
#include "../include/monitor.hpp"
#include <cassert>
#include <cmath>
#include <memory>
#include <mpi.h>
using namespace std;
using namespace otf;

int main(int argc, char* argv[])
{
    // NOTE: test in rank 4 mpi process program: 10 coordinate in each rank
    int rank = -1;
    int size = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    assert(size == 4);  // check the mpi size

    int mockTypesG[40];
    for (auto& mockType : mockTypesG)
    {
        mockType = 2;
    }

    int mockIDsG[40];
    for (auto i = 0; i < 40; ++i)
    {
        mockIDsG[i] = i + 1;
    }

    double mockMassG[40];
    for (auto i = 0U; i < 40; ++i)
    {
        mockMassG[i] = ( double )i + 0.1;
    }

    double mockPotG[40];
    for (auto i = 0U; i < 40; ++i)
    {
        mockPotG[i] = -( double )i - 0.1;
    }

    double mockPosG[120];
    for (auto i = 0U; i < 120; ++i)
    {
        mockPosG[i] = ( double )i + 0.2;
    }

    double mockVelG[120];
    for (auto i = 0U; i < 120; ++i)
    {
        mockVelG[i] = pow(-1, i) * (( double )i + 0.3);
    }

    // mock there are 9, 7, 11, 13 number of particles in the 4 ranks
    int localNums[]  = {9, 7, 11, 13};
    int localNums3[] = {9, 7, 11, 13};
    for (auto& triple : localNums3)
    {
        triple *= 3;
    }

    // local data
    unique_ptr<int[]>    mockTypes(new int[localNums[rank]]());
    unique_ptr<int[]>    mockIDs(new int[localNums[rank]]());
    unique_ptr<double[]> mockMass(new double[localNums[rank]]());
    unique_ptr<double[]> mockPot(new double[localNums[rank]]());
    unique_ptr<double[]> mockPos(new double[localNums[rank] * 3]());
    unique_ptr<double[]> mockVel(new double[localNums[rank] * 3]());

    // MPI scatter
    // calculate the offsets
    int localOffset  = 0;
    int localOffset3 = 0;
    for (int i = 0; i < rank; ++i)
    {
        localOffset += localNums[i];
        localOffset3 += 3 * localNums[i];
    }
    const unique_ptr<int[]> offsets(new int[size]());
    MPI_Allgather(&localOffset, 1, MPI_INT, offsets.get(), 1, MPI_INT,
                  MPI_COMM_WORLD);
    const unique_ptr<int[]> offset3s(new int[size]());
    MPI_Allgather(&localOffset3, 1, MPI_INT, offset3s.get(), 1, MPI_INT,
                  MPI_COMM_WORLD);

    // types
    MPI_Scatterv(mockTypesG, localNums, offsets.get(), MPI_INT, mockTypes.get(),
                 localNums[rank], MPI_INT, 0, MPI_COMM_WORLD);
    // ids
    MPI_Scatterv(mockIDsG, localNums, offsets.get(), MPI_INT, mockIDs.get(),
                 localNums[rank], MPI_INT, 0, MPI_COMM_WORLD);
    // masses
    MPI_Scatterv(mockMassG, localNums, offsets.get(), MPI_DOUBLE,
                 mockMass.get(), localNums[rank], MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    // potentials
    MPI_Scatterv(mockPotG, localNums, offsets.get(), MPI_DOUBLE, mockPot.get(),
                 localNums[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // coordinates
    MPI_Scatterv(mockPosG, localNums3, offset3s.get(), MPI_DOUBLE,
                 mockPos.get(), localNums[rank] * 3, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);
    // velocities
    MPI_Scatterv(mockVelG, localNums3, offset3s.get(), MPI_DOUBLE,
                 mockVel.get(), localNums[rank] * 3, MPI_DOUBLE, 0,
                 MPI_COMM_WORLD);

    double mockTime   = 0;     // mock the time of the simulation
    double mockDeltaT = 0.13;  // mock the time step of the simulation
    double mockDrift  = 1.2;   // mock the drift
    double mockKick   = -1.7;  // mock the kick
    int    maxStep    = 23;    // mock the number of synchronized steps

    monitor otfServer("../validation/orbit_log_test.toml");
    for (auto i = 0; i < maxStep; ++i)
    {
        otfServer.main_analysis_api(
            mockTime, localNums[rank], mockIDs.get(), mockTypes.get(),
            mockMass.get(), mockPot.get(), mockPos.get(), mockVel.get());

        // mock the kick-drift pair
        for (auto j = 0; j < localNums[rank]; ++j)
            for (auto k = 0; k < 3; ++k)
            {
                mockPos[j * 3 + k] += mockDrift;
                mockVel[j * 3 + k] += mockKick;
            }
        // mock the time increment
        mockTime += mockDeltaT;
    }

    MPI_Finalize();
    return 0;
}
