#include "recenter.hpp"
#include "myprompt.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <mpi.h>
using namespace std;

namespace otf {

/**
 * @brief Calculate the center of a system.
 *
 * @param method method used to calculate the center
 * @param partNum particle number
 * @param coordinate coordinates of particles
 * @param mass masses of particles
 * @param potential potentials of particles
 * @param radius enclose radius of the region used for calculation
 * @param previousPos the position of the previous center
 * @return uniqure_ptr to the gotten center of mass
 */
auto recenter::get_center(const recenter_method method,
                          const unsigned&       partNum,
                          const double*         masses,
                          const double*         potentials,
                          const double*         coordinates,
                          const double          radius,
                          const double* previousPos) -> unique_ptr<double[]>
{
    switch (method)
    {
    case recenter_method::COM:
        return center_of_mass(masses, coordinates, partNum, radius,
                              previousPos);
        break;
    case recenter_method::MBP:
        return most_bound_particle(potentials, coordinates, partNum);
        break;
    default:
        ERROR("Get into an unexpected branch!");
        return nullptr;
    }
}

/**
 * @brief Calculate the center of mass in specified range.
 *
 * @param mass masses of particles
 * @param coordinates coordinates of particles
 * @param partNum particle number
 * @param radius enclose radius of the chosen range
 * @param previousPos the position of the previous center
 * @return uniqure_ptr to the gotten center of mass
 */
auto recenter::center_of_mass(const double*   mass,
                              const double*   coordinates,
                              const unsigned& partNum,
                              const double    radius,
                              const double* previousPos) -> unique_ptr<double[]>
{
    // results of the center of mass
    auto   com(make_unique<double[]>(3));
    // summation of mass
    double massSum = 0;
    // summation of mass[i]xcoordinates[i]
    double coordMassSum[3] = {0, 0, 0};

    // indexes for iteration
    static unsigned i = 0;
    static unsigned j = 0;
    // error vector
    static double   error[3] = {0, 0, 0};
    // norm function
    static auto     norm = [](const double vec[3]) -> double {
        return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
    };

    for (i = 0; i < partNum; ++i)
    {
        // get the error
        for (j = 0; j < 3; ++j)
        {
            error[j] = previousPos[j] - coordinates[i * 3 + j];
        }

        // accumulate if the particle locates around the previousPos within some
        // enclosed radius
        if (norm(error) < radius)
        {
            massSum += mass[i];
            for (j = 0; j < 3; ++j)
            {
                coordMassSum[j] += mass[i] * coordinates[i * 3 + j];
            }
        }
    }

    MPI_Allreduce(MPI_IN_PLACE, &massSum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, coordMassSum, 3, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    if (massSum != 0)
    {
        for (j = 0; j < 3; ++j)
        {
            com[j] = coordMassSum[j] / massSum;
        }
    }
    else
    {
        int rank = 0;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_WARN(rank, "Get an Mtot=0 in calculation of CoM, return 0 only.");
    }
    return com;
}

/**
 * @brief Calculate the position of the most bound particle.
 *
 * @param potential potential of particles
 * @param coordinates coordinates of particles
 * @param partNum particle number
 * @return the coordinates of the most bound particle
 */
auto recenter::most_bound_particle(const double*   potential,
                                   const double*   coordinates,
                                   const unsigned& partNum)
    -> std::unique_ptr<double[]>
{
    auto      minPotPosition(make_unique<double[]>(3));
    const int minLocateId =
        min_element(potential, potential + partNum) - potential;
    double min = potential[minLocateId];  // local min
    // global min after reduce
    MPI_Allreduce(MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    // Get the rank number of the minimal potential
    int minLocateRank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &minLocateRank);
    if (min != potential[minLocateId])  // local min!=global min
    {
        minLocateRank = 0;
    }
    else
    {
        for (int i = 0; i < 3; ++i)
        {
            minPotPosition[i] = coordinates[3 * minLocateId + i];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &minLocateRank, 1, MPI_INT, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Bcast(minPotPosition.get(), 3, MPI_DOUBLE, minLocateRank,
              MPI_COMM_WORLD);

    return minPotPosition;
}

}  // namespace otf
