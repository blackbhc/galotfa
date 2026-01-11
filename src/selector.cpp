#include "selector.hpp"
#include "myprompt.hpp"
#include "para.hpp"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iterator>
#include <memory>
#include <mpi.h>
#include <random>
#include <set>
#include <string>
#include <sys/unistd.h>
#include <unistd.h>
#include <utility>
#include <vector>
using namespace std;

namespace otf {

/**
 * @brief Select particle ids with a specified fraction, can be used in any mpi
 * rank.
 *
 * @param raw raw id list
 * @param types pointer to the particle types of simulation data
 * @param sampleTypes vector of the targeted particle types to be sampled
 * @param fraction sampling fraction
 * @return a vector<unsigned> of the selected id list
 */
auto orbit_selector::id_sample(const vector<int>& rawIds,
                               const int*         types,
                               const vector<int>& sampleTypes,
                               double             fraction) -> vector<int>
{
    // check the fraction argument is reasonable
    if (fraction <= 0 or fraction > 1)
    {
        ERROR("Try to select a illegal fraction: [%lf]", fraction);
    }

    auto        partNum = rawIds.size();
    vector<int> filtered(
        partNum);        // exclude the particles without target types
    auto counter = 0UL;  // counter of effective particles in this mpi rank
    for (auto i = 0UL; i < partNum; ++i)
    {
        if (find(sampleTypes.begin(), sampleTypes.end(), types[i])
            != sampleTypes.end())
        {
            filtered[counter++] = rawIds[i];
        }
    }
    filtered.resize(counter);  // remove the tailing rubish values

    if (fraction == 1)  // if it's 100% fraction, directly return the ids
    {
        return filtered;
    }

    // if not 100%, use std::sample to randomly sample the ids
    auto const selectNum =
        ( size_t )(counter * fraction);  // number of data points
    vector<int> res;                     // results
    sample(filtered.begin(), filtered.end(), back_inserter(res), selectNum,
           mt19937{random_device{}()});
    return res;
}

/**
 * @brief Read the ids in the given id list file, can be used in any mpi rank to
 * read the ids into it.
 *
 * @param idFilename filename of the id list file (in ASCII txt format)
 * @return std::vector<int> of the ids.
 */
auto orbit_selector::id_read(const string& idFilename) -> vector<int>
{
    // check the availability of the file
    if (access(idFilename.c_str(), F_OK) != 0)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_ERROR(rank, "Particle ID file not found: [%s]", idFilename.c_str());
        throw "Particle ID File not found!";
    }

    // read the txt file line-by-line
    ifstream    fp;
    vector<int> ids;
    fp.open(idFilename.c_str(), ios::in);
    while (!fp.eof())
    {
        string lineStr;
        getline(fp, lineStr);
        if (!lineStr.empty())
        {
            try
            {
                // reasonable value
                ids.push_back(std::stoi(lineStr));
            }
            catch (...)
            {
                // unreasonable value
                ERROR("Get an unexpected value: %s", lineStr.c_str());
                throw "Get an unexpected value";
            }
        }
    }

    // Remove the possible repeated values
    set<int> tmpSet(ids.begin(), ids.end());
    ids.assign(tmpSet.begin(), tmpSet.end());

    return ids;
}

orbit_selector::orbit_selector(const runtime_para& para) : para(para) { ; }

/**
 * @brief extract the target ids of orbital log based on the specified
 * parameters.
 *
 * @param particleNumber the number of particles, which will be used to
 * determine the array length
 * @param particleID the id of particles
 * @param partType the PartType of particles
 * @return a vector of the extracted ids
 */
auto orbit_selector::extract_target_ids(const unsigned particleNumber,
                                        const int*     particleID,
                                        const int*     partType) const
    -> vector<int>
{
    // get the target id list based on specified parameters
    vector<int> targetIDs;
    if (para.orbit->method
        == otf::orbit::id_selection_method::RANDOM)  // by random sampling
    {
        // restore the raw ids into an vector
        vector<int> rawIds(particleNumber);
        for (auto i = 0U; i < particleNumber; ++i)
        {
            rawIds[i] = particleID[i];
        }
        // random sampling
        auto localTargetIDs = id_sample(
            rawIds, partType, para.orbit->sampleTypes, para.orbit->fraction);

        // gather the target ids in each rank to one vector
        int rank;
        int size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        int localLength =
            localTargetIDs.size();  // the number of ids in local mpi rank
        const unique_ptr<int[]> numInEachRank(
            new int[size]());  // number in each rank
        // collective communication: gather the number of particles in each rank
        MPI_Allgather(&localLength, 1, MPI_INT, numInEachRank.get(), 1, MPI_INT,
                      MPI_COMM_WORLD);

        // get the global total number
        int totalLength = 0;
        for (int i = 0; i < size; ++i)
        {
            totalLength += numInEachRank[i];
        }

        // the offset of the local mpi rank, used for variable length mpi
        // gathering
        int localOffset = 0;
        for (int i = 0; i < rank; ++i)
        {
            localOffset += numInEachRank[i];
        }

        // the array of offset values
        const unique_ptr<int[]> offsets(new int[size]());
        MPI_Allgather(&localOffset, 1, MPI_INT, offsets.get(), 1, MPI_INT,
                      MPI_COMM_WORLD);

        // the global target ids
        vector<int> globalTargetIds(totalLength);
        MPI_Allgatherv(localTargetIDs.data(), localLength, MPI_INT,
                       globalTargetIds.data(), numInEachRank.get(),
                       offsets.get(), MPI_INT, MPI_COMM_WORLD);
        std::swap(globalTargetIds, targetIDs);
    }
    else  // txt file
    {
        // read from a txt file
        targetIDs = id_read(para.orbit->idfile);
    }
    return targetIDs;
}

/**
 * @brief Extracted the data of orbital logs.
 *
 * @param particleNumber number of particles in the local mpi rank
 * @param particleID particle ids
 * @param partType PartTypes of particles
 * @param mass masses of particles
 * @param coordinate coordinates of particles
 * @param velocity velocity of particles
 * @return the extracted data, restore in a dataContainer object
 */
auto orbit_selector::select(const unsigned particleNumber,
                            const int*     particleID,
                            const int*     partType,
                            const double*  mass,
                            const double*  coordinate,
                            const double*  velocity) const
    -> unique_ptr<dataContainer>
{
    if (not para.orbit->m_enable)
    {
        // if the orbital log is not enabled, just ignore the function
        return nullptr;
    }

    static vector<int> targetIDs =
        extract_target_ids(particleNumber, particleID, partType);

    // count of found particles
    unsigned counter = 0;

    // temporary variables restoring extracted data
    vector<double> tmpMass(particleNumber);
    vector<int>    tmpId(particleNumber);
    vector<double> tmpPos(particleNumber * 3);
    vector<double> tmpVel(particleNumber * 3);

    for (auto i = 0U; i < particleNumber; ++i)
    {
        // check whether the id is in the target id list
        if (find(targetIDs.begin(), targetIDs.end(), particleID[i])
            != targetIDs.end())
        {
            tmpMass[counter]        = mass[i];
            tmpId[counter]          = particleID[i];
            tmpPos[counter * 3 + 0] = coordinate[i * 3 + 0];
            tmpPos[counter * 3 + 1] = coordinate[i * 3 + 1];
            tmpPos[counter * 3 + 2] = coordinate[i * 3 + 2];
            tmpVel[counter * 3 + 0] = velocity[i * 3 + 0];
            tmpVel[counter * 3 + 1] = velocity[i * 3 + 1];
            tmpVel[counter * 3 + 2] = velocity[i * 3 + 2];
            // increase the count of particles found
            counter++;
        }
    }

    // remove the rubish value at the end
    tmpMass.resize(counter);
    tmpId.resize(counter);
    tmpPos.resize(counter * 3);
    tmpVel.resize(counter * 3);

    // get the data container
    unique_ptr<dataContainer> container = make_unique<dataContainer>();
    container->count                    = counter;
    container->mass                     = std::move(tmpMass);
    container->id                       = std::move(tmpId);
    container->coordinate               = std::move(tmpPos);
    container->velocity                 = std::move(tmpVel);

    return container;
}

}  // namespace otf
