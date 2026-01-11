/**
 * @file test_orbitalSelection.cpp
 * @brief Test the orbital particle selection.
 */

#define DEBUG 1
#include "../include/myprompt.hpp"
#include "../include/selector.hpp"
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

    int mockTypes[40];
    for (auto& mockType : mockTypes)
    {
        mockType = 2;
    }

    int mockIDs[40];
    for (auto i = 0; i < 40; ++i)
    {
        mockIDs[i] = i + 1;
    }

    double mockMass[40];
    for (auto i = 0U; i < 40; ++i)
    {
        mockMass[i] = ( double )i + 0.1;
    }

    double mockPos[120];
    for (auto i = 0U; i < 120; ++i)
    {
        mockPos[i] = ( double )i + 0.2;
    }

    double mockVel[120];
    for (auto i = 0U; i < 120; ++i)
    {
        mockVel[i] = pow(-1, i) * (( double )i + 0.3);
    }

    runtime_para para("../validation/orbit_select_test.toml");
    if (not para.orbit->m_enable)
    {
        MPI_Finalize();
        return 0;
    }

    auto hasMass = [&mockMass](double value) {
        for (auto& mass : mockMass)
        {
            if (mass == value)
            {
                return true;
            }
        }
        return false;
    };

    auto hasPos = [&mockPos](double* coordPtr) {
        for (int i = 0; i < 40; ++i)
        {
            if (mockPos[3 * i + 0] == coordPtr[0]
                && mockPos[3 * i + 1] == coordPtr[1]
                && mockPos[3 * i + 2] == coordPtr[2])
            {
                return true;
            }
        }
        return false;
    };

    auto hasVel = [&mockVel](double* velPtr) {
        for (int i = 0; i < 40; ++i)
        {
            if (mockVel[3 * i + 0] == velPtr[0]
                && mockVel[3 * i + 1] == velPtr[1]
                && mockVel[3 * i + 2] == velPtr[2])
            {
                return true;
            }
        }
        return false;
    };

    orbit_selector orbitSelector(para);
    auto           getData = orbitSelector.select(
        10, mockIDs + 10 * rank, mockTypes + 10 * rank, mockMass + 10 * rank,
        mockPos + 3 * 10 * rank, mockVel + 3 * 10 * rank);
    if (para.orbit->method == otf::orbit::id_selection_method::RANDOM)
    {
        assert(getData->count == 9);
        for (auto i = 0U; i < getData->count; ++i)
        {
            assert(hasMass(getData->mass[i]));
            assert(hasPos(getData->coordinate.data() + 3 * i));
            assert(hasVel(getData->velocity.data() + 3 * i));
        }
    }
    else
    {
        assert(getData->count == ( unsigned )rank);
        switch (rank)
        {
        case 0:
            for (auto i = 0U; i < getData->count; ++i)
            {
                // by design, the programe should never reach here
                assert(false);
            }
            break;
        case 1:
            for (auto i = 0U; i < getData->count; ++i)
            {
                // 10
                assert(getData->mass[i] == mockMass[10 * rank + rank - 1 + i]);
                for (int j = 0; j < 3; ++j)
                {
                    assert(getData->coordinate[3 * i + j]
                           == mockPos[3 * (10 * rank + rank - 1 + i) + j]);
                }
                for (int j = 0; j < 3; ++j)
                {
                    assert(getData->velocity[3 * i + j]
                           == mockVel[3 * (10 * rank + rank - 1 + i) + j]);
                }
            }
            break;
        case 2:
            for (auto i = 0U; i < getData->count; ++i)
            {
                // 21-22
                assert(getData->mass[i] == mockMass[10 * rank + rank - 1 + i]);
                for (int j = 0; j < 3; ++j)
                {
                    assert(getData->coordinate[3 * i + j]
                           == mockPos[3 * (10 * rank + rank - 1 + i) + j]);
                }
                for (int j = 0; j < 3; ++j)
                {
                    assert(getData->velocity[3 * i + j]
                           == mockVel[3 * (10 * rank + rank - 1 + i) + j]);
                }
            }
            break;
        case 3:
            for (auto i = 0U; i < getData->count; ++i)
            {
                // 33-35
                assert(getData->mass[i] == mockMass[10 * rank + rank + i]);
                for (int j = 0; j < 3; ++j)
                {
                    assert(getData->coordinate[3 * i + j]
                           == mockPos[3 * (10 * rank + rank + i) + j]);
                }
                for (int j = 0; j < 3; ++j)
                {
                    assert(getData->velocity[3 * i + j]
                           == mockVel[3 * (10 * rank + rank + i) + j]);
                }
            }
            break;
        }
    }

    MPI_Finalize();
    return 0;
}
