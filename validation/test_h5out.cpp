/**
 * @file test_h5out.cpp
 * @brief Test the hdf5 output class.
 */

#define DEBUG 1
#include "../include/h5out.hpp"
#include "H5Tpublic.h"
#include <hdf5.h>

int main()
{
    int testData[30 * 3];
    int count = 0;
    for (int i = 0; i < 30; ++i)
        for (int j = 0; j < 3; ++j)
            testData[i * 3 + j] = count++;


    h5_out organizer("./test_log/", "galotfa.hdf5");
    organizer.create_dataset_in_group("A2", "component1", {3}, H5T_NATIVE_INT);
    for (int i = 0; i < 30; ++i)
    {
        auto returnCode =
            organizer.flush_single_block("component1", "A2", testData + i * 3);
        if (returnCode != 0)
            return returnCode;
    }

    double mockImage[4 * 4];
    for (int i = 0; i < 16; ++i)
    {
        mockImage[i] = i + 0.314;
    }
    organizer.create_dataset_in_group("Image", "component1", {4, 4},
                                      H5T_NATIVE_DOUBLE_g);
    for (int i = 0; i < 4; ++i)
    {
        auto returnCode =
            organizer.flush_single_block("component1", "Image", mockImage);
        if (returnCode != 0)
            return returnCode;
    }

    double mockOrbit[3 * 7];
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 7; ++j)
            mockOrbit[i * 7 + j] = i * 7 + j + 0.314;
    }

    organizer.create_dataset_in_group("particle-x", "Orbit", {7},
                                      H5T_NATIVE_DOUBLE_g);
    for (int i = 0; i < 3; ++i)
    {
        auto returnCode = organizer.flush_single_block("Orbit", "particle-x",
                                                       mockOrbit + i * 7);
        if (returnCode != 0)
            return returnCode;
    }

    return 0;
}
