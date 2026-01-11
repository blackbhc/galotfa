/**
 * @file test_coordTrans.cpp
 * @brief Test the coordinate transformation class.
 */

#define DEBUG 1
#include "../include/coordinate.hpp"
#include <iostream>
#include <vector>
using namespace std;

int main()
{
    coordinate_transformer transfer;
    vector<double*>        cartesians;
    double                 car1[3] = {1, 0, 0};
    double                 car2[3] = {0, 1, 0};
    double                 car3[3] = {0, 0, 1};
    double                 car4[3] = {1, 1, 0};
    double                 car5[3] = {0, 1, 1};
    double                 car6[3] = {1, 0, 1};
    double                 car7[3] = {1, 2, 3};
    cartesians.push_back(car1);
    cartesians.push_back(car2);
    cartesians.push_back(car3);
    cartesians.push_back(car4);
    cartesians.push_back(car5);
    cartesians.push_back(car6);
    cartesians.push_back(car7);

    cout << "--------------" << endl;
    cout << "Original cartesian coordinates:" << endl;
    for (auto& array : cartesians)
    {
        cout << array[0] << " " << array[1] << " " << array[2] << " " << endl;
    }

    cout << "--------------" << endl;
    cout << "Spherical coordinates:" << endl;
    for (auto& array : cartesians)
    {
        transfer.car2sph(array);
    }
    for (auto& array : cartesians)
    {
        cout << array[0] << " " << array[1] << " " << array[2] << " " << endl;
    }

    cout << "--------------" << endl;
    cout << "Spherical to cylindrical coordinates:" << endl;
    for (auto& array : cartesians)
    {
        transfer.sph2cyl(array);
    }
    for (auto& array : cartesians)
    {
        cout << array[0] << " " << array[1] << " " << array[2] << " " << endl;
    }

    cout << "--------------" << endl;
    cout << "Restored spherical coordinates:" << endl;
    for (auto& array : cartesians)
    {
        transfer.cyl2sph(array);
    }
    for (auto& array : cartesians)
    {
        cout << array[0] << " " << array[1] << " " << array[2] << " " << endl;
    }

    return 0;
}
