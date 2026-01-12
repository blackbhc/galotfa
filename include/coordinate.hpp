/**
 * @file coordinate.hpp
 * @brief Some functions of coordinate transformation.
 */

#pragma once

#include <cstdint>
enum class coordate_type : std::uint8_t
{
    CARTESIAN = 0,
    SPHERICAL,
    CYLINDRICAL
};

/**
 * @class coordinate_transformer
 * @brief transform the coordinates from a given type to another.
 *
 */
class coordinate_transformer
{
public:
    static void transform(const unsigned&      num,
                          double*              data,
                          const coordate_type& from,
                          const coordate_type& to);

#ifdef DEBUG

#else
private:
#endif
    static void car2sph(double data[3]);
    static void sph2car(double data[3]);
    static void car2cyl(double data[3]);
    static void cyl2car(double data[3]);
    static void sph2cyl(double data[3]);
    static void cyl2sph(double data[3]);
};
