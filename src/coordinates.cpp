#include "coordinate.hpp"
#include "myprompt.hpp"
#include <cmath>
#include <numbers>
using namespace std;

/**
 * @brief Transform the coordinate from cartesian to spherical.
 *
 * @param data double[3] for x, y, z to r, theta, phi
 */
void coordinate_transformer::car2sph( double data[ 3 ] )
{
    static double r = 0;
    // theta for latitude, phi for longitude
    static double theta = 0;
    static double phi   = 0;
    r = sqrt( data[ 0 ] * data[ 0 ] + data[ 1 ] * data[ 1 ] + data[ 2 ] * data[ 2 ] );

    if ( r != 0 )
    {
        theta = acos( data[ 2 ] / r );
    }
    else if ( data[ 2 ] >= 0 )
    {
        theta = 0;
    }
    else
    {
        theta = -numbers::pi;
    }

    phi = atan2( data[ 1 ], data[ 0 ] );

    data[ 0 ] = r;
    data[ 1 ] = theta;
    data[ 2 ] = phi;
}

/**
 * @brief Transform the coordinate from spherical to cartesian.
 *
 * @param data double[3] for r, theta (latitude), phi (longitude) to x, y, z
 */
void coordinate_transformer::sph2car( double data[ 3 ] )
{
    static double x = 0;
    static double y = 0;
    static double z = 0;
    x               = data[ 0 ] * sin( data[ 1 ] ) * cos( data[ 2 ] );
    y               = data[ 0 ] * sin( data[ 1 ] ) * sin( data[ 2 ] );
    z               = data[ 0 ] * cos( data[ 1 ] );
    data[ 0 ]       = x;
    data[ 1 ]       = y;
    data[ 2 ]       = z;
}

/**
 * @brief Transform the coordinate from cartesian to cylindrical.
 *
 * @param data double[3] for x, y, z to R, phi, z
 */
void coordinate_transformer::car2cyl( double data[ 3 ] )
{
    static double R   = 0;
    static double phi = 0;
    R                 = sqrt( data[ 0 ] * data[ 0 ] + data[ 1 ] * data[ 1 ] );
    phi               = atan2( data[ 1 ], data[ 0 ] );

    data[ 0 ] = R;
    data[ 1 ] = phi;
}

/**
 * @brief Transform the coordinate from cylindrical to cartesian.
 *
 * @param data double[3] for R, phi, z to x, y, z
 */
void coordinate_transformer::cyl2car( double data[ 3 ] )
{
    static double x = 0;
    static double y = 0;
    x               = data[ 0 ] * cos( data[ 1 ] );
    y               = data[ 0 ] * sin( data[ 1 ] );

    data[ 0 ] = x;
    data[ 1 ] = y;
}

/**
 * @brief Transform the coordinate from cartesian to cylindrical.
 *
 * @param data double[3] for r, theta, phi to R, phi, z
 */
void coordinate_transformer::sph2cyl( double data[ 3 ] )
{
    static double R   = 0;
    static double phi = 0;
    static double z   = 0;
    R                 = data[ 0 ] * sin( data[ 1 ] );
    z                 = data[ 0 ] * cos( data[ 1 ] );
    phi               = data[ 2 ];
    data[ 0 ]         = R;
    data[ 1 ]         = phi;
    data[ 2 ]         = z;
}

/**
 * @brief Transform the coordinate from cylindrical to cartesian.
 *
 * @param data double[3] for R, phi, z to r, theta, phi
 */
void coordinate_transformer::cyl2sph( double data[ 3 ] )
{
    static double r     = 0;
    static double theta = 0;
    static double phi   = 0;  // theta for latitude, phi for longitude
    r                   = sqrt( data[ 0 ] * data[ 0 ] + data[ 2 ] * data[ 2 ] );
    theta               = atan( data[ 0 ] / data[ 2 ] );
    if ( theta < 0 )
    {
        theta += numbers::pi;
    }
    phi       = data[ 1 ];
    data[ 0 ] = r;
    data[ 1 ] = theta;
    data[ 2 ] = phi;
}

/**
 * @brief Transform the array of coordinates into other frame.
 *
 * @param num array length (array is <num>x3)
 * @param data the pointer to the array
 * @param from original frame type
 * @param to target frame type
 */
void coordinate_transformer::transform( const unsigned& num, double* data,
                                        const coordate_type& from, const coordate_type& to )
{
    if ( from == to )
    {
        WARN( "Try to perform identical transformation!" );
        return;
    }

    if ( from == coordate_type::CARTESIAN and to == coordate_type::SPHERICAL )
    {
        for ( auto i = 0U; i < num; ++i )
        {
            car2sph( data + static_cast< unsigned long >( 3 * i ) );
        }
        return;
    }

    if ( from == coordate_type::CARTESIAN and to == coordate_type::CYLINDRICAL )
    {
        for ( auto i = 0U; i < num; ++i )
        {
            car2cyl( data + static_cast< unsigned long >( 3 * i ) );
        }
        return;
    }

    if ( from == coordate_type::SPHERICAL and to == coordate_type::CARTESIAN )
    {
        for ( auto i = 0U; i < num; ++i )
        {
            sph2car( data + static_cast< unsigned long >( 3 * i ) );
        }
        return;
    }

    if ( from == coordate_type::SPHERICAL and to == coordate_type::CYLINDRICAL )
    {
        for ( auto i = 0U; i < num; ++i )
        {
            sph2cyl( data + static_cast< unsigned long >( 3 * i ) );
        }
        return;
    }

    if ( from == coordate_type::CYLINDRICAL and to == coordate_type::CARTESIAN )
    {
        for ( auto i = 0U; i < num; ++i )
        {
            cyl2car( data + static_cast< unsigned long >( 3 * i ) );
        }
        return;
    }

    if ( from == coordate_type::CYLINDRICAL and to == coordate_type::SPHERICAL )
    {
        for ( auto i = 0U; i < num; ++i )
        {
            cyl2sph( data + static_cast< unsigned long >( 3 * i ) );
        }
        return;
    }
}
