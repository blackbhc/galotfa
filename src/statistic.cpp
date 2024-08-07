/**
 * @file
 * @brief This file includes a class as a wrapper for statistic functions. At now, mainly the 2D
 * evenly binning statistics for some methods.
 */

#include "../include/statistic.hpp"
#include "../include/myprompt.hpp"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <mpi.h>
using namespace std;

/**
 * @brief 2D binning statistics with chosen method, support count, sum, mean and standard deviation.
 *
 * @param xData: pointing to the first coordinates
 * @param xLowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param xUpperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param xBinNum: binnum of the first coordinate
 * @param yData: pointing to the second coordinates
 * @param yLowerBound: the lower inclusive limit of the second coordinates to be analyzed
 * @param yUpperBound: the upper exclusive limit of the second coordinates to be analyzed
 * @param yBinNum: binnum of the second coordinate
 * @param method: statistic method
 * @param dataNum: number of data points to be analyzed
 * @param data: pointing to target data points
 * @return a unique_ptr pointing to the 1D array of the 2D resutls, in row-major order.
 */
auto statistic::bin2d( int mpiRank, const double* xData, const double xLowerBound,
                       const double xUpperBound, const unsigned long& xBinNum, const double* yData,
                       const double yLowerBound, const double yUpperBound,
                       const unsigned long& yBinNum, statistic_method method,
                       const unsigned long& dataNum, const double* data ) -> unique_ptr< double[] >
{
    switch ( method )
    {
    case statistic_method::COUNT: {
        return bin2dcount( mpiRank, xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound,
                           yUpperBound, yBinNum, dataNum );
    }
    case statistic_method::SUM: {
        return bin2dsum( mpiRank, xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound,
                         yUpperBound, yBinNum, dataNum, data );
    }
    case statistic_method::MEAN: {
        return bin2dmean( mpiRank, xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound,
                          yUpperBound, yBinNum, dataNum, data );
    }
    case statistic_method::STD: {
        return bin2dstd( mpiRank, xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound,
                         yUpperBound, yBinNum, dataNum, data );
    }
    default: {
        ERROR( "Get an unsupported statistic method!" );
        return nullptr;
    }
    }

    return nullptr;
}

/**
 * @brief 2D binning statistics for count
 *
 * @param xData: pointing to the first coordinates
 * @param xLowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param xUpperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param xBinNum: binnum of the first coordinate
 * @param yData: pointing to the second coordinates
 * @param yLowerBound: the lower inclusive limit of the second coordinates to be analyzed
 * @param yUpperBound: the upper exclusive limit of the second coordinates to be analyzed
 * @param yBinNum: binnum of the second coordinate
 * @param dataNum: number of data points to be analyzed
 * @return a unique_ptr pointing to the 1D array of the 2D resutls, in row-major order.
 */
auto statistic::bin2dcount( int mpiRank, const double* xData, const double xLowerBound,
                            const double xUpperBound, const unsigned long& xBinNum,
                            const double* yData, const double yLowerBound, const double yUpperBound,
                            const unsigned long& yBinNum,
                            const unsigned long& dataNum ) -> unique_ptr< double[] >

{
    unsigned long                      idx = 0;
    unsigned long                      idy = 0;
    unique_ptr< double[] >             statisticResutls( new double[ xBinNum * yBinNum ]() );
    const unique_ptr< unsigned int[] > count( new unsigned int[ xBinNum * yBinNum ]() );
    unique_ptr< unsigned int[] >       countRecv = nullptr;

    if ( mpiRank == 0 )
    {
        countRecv = make_unique< unsigned int[] >( xBinNum * yBinNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] >= xLowerBound and xData[ i ] < xUpperBound and yData[ i ] >= yLowerBound
             and yData[ i ] < yUpperBound )
        {
            idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
            idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
            ++count[ idx * yBinNum + idy ];
        }
    }

    MPI_Reduce( count.get(), countRecv.get(), xBinNum * yBinNum, MPI_UNSIGNED, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( mpiRank == 0 )  // effectively update the results in the root process
    {
        for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
        {
            statisticResutls[ i ] = ( double )countRecv[ i ];
        }
    }
    return statisticResutls;
}

/**
 * @brief 2D binning statistics for summation
 *
 * @param xData: pointing to the first coordinates
 * @param xLowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param xUpperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param xBinNum: binnum of the first coordinate
 * @param yData: pointing to the second coordinates
 * @param yLowerBound: the lower inclusive limit of the second coordinates to be analyzed
 * @param yUpperBound: the upper exclusive limit of the second coordinates to be analyzed
 * @param yBinNum: binnum of the second coordinate
 * @param dataNum: number of data points to be analyzed
 * @param data: pointing to target data points
 * @return a unique_ptr pointing to the 1D array of the 2D resutls, in row-major order.
 */
auto statistic::bin2dsum( int mpiRank, const double* xData, const double xLowerBound,
                          const double xUpperBound, const unsigned long& xBinNum,
                          const double* yData, const double yLowerBound, const double yUpperBound,
                          const unsigned long& yBinNum, const unsigned long& dataNum,
                          const double* data ) -> unique_ptr< double[] >

{
    unsigned long                idx = 0;
    unsigned long                idy = 0;
    unique_ptr< double[] >       statisticResutls( new double[ xBinNum * yBinNum ]() );
    const unique_ptr< double[] > sum( new double[ xBinNum * yBinNum ]() );
    unique_ptr< double[] >       sumRecv = nullptr;

    if ( mpiRank == 0 )
    {
        sumRecv = make_unique< double[] >( xBinNum * yBinNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] >= xLowerBound and xData[ i ] < xUpperBound and yData[ i ] >= yLowerBound
             and yData[ i ] < yUpperBound )
        {
            idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
            idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
            sum[ idx * yBinNum + idy ] += data[ i ];
        }
    }

    MPI_Reduce( sum.get(), sumRecv.get(), xBinNum * yBinNum, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( mpiRank == 0 )  // effectively update the results in the root process
    {
        for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
        {
            statisticResutls[ i ] = sumRecv[ i ];
        }
    }
    return statisticResutls;
}

/**
 * @brief 2D binning statistics for mean values
 *
 * @param xData: pointing to the first coordinates
 * @param xLowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param xUpperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param xBinNum: binnum of the first coordinate
 * @param yData: pointing to the second coordinates
 * @param yLowerBound: the lower inclusive limit of the second coordinates to be analyzed
 * @param yUpperBound: the upper exclusive limit of the second coordinates to be analyzed
 * @param yBinNum: binnum of the second coordinate
 * @param dataNum: number of data points to be analyzed
 * @param data: pointing to target data points
 * @return a unique_ptr pointing to the 1D array of the 2D resutls, in row-major order.
 */
auto statistic::bin2dmean( int mpiRank, const double* xData, const double xLowerBound,
                           const double xUpperBound, const unsigned long& xBinNum,
                           const double* yData, const double yLowerBound, const double yUpperBound,
                           const unsigned long& yBinNum, const unsigned long& dataNum,
                           const double* data ) -> unique_ptr< double[] >
{
    unsigned long                      idx = 0;
    unsigned long                      idy = 0;
    unique_ptr< double[] >             statisticResutls( new double[ xBinNum * yBinNum ]() );
    const unique_ptr< double[] >       sum( new double[ xBinNum * yBinNum ]() );
    const unique_ptr< unsigned int[] > count( new unsigned int[ xBinNum * yBinNum ]() );
    unique_ptr< double[] >             sumRecv( new double[ xBinNum * yBinNum ]() );
    unique_ptr< unsigned int[] >       countRecv( new unsigned int[ xBinNum * yBinNum ]() );

    if ( mpiRank == 0 )
    {
        sumRecv   = make_unique< double[] >( xBinNum * yBinNum );
        countRecv = make_unique< unsigned int[] >( xBinNum * yBinNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] >= xLowerBound and xData[ i ] < xUpperBound and yData[ i ] >= yLowerBound
             and yData[ i ] < yUpperBound )
        {
            idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
            idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
            ++count[ idx * yBinNum + idy ];
            sum[ idx * yBinNum + idy ] += data[ i ];
        }
    }

    MPI_Reduce( count.get(), countRecv.get(), xBinNum * yBinNum, MPI_UNSIGNED, MPI_SUM, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( sum.get(), sumRecv.get(), xBinNum * yBinNum, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( mpiRank == 0 )
    {
        for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
        {
            if ( count[ i ] != 0 )
            {
                statisticResutls[ i ] = sumRecv[ i ] / countRecv[ i ];
            }
            else
            {
                statisticResutls[ i ] = nan( "" );
            }
        }
    }

    return statisticResutls;
}

/**
 * @brief 2D binning statistics for standard deviation, without Bessel correction.
 *
 * @param xData: pointing to the first coordinates
 * @param xLowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param xUpperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param xBinNum: binnum of the first coordinate
 * @param yData: pointing to the second coordinates
 * @param yLowerBound: the lower inclusive limit of the second coordinates to be analyzed
 * @param yUpperBound: the upper exclusive limit of the second coordinates to be analyzed
 * @param yBinNum: binnum of the second coordinate
 * @param dataNum: number of data points to be analyzed
 * @param data: pointing to target data points
 * @return a unique_ptr pointing to the 1D array of the 2D resutls, in row-major order.
 */
auto statistic::bin2dstd( int mpiRank, const double* xData, const double xLowerBound,
                          const double xUpperBound, const unsigned long& xBinNum,
                          const double* yData, const double yLowerBound, const double yUpperBound,
                          const unsigned long& yBinNum, const unsigned long& dataNum,
                          const double* data ) -> unique_ptr< double[] >

{
    unsigned long          idx = 0;
    unsigned long          idy = 0;
    unique_ptr< double[] > statisticResutls( new double[ xBinNum * yBinNum ]() );

    const unique_ptr< double[] >       sum( new double[ xBinNum * yBinNum ]() );
    const unique_ptr< unsigned int[] > count( new unsigned int[ xBinNum * yBinNum ]() );
    unique_ptr< double[] >             sumRecv( new double[ xBinNum * yBinNum ]() );
    unique_ptr< unsigned int[] >       countRecv( new unsigned int[ xBinNum * yBinNum ]() );

    if ( mpiRank == 0 )
    {
        sumRecv   = make_unique< double[] >( xBinNum * yBinNum );
        countRecv = make_unique< unsigned int[] >( xBinNum * yBinNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] >= xLowerBound and xData[ i ] < xUpperBound and yData[ i ] >= yLowerBound
             and yData[ i ] < yUpperBound )
        {
            idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
            idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
            ++count[ idx * yBinNum + idy ];
            sum[ idx * yBinNum + idy ] += data[ i ];
        }
    }

    MPI_Reduce( count.get(), countRecv.get(), xBinNum * yBinNum, MPI_UNSIGNED, MPI_SUM, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( sum.get(), sumRecv.get(), xBinNum * yBinNum, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( mpiRank == 0 )
    {
        for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
        {
            if ( count[ i ] != 0 )
            {
                statisticResutls[ i ] = sumRecv[ i ] / countRecv[ i ];
            }
            else
            {
                statisticResutls[ i ] = nan( "" );
            }
        }
    }

    MPI_Bcast( statisticResutls.get(), xBinNum * yBinNum, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    std::memset( sum.get(), 0, sizeof( double ) * xBinNum * yBinNum );  // reset the sum to 0
    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] >= xLowerBound and xData[ i ] < xUpperBound and yData[ i ] >= yLowerBound
             and yData[ i ] < yUpperBound )
        {
            idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
            idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
            sum[ idx * yBinNum + idy ] += ( data[ i ] - statisticResutls[ idx * yBinNum + idy ] )
                                          * ( data[ i ] - statisticResutls[ idx * yBinNum + idy ] );
        }
    }

    MPI_Reduce( sum.get(), sumRecv.get(), xBinNum * yBinNum, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( mpiRank == 0 )
    {
        for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
        {
            if ( count[ i ] != 0 )
            {
                statisticResutls[ i ] = sqrt( sumRecv[ i ] / countRecv[ i ] );
            }
        }
    }

    return statisticResutls;
}

/**
 * @brief function that determines the bin a data point should be located, for evenly distributed
 * bins only and there is no boundary check!
 *
 * @param lowerBound: the lower boundary of the data range.
 * @param upperBound: the upper boundary of the data range.
 * @param binNum: the number of bins.
 * @param value: the value of the data point.
 * @return
 */
auto statistic::find_index( double lowerBound, double upperBound, const unsigned long& binNum,
                            double value ) -> unsigned long
{
    return ( value - lowerBound ) / ( upperBound - lowerBound ) * binNum;
}
