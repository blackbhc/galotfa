/**
 * @file statistic.hpp
 * @brief This file includes a class as a wrapper for statistic functions. At
 * now, mainly the 1D/2D evenly binning statistics for limited methods.
 */

#pragma once

#include <cstdint>
#include <memory>
enum class statistic_method : std::uint8_t
{
    COUNT = 0,
    MEAN,
    STD,
    SUM
};

/**
 * @class statistic
 * @brief Wrapper class of the internal statistic system.
 *
 */
class statistic
{
public:
    static auto bin2d(int              mpiRank,
                      const double*    xData,
                      double           xLowerBound,
                      double           xUpperBound,
                      unsigned long    xBinNum,
                      const double*    yData,
                      double           yLowerBound,
                      double           yUpperBound,
                      unsigned long    yBinNum,
                      statistic_method method,
                      unsigned long    dataNum,
                      const double*    data = nullptr)
        -> std::unique_ptr<double[]>;
    static auto bin1d(int              mpiRank,
                      const double*    coord,
                      double           lowerBound,
                      double           upperBound,
                      unsigned long    binNum,
                      statistic_method method,
                      unsigned long    dataNum,
                      const double*    data = nullptr)
        -> std::unique_ptr<double[]>;

#ifdef DEBUG

#else
private:
#endif
    // function that determines the bin a data point should be located (evenly
    // distributed bins)
    static auto find_index(double        lowerBound,
                           double        upperBound,
                           unsigned long binNum,
                           double        value) -> unsigned long;
    static auto bin2dcount(int           mpiRank,
                           const double* xData,
                           double        xLowerBound,
                           double        xUpperBound,
                           unsigned long xBinNum,
                           const double* yData,
                           double        yLowerBound,
                           double        yUpperBound,
                           unsigned long yBinNum,
                           unsigned long dataNum) -> std::unique_ptr<double[]>;
    static auto bin2dsum(int           mpiRank,
                         const double* xData,
                         double        xLowerBound,
                         double        xUpperBound,
                         unsigned long xBinNum,
                         const double* yData,
                         double        yLowerBound,
                         double        yUpperBound,
                         unsigned long yBinNum,
                         unsigned long dataNum,
                         const double* data) -> std::unique_ptr<double[]>;
    static auto bin2dmean(int           mpiRank,
                          const double* xData,
                          double        xLowerBound,
                          double        xUpperBound,
                          unsigned long xBinNum,
                          const double* yData,
                          double        yLowerBound,
                          double        yUpperBound,
                          unsigned long yBinNum,
                          unsigned long dataNum,
                          const double* data) -> std::unique_ptr<double[]>;
    static auto bin2dstd(const double* xData,
                         double        xLowerBound,
                         double        xUpperBound,
                         unsigned long xBinNum,
                         const double* yData,
                         double        yLowerBound,
                         double        yUpperBound,
                         unsigned long yBinNum,
                         unsigned long dataNum,
                         const double* data) -> std::unique_ptr<double[]>;
    static auto bin1dcount(int           mpiRank,
                           const double* coord,
                           double        lowerBound,
                           double        upperBound,
                           unsigned long binNum,
                           unsigned long dataNum) -> std::unique_ptr<double[]>;
    static auto bin1dsum(int           mpiRank,
                         const double* coord,
                         double        lowerBound,
                         double        upperBound,
                         unsigned long binNum,
                         unsigned long dataNum,
                         const double* data = nullptr)
        -> std::unique_ptr<double[]>;
    static auto bin1dmean(int           mpiRank,
                          const double* coord,
                          double        lowerBound,
                          double        upperBound,
                          unsigned long binNum,
                          unsigned long dataNum,
                          const double* data = nullptr)
        -> std::unique_ptr<double[]>;
    static auto bin1dstd(const double* coord,
                         double        lowerBound,
                         double        upperBound,
                         unsigned long binNum,
                         unsigned long dataNum,
                         const double* data = nullptr)
        -> std::unique_ptr<double[]>;
};
