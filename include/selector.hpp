/**
 * @file selector.hpp
 * @brief ID selector and reader.
 */

#ifndef SELECTOR_HEADER
#define SELECTOR_HEADER
#include "para.hpp"
#include <memory>
#include <string>
#include <vector>
namespace otf {

/**
 * @class dataContainer
 * @brief Collection of points to the data for component-analysis and orbital
 * log, aiming for convenience of data selection.
 *
 */
struct dataContainer
{
public:
    unsigned            count = 0;
    std::vector<int>    id;
    std::vector<double> mass;
    std::vector<double> coordinate;
    std::vector<double> velocity;
};

/**
 * @class orbit_selector
 * @brief Data selection for orbital log.
 *
 */
class orbit_selector
{
public:
    orbit_selector(const runtime_para& para);
    auto select(unsigned      particleNumber,
                const int*    particleID,
                const int*    partType,
                const double* mass,
                const double* coordinate,
                const double* velocity) const -> std::unique_ptr<dataContainer>;

#ifdef DEBUG

#else
private:
#endif
    const runtime_para& para;
    static auto         id_sample(const std::vector<int>& rawIds,
                                  const int*              types,
                                  const std::vector<int>& sampleTypes,
                                  double                  fraction) -> std::vector<int>;
    static auto id_read(const std::string& idFilename) -> std::vector<int>;
    auto        extract_target_ids(unsigned   particleNumber,
                                   const int* particleID,
                                   const int* partType) const -> std::vector<int>;
};

}  // namespace otf
#endif
