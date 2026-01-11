/**
 * @file recenter.hpp
 * @brief Some utilities for recenter of N-body system.
 */

#ifndef RECENTER_HEADER
#define RECENTER_HEADER
#include <cstdint>
#include <memory>

namespace otf {

enum class recenter_method : std::uint8_t
{
    COM = 0,
    MBP = 1
};

/**
 * @class recenter
 * @brief Wrapper class of the recenter APIs.
 *
 */
class recenter
{
public:
    static auto get_center(recenter_method method,
                           const unsigned& partNum,
                           const double*   masses,
                           const double*   potentials,
                           const double*   coordinates,
                           double          radius,
                           const double*   previousPos = nullptr)
        -> std::unique_ptr<double[]>;

#ifdef DEBUG

#else
private:
#endif
    static auto center_of_mass(const double*   mass,
                               const double*   coordinates,
                               const unsigned& partNum,
                               double          radius,
                               const double*   previousPos = nullptr)
        -> std::unique_ptr<double[]>;
    static auto most_bound_particle(const double*   potentials,
                                    const double*   coordinates,
                                    const unsigned& partNum)
        -> std::unique_ptr<double[]>;
};

}  // namespace otf
#endif
