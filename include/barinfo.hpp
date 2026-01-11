/**
 * @file barinfo.hpp
 * @brief Wrapper class of bar quantification functions.
 */

#ifndef BARINFO_HEADER
#define BARINFO_HEADER

namespace otf {

/**
 * @class bar_info
 * @brief A0, A2, Sbar, Sbuckle, bar ellipticity (to be implemented)
 *
 */
class bar_info
{
public:
    struct A2info
    {
        double amplitude;  // amplitude of the m=2 Fourier mode
        double phase;      // phase angle of the m=2 Fourier mode
    };
    static auto A0(unsigned partNum, const double* masses) -> double;
    static auto A2(unsigned partNum, const double* masses, const double* phis)
        -> double;
    static auto bar_angle(unsigned      partNum,
                          const double* masses,
                          const double* phis) -> double;
    static auto Sbuckle(unsigned      partNum,
                        const double* masses,
                        const double* phis,
                        const double* zeds) -> double;
};

}  // namespace otf
#endif
