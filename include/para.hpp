/**
 * @file para.hpp
 * @brief The parameter parser and its utils.
 */

#pragma once

#include "recenter.hpp"
#include "toml.hpp"
#include <cstdint>
#include <memory>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace otf {

constexpr auto vecDim = 3;

/**
 * @class recenter_para
 * @brief The parameters used for recenter.
 *
 */
struct RecenterPara
{
    bool   enable;
    double radius;                // enclose radius used for coordinate recenter
    double initialGuess[vecDim];  // initial guess of the coordinate center
    otf::recenter_method method;  // recenter method
};

/**
 * @class align_para
 * @brief The parameters used for alignment.
 *
 */
struct AlignPara
{
    bool   enable;
    double radius;  // enclosing radius of the inertia tensor calculation
};

/**
 * @class image_para
 * @brief The parameters used for image calculation.
 *
 */
struct ImagePara
{
    bool     enable;
    double   halfLength;  // half length of the box size to be plotted
    unsigned binNum;      // binnum of the image
};

/**
 * @class BarPara
 * @brief The parameters used for bar info calculation, Sbar et al.
 *
 */
struct BarPara
{
    bool   enable;
    double rmin;
    double rmax;
};

/**
 * @class a2_profile_para
 * @brief The parameters used for A2 radial profile calculation.
 *
 */
struct a2_profile_para
{
    bool     enable;
    double   rmin;
    double   rmax;
    unsigned binNum;
};

/**
 * @class orbit_recenter_para
 * @brief The parameters used for recenter in orbital log.
 *
 */
struct OrbitRecenterPara : RecenterPara
{
    // the anchor type of particles used for recenter
    std::vector<unsigned> anchorIds;
};

enum class coordinate_frame : std::uint8_t
{
    CYLINDRICAL = 0,
    SPHERICAL,
    CARTESIAN
};

/**
 * @class component
 * @brief The wrapper of parameters used for each component.
 *
 */
struct Component
{
    Component(std::string_view& compName, toml::table& compNodeTable);
    std::string           compName;  // name of the component
    std::vector<unsigned> types;     // particle types in this component
    int                   period;    // analysis period
    RecenterPara          recenter;  // parameter of coordinate recenter
    coordinate_frame      frame;     // coordinate frame type
    AlignPara       align;  // whether align coordinates with the inertia tensor
    ImagePara       image;  // parameter of the spatial image part
    BarPara         sBar;   // bar strength parameter
    BarPara         barAngle;   // bar angle parameter
    BarPara         sBuckle;    // buckling strength parameter
    a2_profile_para A2profile;  // A2(R) profile parameter
};

/**
 * @class orbit
 * @brief The wrapper of parameter blocks used for orbital log.
 *
 */
class orbit
{
public:
    explicit orbit(toml::table& orbitNodeTable);
    // method for id log: TXTFILE to use a text file of id list, and RANDOM for
    // random selection according to specified parameters.
    enum class id_selection_method : std::uint8_t
    {
        TXTFILE = 0,
        RANDOM
    };
    [[nodiscard]] auto enable() const { return m_enable; }
    [[nodiscard]] auto period() const { return m_period; }
    [[nodiscard]] auto method() const { return m_method; }
    [[nodiscard]] auto idfile() const { return m_idfile; }
    [[nodiscard]] auto fraction() const { return m_fraction; }
    [[nodiscard]] auto sampleTypes() const -> const auto&
    {
        return m_sampleTypes;
    }
    // TODO: move this to private
    OrbitRecenterPara recenter;  // whether recenter the coordinate of orbits


private:
    bool                m_enable;  // enable orbital log
    int                 m_period;  // log period
    id_selection_method m_method;  // id determination method
    std::string         m_idfile;  // if method is txt file, give the file name
    double m_fraction = -1;  // if method is random sample, give the fraction
    std::vector<int> m_sampleTypes;  // particle types to be sampled
};

/**
 * @class runtime_para
 * @brief Container of the runtime parameter, designed to be work in each mpi
 * rank.
 *
 */
class RuntimePara
{
public:
    explicit RuntimePara(const std::string_view& tomlParaFile);
    [[nodiscard]] auto enableOtf() const { return m_enableOtf; }
    [[nodiscard]] auto outputDir() const { return m_outputDir; }
    [[nodiscard]] auto filename() const { return m_filename; }
    [[nodiscard]] auto maxIter() const { return m_maxIter; }
    [[nodiscard]] auto epsilon() const { return m_epsilon; }

    // TODO: make the following two members private
    // hash map of parameter for each component
    std::unordered_map<std::string, std::unique_ptr<otf::Component>> comps;
    // parameter pointer of orbital logs
    std::unique_ptr<otf::orbit>                                      orbit;

private:
    bool        m_enableOtf;   // whether enable on-the-fly analysis
    std::string m_outputDir;   // output directory of the logs
    std::string m_filename;    // prefix of the log file
    int         m_maxIter{0};  // specify the maximal iteration times
    // the equal threshold of floating-point numbers
    double      m_epsilon{0};
};

}  // namespace otf
