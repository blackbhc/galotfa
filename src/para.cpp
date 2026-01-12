#include "para.hpp"
#include "myprompt.hpp"
#include "recenter.hpp"
#include "toml.hpp"
#include <cassert>
#include <memory>
#include <mpi.h>
#include <string>
#include <string_view>
#include <unistd.h>
using namespace std;

namespace otf {

runtime_para::runtime_para(const std::string_view& tomlParaFile)
{
    if (access(tomlParaFile.data(), F_OK) != 0)
    {
        ERROR("The toml file [%s] not found!", tomlParaFile.data());
        throw;
    }

    toml::table paraTable = toml::parse_file(tomlParaFile);

    // check whether enable the on-the-fly analysis
    m_enableOtf = *paraTable["global"]["enable"].value<bool>();
    if (not enableOtf())
    {
        INFO("The orbital log is not enabled");
        return;
    }

    // other global parameters
    const string tmpFileName(
        *paraTable["global"]["filename"].value<string_view>());
    const string tmpDir(*paraTable["global"]["outdir"].value<string_view>());
    // fileName                              = std::move( tmpFileName );
    m_outputDir                       = tmpDir;
    m_filename                        = tmpFileName;
    constexpr unsigned defaultMaxIter = 25;
    constexpr double   defaultEpsilon =
        1e-8;  // floating-point number equal threshold
    m_maxIter = static_cast<int>(
        paraTable["global"]["maxiter"].value_or(defaultMaxIter));
    m_epsilon = paraTable["global"]["outdir"].value_or(defaultEpsilon);
    if (not(m_maxIter > 0))
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_ERROR(rank, "maxIter must be positive!");
        throw;
    }

    if (not(m_epsilon > 0))
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_ERROR(rank, "epsilon must be positive!");
        throw;
    };

    // orbital log parameters
    orbit = make_unique<otf::orbit>(*paraTable["orbit"].as_table());

    // parameters of each components
    paraTable.for_each([this](auto& key, auto& value) {
        if constexpr (toml::is_key<decltype(key)>
                      and toml::is_table<decltype(value)>)
        {
            auto key_view = string_view(key);
            if (key_view.substr(0, 9) == "component")
            {
                this->comps[string(key_view)] =
                    make_unique<otf::component>(key_view, *value.as_table());
            }
        }
    });

    // remove useless components: a component without effective analysis
    vector<string> toRemove;  // vector of component names
    for (auto& comp : comps)
    {
        bool effective;
        // if at least one information is enabled, it's an effective component
        effective = comp.second->recenter.enable or comp.second->align.enable
                    or comp.second->sBar.enable or comp.second->barAngle.enable
                    or comp.second->sBuckle.enable or comp.second->image.enable
                    or comp.second->A2profile.enable;

        if (not effective)
        {
            toRemove.push_back(comp.first);
        }
    }
    // remove uneffective components
    for (auto& str : toRemove)
    {
        comps.erase(str);
    }

    // NOTE: if there is no any component and orbital logs are enables, then
    // toggle off the on-the-fly analysis
    if (comps.size() == 0 and (not orbit->enable()))
    {
        m_enableOtf = false;
    }
}

component::component(string_view& compName, toml::table& compNodeTable)
{
    this->compName = compName;
    // particle types in this component
    auto typeIDs = compNodeTable["types"];
    if (toml::array* arr = typeIDs.as_array())
    {
        // visitation with for_each() helps deal with heterogeneous data
        arr->for_each([this](auto&& element) {
            if constexpr (toml::is_number<decltype(element)>)
            {
                types.push_back(*element);
            }
        });
    }

    // period
    period = *compNodeTable["period"].value<int>();
    if (not(period > 0))
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_ERROR(rank, "period must be positive!");
        throw;
    };

    // recenter parameters
    recenter.enable = *compNodeTable["recenter"]["enable"].value<bool>();
    if (recenter.enable)
    {
        recenter.radius = *compNodeTable["recenter"]["radius"].value<double>();
        const string_view str =
            *compNodeTable["recenter"]["method"].value<string_view>();
        if (str == "com")
        {
            recenter.method = recenter_method::COM;
        }
        else if (str == "mbp")
        {
            recenter.method = recenter_method::MBP;
        }
        else
        {
            ERROR("Get an unknown value for [recenter method] of [%s]: [%s]",
                  compName.data(), str.data());
            ERROR("Must be 'com' (for center of mass) or 'mbp' (for most bound "
                  "particle).");
            exit(-1);
        }
        auto iguess = compNodeTable["recenter"]["iguess"];
        for (auto i = 0; i < 3; ++i)
        {
            recenter.initialGuess[i] = *iguess[i].value<double>();
        }
    }

    // NOTE: the frame parameter is unused in the current version
    //
    // // frame
    // auto str = *compNodeTable[ "frame" ].value< string_view >();
    // if ( str == "cyl" )
    // {
    //     frame = coordinate_frame::CYLINDRICAL;
    // }
    // else if ( str == "sph" )
    // {
    //     frame = coordinate_frame::SPHERICAL;
    // }
    // else if ( str == "car" )
    // {
    //     frame = coordinate_frame::CARTESIAN;
    // }
    // else
    // {
    //     ERROR( "Get an unknown value for [coordinate frame] of [%s]: [%s]",
    //            compName.data(), str.data() );
    //     ERROR( "Must be one of 'cyl' (for cylindrical), 'car' (for
    //     Cartesian), or 'sph' (for "
    //            "spherical)." );
    //     exit( -1 );
    // }

    // align
    align.enable = *compNodeTable["align"]["enable"].value<bool>();
    if (align.enable)
    {
        align.radius = *compNodeTable["align"]["radius"].value<double>();
    }

    // image
    image.enable = *compNodeTable["image"]["enable"].value<bool>();
    if (image.enable)
    {
        image.halfLength =
            *compNodeTable["image"]["halflength"].value<double>();
        image.binNum = *compNodeTable["image"]["binnum"].value<unsigned>();
    }

    // bar info parameters
    // sBar
    sBar.enable = *compNodeTable["A2"]["enable"].value<bool>();
    if (sBar.enable)
    {
        sBar.rmin = *compNodeTable["A2"]["rmin"].value<double>();
        sBar.rmax = *compNodeTable["A2"]["rmax"].value<double>();
        if (not(sBar.rmin >= 0 and sBar.rmin < sBar.rmax))
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_ERROR(rank,
                      "The radial range for bar strength calculation of [%s] "
                      "is illegal.",
                      compName.data());
            throw;
        };
    }
    // bar angle
    barAngle.enable = *compNodeTable["barangle"]["enable"].value<bool>();
    if (barAngle.enable)
    {
        barAngle.rmin = *compNodeTable["barangle"]["rmin"].value<double>();
        barAngle.rmax = *compNodeTable["barangle"]["rmax"].value<double>();
        if (not(barAngle.rmin >= 0 and barAngle.rmin < barAngle.rmax))
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_ERROR(rank,
                      "The radial range for bar angle calculation of [%s] is "
                      "illegal.",
                      compName.data());
            throw;
        };
    }
    // buckling strength
    sBuckle.enable = *compNodeTable["buckle"]["enable"].value<bool>();
    if (sBuckle.enable)
    {
        sBuckle.rmin = *compNodeTable["buckle"]["rmin"].value<double>();
        sBuckle.rmax = *compNodeTable["buckle"]["rmax"].value<double>();
        if (not(sBuckle.rmin >= 0 and sBuckle.rmin < sBuckle.rmax))
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_ERROR(rank,
                      "The radial range for buckling strength calculation of "
                      "[%s] is illegal.",
                      compName.data());
            throw;
        };
    }
    // buckling strength
    A2profile.enable = *compNodeTable["A2profile"]["enable"].value<bool>();
    if (A2profile.enable)
    {
        A2profile.rmin = *compNodeTable["A2profile"]["rmin"].value<double>();
        A2profile.rmax = *compNodeTable["A2profile"]["rmax"].value<double>();
        A2profile.binNum =
            *compNodeTable["A2profile"]["binnum"].value<unsigned>();
        if (not(A2profile.rmin >= 0 and A2profile.rmin < A2profile.rmax
                and A2profile.binNum > 0))
        {
            int rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_ERROR(rank,
                      "The parameters for radial A2 profile calculation of "
                      "[%s] is illegal.",
                      compName.data());
            throw;
        };
    }
}

orbit::orbit(toml::table& orbitNode)
{
    // whether enable orbital logs
    m_enable = *orbitNode["enable"].value<bool>();

    if (not m_enable)
    {
        return;
    }

    // period
    m_period = *orbitNode["period"].value<int>();
    if (not(m_period > 0))
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_ERROR(rank, "period must be positive!");
        throw;
    };

    // particle types to be logged
    auto typeIDs = orbitNode["logtypes"];
    if (toml::array* arr = typeIDs.as_array())
    {
        // visitation with for_each() helps deal with heterogeneous data
        arr->for_each([this](auto&& el) {
            if constexpr (toml::is_number<decltype(el)>)
            {
                m_sampleTypes.push_back(*el);
            }
        });
    }

    auto str = *orbitNode["method"].value<string_view>();
    if (str == "txtfile")
    {
        m_method = id_selection_method::TXTFILE;
    }
    else if (str == "random")
    {
        m_method = id_selection_method::RANDOM;
    }
    else
    {
        ERROR("Get an unknown value for [orbital log method]: [%s]",
              str.data());
        ERROR("Must be one 'txtfile' (for a text file of id list) or 'random' "
              "(for random "
              "selection).");
        exit(-1);
    }

    if (m_method == id_selection_method::RANDOM)
    {
        // random selection
        m_fraction = *orbitNode["fraction"].value<double>();
        auto frac  = fraction();
        if (frac <= 0 or frac > 1)
        {
            throw std::runtime_error(
                "The fraction for random selection must be in (0, 1]");
        }
    }
    else
    {
        // By a text file
        const string tmpIdFileName(*orbitNode["idfile"].value<string_view>());
        m_idfile = tmpIdFileName;
    }

    // recenter parameters
    recenter.enable = *orbitNode["recenter"]["enable"].value<bool>();
    if (recenter.enable)
    {
        // particle types to be used as anchors of recenter
        auto typeIDs = orbitNode["recenter"]["anchorids"];
        if (toml::array* arr = typeIDs.as_array())
        {
            // visitation with for_each() helps deal with heterogeneous data
            arr->for_each([this](auto&& el) {
                if constexpr (toml::is_number<decltype(el)>)
                {
                    recenter.anchorIds.push_back(*el);
                }
            });
        }
        recenter.radius = *orbitNode["recenter"]["radius"].value<double>();
        str             = *orbitNode["recenter"]["method"].value<string_view>();
        if (str == "com")
        {
            recenter.method = recenter_method::COM;
        }
        else if (str == "mbp")
        {
            recenter.method = recenter_method::MBP;
        }
        else
        {
            ERROR("Get an unknown value for [recenter method] in orbital logs: "
                  "[%s]",
                  str.data());
            ERROR("Must be 'com' (for center of mass) or 'mbp' (for most bound "
                  "particle).");
            exit(-1);
        }
        auto iguess = orbitNode["recenter"]["iguess"];
        for (auto i = 0; i < 3; ++i)
        {
            recenter.initialGuess[i] = *iguess[i].value<double>();
        }
    }
}

}  // namespace otf
