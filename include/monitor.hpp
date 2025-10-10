/**
 * @file monitor.hpp
 * @brief The monitor/server of on-the-fly analysis.
 */

#ifndef MONITOR_HEADER
#define MONITOR_HEADER
#include "../include/h5out.hpp"
#include "../include/para.hpp"
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace otf {

/**
 * @class monitor
 * @brief The main server of the on-the-fly analysis, which organize the other modules to work
 * together
 *
 */
class monitor
{
public:
    monitor( const std::string_view& tomlParaFile );  // initialize with the toml file name
    ~monitor();
    void main_analysis_api( double time, unsigned particleNumber, const int* id,
                            const int* partTypes, const double* masses, const double* potentials,
                            const double* coordinates, const double* velocities );

#ifdef DEBUG

#else
private:
#endif
    int          mpiRank;
    int          mpiSize;
    bool         isRootRank;
    unsigned     stepCounter;             // counter of the synchronized time step
    bool         mpiInitialzedByMonitor;  // whether the MPI_init is called by the monitor object
    runtime_para para;                    // ptr to the runtime paramter
    std::vector< std::string > orbitDatasetNames;
    static constexpr unsigned  orbitPointDim = 7;
    using orbitPoint                         = struct
    {
        int    particleID;
        double data[ orbitPointDim ];
    };

    // the container of data for a single component
    using compDataContainer = struct compDataStruct
    {
        unsigned                    partNum     = 0;        // number of particles in this component
        std::unique_ptr< double[] > masses      = nullptr;  // masses of particles
        std::unique_ptr< double[] > potentials  = nullptr;  // potentials of particles
        std::unique_ptr< double[] > coordinates = nullptr;  // coordinates of particles
        std::unique_ptr< double[] > velocities  = nullptr;  // velocities of particles
    };

    // the container of analysis results for a single component
    using compResContainer = struct compResStruct
    {
        double                      center[ 3 ] = { 0, 0, 0 };  // center of the component
        double                      sBar        = 0;            // bar strength parameter
        double                      barAngle    = 0;            // bar angle
        double                      sBuckle     = 0;            // buckling strength
        unsigned                    imageBinNum = 0;            // image matrix rank
        std::unique_ptr< double[] > imageXY     = nullptr;      // image matrix x-y
        std::unique_ptr< double[] > imageXZ     = nullptr;      // image matrix x-z
        std::unique_ptr< double[] > imageYZ     = nullptr;      // image matrix y-z
        // For radial A2 profile
        std::unique_ptr< double[] > A2Re = nullptr;  // real parts of the radial A2 profile
        std::unique_ptr< double[] > A2Im = nullptr;  // imaginary parts of the radial A2 profile
    };

    // extract the data used for orbital log
    auto id_data_process( double time, unsigned particleNumber, const int* particleIDs,
                          const int* particleTypes, const double* masses, const double* coordinates,
                          const double* velocities ) const -> std::vector< monitor::orbitPoint >;
    // NOTE: API of orbital log
    void orbital_log( double time, unsigned particleNumber, const int* ids, const int* partTypes,
                      const double* masses, const double* coordinates, const double* velocities );
    // extract the data of a single component
    static auto component_data_extract( unsigned particleNumber, const int* partTypes,
                                        const double* masses, const double* potentials,
                                        const double* coordinates, const double* velocities,
                                        std::unique_ptr< otf::component >& comp )
        -> monitor::compDataContainer;
    // analyze the data of a single component
    auto component_data_analyze( monitor::compDataContainer&        dataContainer,
                                 std::unique_ptr< otf::component >& comp ) const
        -> monitor::compResContainer;
    // NOTE: API to analyze the data of a single component
    void component_analysis( double time, unsigned particleNumber, const int* partTypes,
                             const double* masses, const double* potentials,
                             const double* coordinates, const double* velocities,
                             std::unique_ptr< otf::component >& comp ) const;

    // NOTE: APIs used in component analysis

    // recenter the coordinates
    static void recenter_coordinate( monitor::compDataContainer&        dataContainer,
                                     std::unique_ptr< otf::component >& comp,
                                     compResContainer&                  res );
    // align the z axis to the direction of total angular momentum
    static void align_angular_momentum( monitor::compDataContainer&        dataContainer,
                                        std::unique_ptr< otf::component >& comp );
    // align the coordinates to the eigenvalues of the inertia tensor
    // NOTE: This is the legacy align_coordinate API (replaced by align_angular_momentum in current
    // version), and the code is remained for possible future usage on bar alignment.
    static void align_inertia_tensor( monitor::compDataContainer&        dataContainer,
                                      std::unique_ptr< otf::component >& comp );
    // bar info calculation
    static void bar_info( monitor::compDataContainer&        dataContainer,
                          std::unique_ptr< otf::component >& comp, compResContainer& res );
    // radial A2 profile calculation
    void a2_profile( monitor::compDataContainer&        dataContainer,
                     std::unique_ptr< otf::component >& comp, compResContainer& res ) const;
    // image calculation
    void image( monitor::compDataContainer& dataContainer, std::unique_ptr< otf::component >& comp,
                compResContainer& res ) const;
    // the smart pointer to the HDF5 file organizer, it is more memory efficient in the none-root
    // rank
    std::unique_ptr< h5_out > h5Organizer;
};

}  // namespace otf
#endif
