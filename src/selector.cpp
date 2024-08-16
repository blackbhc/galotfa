#include "../include/selector.hpp"
#include "../include/myprompt.hpp"
#include "../include/para.hpp"
#include "mpi.h"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iterator>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <sys/unistd.h>
#include <unistd.h>
#include <vector>
using namespace std;

namespace otf {

/**
 * @brief Select particle ids with a specified fraction, can be used in any mpi rank.
 *
 * @param raw raw id list
 * @param fraction fraction
 * @return a vector<unsigned int> of the selected id list
 */
auto orbit_selector::id_select( const vector< unsigned int >& raw,
                                const double                  fraction ) -> vector< unsigned int >
{
    if ( fraction <= 0 or fraction > 1 )
    {
        ERROR( "Try to select a illegal fraction: [%lf]", fraction );
    }
    else if ( fraction == 1 )
    {
        return raw;
    }

    auto const             selectNum = ( size_t )( raw.size() * fraction );
    vector< unsigned int > res;
    sample( raw.begin(), raw.end(), back_inserter( res ), selectNum, mt19937{ random_device{}() } );
    return res;
}

/**
 * @brief Read the ids in the given id list file, can be used in any mpi rank to read the ids into
 * it.
 *
 * @param idFilename filename of the id list file (in ASCII txt format)
 * @return std::vector<int> of the ids.
 */
auto orbit_selector::id_read( const string& idFilename ) -> vector< unsigned int >
{
    // check the availability of the file
    if ( access( idFilename.c_str(), F_OK ) != 0 )
    {
        int rank;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        MPI_ERROR( rank, "Particle ID file not found: [%s]!", idFilename.c_str() );
        exit( -1 );
    }

    ifstream               fp;
    vector< unsigned int > ids;
    fp.open( idFilename.c_str(), ios::in );
    while ( !fp.eof() )
    {
        string lineStr;
        getline( fp, lineStr );
        if ( !lineStr.empty() )
        {
            try
            {
                ids.push_back( std::stoi( lineStr ) );
            }
            catch ( ... )
            {
                ERROR( "Get an unexpected value: %s", lineStr.c_str() );
                exit( -1 );
            }
        }
    }

    // Remove the possible repeated values
    set< int > tmpSet( ids.begin(), ids.end() );
    ids.assign( tmpSet.begin(), tmpSet.end() );

    return ids;
}

orbit_selector::orbit_selector( unique_ptr< runtime_para >& para ) : para( para )
{
    ;
}

auto orbit_selector::select( const unsigned int particleNumber, const unsigned int* particleIDs,
                             const unsigned int* particleTypes, const double* masses,
                             const double* coordiantes,
                             const double* velocities ) -> const std::unique_ptr< dataContainer >
{
    if ( not para->orbit->enable )
    {
        return nullptr;
    }

    unsigned int counter = 0;

    vector< double > tmpMass( particleNumber );
    vector< double > tmpPos( particleNumber );
    vector< double > tmpVel( particleNumber );

    ( void )particleIDs;
    ( void )particleTypes;
    ( void )masses;
    ( void )coordiantes;
    ( void )velocities;
    for ( auto i = 0U; i < particleNumber; ++i )
    {
        counter++;
    }

    tmpMass.resize( counter );
    tmpPos.resize( counter );
    tmpVel.resize( counter );

    // get the data container
    unique_ptr< dataContainer > container;
    container->count      = counter;
    container->mass       = std::move( tmpMass );
    container->coordinate = std::move( tmpPos );
    container->velocity   = std::move( tmpVel );

    return container;
}

}  // namespace otf
