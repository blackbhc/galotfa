// The file define the default value of the parameters
#ifndef GALOTFA_DEFAULT_H
#define GALOTFA_DEFAULT_H
#include "ini_parser.h"
#include <string>
#include <vector>
using std::vector;
namespace galotfa {
struct para
{
    // the switchs of each part
    bool glb_switch_on = true, md_switch_on = false, ptc_switch_on = false, orb_switch_on = false,
         grp_switch_on = false, post_switch_on = false;

    // the file/path of each part
    std::string glb_output_dir = "./otfoutput", md_filename = "model", ptc_filename = "particle",
                orb_filename = "orbit", grp_filename = "group", post_filename = "post";

    // the analysis periods
    int md_period = 10, ptc_period = 10000, orb_period = 1, grp_period = 10;

    // other global parameters
    std::string glb_convergence_type = "absolute", glb_sim_type = "galaxy";
    int         glb_pot_tracer = -10086, glb_max_iter = 25;
    double      glb_convergence_threshold = 0.001, glb_equal_threshold = 1e-10;

    // pre section parameters
    bool          pre_recenter    = true;
    double        pre_region_size = 20, pre_axis_ratio = 1;
    std::string   pre_region_shape = "cylinder", pre_recenter_method = "density";
    vector< int > pre_recenter_anchors;

    // other model section parameters
    bool md_image = false, md_bar_major_axis = false, md_bar_radius = false, md_sbar = false,
         md_sbuckle = false, md_inertia_tensor = false, md_align_bar = true,
         md_dispersion_tensor           = false;
    bool                    md_multiple = false;
    vector< int >           md_particle_types;
    vector< std::string >   md_classification;
    vector< vector< int > > md_target_sets;
    int                     md_image_bins = 100, md_rbins = 20;
    double md_region_size = 20, md_axis_ratio = 1, md_bar_threshold = 0.15, md_rmin = 0.0,
           md_rmax = 0.0, md_deg = 3, md_percentage = 70;
    std::string           md_region_shape = "cylinder";
    vector< int >         md_an;  // An, lowercase for ini key
    vector< std::string > md_colors;

    // other particle section parameters
    bool ptc_circularity = false, ptc_circularity_3d = false, ptc_rg = false, ptc_freq = false;
    vector< int > ptc_particle_types;

    // other orbit section parameters
    std::string orb_idfile;

    // other group section parameters
    bool                  grp_vmg = false, grp_rmg = false, grp_ellipticity = false;
    vector< std::string > grp_group_types;

    // other post section parameters
    bool post_pattern_speed = false;

    // constructor and destructor
    para( ini_parser& parser );
    int        check( void );  // the function to check the dependencies between the parameters
    inline int target_sets_parser();  // parse the multiple target sets
#ifdef debug_parameter
    int test_print();
#endif
};
}  // namespace galotfa
#endif
