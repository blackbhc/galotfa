/** @mainpage
 *
 * @author Bin-Hui Chen (GitHub: https://github.com/blackbhc)
 *
 * @section over Overview
 * This page explains the code structure of the galotfa project.
 *
 * - The C-style API be called at the synchronized time steps of N-body
 * simulations.
 * -# include/galotfa.h
 * - The server of the on-the-fly analysis eigen, which calls the orbital
 * analysis and component analysis modules implemented in other files to execute
 * the analysis.
 * -# include/monitor.hpp
 * - The hdf5 output organizer.
 * -# include/h5out.hpp Organize the hdf5 file output.
 * - The runtime parameters organizer.
 * -# include/toml.hpp The toml parser taken from toml++
 * (https://marzer.github.io/tomlplusplus/#mainpage-example)
 * -# include/para.hpp Organize the toml file parsing and restoring of runtime
 * parameters.
 * - The orbital analysis module.
 * -# include/selector.hpp Select the particles to be logged.
 * - The component analysis module.
 * -# include/recenter.hpp Recenter the coordinates.
 * -# include/barinfo.hpp Calculate the A2, Sbuckle and other bar informations.
 * -# include/statistic.hpp Calculate the 1D or 2D binning statistics.
 * -# include/eigen.hpp Calculate the eigenvalues and eigenvectors of a given
 * matrix.
 *
 * <hr>
 * @todo Bar length calculation.
 * @todo Toomre Q parameter calculation.
 */
