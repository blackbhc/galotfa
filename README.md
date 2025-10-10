# Overview

`galotfa`: galactic on-the-fly analysis, a library for on-the-fly analysis of
galactic simulations.

At present, it is mainly focused on disk or bar galaxies. This library is implemented in C++ and offers a C API. It uses `HDF5` for data storage and `MPI` for parallelization.

---

## Dependencies

- A C++ compiler that supports C++20.

- [`HDF5`](https://www.hdfgroup.org/solutions/hdf5/)

- [`gsl`](https://www.gnu.org/software/gsl/)

- A `MPI` library (e.g. [`OpenMPI`](https://www.open-mpi.org/), [`MPICH`](https://www.mpich.org/)).

- [`CMake`](https://cmake.org/) >= 3.12

- [`toml++`](https://marzer.github.io/tomlplusplus/#mainpage-example) for
  `toml` parser, open source with MIT license, already included in this project.

---

## Installation

The installation steps are straightforward if you have installed all the dependent software and libraries. Just follow the below steps:

1. Clone the repository: `git clone https://github.com/blackbhc/galotfa.git`
2. Get into the directory: `cd galotfa`
3. Configuration of the building:

   ```shell
   cmake -S . -DCMAKE_BUILD_TYPE=Release -B build -DCMAKE_INSTALL_PREFIX=<dir you want to install with>
   ```

   The last `CMAKE_INSTALL_PREFIX` can be ignored, if you want to install `galotfa` under `~/.local/`.

4. Compile and install the library: `cmake --build build --target install`, then everything done.

5. If you encounter any error during the above steps, just run `rm -rf ./build/` and go back to Step 3 to reinstall the library.

6. If you want to uninstall the `galotfa`, you can run `cat ./build/install_manifest.txt` to see all installed files and remove them manually.

---

## Usage

1. Modify your simulation program:
   - Use my `gadget4` fork with `galotfa`

     `git clone https://github.com/blackbhc/gadget4-fork.git ./gadget4`

     And then compile it as a conventional `gadget4`, except open the `GALOTFA_ON` option in
     the configuration file `Config.sh`.

     <font color="red">Note: the compilation of the simulation programs is independent
     of the `galotfa` installation.</font>

   - Or from scratch

     To use `galotfa` in your N-body simulation, you need to modify your simulation program.
     Call the

     ```C
     OnTheFly_Analysis_Nbody( const double currentTime, const unsigned particleNumber,
                              const int* particleIDs, const int* particleTypes,
                              const double* masses, const double* potentials,
                              const double* coordinates, const double* velocities )
     ```

     API in the main loop (synchronized time steps) of your simulation program, and offer the
     require quantities. And then recompile your program to enable `galotfa`.

2. Setup the runtime parameters of `galotfa`.

   You can copy `./examples/galotfa.toml` to the working directory of your simulation, modify
   the parameters in it according to your requirement.
   Note that you need to put the `galotfa.toml` file under the same directory of the
   simulation program (for example `Gadget4`).

3. Run the simulation, and the on-the-fly analysis results is restored in a file specified
   by `galotfa.toml`.

---

## Examples

The structure of the `HDF5` output file produced by `galotfa`:

- `component1`
  - Time
  - A2
  - BarAngle
  - ImageXY
  - ...

- `Orbit`
  - Particle-X

    Restore the orbit of Particle-X as datasets, where X is the id of the particle. Each
    dataset has N x 7 numbers: N is the number of logged steps, 7 for 1 time, 3 coordinates,
    and 3 velocities in Cartesian coordinate frame.

---

## Document for developer

If you want to learn about the programming details, you can use `doxygen` to create the code document.
