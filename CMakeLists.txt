# GLOBAL PART
# minimal version of cmake
cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
# project name, language
project(galotfa LANGUAGES CXX)
# standard: C++20
set(CMAKE_CXX_STANDARD 20)
# no compiler dependent extensions
set(CMAKE_CXX_EXTENSIONS OFF)
# check the compiler support the specified standard
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Build type, default: Debug
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()
# export the compilation commands
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# COMPILATION FLAGS PART
# customized debug flags: default one only includes -g
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -Wextra -Wpedantic -Werror -fPIC")
# print fundamental flags: -g, -O, and whether -DNDEBUG
message(STATUS "C++ flags, Debug configuration: ${CMAKE_CXX_FLAGS_DEBUG}")
message(STATUS "C++ flags, Release configuration: ${CMAKE_CXX_FLAGS_RELEASE}")
message(
    STATUS
    "C++ flags, Release configuration with Debug info: ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}"
)
message(
    STATUS
    "C++ flags, minimal Release configuration: ${CMAKE_CXX_FLAGS_MINSIZEREL}"
)

# TARGET PART
add_library(galotfa SHARED include/galotfa.hpp src/galotfa.cpp)
target_link_libraries(galotfa PRIVATE gsl gslcblas hdf5)
# Library with the same name
#[[ add_library(message-shared
  SHARED
    $<TARGET_OBJECTS:message-objs>
 )

set_target_properties(message-shared
 PROPERTIES
  OUTPUT_NAME "message"
 )

add_library(message-static
 STATIC
  $<TARGET_OBJECTS:message-objs>
 )

set_target_properties(message-static
 PROPERTIES
  OUTPUT_NAME "message"
 ) ]]

# UNIT TESTS
# C++ programes
enable_testing()
add_executable(
    coordinate
    ./validation/test_coordTrans.cpp
    ./src/coordinates.cpp
)
target_link_options(
    coordinate
    PRIVATE
        -fsanitize=address
        # -fsanitize=leak
        -fsanitize=undefined
)
add_test(NAME test_coordinate COMMAND $<TARGET_FILE:coordinate>)

add_executable(bin2d ./validation/test_bin2d.cpp ./src/statistic.cpp)
target_link_libraries(bin2d PRIVATE gsl gslcblas)
find_package(MPI REQUIRED) # first: find the package
target_link_libraries(bin2d PUBLIC MPI::MPI_CXX)
target_link_options(
    bin2d
    PRIVATE
        -fsanitize=address
        # -fsanitize=leak
        -fsanitize=undefined
)
add_test(NAME test_bin2d COMMAND mpirun -np 4 $<TARGET_FILE:bin2d>)

add_executable(eigen ./validation/test_eigen.cpp ./src/eigen.cpp)
target_link_libraries(eigen PRIVATE gsl gslcblas)
target_link_options(
    eigen
    PRIVATE
        -fsanitize=address
        # -fsanitize=leak
        -fsanitize=undefined
)
add_test(NAME test_eigen COMMAND $<TARGET_FILE:eigen>)

add_executable(h5out ./validation/test_h5out.cpp ./src/h5out.cpp)
target_link_libraries(h5out PRIVATE hdf5)
target_link_options(
    h5out
    PRIVATE
        -fsanitize=address
        # -fsanitize=leak
        -fsanitize=undefined
)
add_test(NAME test_h5out COMMAND $<TARGET_FILE:h5out>)

# INSTALL PART
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "~/.local" CACHE PATH "install prefix" FORCE)
endif()
include(GNUInstallDirs)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR}
)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_BINDIR}
)
# Offer the user the choice of overriding the installation directories
set(INSTALL_LIBDIR
    ${CMAKE_INSTALL_LIBDIR}
    CACHE PATH
    "Installation directory for libraries"
)
set(INSTALL_BINDIR
    ${CMAKE_INSTALL_BINDIR}
    CACHE PATH
    "Installation directory for executables"
)
set(INSTALL_INCLUDEDIR
    ${CMAKE_INSTALL_INCLUDEDIR}
    CACHE PATH
    "Installation directory for header files"
)
set(DEF_INSTALL_CMAKEDIR CMake)
set(INSTALL_CMAKEDIR
    ${DEF_INSTALL_CMAKEDIR}
    CACHE PATH
    "Installation directory for CMake files"
)
# Report to user
foreach(p LIB BIN INCLUDE CMAKE)
    file(TO_NATIVE_PATH ${CMAKE_INSTALL_PREFIX}/${INSTALL_${p}DIR} _path)
    message(STATUS "Installing ${p} components to ${_path}")
    unset(_path)
endforeach()
install(
    TARGETS galotfa
    ARCHIVE DESTINATION ${INSTALL_LIBDIR} COMPONENT lib
    RUNTIME DESTINATION ${INSTALL_BINDIR} COMPONENT bin
    LIBRARY DESTINATION ${INSTALL_LIBDIR} COMPONENT lib
    PUBLIC_HEADER DESTINATION ${INSTALL_INCLUDEDIR}/message COMPONENT dev
)
# Offer the user the choice of overriding the installation directories
set(INSTALL_LIBDIR
    ${CMAKE_INSTALL_LIBDIR}
    CACHE PATH
    "Installation directory for libraries"
)
set(INSTALL_BINDIR
    ${CMAKE_INSTALL_BINDIR}
    CACHE PATH
    "Installation directory for executables"
)
set(INSTALL_INCLUDEDIR
    ${CMAKE_INSTALL_INCLUDEDIR}
    CACHE PATH
    "Installation directory for header files"
)

#[[
Others
1. Run shell commands: see the blogs at
https://blog.csdn.net/guotianqing/article/details/121661067

2. MPI programe
find_package(MPI REQUIRED) # first: find the package
add_executable(hello-mpi hello-mpi.cpp)
target_link_libraries(hello-mpi
  PUBLIC
 	  MPI::MPI_CXX
  )

3. Find all files (non-recursive) under a directory:
file(GLOB SOURCE_FILES ./src/*.cpp)
message(
  STATUS
    "All files under src: ${SOURCE_FILES}")

4. Unit test with python script:
find_package (Python COMPONENTS Interpreter Development)
message(
  STATUS
    "Find python3 interpreter: ${Python_EXECUTABLE}")
add_test(
  NAME python_test
  COMMAND ${Python_EXECUTABLE} ../validation/test_hello.py
  )
]]
