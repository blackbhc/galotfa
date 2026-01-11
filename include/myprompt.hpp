/**
 * @file myprompt.hpp
 * @brief Define some utility functions for convenient print.
 */

#ifndef MY_UTILS_HEADER
#define MY_UTILS_HEADER
#include <cstdio>
#include <mpi.h>

// print with a line break at the end
#define myprint(...)              \
    {                             \
        std::printf(__VA_ARGS__); \
        std::printf("\n");        \
    }

// similar but for fprint
#define myfprint(file_ptr, ...)              \
    {                                        \
        std::fprintf(file_ptr, __VA_ARGS__); \
        std::fprintf(file_ptr, "\n");        \
    }
// similar but used for mpi case
#define mpi_print(rank, ...)          \
    {                                 \
        if (rank == 0)                \
        {                             \
            std::printf(__VA_ARGS__); \
            std::printf("\n");        \
        }                             \
    }
#define mpi_fprint(rank, file_prt, ...)          \
    {                                            \
        if (rank == 0)                           \
        {                                        \
            std::fprintf(file_ptr, __VA_ARGS__); \
            std::fprintf(file_prt, "\n");        \
        }                                        \
    }

#define test_reach myprint("Call from file [%s] line[%d].", __FILE__, __LINE__);
#define mpi_test_reach                                                      \
    {                                                                       \
        int rank;                                                           \
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);                               \
        myprint("Rank [%d]: Call from file [%s] line[%d].", rank, __FILE__, \
                __LINE__);                                                  \
    }

// for info, warning, and error.
#define INFO(...) myprint("[INFO]: " __VA_ARGS__);
#define WARN(...)                                                              \
    {                                                                          \
        myfprint(stderr, "Call from file [%s] line[%d].", __FILE__, __LINE__); \
        myfprint(stderr, "[WARN]: " __VA_ARGS__);                              \
    }
#define ERROR(...)                                                             \
    {                                                                          \
        myfprint(stderr, "Call from file [%s] line[%d].", __FILE__, __LINE__); \
        myfprint(stderr, "[ERROR]: " __VA_ARGS__);                             \
    }
#define MPI_INFO(rank, ...) mpi_print(rank, "[INFO]: " __VA_ARGS__);
#define MPI_WARN(rank, ...)                                                \
    if (rank == 0)                                                         \
    {                                                                      \
        std::fprintf(stderr, "Rank [%d]: Call from file [%s] line[%d].\n", \
                     rank, __FILE__, __LINE__);                            \
        std::fprintf(stderr, "[WARN]: " __VA_ARGS__);                      \
        std::fprintf(stderr, "\n");                                        \
    }
#define MPI_ERROR(rank, ...)                                               \
    if (rank == 0)                                                         \
    {                                                                      \
        std::fprintf(stderr, "Rank [%d]: Call from file [%s] line[%d].\n", \
                     rank, __FILE__, __LINE__);                            \
        std::fprintf(stderr, "[ERROR]: " __VA_ARGS__);                     \
        std::fprintf(stderr, "\n");                                        \
    }
#endif
