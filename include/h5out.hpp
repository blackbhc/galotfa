/**
 * @file h5out.hpp
 * @brief Organizer of the hdf5 output.
 */

#ifndef MY_H5_OUTPUT_HEADER
#define MY_H5_OUTPUT_HEADER
#include "H5Ipublic.h"
#include "H5public.h"
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#ifdef DEBUG
#define CHUCK_SIZE 7
#else
#define CHUCK_SIZE 1024
#endif
/**
 * @class dataset_handle
 * @brief The class that handles the basic operations of a dataset, including
 * memory space, chuck size, dataset compression, data flush etc.
 *
 */
class dataset_handle
{
public:
    dataset_handle(hid_t&                       parent,
                   const std::string&           datasetName,
                   const std::vector<unsigned>& sizeInEachDim,
                   hid_t&                       dataType);
    ~dataset_handle();
    auto flush_single_block(const void* dataBuffer) -> int;

#ifdef DEBUG

#else
private:
#endif
    std::string datasetName;
    hid_t       singleBlockSpace =
        H5I_INVALID_HID;  // memory space id for single step
    hid_t                      dataset  = H5I_INVALID_HID;
    hid_t                      property = H5I_INVALID_HID;
    hid_t                      dataType = H5I_INVALID_HID;
    std::vector<unsigned>      sizeInEachDim;
    std::unique_ptr<hsize_t[]> fileSize;
    std::unique_ptr<hsize_t[]> countOfSingleBlock;
    std::unique_ptr<hsize_t[]> offset;
    // NOTE: curIndex always points to the current to be logged index
    unsigned long long         curIndex = 0;
};

/**
 * @class h5_out
 * @brief Handle the data output, which should be used only in the main process.
 *
 */
class h5_out
{
public:
    h5_out(const std::string& dir, const std::string& filename);
    ~h5_out();
    auto create_dataset_in_group(const std::string&           datasetName,
                                 const std::string&           groupName,
                                 const std::vector<unsigned>& sizeInEachDim,
                                 hid_t                        dataType) -> int;
    auto flush_single_block(const std::string& groupName,
                            const std::string& datasetName,
                            const void*        dataBuffer) -> int;

#ifdef DEBUG
public:
#else
private:
#endif
    hid_t                                  m_file;
    std::unordered_map<std::string, hid_t> m_groups;  // group handles

    // dataset handles under a parent node
    // relation: group -> dataset name -> pointer of a handle
    std::unordered_map<
        hid_t,
        std::unordered_map<std::string, std::unique_ptr<dataset_handle>>>
                m_datasetPtrs;
    auto        create_group_if_necessary(const std::string& groupName) -> int;
    auto        ensure_dataset_empty(const std::string& groupName,
                                     const std::string& datasetName) -> int;
    std::string m_filename;
};

void Backup_Old_Logs_If_Necessary(const std::string& dir,
                                  const std::string& filename);
#endif
