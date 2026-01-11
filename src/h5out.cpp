#include "h5out.hpp"
#include "myprompt.hpp"
#include <H5Dpublic.h>
#include <H5Fpublic.h>
#include <H5Gpublic.h>
#include <H5Ipublic.h>
#include <H5Ppublic.h>
#include <H5Spublic.h>
#include <H5public.h>
#include <algorithm>
#include <cstdio>
#include <memory>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <utility>
#include <vector>
using namespace std;

/**
 * @brief Make backups of previous on-the-fly logs under the same directory.
 *
 * @param dir path the output, for which do not add a "/" at the end
 * @param filename hdf5 filename of the logs
 */
void Backup_Old_Logs_If_Necessary(const string& dir, const string& filename)
{
    string const path_to_file = dir + "/" + filename;
    // check whether the file exists, if so, make a backup for it.
    if (access(path_to_file.c_str(), F_OK) == 0)
    {
        // get the minimal prefix that avoids overwrite
        int    min_prefix = 1;
        string safe_path_to_file =
            dir + "/" + "bak" + to_string(min_prefix) + "-" + filename;
        while (access(safe_path_to_file.c_str(), F_OK) == 0)
        {
            ++min_prefix;
            safe_path_to_file = dir + "/bak" + to_string(min_prefix) + "-";
            safe_path_to_file += filename;
        }
        // backup the file
        rename(path_to_file.c_str(), safe_path_to_file.c_str());
    }
}

/**
 * @brief Create the given directory if it doesn't exist.
 *
 * @param dir
 */
void Create_Dir_If_Necessary(const string& dir)
{
    // check whether the directory exists, if not, create it.
    if (access(dir.c_str(), F_OK) != 0)
    {
        mkdir(dir.c_str(), 0755);
    }
}

h5_out::h5_out(const string& dir, const string& filename)
{
    Create_Dir_If_Necessary(dir);
    Backup_Old_Logs_If_Necessary(dir, filename);
    this->filename = dir + "/" + filename;
    this->file = H5Fcreate(this->filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                           H5P_DEFAULT);
}

h5_out::~h5_out()
{
    datasetPtrs.clear();

    // Second: close the groups
    for (auto& group : groups)
    {
        H5Gclose(group.second);
    }
    groups.clear();

    // Last: close the file
    H5Fclose(file);
}

/**
 * @brief Create the dataset in a group, based on the give parameters. Note that
 * the function can not create a dataset in the root group!
 *
 * @param datasetName string to the name of the dataset, including no "/"
 * @param groupName similar to datasetName but for the group name.
 * @param sizeInEachDim the length of each dimension in the dataset, except the
 * first extensible one.
 * @param dataType
 * @return
 */
auto h5_out::create_dataset_in_group(const string&           datasetName,
                                     const string&           groupName,
                                     const vector<unsigned>& sizeInEachDim,
                                     hid_t                   dataType) -> int
{
    // ensure the parent group exists
    int const returnCode = create_group_if_necessary(groupName);
    // get the group id
    hid_t     groupId = groups[groupName];
    if (returnCode != 0)
    {
        return returnCode;
    }

    // create the group based on the specified info
    int const setExistsInGroup = ensure_dataset_empty(groupName, datasetName);
    if (setExistsInGroup != 0)
    {
        WARN("Try to create an existing dataset [%s] in group [%s]!",
             datasetName.c_str(), groupName.c_str());
        return -1;
    }
    unique_ptr<dataset_handle> ptrToHandle = make_unique<dataset_handle>(
        groupId, datasetName, sizeInEachDim, dataType);
    datasetPtrs[groupId][datasetName] = std::move(ptrToHandle);

    return 0;
}

/**
 * @brief Check whether the dataset already exists.
 *
 * @param groupName
 * @param datasetName
 * @return
 */
auto h5_out::ensure_dataset_empty(const string& groupName,
                                  const string& datasetName) -> int
{
    hid_t const groupId  = groups[groupName];
    bool const  nonEmpty = any_of(
        datasetPtrs[groupId].begin(), datasetPtrs[groupId].end(),
        [datasetName](const auto& pair) { return pair.first == datasetName; });
    if (nonEmpty)
    {
        return 1;
    }
    return 0;
}

/**
 * @brief Create the group if necessary.
 *
 * @param groupName
 * @return
 */
auto h5_out::create_group_if_necessary(const string& groupName) -> int
{
    // check whether the group already exists, if exists then return
    if (static_cast<unsigned>(groups.contains(groupName)) != 0U)
    {
        return 0;
    }

    hid_t const group = H5Gcreate2(file, groupName.c_str(), H5P_DEFAULT,
                                   H5P_DEFAULT, H5P_DEFAULT);
    if (group == H5I_INVALID_HID)  // if creation failed
    {
        ERROR("File [%s]: group [%s] creatation failed.", filename.c_str(),
              groupName.c_str());
        return -1;
    }

    groups[groupName] = group;


    return 0;
}

dataset_handle::dataset_handle(hid_t&                  parent,
                               const string&           datasetName,
                               const vector<unsigned>& sizeInEachDim,
                               hid_t&                  dataType)
    : datasetName(datasetName), dataType(dataType), sizeInEachDim(sizeInEachDim)
{
    // set the chunk size and compression at here
    auto                        rank = ( int )sizeInEachDim.size() + 1;
    unique_ptr<hsize_t[]> const chunk(new hsize_t[rank]());
    fileSize           = make_unique<hsize_t[]>(rank);
    countOfSingleBlock = make_unique<hsize_t[]>(rank);
    offset             = make_unique<hsize_t[]>(rank);
    unique_ptr<hsize_t[]> const maxDims(new hsize_t[rank]());
    fileSize[0]           = 0;
    countOfSingleBlock[0] = 1;
    for (auto i = 0; i < rank; ++i)
    {
        offset[i] = 0;
    }
    chunk[0]   = CHUCK_SIZE;
    maxDims[0] = H5S_UNLIMITED;

    for (auto i = 1; i < rank; ++i)
    {
        countOfSingleBlock[i] = fileSize[i] = chunk[i] = maxDims[i] =
            sizeInEachDim[i - 1];
    }

    // create property list and set chunk and compression
    property = H5Pcreate(H5P_DATASET_CREATE);
    // set chunk
    herr_t status = H5Pset_chunk(property, rank, chunk.get());
    if (status < 0)
    {
        ERROR("Failed to set chuck!");
        return;
    }
    // set compression
    status = H5Pset_deflate(property, 6);
    if (status < 0)
    {
        ERROR("Failed to set compression!");
        return;
    }

    // create the memory space for single-step
    fileSize[0]      = 1;
    singleBlockSpace = H5Screate_simple(rank, fileSize.get(), maxDims.get());

    // create the zero-size dataspace for dataset creation
    fileSize[0] = 0;
    hid_t const zeroSizeSpace =
        H5Screate_simple(rank, fileSize.get(), maxDims.get());

    // create the dataset
    dataset = H5Dcreate2(parent, datasetName.c_str(), dataType, zeroSizeSpace,
                         H5P_DEFAULT, property, H5P_DEFAULT);
}

/**
 * @brief A wrapper to log the analyses results in single-step to the log file,
 * without boundary check.
 *
 * @param datasetName
 * @param groupName
 * @param dataBuffer pointer to the data buffer.
 * @return
 */
auto h5_out::flush_single_block(const string& groupName,
                                const string& datasetName,
                                const void*   dataBuffer) -> int
{
    if (static_cast<unsigned>(groups.contains(groupName)) == 0U)
    {
        ERROR("Flushing to a non existent group [%s]!", groupName.c_str())
        return -1;
    }

    hid_t const group = groups[groupName];
    if (static_cast<unsigned>(datasetPtrs[group].contains(datasetName)) == 0U)
    {
        ERROR("Flushing to a non existent dataset [%s] in group [%s]!",
              datasetName.c_str(), groupName.c_str())
        return -1;
    }

    datasetPtrs[group][datasetName]->flush_single_block(dataBuffer);
    return 0;
}

dataset_handle::~dataset_handle()
{
    // remove the additional garbage values in the dataset
    fileSize[0]         = curIndex;
    herr_t const status = H5Dset_extent(dataset, fileSize.get());
    ( void )status;
    H5Dflush(dataset);

    // flush the data buffers and close the dataset
    if (dataset != H5I_INVALID_HID)
    {
        H5Dflush(dataset);
        H5Dclose(dataset);
    }
    // close the used objects
    if (property != H5I_INVALID_HID)
    {
        H5Pclose(property);
    }
    if (singleBlockSpace != H5I_INVALID_HID)
    {
        H5Sclose(singleBlockSpace);
    }
}

/**
 * @brief A wrapper to log the analyses results to the log file, without
 * boundary check.
 *
 * @param dataBuffer pointer to the data buffer, which includes the data in a
 * single-step.
 * @return
 */
auto dataset_handle::flush_single_block(const void* dataBuffer) -> int
{
    // static variables
    static herr_t status    = -1;
    static hid_t  fileSpace = H5I_INVALID_HID;

    // extend the dataspace at the begin or any time reach the end of a chuck
    if (curIndex % CHUCK_SIZE == 0)
    {
        // calculate the number of used chucks
        static unsigned long long chuckNum = 0;
        chuckNum                           = curIndex / CHUCK_SIZE;

        // extend the file size
        fileSize[0] = (chuckNum + 1) * CHUCK_SIZE;
        status      = H5Dextend(dataset, fileSize.get());
        H5Dflush(dataset);

        // check whether succeed or not
        if (status < 0)
        {
            ERROR("Extention of the dataset [%s] failed, there may be no "
                  "enough memory!",
                  datasetName.c_str());
            return -1;
        }
    }

    // write data to a single block
    fileSpace = H5Dget_space(dataset);  // get file space
    offset[0] = curIndex;               // get offset
    // select the subset to be written
    status = H5Sselect_hyperslab(fileSpace, H5S_SELECT_SET, offset.get(),
                                 nullptr, countOfSingleBlock.get(), nullptr);
    status = H5Dwrite(dataset, dataType, singleBlockSpace, fileSpace,
                      H5P_DEFAULT, dataBuffer);
    H5Sclose(fileSpace);

    ++curIndex;  // update the currect index in the file space

    if (status < 0)
    {
        ERROR("Dataset write failed!");
        return -1;
    }

    return 0;
}
