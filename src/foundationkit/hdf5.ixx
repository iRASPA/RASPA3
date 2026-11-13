module;

//#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(__NT__)
//#define H5_BUILT_AS_DYNAMIC_LIB 0
//#endif

#include <H5Cpp.h>

export module hdf5;

import std;

/*!
 * \brief HDF5 writer class.
 *
 *  Takes care of writing data to hdf5 file. The HDF5 file format can be regarded similar to a linux file system.
 * It has a main group "/" and a set of subgroups "/foo". Each group (including main) can hold attributes, which
 * we refer to as metadata here and datasets. Attributes are key-value pairs and for example the main group "/" has
 * an attribute "RASPA version": "x.x.x". We define a dataset as a vector of data and can hold any type of data.
 * A dataset is first created by giving a name, size and metadata. The metadata given as a vector of key-value pairs
 * holds important information about the dataset and a general guide is to at least include dimensions, for example
 * {"dimensions", "(N, M)"} for a NxM matrix.
 *
 */

export class HDF5Writer
{
 public:
  HDF5Writer(const std::string& filename) : file(filename.c_str(), H5F_ACC_TRUNC) {}
  void createGroup(const std::string& groupName) { file.createGroup(groupName.c_str()); }

  // It is not yet possible with H5Cpp.h to add dimension scales. This should be added in future versions,
  // but for now we should document the dimension names to datasets with multiple dimensions
  // void createDimensionScale(const std::string& groupName, const std::string& scaleName, hsize_t size);

  template <typename T>
  void createDataset(const std::string& groupName, const std::string& datasetName,
                     const std::vector<size_t>& dimensions,
                     const std::vector<std::pair<std::string, std::string>>& metadata)
  {
    H5::Group group = file.openGroup(groupName.c_str());
    std::vector<hsize_t> vec_hsize_t(dimensions.begin(), dimensions.end());
    H5::DataSpace dataspace(static_cast<int>(vec_hsize_t.size()), vec_hsize_t.data());
    H5::DataSet dataset = group.createDataSet(datasetName.c_str(), getH5Type<T>(), dataspace);

    // metadata
    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
    H5::DataSpace attrSpace(H5S_SCALAR);

    for (auto pair : metadata)
    {
      H5::Attribute about = dataset.createAttribute(pair.first, strType, attrSpace);
      about.write(strType, pair.second);
    }
  }

  void createStringDataset(const std::string& groupName, const std::string& datasetName,
                           const std::vector<size_t>& dimensions, size_t maxLength)
  {
    H5::Group group = file.openGroup(groupName.c_str());
    std::vector<hsize_t> vec_hsize_t(dimensions.begin(), dimensions.end());
    H5::DataSpace dataspace(static_cast<int>(vec_hsize_t.size()), vec_hsize_t.data());
    H5::StrType strtype(H5::PredType::C_S1, maxLength);
    H5::DataSet dataset = group.createDataSet(datasetName.c_str(), strtype, dataspace);
  }

  template <typename T>
  void writeVector(const std::string& groupName, const std::string& datasetName, const std::vector<T>& data);

  template <typename T>
  void writeSingleValue(const std::string& groupName, const std::string& datasetName,
                        const std::vector<hsize_t>& indices, T value)
  {
    H5::Group group = file.openGroup(groupName.c_str());
    H5::DataSet dataset = group.openDataSet(datasetName.c_str());
    H5::DataSpace dataspace = dataset.getSpace();

    // one element: 'indices' holds its coordinate, one entry per dimension of the dataset (the selection reads
    // that many coordinates from 'indices', so fewer would be read past its end)
    if (indices.size() != static_cast<std::size_t>(dataspace.getSimpleExtentNdims()))
    {
      throw std::invalid_argument(
          std::format("HDF5Writer: dataset '{}/{}' has {} dimensions, but {} indices were given", groupName,
                      datasetName, dataspace.getSimpleExtentNdims(), indices.size()));
    }
    const hsize_t oneElement = 1;
    H5::DataSpace memspace(1, &oneElement, NULL);
    dataspace.selectElements(H5S_SELECT_SET, 1, indices.data());
    if constexpr (std::is_same_v<T, bool>)
    {
      // the memory type is NATIVE_HBOOL: the buffer is an hbool_t, whatever its size in this HDF5 build
      const hbool_t converted = static_cast<hbool_t>(value);
      dataset.write(&converted, getH5Type<T>(), memspace, dataspace);
    }
    else
    {
      dataset.write(&value, getH5Type<T>(), memspace, dataspace);
    }
  }

  void writeMetaInfo(const std::string& groupName, const std::string& infoName, const std::string& infoValue)
  {
    // Create a string attribute at the root of the file
    H5::Group group = file.openGroup(groupName.c_str());
    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
    H5::DataSpace attrSpace(H5S_SCALAR);
    H5::Attribute attr = group.createAttribute(infoName.c_str(), strType, attrSpace);
    attr.write(strType, infoValue);
  }

  void writeMetaInfo(const std::string& groupName, const std::string& infoName, const double& infoValue)
  {
    // Create a string attribute at the root of the file
    H5::Group group = file.openGroup(groupName.c_str());
    H5::PredType datatype = getH5Type<double>();
    H5::DataSpace attrSpace(H5S_SCALAR);
    H5::Attribute attr = group.createAttribute(infoName.c_str(), datatype, attrSpace);
    attr.write(datatype, &infoValue);
  }

  void writeMetaInfo(const std::string& groupName, const std::string& infoName, const size_t& infoValue)
  {
    // Create a string attribute at the root of the file
    H5::Group group = file.openGroup(groupName.c_str());
    H5::PredType datatype = getH5Type<size_t>();
    H5::DataSpace attrSpace(H5S_SCALAR);
    H5::Attribute attr = group.createAttribute(infoName.c_str(), datatype, attrSpace);
    attr.write(datatype, &infoValue);
  }

 private:
  H5::H5File file;

  // a whole-dataset write reads as many elements from the buffer as the dataset holds: a shorter vector would be
  // read past its end
  void checkElementCount(const std::string& groupName, const std::string& datasetName, const H5::DataSpace& dataspace,
                         std::size_t count)
  {
    const std::size_t points = static_cast<std::size_t>(dataspace.getSimpleExtentNpoints());
    if (count != points)
    {
      throw std::invalid_argument(std::format("HDF5Writer: dataset '{}/{}' holds {} elements, but {} were given",
                                              groupName, datasetName, points, count));
    }
  }

  template <typename T>
  H5::PredType getH5Type()
  {
    if (std::is_same<T, int>::value)
    {
      return H5::PredType::NATIVE_INT;
    }
    else if (std::is_same<T, float>::value)
    {
      return H5::PredType::NATIVE_FLOAT;
    }
    else if (std::is_same<T, double>::value)
    {
      return H5::PredType::NATIVE_DOUBLE;
    }
    else if (std::is_same<T, bool>::value)
    {
      return H5::PredType::NATIVE_HBOOL;
    }
    else if (std::is_same<T, size_t>::value)
    {
      static_assert(sizeof(std::size_t) == 8 || sizeof(std::size_t) == 4);
      return sizeof(std::size_t) == 8 ? H5::PredType::NATIVE_UINT64 : H5::PredType::NATIVE_UINT32;
    }
    else
    {
      throw std::runtime_error("HDF5Writer: unsupported data type");
    }
  }
};

template <typename T>
void HDF5Writer::writeVector(const std::string& groupName, const std::string& datasetName, const std::vector<T>& data)
{
  H5::Group group = file.openGroup(groupName.c_str());
  H5::DataSet dataset = group.openDataSet(datasetName.c_str());
  H5::DataSpace dataspace = dataset.getSpace();
  H5::PredType datatype = getH5Type<T>();
  checkElementCount(groupName, datasetName, dataspace, data.size());
  dataset.write(data.data(), datatype);
}

template <>
void HDF5Writer::writeVector<bool>(const std::string& groupName, const std::string& datasetName, const std::vector<bool>& data)
{
  H5::Group group = file.openGroup(groupName.c_str());
  H5::DataSet dataset = group.openDataSet(datasetName.c_str());
  H5::DataSpace dataspace = dataset.getSpace();
  H5::PredType datatype = getH5Type<bool>();
  checkElementCount(groupName, datasetName, dataspace, data.size());

  // the memory type is NATIVE_HBOOL, so the buffer holds hbool_t, whatever its size in this HDF5 build (an array,
  // since hbool_t may be bool, whose std::vector has no contiguous buffer)
  std::unique_ptr<hbool_t[]> converted = std::make_unique_for_overwrite<hbool_t[]>(data.size());
  std::transform(data.begin(), data.end(), converted.get(), [](bool b) { return static_cast<hbool_t>(b); });
  dataset.write(converted.get(), datatype);
}

template <>
void HDF5Writer::writeVector<std::string>(const std::string& groupName, const std::string& datasetName,
                              const std::vector<std::string>& data)
{
  H5::Group group = file.openGroup(groupName.c_str());
  H5::DataSet dataset = group.openDataSet(datasetName.c_str());
  H5::DataSpace dataspace = dataset.getSpace();
  H5::StrType strtype = dataset.getStrType();
  size_t maxLength = strtype.getSize();
  checkElementCount(groupName, datasetName, dataspace, data.size());

  std::vector<char> buffer(data.size() * maxLength, '\0');
  for (size_t i = 0; i < data.size(); ++i)
  {
    std::strncpy(&buffer[i * maxLength], data[i].c_str(), maxLength);
  }

  dataset.write(buffer.data(), strtype);
}
