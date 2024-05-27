module;

#ifdef USE_LEGACY_HEADERS
#include <iostream>
#include <vector>

#endif

#include "H5Cpp.h"
export module hdf5;

#ifndef USE_LEGACY_HEADERS
import iostream;
import vector;
#endif

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
  HDF5Writer(const std::string& filename) : file(filename, H5F_ACC_TRUNC) {}
  void createGroup(const std::string& groupName) { file.createGroup(groupName); }

  // It is not yet possible with H5Cpp.h to add dimension scales. This should be added in future versions,
  // but for now we should document the dimension names to datasets with multiple dimensions
  // void createDimensionScale(const std::string& groupName, const std::string& scaleName, hsize_t size);

  template <typename T>
  void createDataset(const std::string& groupName, const std::string& datasetName,
                     const std::vector<size_t>& dimensions,
                     const std::vector<std::pair<std::string, std::string>>& metadata)
  {
    H5::Group group = file.openGroup(groupName);
    std::vector<hsize_t> vec_hsize_t(dimensions.begin(), dimensions.end());
    H5::DataSpace dataspace(static_cast<int>(vec_hsize_t.size()), vec_hsize_t.data());
    H5::DataSet dataset = group.createDataSet(datasetName, getH5Type<T>(), dataspace);

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
    H5::Group group = file.openGroup(groupName);
    std::vector<hsize_t> vec_hsize_t(dimensions.begin(), dimensions.end());
    H5::DataSpace dataspace(static_cast<int>(vec_hsize_t.size()), vec_hsize_t.data());
    H5::StrType strtype(H5::PredType::C_S1, maxLength);
    H5::DataSet dataset = group.createDataSet(datasetName, strtype, dataspace);
  }

  template <typename T>
  void writeVector(const std::string& groupName, const std::string& datasetName, const std::vector<T>& data)
  {
    H5::Group group = file.openGroup(groupName);
    H5::DataSet dataset = group.openDataSet(datasetName);
    H5::DataSpace dataspace = dataset.getSpace();
    H5::PredType datatype = getH5Type<T>();
    dataset.write(data.data(), datatype);
  }

  template <>
  void writeVector<bool>(const std::string& groupName, const std::string& datasetName, const std::vector<bool>& data)
  {
    H5::Group group = file.openGroup(groupName);
    H5::DataSet dataset = group.openDataSet(datasetName);
    H5::DataSpace dataspace = dataset.getSpace();
    H5::PredType datatype = getH5Type<bool>();

    std::vector<char> converted(data.size());
    std::transform(data.begin(), data.end(), converted.begin(), [](bool b) { return static_cast<char>(b); });
    dataset.write(converted.data(), datatype);
  }

  template <>
  void writeVector<std::string>(const std::string& groupName, const std::string& datasetName,
                                const std::vector<std::string>& data)
  {
    H5::Group group = file.openGroup(groupName);
    H5::DataSet dataset = group.openDataSet(datasetName);
    H5::DataSpace dataspace = dataset.getSpace();
    H5::StrType strtype = dataset.getStrType();
    size_t maxLength = strtype.getSize();

    std::vector<char> buffer(data.size() * maxLength, '\0');
    for (size_t i = 0; i < data.size(); ++i)
    {
      std::strncpy(&buffer[i * maxLength], data[i].c_str(), maxLength);
    }

    dataset.write(buffer.data(), strtype);
  }

  template <typename T>
  void writeSingleValue(const std::string& groupName, const std::string& datasetName,
                        const std::vector<hsize_t>& indices, T value)
  {
    H5::Group group = file.openGroup(groupName);
    H5::DataSet dataset = group.openDataSet(datasetName);
    H5::DataSpace dataspace = dataset.getSpace();

    H5::DataSpace memspace(1, &indices[indices.size() - 1], NULL);
    dataspace.selectElements(H5S_SELECT_SET, indices.size(), indices.data());
    dataset.write(&value, getH5Type<T>(), memspace, dataspace);
  }

  void writeMetaInfo(const std::string& groupName, const std::string& infoName, const std::string& infoValue)
  {
    // Create a string attribute at the root of the file
    H5::Group group = file.openGroup(groupName);
    H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
    H5::DataSpace attrSpace(H5S_SCALAR);
    H5::Attribute attr = group.createAttribute(infoName, strType, attrSpace);
    attr.write(strType, infoValue);
  }

  void writeMetaInfo(const std::string& groupName, const std::string& infoName, const double& infoValue)
  {
    // Create a string attribute at the root of the file
    H5::Group group = file.openGroup(groupName);
    H5::PredType datatype = getH5Type<double>();
    H5::DataSpace attrSpace(H5S_SCALAR);
    H5::Attribute attr = group.createAttribute(infoName, datatype, attrSpace);
    attr.write(datatype, &infoValue);
  }

  void writeMetaInfo(const std::string& groupName, const std::string& infoName, const size_t& infoValue)
  {
    // Create a string attribute at the root of the file
    H5::Group group = file.openGroup(groupName);
    H5::PredType datatype = getH5Type<size_t>();
    H5::DataSpace attrSpace(H5S_SCALAR);
    H5::Attribute attr = group.createAttribute(infoName, datatype, attrSpace);
    attr.write(datatype, &infoValue);
  }

 private:
  H5::H5File file;

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
      return H5::PredType::NATIVE_UINT;
    }
    else
    {
      throw std::runtime_error("Unsupported data type");
    }
  }
};
