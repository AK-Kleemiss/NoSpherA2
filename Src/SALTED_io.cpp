#include "SALTED_io.h"
#include <filesystem>

template <typename T>
void readHDF5Data(H5::DataSet &dataset, std::vector<T> &data)
{
    // Select correct datatype depending on T
    H5::PredType type = H5::PredType::PREDTYPE_CONST;
    if constexpr (std::is_same_v<T, double>)
    {
        type = H5::PredType::NATIVE_DOUBLE;
    }
    else if constexpr (std::is_same_v<T, float>)
    {
        type = H5::PredType::NATIVE_FLOAT;
    }
    else if constexpr (std::is_same_v<T, int>)
    {
        type = H5::PredType::NATIVE_INT;
    }
    else if constexpr (std::is_same_v<T, int64_t>)
    {
        type = H5::PredType::NATIVE_INT64;
    }
    else
    {
        throw std::runtime_error("Unsupported datatype");
    }
    dataset.read(data.data(), type);
}
template void readHDF5Data(H5::DataSet &dataset, std::vector<double> &data);
template void readHDF5Data(H5::DataSet &dataset, std::vector<float> &data);
template void readHDF5Data(H5::DataSet &dataset, std::vector<int> &data);
template void readHDF5Data(H5::DataSet &dataset, std::vector<int64_t> &data);

template <typename T>
std::vector<T> readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t> &dims_out)
{

    try
    {
        H5::Exception::dontPrint();
        H5::DataSet dataset = file.openDataSet(dataset_name);
        H5::DataSpace dataspace = dataset.getSpace();

        const int rank = dataspace.getSimpleExtentNdims();
        dims_out.resize(rank);
        dataspace.getSimpleExtentDims(dims_out.data(), NULL);

        size_t totalSize = 1;
        for (const auto &dim : dims_out)
        {
            totalSize *= dim;
        }

        std::vector<T> flatData(totalSize);
        readHDF5Data(dataset, flatData);
        dataset.close();

        return flatData;
    }
    // catch failure caused by the H5File operations
    catch (H5::FileIException error)
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSet operations
    catch (H5::DataSetIException error)
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSpace operations
    catch (H5::DataSpaceIException error)
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSpace operations
    catch (H5::DataTypeIException error)
    {
        error.printErrorStack();
    }
    catch (H5::Exception error)
    {
        error.printErrorStack();
    }
    return {};
}
template std::vector<double> readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t> &dims_out);
template std::vector<float> readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t> &dims_out);
template std::vector<int> readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t> &dims_out);
template std::vector<int64_t> readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t> &dims_out);
template std::vector<unsigned long long> readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t> &dims_out);

template <typename T>
Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent>> readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t>& dims_out)
{

    try
    {
        H5::Exception::dontPrint();
        H5::DataSet dataset = file.openDataSet(dataset_name);
        H5::DataSpace dataspace = dataset.getSpace();

        const int rank = dataspace.getSimpleExtentNdims();
        dims_out.resize(rank);
        dataspace.getSimpleExtentDims(dims_out.data(), NULL);

        size_t totalSize = 1;
        for (const auto& dim : dims_out)
        {
            totalSize *= dim;
        }

        std::vector<T> flatData(totalSize);
        readHDF5Data(dataset, flatData);
        dataset.close();

        return Kokkos::mdspan<T, Kokkos::extents<unsigned long long, std::dynamic_extent>>(flatData.data(), (unsigned long long) totalSize);
    }
    // catch failure caused by the H5File operations
    catch (H5::FileIException error)
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSet operations
    catch (H5::DataSetIException error)
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSpace operations
    catch (H5::DataSpaceIException error)
    {
        error.printErrorStack();
    }
    // catch failure caused by the DataSpace operations
    catch (H5::DataTypeIException error)
    {
        error.printErrorStack();
    }
    catch (H5::Exception error)
    {
        error.printErrorStack();
    }
    return {};
}
template dMatrix1 readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t>& dims_out);
template iMatrix1 readHDF5(H5::H5File file, std::string dataset_name, std::vector<hsize_t>& dims_out);


std::filesystem::path find_first_h5_file(const std::filesystem::path &directory_path)
{
    try
    {
        // Iterate through the directory
        for (const auto &entry : std::filesystem::directory_iterator(directory_path))
        {
            // Check if the entry is a regular file and has a .h5 extension
            if (entry.is_regular_file() && entry.path().extension() == ".h5")
            {
                return entry.path().filename().string(); // Return the filename
            }
        }
    }
    catch (const std::filesystem::filesystem_error &e)
    {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cerr << "General error: " << e.what() << std::endl;
    }

    return std::filesystem::path(); // Return an empty path if no .h5 file is found
}

template <typename Scalar>
void read_npy(std::filesystem::path &filename, std::vector<Scalar> &data)
{
    std::vector<unsigned long> shape{};
    bool fortran_order;
    npy::LoadArrayFromNumpy(filename, shape, fortran_order, data);
}
template void read_npy(std::filesystem::path &filename, std::vector<double> &data);
template void read_npy(std::filesystem::path &filename, std::vector<float> &data);

template <typename T>
std::unordered_map<int, std::vector<T>> read_fps(std::filesystem::path &filename, int lmax_max)
{
    std::unordered_map<int, std::vector<T>> vfps{};
    std::filesystem::path p;
    for (int lam = 0; lam < lmax_max + 1; lam++)
    {
        p = (filename.string() + std::to_string(lam) + ".npy");
        read_npy(p, vfps[lam]);
    }
    return vfps;
}
template std::unordered_map<int, std::vector<double>> read_fps(std::filesystem::path &filename, int lmax_max);
template std::unordered_map<int, std::vector<int64_t>> read_fps(std::filesystem::path &filename, int lmax_max);

template <class T>
std::vector<T> readVectorFromFile(const std::filesystem::path &filename)
{
    std::vector<T> result;
    std::ifstream file(filename);
    std::string line;

    err_checkf(file.is_open(), "Could not open file: " + filename.string(), std::cout);

    while (std::getline(file, line))
    {
        try
        {
            T value;
            if constexpr (std::is_same_v<T, double>)
            {
                value = std::stod(line);
            }
            else if constexpr (std::is_same_v<T, int>)
            {
                value = std::stoi(line);
            }
            result.push_back(value);
        }
        catch (const std::invalid_argument &e)
        {
            // Handle the case where the line cannot be converted to double
            // For now, we'll just skip the line
            std::string message = "Could not convert line to double: " + line + ": " + e.what();
            continue;
        }
    }

    return result;
}
template std::vector<double> readVectorFromFile(const std::filesystem::path &filename);
template std::vector<int> readVectorFromFile(const std::filesystem::path &filename);

// ----------------- Functions to populate the Config struct -----------------

// From txt file
void Config::populateFromFile(const std::filesystem::path &filename)
{
    err_checkf(std::filesystem::exists(filename), "Couldn't open or find " + filename.string() + ", leaving", std::cout);
    std::ifstream file(filename);
    std::string line;
    while (std::getline(file, line))
    {
        std::istringstream iss(line);
        std::string key;
        if (std::getline(iss, key, '='))
        {
            std::string value;
            if (std::getline(iss, value))
            {
                value = trim(value);
                key = trim(key);

                if (key == "average")
                    this->average = (value == "True");
                else if (key == "field")
                    this->field = (value == "True");
                else if (key == "sparsify")
                    this->sparsify = (value == "True");
                else if (key == "ncut")
                    this->ncut = std::stoi(value);
                else if (key == "species")
                    this->species = parseVector(value);
                else if (key == "rcut1")
                    this->rcut1 = std::stod(value);
                else if (key == "rcut2")
                    this->rcut2 = std::stod(value);
                else if (key == "nang1")
                    this->nang1 = std::stoi(value);
                else if (key == "nang2")
                    this->nang2 = std::stoi(value);
                else if (key == "nrad1")
                    this->nrad1 = std::stoi(value);
                else if (key == "nrad2")
                    this->nrad2 = std::stoi(value);
                else if (key == "sig1")
                    this->sig1 = std::stod(value);
                else if (key == "sig2")
                    this->sig2 = std::stod(value);
                else if (key == "neighspe1")
                    this->neighspe1 = parseVector(value);
                else if (key == "neighspe2")
                    this->neighspe2 = parseVector(value);
                else if (key == "zeta")
                    this->zeta = std::stod(value);
                else if (key == "Menv")
                    this->Menv = std::stoi(value);
                else if (key == "Ntrain")
                    this->Ntrain = std::stoi(value);
                else if (key == "trainfrac")
                    this->trainfrac = std::stod(value);
                else if (key == "dfbasis")
                    this->dfbasis = value;
            }
        }
    }
    this->nspe1 = static_cast<int>(this->neighspe1.size());
    this->nspe2 = static_cast<int>(this->neighspe2.size());
    this->from_h5 = false;
}
std::vector<std::string> Config::parseVector(const std::string &value)
{
    std::vector<std::string> result;

    // Find the start and end of the list in the string
    size_t start = value.find("[");
    size_t end = value.find("]");
    if (start == std::string::npos || end == std::string::npos || start >= end)
    {
        return result; // Return empty vector if no valid list is found
    }

    // Extract the list content
    std::string listContent = value.substr(start + 1, end - start - 1);

    // Remove spaces, quotes, and other unnecessary characters
    listContent.erase(std::remove_if(listContent.begin(), listContent.end(), [](char c)
                                     { return std::isspace(c) || c == '"' || c == '\''; }),
                      listContent.end());

    // Split the cleaned string into items
    std::istringstream stream(listContent);
    std::string item;
    while (std::getline(stream, item, ','))
    {
        result.push_back(item);
    }

    return result;
}

// Included in h5 file
void Config::populateFromFile(const H5::H5File file)
{
    err_checkf(!file.attrExists("input"), "Group 'input' not found in h5-File!", std::cout);
    H5::Group input_grp = file.openGroup("input");

    std::vector<std::string> *dataset_names = new std::vector<std::string>();
    // Get all datasets in the group
    herr_t status;
    status = H5Literate(input_grp.getId(), H5_INDEX_NAME, H5_ITER_NATIVE, NULL, &Config::op_func, dataset_names);
    err_checkf((status == 0), "Error reading h5 input file", std::cout);

    for (int i = 0; i < dataset_names->size(); i++)
    {
        std::string dataset_name = dataset_names->at(i);
        H5::DataSet dataSet = input_grp.openDataSet(dataset_name);
        // read type of dataset
        H5T_class_t type_class = dataSet.getTypeClass();
        // call the correct handler
        if (this->handlers.find(type_class) != this->handlers.end())
        {
            this->handlers[type_class](dataset_name, dataSet);
        }
    }
    this->nspe1 = static_cast<int>(this->neighspe1.size());
    this->nspe2 = static_cast<int>(this->neighspe2.size());
    this->from_h5 = true;
}

void Config::populate_config(const std::filesystem::path &dataset_name, const int &data)
{
    if (dataset_name == "average" || dataset_name == "averages")
        this->average = data != 0; // Assuming boolean stored as integer
    else if (dataset_name == "field")
        this->field = data != 0; // Assuming boolean stored as integer
    else if (dataset_name == "sparsify")
        this->sparsify = data != 0; // Assuming boolean stored as integer
    else if (dataset_name == "ncut")
        this->ncut = data;
    else if (dataset_name == "nang1")
        this->nang1 = data;
    else if (dataset_name == "nang2")
        this->nang2 = data;
    else if (dataset_name == "nrad1")
        this->nrad1 = data;
    else if (dataset_name == "nrad2")
        this->nrad2 = data;
    else if (dataset_name == "Menv")
        this->Menv = data;
    else if (dataset_name == "Ntrain")
        this->Ntrain = data;
    else
        std::cout << "Unknown dataset name: " << dataset_name << std::endl;
}
void Config::populate_config(const std::filesystem::path &dataset_name, const float &data)
{
    if (dataset_name == "rcut1")
        this->rcut1 = data;
    else if (dataset_name == "rcut2")
        this->rcut2 = data;
    else if (dataset_name == "sig1")
        this->sig1 = data;
    else if (dataset_name == "sig2")
        this->sig2 = data;
    else if (dataset_name == "zeta")
        this->zeta = data;
    else if (dataset_name == "trainfrac")
        this->trainfrac = data;
    else
        std::cout << "Unknown dataset name: " << dataset_name << std::endl;
}
void Config::populate_config(const std::filesystem::path &dataset_name, const double &data)
{
    if (dataset_name == "rcut1")
        this->rcut1 = data;
    else if (dataset_name == "rcut2")
        this->rcut2 = data;
    else if (dataset_name == "sig1")
        this->sig1 = data;
    else if (dataset_name == "sig2")
        this->sig2 = data;
    else if (dataset_name == "zeta")
        this->zeta = data;
    else if (dataset_name == "trainfrac")
        this->trainfrac = data;
    else
        std::cout << "Unknown dataset name: " << dataset_name << std::endl;
}
void Config::populate_config(const std::filesystem::path &dataset_name, const std::string &data)
{
    if (dataset_name == "dfbasis")
        this->dfbasis = data;
    else
        std::cout << "Unknown dataset name: " << dataset_name << std::endl;
}
void Config::populate_config(const std::filesystem::path &dataset_name, const std::vector<std::string> &data)
{
    if (dataset_name == "species")
        this->species = data;
    else if (dataset_name == "neighspe1")
        this->neighspe1 = data;
    else if (dataset_name == "neighspe2")
        this->neighspe2 = data;
    else
        std::cout << "Unknown dataset name: " << dataset_name << std::endl;
}

void Config::handle_int_dataset(const std::filesystem::path &dataset_name, H5::DataSet &dataSet)
{
    int data{};
    dataSet.read(&data, H5::PredType::NATIVE_INT);
    populate_config(dataset_name, data);
}
void Config::handle_float_dataset(const std::filesystem::path &dataset_name, H5::DataSet &dataSet)
{
    size_t byteSize = dataSet.getFloatType().getSize();
    if (byteSize == 4)
    {
        float data{};
        dataSet.read(&data, H5::PredType::NATIVE_FLOAT);
        populate_config(dataset_name, data);
    }
    else if (byteSize == 8)
    {
        double data{};
        dataSet.read(&data, H5::PredType::NATIVE_DOUBLE);
        populate_config(dataset_name, data);
    }
}
void Config::handle_string_dataset(const std::filesystem::path &dataset_name, H5::DataSet &dataSet)
{
    H5::DataSpace dataspace = dataSet.getSpace();
    hsize_t dims_out[2];
    int ndims = dataspace.getSimpleExtentDims(dims_out, NULL);
    if (ndims == 0)
    {
        std::string data{};
        dataSet.read(data, dataSet.getStrType(), H5S_SCALAR);
        populate_config(dataset_name, data);
    }
    else
    {
        H5::StrType str_type(H5::PredType::C_S1, H5T_VARIABLE);
        std::vector<char *> data{};
        data.resize(dims_out[0]);
        dataSet.read(data.data(), dataSet.getStrType());
        std::vector<std::string> data_str;
        for (int i = 0; i < data.size(); i++)
        {
            data_str.push_back(data[i]);
        }
        populate_config(dataset_name, data_str);
    }
}

herr_t Config::op_func(hid_t loc_id, const char *name, const H5L_info_t *info, void *operator_data)
{
    // Get the type of the object and read the name into the vector
    herr_t status;
    H5O_info_t infobuf;
#if H5_VERSION_GE(1, 12, 0) && !defined(H5_USE_110_API) && !defined(H5_USE_18_API) && !defined(H5_USE_16_API)
    status = H5Oget_info_by_name(loc_id, name, &infobuf, H5O_INFO_ALL, H5P_DEFAULT);
#else
    status = H5Oget_info_by_name(loc_id, name, &infobuf, H5P_DEFAULT);
#endif
    std::vector<std::string> *data = (std::vector<std::string> *)operator_data;
    data->push_back(name);
    return 0;
    (void)info;
}