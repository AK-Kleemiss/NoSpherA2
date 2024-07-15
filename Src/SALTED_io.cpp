#include "SALTED_io.h"

using namespace std;

#if has_RAS
vec2 readHDF5(H5::H5File file, string dataset_name)
{
    /*
     * Try block to detect exceptions raised by any of the calls inside it
     */
    try
    {
        /*
         * Turn off the auto-printing when failure occurs so that we can
         * handle the errors appropriately
         */
        H5::Exception::dontPrint();

        /*
         * Open the specified dataset in the file.
         */
        H5::DataSet dataset = file.openDataSet(dataset_name);

        /*
         * Get dataspace of the dataset.
         */
        H5::DataSpace dataspace = dataset.getSpace();

        /*
         * Get the number of dimensions in the dataspace.
         */
        const int rank = dataspace.getSimpleExtentNdims();

        /*
         * Get the dimension size of each dimension in the dataspace and
         * display them.
         */
        vector<hsize_t> dims_out(rank);
        (void)dataspace.getSimpleExtentDims(dims_out.data(), NULL);

        vector<double> data(dims_out[0] * dims_out[1]);
        /*
         * Read data from hyperslab in the file into the hyperslab in
         * memory and display the data.
         */
        dataset.read(data.data(), H5::PredType::NATIVE_DOUBLE);

        // Reorder the results to match the original matrix
        vec2 retData(dims_out[0], vector<double>(dims_out[1]));
        for (int i = 0; i < dims_out[0]; i++)
        {
            for (int j = 0; j < dims_out[1]; j++)
            {
                retData[i][j] = data[i * dims_out[1] + j];
            }
        }
        /*
         * Close the dataset and file
         */
        dataset.close();

        // successfully terminated
        return retData;
    } // end of try block

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
#endif

template <typename Scalar>
void read_npy(string filename, vector<Scalar> &data)
{
    vector<unsigned long> shape{};
    bool fortran_order;
    npy::LoadArrayFromNumpy(filename, shape, fortran_order, data);
}
template void read_npy(string filename, vector<double> &data);
template void read_npy(string filename, vector<float> &data);

template <typename T>
unordered_map<int, vector<T>> read_fps(string filename, int lmax_max)
{
    unordered_map<int, vector<T>> vfps{};
    vector<T> data{};
    for (int lam = 0; lam < lmax_max + 1; lam++)
    {
        data.clear();
        read_npy(filename + to_string(lam) + ".npy", data);
        vfps[lam] = data;
    }
    return vfps;
}
template unordered_map<int, vector<double>> read_fps(string filename, int lmax_max);
template unordered_map<int, vector<int64_t>> read_fps(string filename, int lmax_max);

template <class T>
vector<T> readVectorFromFile(const string &filename)
{
    std::vector<T> result;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open())
    {
        throw std::runtime_error("Could not open file: " + filename);
    }

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
            string message = "Could not convert line to double: " + line + ": " + e.what();
            continue;
        }
    }

    return result;
}
template vector<double> readVectorFromFile(const string &filename);
template vector<int> readVectorFromFile(const string &filename);

void Config::populateFromFile(const std::string &filename)
{
    err_checkf(exists(filename), "Couldn't open or find " + filename + ", leaving", std::cout);
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
                    average = (value == "True");
                else if (key == "field")
                    field = (value == "True");
                else if (key == "sparsify")
                    sparsify = (value == "True");
                else if (key == "ncut")
                    ncut = std::stoi(value);
                else if (key == "species")
                    species = parseVector(value);
                else if (key == "rcut1")
                    rcut1 = std::stod(value);
                else if (key == "rcut2")
                    rcut2 = std::stod(value);
                else if (key == "nang1")
                    nang1 = std::stoi(value);
                else if (key == "nang2")
                    nang2 = std::stoi(value);
                else if (key == "nrad1")
                    nrad1 = std::stoi(value);
                else if (key == "nrad2")
                    nrad2 = std::stoi(value);
                else if (key == "sig1")
                    sig1 = std::stod(value);
                else if (key == "sig2")
                    sig2 = std::stod(value);
                else if (key == "neighspe1")
                    neighspe1 = parseVector(value);
                else if (key == "neighspe2")
                    neighspe2 = parseVector(value);
                else if (key == "zeta")
                    zeta = std::stod(value);
                else if (key == "Menv")
                    Menv = std::stoi(value);
                else if (key == "Ntrain")
                    Ntrain = std::stoi(value);
                else if (key == "trainfrac")
                    trainfrac = std::stof(value);
            }
        }
    }
    this->nspe1 = static_cast<int>(this->neighspe1.size());
    this->nspe2 = static_cast<int>(this->neighspe2.size());
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