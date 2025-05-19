#include "pch.h"
#include "SALTED_io.h"
#include <filesystem>
#include <iostream>
#include "nos_math.h"

std::filesystem::path find_first_salted_file(const std::filesystem::path &directory_path)
{
    try
    {
        // Iterate through the directory
        for (const auto &entry : std::filesystem::directory_iterator(directory_path))
        {
            // Check if the entry is a regular file and has a .salted extension
            if (entry.is_regular_file() && entry.path().extension() == ".salted")
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

    return std::filesystem::path(); // Return an empty path if no .salted file is found
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
    this->from_binary = false;
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


std::string SALTED_BINARY_FILE::read_string_remove_NULL(const int lengh) {
    std::vector<char> string_out(lengh, '\0');
    file.read(string_out.data(), lengh);
    //Remove null characters from the string
    string_out.erase(std::remove(string_out.begin(), string_out.end(), '\0'), string_out.end());
    return std::string(string_out.begin(), string_out.end());
}

void SALTED_BINARY_FILE::open_file() {
    //Check if file exists and if it is already open
    err_checkf(std::filesystem::exists(filepath), "Couldn't open or find " + filepath.string() + ", leaving", std::cout);
    err_checkf(!file.is_open(), "File already open, leaving", std::cout);
    file.open(filepath, std::ios::in | std::ios::binary);
    err_checkf(file.is_open(), "Couldn't open file: " + filepath.string(), std::cout);
}

bool SALTED_BINARY_FILE::read_header() {
    // Read and verify magic number
    file.seekg(0, std::ios::beg);
    char magic[HEADER_SIZE];
    file.read(magic, HEADER_SIZE);
    if (std::string(magic, HEADER_SIZE) != MAGIC_NUMBER) {
        std::cerr << "Invalid file format!" << std::endl;
        return false;
    }
    // Read version
    file.read((char*)&version, sizeof(int));
    if (debug) std::cout << "File Version: " << version << std::endl;

    //Read number of blocks
    file.read((char*)&numBlocks, sizeof(int));
    if (debug) std::cout << "Number of blocks: " << numBlocks << std::endl;

    //Now follows (Chunkname (5b str), location (4b int)) * numBlocks
    // Chunknames that are not 5 bytes long are padded with ' '
    for (int i = 0; i < numBlocks; i++) {
        std::string chunkname = read_string_remove_NULL(5);
        int location;
        file.read((char*)&location, sizeof(int));
        table_of_contents[chunkname] = location;
        if (debug) std::cout << "Chunk: " << chunkname << " at location: " << location << std::endl;
    }


    header_end = file.tellg();

    return true;
}

//Small macro to read a block of data
//Skips the first 9 bytes of the block with info about the key (5 bytes str) +  datatype (4 bytes int32) 
//Reads the data into the container
//Increments the block_id
#define READ_BLOCK(container, size) file.seekg(9, std::ios::cur); file.read((char*)container, size);
#define READ_BLOCK_STRING(container) file.seekg(9, std::ios::cur); file.read((char*)&string_size, 4); container.resize(string_size, '\0');  file.read((char*)container.data(), string_size);

void SALTED_BINARY_FILE::populate_config(Config &config) {
    err_checkf(header_end != -1, "Header not read yet! Aborting", std::cout);
    file.seekg(table_of_contents["CONFG"], std::ios::beg);

    // Read config data
    READ_BLOCK(&config.average, 1); //Bool data
    READ_BLOCK(&config.field, 1); 
    READ_BLOCK(&config.sparsify, 1); 
    READ_BLOCK(&config.ncut, 4); //Int data
    READ_BLOCK(&config.nang1, 4); 
    READ_BLOCK(&config.nang2, 4);
    READ_BLOCK(&config.nrad1, 4); 
    READ_BLOCK(&config.nrad2, 4); 
    READ_BLOCK(&config.Menv, 4); 
    READ_BLOCK(&config.Ntrain, 4);
    READ_BLOCK(&config.rcut1, 8);  //Double data
    READ_BLOCK(&config.rcut2, 8);
    READ_BLOCK(&config.sig1, 8);
    READ_BLOCK(&config.sig2, 8);
    READ_BLOCK(&config.zeta, 8);
    READ_BLOCK(&config.trainfrac, 8);

    int string_size;
    std::vector<char> string_out;
    READ_BLOCK_STRING(string_out);
    config.species = split_string<std::string>(std::string(string_out.begin(), string_out.end()), " ");

    READ_BLOCK_STRING(string_out);
    config.neighspe1 = split_string<std::string>(std::string(string_out.begin(), string_out.end()), " ");
    config.nspe1 = (int)config.neighspe1.size();

    READ_BLOCK_STRING(string_out);
    config.neighspe2 = split_string<std::string>(std::string(string_out.begin(), string_out.end()), " ");
    config.nspe2 = (int)config.neighspe2.size();

    READ_BLOCK_STRING(string_out);
    config.dfbasis = std::string(string_out.begin(), string_out.end());
}


template <typename T>
void SALTED_BINARY_FILE::read_dataset(std::vector<T>& data, std::vector<size_t>& dims) {
    int ndims;
    file.read((char*)&ndims, 4);
    dims.resize(ndims, 0);
    for (int i = 0; i < ndims; i++) {
        file.read((char*)&dims[i], 4);
    }
    int size = static_cast<int>(std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>()));
    data.resize(size);
    file.read((char*)data.data(), size * sizeof(T));
}

template <typename T>
T SALTED_BINARY_FILE::read_generic_blocks(const std::string& key, std::function<void(T&, int)> process_block) {
    err_checkf(header_end != -1, "Header not read yet! Aborting", std::cout);
    file.seekg(table_of_contents[key], std::ios::beg);
    file.seekg(4, std::ios::cur); // Skip the datatype

    int n_blocks;
    file.read(reinterpret_cast<char*>(&n_blocks), sizeof(n_blocks));

    T result(n_blocks);
    for (int i = 0; i < n_blocks; i++) {
        process_block(result, i);
    }
    return result;
}

std::unordered_map<int, std::vector<int64_t>> SALTED_BINARY_FILE::read_fps() {
    return read_generic_blocks<std::unordered_map<int, std::vector<int64_t>>>("FPS",
        [this](std::unordered_map<int, std::vector<int64_t>>& fps, int i) {
            std::vector<size_t> dims;
            std::vector<int64_t> data;
            read_dataset(data, dims);
            fps[i] = data;
        }
    );
}

std::unordered_map<std::string, vec> SALTED_BINARY_FILE::read_averages() {
    return read_generic_blocks<std::unordered_map<std::string, vec>>("AVERG",
        [this](std::unordered_map<std::string, vec>& averages, int i) {
            std::string element = read_string_remove_NULL(5);
            std::vector<size_t> dims;
            vec data;
            read_dataset(data, dims);
            averages[element] = data;
        }
    );
}

std::unordered_map<int, vec> SALTED_BINARY_FILE::read_wigners() {
    return read_generic_blocks<std::unordered_map<int, vec>>("WIG",
        [this](std::unordered_map<int, vec>& wigners, int i) {
            std::vector<size_t> dims;
            vec data;
            read_dataset(data, dims);
            wigners[i] = data;
        }
    );
}

vec SALTED_BINARY_FILE::read_weights() {
    vec weights;
    read_generic_blocks<std::vector<vec>>("WEIGH",
        [this, &weights](std::vector<vec>&, int i) {
            std::vector<size_t> dims;
            vec data;
            read_dataset(data, dims);
            weights = data;
        }
    );
    return weights;
}

std::unordered_map<std::string, dMatrix2> SALTED_BINARY_FILE::read_projectors() {
    return read_lambda_based_data("PROJ");
}

std::unordered_map<std::string, dMatrix2> SALTED_BINARY_FILE::read_features() {
    return read_lambda_based_data("FEATS");
}

std::unordered_map<std::string, dMatrix2> SALTED_BINARY_FILE::read_lambda_based_data(const std::string& key) {
    return read_generic_blocks<std::unordered_map<std::string, dMatrix2>>(key,
        [this](std::unordered_map<std::string, dMatrix2>& container, int i) {
            std::string element = read_string_remove_NULL(5);
            int nlambda;
            file.read(reinterpret_cast<char*>(&nlambda), sizeof(nlambda));
            for (int lam = 0; lam < nlambda; lam++) {
                std::vector<size_t> dims;
                vec data;
                read_dataset(data, dims);
                //container.emplace(std::piecewise_construct, std::forward_as_tuple(element + std::to_string(lam)), std::forward_as_tuple(dims[0], dims[1]));
                //dMatrix2 A = reshape<dMatrix2>(data, Shape2D{ dims[0], dims[1] });
                //std::copy(data.begin(), data.end(), container[element + std::to_string(lam)].data());
                container[element + std::to_string(lam)] = reshape<dMatrix2>(data, Shape2D{ dims[0], dims[1] });
            }
        }
    );
}
