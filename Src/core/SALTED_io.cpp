#include "pch.h"
#include "SALTED_io.h"
#include <filesystem>
#include <iostream>
#include "nos_math.h"
#include "basis_set.h"

namespace {
constexpr int kMaxSaltedBlocks = 10000;
constexpr int kMaxSaltedStringLength = 1 << 20;
constexpr int kMaxSaltedDatasetDims = 8;
constexpr size_t kMaxSaltedDatasetValues = static_cast<size_t>(1) << 31;

template <typename T>
void read_exact(std::ifstream& file, T& value, const std::string& context)
{
    file.read(reinterpret_cast<char*>(&value), sizeof(T));
    err_checkf(static_cast<bool>(file), "Error reading SALTED binary " + context, std::cout);
}

void read_exact_bytes(std::ifstream& file, void* data, const std::streamsize size, const std::string& context)
{
    file.read(reinterpret_cast<char*>(data), size);
    err_checkf(static_cast<bool>(file), "Error reading SALTED binary " + context, std::cout);
}

void skip_exact(std::ifstream& file, const std::streamoff offset, const std::string& context)
{
    file.seekg(offset, std::ios::cur);
    err_checkf(static_cast<bool>(file), "Error seeking SALTED binary " + context, std::cout);
}
}

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


std::string SALTED_BINARY_FILE::read_string_remove_NULL(const int length) {
    // Validate input length
    if (length <= 0 || length > 1000) { // Reasonable upper limit for string length
        std::cerr << "Invalid string length: " << length << std::endl;
        return std::string();
    }
    
    std::vector<char> string_out(length, '\0');
    file.read(string_out.data(), length);
    
    // Check if the read operation succeeded
    if (!file.good()) {
        std::cerr << "Error reading string data from file!" << std::endl;
        return std::string();
    }
    
    //Remove null characters from the string
    string_out.erase(std::remove(string_out.begin(), string_out.end(), '\0'), string_out.end());
    return trim(std::string(string_out.begin(), string_out.end()));
}

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#endif

// A raw positioned read of the model blocks. std::ifstream costs about 2.7x the
// same bytes read this way; see the note in the header for why the file is not
// mapped instead.
void SALTED_BINARY_FILE::open_raw() {
#ifdef _WIN32
    if (raw_handle_) return;
    HANDLE fh = CreateFileW(filepath.wstring().c_str(), GENERIC_READ, FILE_SHARE_READ,
        NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
    if (fh != INVALID_HANDLE_VALUE) raw_handle_ = fh;
#else
    if (raw_fd_ >= 0) return;
    raw_fd_ = ::open(filepath.string().c_str(), O_RDONLY);
#endif
}

void SALTED_BINARY_FILE::close_raw() {
#ifdef _WIN32
    if (raw_handle_) { CloseHandle(static_cast<HANDLE>(raw_handle_)); raw_handle_ = nullptr; }
#else
    if (raw_fd_ >= 0) { ::close(raw_fd_); raw_fd_ = -1; }
#endif
}

// false means "could not", and every caller falls back to the stream.
bool SALTED_BINARY_FILE::read_at(std::streamoff offset, void* dest, std::size_t bytes) {
    char* out = static_cast<char*>(dest);
#ifdef _WIN32
    if (!raw_handle_) return false;
    std::size_t done = 0;
    while (done < bytes)
    {
        const DWORD want = static_cast<DWORD>(std::min<std::size_t>(bytes - done, 1u << 30));
        OVERLAPPED ov{};
        const unsigned long long at = static_cast<unsigned long long>(offset) + done;
        ov.Offset = static_cast<DWORD>(at & 0xFFFFFFFFull);
        ov.OffsetHigh = static_cast<DWORD>(at >> 32);
        DWORD got = 0;
        if (!ReadFile(static_cast<HANDLE>(raw_handle_), out + done, want, &got, &ov) || got == 0)
            return false;
        done += got;
    }
    return true;
#else
    if (raw_fd_ < 0) return false;
    std::size_t done = 0;
    while (done < bytes)
    {
        const ssize_t got = ::pread(raw_fd_, out + done, bytes - done,
            static_cast<off_t>(offset) + static_cast<off_t>(done));
        if (got <= 0) return false;
        done += static_cast<std::size_t>(got);
    }
    return true;
#endif
}

void SALTED_BINARY_FILE::open_file() {
    //Check if file exists and if it is already open
    err_checkf(std::filesystem::exists(filepath), "Couldn't open or find " + filepath.string() + ", leaving", std::cout);
    err_checkf(std::filesystem::is_regular_file(filepath), "Refusing to open non-regular file: " + filepath.string(), std::cout);
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
    if (!file.good()) {
        std::cerr << "Error reading number of blocks from file!" << std::endl;
        return false;
    }
    
    // Validate numBlocks to prevent excessive memory usage or infinite loops
    if (numBlocks < 0 || numBlocks > kMaxSaltedBlocks) { // Reasonable upper limit
        std::cerr << "Invalid number of blocks: " << numBlocks << std::endl;
        return false;
    }
    
    if (debug) std::cout << "Number of blocks: " << numBlocks << std::endl;

    //Now follows (Chunkname (5b str), location (4b int)) * numBlocks
    // Chunknames that are not 5 bytes long are padded with ' '
    for (int i = 0; i < numBlocks; i++) {
        std::string chunkname = read_string_remove_NULL(5);
        if (!file.good()) {
            std::cerr << "Error reading chunk name at block " << i << std::endl;
            return false;
        }
        
        int location;
        file.read((char*)&location, sizeof(int));
        if (!file.good()) {
            std::cerr << "Error reading location at block " << i << std::endl;
            return false;
        }
        
        table_of_contents[chunkname] = location;
        if (debug) std::cout << "Chunk: " << chunkname << " at location: " << location << std::endl;
    }


    header_end = file.tellg();

    return true;
}

void SALTED_BINARY_FILE::populate_config(Config &config) {
    err_checkf(header_end != -1, "Header not read yet! Aborting", std::cout);
    err_checkf(table_of_contents.find("CONFG") != table_of_contents.end(),
        "SALTED binary file does not contain required CONFG block: " + filepath.string(), std::cout);
    file.seekg(table_of_contents["CONFG"], std::ios::beg);

    // Read config data
    auto read_config_block = [this](void* target, const std::streamsize size, const std::string& name) {
        skip_exact(file, 9, "config field header for " + name);
        read_exact_bytes(file, target, size, "config field " + name);
    };
    read_config_block(&config.average, 1, "average"); //Bool data
    read_config_block(&config.field, 1, "field");
    read_config_block(&config.sparsify, 1, "sparsify");
    read_config_block(&config.ncut, 4, "ncut"); //Int data
    read_config_block(&config.nang1, 4, "nang1");
    read_config_block(&config.nang2, 4, "nang2");
    read_config_block(&config.nrad1, 4, "nrad1");
    read_config_block(&config.nrad2, 4, "nrad2");
    read_config_block(&config.Menv, 4, "Menv");
    read_config_block(&config.Ntrain, 4, "Ntrain");
    read_config_block(&config.rcut1, 8, "rcut1");  //Double data
    read_config_block(&config.rcut2, 8, "rcut2");
    read_config_block(&config.sig1, 8, "sig1");
    read_config_block(&config.sig2, 8, "sig2");
    read_config_block(&config.zeta, 8, "zeta");
    read_config_block(&config.trainfrac, 8, "trainfrac");

    int string_size;
    std::vector<char> string_out;
    auto read_config_string = [this, &string_size, &string_out](const std::string& name) {
        skip_exact(file, 9, "config string header for " + name);
        read_exact(file, string_size, "config string length for " + name);
        err_checkf(string_size >= 0 && string_size <= kMaxSaltedStringLength,
            "Invalid SALTED config string length for " + name + ": " + std::to_string(string_size), std::cout);
        string_out.assign(static_cast<size_t>(string_size), '\0');
        if (string_size > 0)
            read_exact_bytes(file, string_out.data(), string_size, "config string " + name);
    };

    read_config_string("species");
    config.species = split_string<std::string>(std::string(string_out.begin(), string_out.end()), " ");

    read_config_string("neighspe1");
    config.neighspe1 = split_string<std::string>(std::string(string_out.begin(), string_out.end()), " ");
    config.nspe1 = (int)config.neighspe1.size();

    read_config_string("neighspe2");
    config.neighspe2 = split_string<std::string>(std::string(string_out.begin(), string_out.end()), " ");
    config.nspe2 = (int)config.neighspe2.size();

    read_config_string("dfbasis");
    config.dfbasis = std::string(string_out.begin(), string_out.end());
}


template <typename T>
void SALTED_BINARY_FILE::read_dataset(std::vector<T>& data, std::vector<size_t>& dims) {
    int ndims;
    read_exact(file, ndims, "dataset dimension count");
    err_checkf(ndims >= 0 && ndims <= kMaxSaltedDatasetDims,
        "Invalid SALTED dataset dimension count: " + std::to_string(ndims), std::cout);
    dims.assign(static_cast<size_t>(ndims), 0);
    for (int i = 0; i < ndims; i++) {
        uint32_t dim = 0;
        read_exact(file, dim, "dataset dimension " + std::to_string(i));
        dims[i] = static_cast<size_t>(dim);
    }
    size_t size = 1;
    for (const size_t dim : dims) {
        err_checkf(dim <= kMaxSaltedDatasetValues && size <= kMaxSaltedDatasetValues / std::max<size_t>(dim, 1),
            "SALTED dataset is unreasonably large or dimension product overflows", std::cout);
        size *= dim;
    }
    err_checkf(size <= kMaxSaltedDatasetValues,
        "SALTED dataset is unreasonably large: " + std::to_string(size) + " values", std::cout);
    data.resize(size);
    if (size > 0)
        read_exact_bytes(file, data.data(), static_cast<std::streamsize>(size * sizeof(T)), "dataset payload");
}

// Same walk as read_dataset, without the payload. Used for species the structure
// does not contain: their shape is still needed to index the flat weight vector,
// their numbers never are.
void SALTED_BINARY_FILE::skip_dataset(std::vector<size_t>& dims, const size_t element_size,
                                      std::streamoff* payload_offset) {
    int ndims;
    read_exact(file, ndims, "dataset dimension count");
    err_checkf(ndims >= 0 && ndims <= kMaxSaltedDatasetDims,
        "Invalid SALTED dataset dimension count: " + std::to_string(ndims), std::cout);
    dims.assign(static_cast<size_t>(ndims), 0);
    for (int i = 0; i < ndims; i++) {
        uint32_t dim = 0;
        read_exact(file, dim, "dataset dimension " + std::to_string(i));
        dims[i] = static_cast<size_t>(dim);
    }
    size_t size = 1;
    for (const size_t dim : dims) {
        err_checkf(dim <= kMaxSaltedDatasetValues && size <= kMaxSaltedDatasetValues / std::max<size_t>(dim, 1),
            "SALTED dataset is unreasonably large or dimension product overflows", std::cout);
        size *= dim;
    }
    if (payload_offset) *payload_offset = static_cast<std::streamoff>(file.tellg());
    file.seekg(static_cast<std::streamoff>(size * element_size), std::ios::cur);
    err_checkf(file.good(), "Could not step over a SALTED dataset payload", std::cout);
}

template <typename T>
T SALTED_BINARY_FILE::read_generic_blocks(const std::string& key, std::function<void(T&, int)> process_block) {
    err_checkf(header_end != -1, "Header not read yet! Aborting", std::cout);
    const auto block = table_of_contents.find(key);
    err_checkf(block != table_of_contents.end(),
        "SALTED binary file does not contain required " + key + " block: " + filepath.string(), std::cout);
    file.seekg(block->second, std::ios::beg);
    err_checkf(static_cast<bool>(file), "Error seeking SALTED binary block " + key, std::cout);
    skip_exact(file, 4, "block datatype for " + key); // Skip the datatype

    int n_blocks;
    read_exact(file, n_blocks, "block count for " + key);
    err_checkf(n_blocks >= 0 && n_blocks <= kMaxSaltedBlocks,
        "Invalid SALTED block count for " + key + ": " + std::to_string(n_blocks), std::cout);
    
    T result;
    if constexpr (std::is_same_v<T, std::shared_ptr<BasisSet>>) {
        result = std::make_shared<BasisSet>();
    }
    else {
        result = T(n_blocks);
    }
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

std::unordered_map<std::string, dMatrix2> SALTED_BINARY_FILE::read_projectors(
    const std::unordered_set<std::string>* wanted,
    std::unordered_map<std::string, std::array<size_t, 2>>* dims_out) {
    return read_lambda_based_data("PROJ", wanted, dims_out);
}

std::unordered_map<std::string, dMatrix2> SALTED_BINARY_FILE::read_features(
    const std::unordered_set<std::string>* wanted) {
    return read_lambda_based_data("FEATS", wanted);
}

// One pass over a block that reads no numbers: for every (species, lambda) it
// records where the matrix starts and how big it is, then steps over it. That
// index is what lets the prediction fetch a lambda when it needs it.
std::unordered_map<std::string, SALTED_BINARY_FILE::block_ref>
SALTED_BINARY_FILE::index_lambda_based_data(const std::string& key) {
    open_raw();   // the blocks this indexes are read through it
    return read_generic_blocks<std::unordered_map<std::string, block_ref>>(key,
        [this, &key](std::unordered_map<std::string, block_ref>& container, int i) {
            std::string element = read_string_remove_NULL(5);
            int nlambda;
            read_exact(file, nlambda, "lambda count for " + element);
            err_checkf(nlambda >= 0 && nlambda <= kMaxSaltedBlocks,
                "Invalid SALTED lambda count for " + element + ": " + std::to_string(nlambda), std::cout);
            for (int lam = 0; lam < nlambda; lam++) {
                std::vector<size_t> dims;
                std::streamoff payload = 0;
                skip_dataset(dims, sizeof(double), &payload);
                err_checkf(dims.size() == 2,
                    "Expected 2D SALTED " + key + " dataset for " + element + " lambda " + std::to_string(lam), std::cout);
                container[element + std::to_string(lam)] = block_ref{ payload, dims[0], dims[1] };
            }
        }
    );
}

dMatrix2 SALTED_BINARY_FILE::load_block(const block_ref& ref) {
    if (ref.rows == 0 || ref.cols == 0) return dMatrix2{};
    // Straight into the matrix. Reading into a vec and reshaping would allocate
    // the whole block twice and memcpy between them, and this runs over the whole
    // 800 MB model once per part - reshape() copies, it does not adopt.
    using ext_t = typename dMatrix2::extents_type;
    dMatrix2 out(ext_t(ref.rows, ref.cols));
    const std::size_t bytes = ref.rows * ref.cols * sizeof(double);
    if (!read_at(ref.offset, out.data(), bytes))
    {
        file.clear();
        file.seekg(ref.offset, std::ios::beg);
        read_exact_bytes(file, out.data(), static_cast<std::streamsize>(bytes),
            "lazily loaded dataset");
    }
    return out;
}

std::unordered_map<std::string, dMatrix2> SALTED_BINARY_FILE::read_lambda_based_data(
    const std::string& key,
    const std::unordered_set<std::string>* wanted,
    std::unordered_map<std::string, std::array<size_t, 2>>* dims_out) {
    return read_generic_blocks<std::unordered_map<std::string, dMatrix2>>(key,
        [this, &key, wanted, dims_out](std::unordered_map<std::string, dMatrix2>& container, int i) {
            std::string element = read_string_remove_NULL(5);
            int nlambda;
            read_exact(file, nlambda, "lambda count for " + element);
            err_checkf(nlambda >= 0 && nlambda <= kMaxSaltedBlocks,
                "Invalid SALTED lambda count for " + element + ": " + std::to_string(nlambda), std::cout);
            const bool load = (wanted == nullptr) || (wanted->find(element) != wanted->end());
            for (int lam = 0; lam < nlambda; lam++) {
                std::vector<size_t> dims;
                if (load) {
                    vec data;
                    read_dataset(data, dims);
                    err_checkf(dims.size() == 2,
                        "Expected 2D SALTED " + key + " dataset for " + element + " lambda " + std::to_string(lam), std::cout);
                    container[element + std::to_string(lam)] = reshape<dMatrix2>(data, Shape2D{ dims[0], dims[1] });
                }
                else {
                    skip_dataset(dims, sizeof(double));
                    err_checkf(dims.size() == 2,
                        "Expected 2D SALTED " + key + " dataset for " + element + " lambda " + std::to_string(lam), std::cout);
                }
                if (dims_out) (*dims_out)[element + std::to_string(lam)] = { dims[0], dims[1] };
            }
        }
    );
}


std::shared_ptr<BasisSet> SALTED_BINARY_FILE::read_basis_set() {
    return read_generic_blocks<std::shared_ptr<BasisSet>>("BASIS",
        [this](std::shared_ptr<BasisSet>& bs, int i) {
            int atomic_nr;
            read_exact(file, atomic_nr, "basis atomic number");
            err_checkf(atomic_nr >= 1 && atomic_nr <= 118,
                "Invalid atomic number in SALTED basis block: " + std::to_string(atomic_nr), std::cout);
            std::cout << "Reading basis set for atomic number: " << atomic_nr << std::endl;
            std::vector<size_t> dims;
            ivec contractions;
            read_dataset(contractions, dims);
            ivec angular_momenta_per_shell;
            read_dataset(angular_momenta_per_shell, dims);
            vec exponents_per_shell;
            read_dataset(exponents_per_shell, dims);
            vec coefficients;
            read_dataset(coefficients, dims);
            size_t primitive_index = 0;
            err_checkf(!dims.empty(), "Invalid empty SALTED basis dimensions", std::cout);
            bs->set_count_for_element(atomic_nr - 1, dims[0]);
            for (int contraction = 0; contraction < contractions.size(); contraction++) {
                err_checkf(contraction < angular_momenta_per_shell.size(),
                    "SALTED basis angular-momentum array shorter than contraction array", std::cout);
                int angular_momentum = angular_momenta_per_shell[contraction];
                for (int func = 0; func < contractions[contraction]; func++, primitive_index++) {
                    err_checkf(primitive_index < angular_momenta_per_shell.size()
                        && primitive_index < exponents_per_shell.size()
                        && primitive_index < coefficients.size(),
                        "SALTED basis primitive arrays have inconsistent sizes", std::cout);
                    bs->add_owned_primitive({ 1, angular_momenta_per_shell[primitive_index], exponents_per_shell[primitive_index], coefficients[primitive_index], contraction });
                }
            }
        }
    );
}
