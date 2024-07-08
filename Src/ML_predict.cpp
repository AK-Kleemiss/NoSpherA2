#include "ML_predict.h"
using namespace std;

// Define types for simplicity
using cvec4 = vector<vector<cvec2>>;

#if defined(_WIN32) || defined(__RASCALINE__)
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
void read_npy(string filename, std::vector<Scalar> &data)
{
    vector<unsigned long> shape{};
    bool fortran_order;
    npy::LoadArrayFromNumpy(filename, shape, fortran_order, data);
}

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

template <class T>
std::vector<T> readVectorFromFile(const std::string &filename)
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
            double value = std::stod(line);
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

// Reshape
struct Shape2D
{
    int rows;
    int cols;
};

struct Shape3D
{
    int depth;
    int rows;
    int cols;
};
// To_2D
template <typename T>
vector<vector<T>> reshape(vector<T> flatVec, Shape2D sizes)
{
    vector<vector<T>> reshapedVec(sizes.rows, vector<T>(sizes.cols));
    for (int i = 0; i < sizes.rows; ++i)
    {
        for (int j = 0; j < sizes.cols; ++j)
        {
            reshapedVec[i][j] = flatVec[i * sizes.cols + j];
        }
    }
    return reshapedVec;
}

// To_3D
template <typename T>
vector<vector<vector<T>>> reshape(vector<T> flatVec, Shape3D sizes)
{
    vector<vector<vector<T>>> reshapedVec(sizes.depth, vector<vector<T>>(sizes.rows, vector<T>(sizes.cols)));
    for (int i = 0; i < sizes.depth; ++i)
    {
        for (int j = 0; j < sizes.rows; ++j)
        {
            for (int k = 0; k < sizes.cols; ++k)
            {
                reshapedVec[i][j][k] = flatVec[i * sizes.rows * sizes.cols + j * sizes.cols + k];
            }
        }
    }
    return reshapedVec;
}

// Flatten Vectors 2D
template <typename T>
vector<T> flatten(const vector<vector<T>> &vec2D)
{
    vector<T> flatVec;
    size_t totalSize = 0;
    for (const auto& row : vec2D)
    {
        totalSize += row.size();
    }
    flatVec.reserve(totalSize);
    for (const vector<T> &row : vec2D)
    {
        // Use std::copy to copy the entire row at once
        flatVec.insert(flatVec.end(), row.begin(), row.end());
    }
    return flatVec;
}

// Flatten Vectors 3D
template <typename T>
vector<T> flatten(const vector<vector<vector<T>>> &vec3D)
{
    vector<T> flatVec;
    size_t totalSize = 0;
    for (const auto& row : vec3D)
    {
      for (const auto& innerRow : row)
      {
        totalSize += innerRow.size();
      }
    }
    flatVec.reserve(totalSize);
    for (const vector<vector<T>> &row : vec3D)
    {
        for (const vector<T> &innerRow : row)
        {
          // Use std::copy to copy the entire innerRow at once
          flatVec.insert(flatVec.end(), innerRow.begin(), innerRow.end());
        }
    }
    return flatVec;
}

// SLICE Operation
std::vector<double> slice(const std::vector<double> &vec, size_t start, size_t length)
{
    if (start + length > vec.size())
    {
        throw std::out_of_range("Slice range is out of bounds.");
    }

    std::vector<double> result;
    for (size_t i = start; i < start + length; ++i)
    {
        result.push_back(vec[i]);
    }

    return result;
}

// Matrix multiplication
// 2D x 2D MATRIX MULTIPLICATION
template <typename T>
vector<vector<T>> dot(const vector<vector<T>> &mat1, const vector<vector<T>> &mat2)
{
    // if either of the matrices is empty, return a empty matrix
    if (mat1.empty() || mat2.empty())
    {
        return {};
    }

    size_t rows1 = mat1.size();
    size_t cols1 = mat1[0].size();
    size_t rows2 = mat2.size();
    size_t cols2 = mat2[0].size();

    // Check if matrix multiplication is possible
    if (cols1 != rows2)
    {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    }

    vector<vector<T>> result(rows1, vector<T>(cols2, 0.0));
    const long long int totalIterations = static_cast<long long int>(rows1 * cols2 * cols1);
    size_t total_size = rows1 * cols2;
#pragma omp parallel
    {
        vector<T> local_result(total_size, 0.0);
        int i, j, k, flatIndex;

#pragma omp for schedule(static) private(i, j, k, flatIndex) nowait
        for (long long int n = 0; n < totalIterations; ++n)
        {
            i = static_cast<int>(n / (cols2 * cols1));
            j = static_cast<int>((n / cols1) % cols2);
            k = static_cast<int>(n % cols1);
            flatIndex = i * cols2 + j;
            local_result[flatIndex] += mat1[i][k] * mat2[k][j];
        }

#pragma omp critical
        {
            for (i = 0; i < rows1; ++i)
            {
                for (j = 0; j < cols2; ++j)
                {
                    flatIndex = i * cols2 + j;
                    result[i][j] += local_result[flatIndex];
                }
            }
        }
    }

    return result;
}

// 2D x 1D MATRIX MULTIPLICATION
template <typename T>
vector<T> dot(const vector<vector<T>> &mat, const vector<T> &vec)
{
    int mat_rows = static_cast<int>(mat.size());
    int mat_cols = static_cast<int>(mat[0].size());
    int vec_size = static_cast<int>(vec.size());

    // Check if matrix multiplication is possible
    if (mat_cols != vec_size)
    {
        throw std::invalid_argument("Matrix dimensions do not match for multiplication");
    }

    vector<T> result(mat_rows, 0.0);
    int totalIterations = mat_rows * mat_cols;
#pragma omp parallel
    {
        vector<T> local_result(mat_rows, 0.0);
#pragma omp for schedule(static)
        for (int n = 0; n < totalIterations; ++n)
        {
            int i = n / mat_cols;
            int j = n % mat_cols;
            local_result[i] += mat[i][j] * vec[j];
        }
#pragma omp critical
        {
            for (int i = 0; i < mat_rows; ++i)
            {
                result[i] += local_result[i];
            }
        }
    }

    return result;
}

// TRANSPOSES
// 3D MATRIX
template <typename T>
vector<vector<vector<T>>> transpose(const vector<vector<vector<T>>> &originalVec)
{
    if (originalVec.empty() || originalVec[0].empty() || originalVec[0][0].empty())
    {
        return {}; // Return an empty vector if the original vector is empty or not 3D
    }

    int newDim1 = originalVec[0][0].size(); // New first dimension is the old third dimension
    int newDim2 = originalVec.size();       // New second dimension is the old first dimension
    int newDim3 = originalVec[0].size();    // New third dimension is the old second dimension

    vector<vector<vector<T>>> transposedVec(newDim1, vector<vector<T>>(newDim2, std::vector<T>(newDim3)));

    for (int i = 0; i < originalVec.size(); ++i)
    {
        for (int j = 0; j < originalVec[i].size(); ++j)
        {
            for (int k = 0; k < originalVec[i][j].size(); ++k)
            {
                transposedVec[k][i][j] = originalVec[i][j][k];
            }
        }
    }

    return transposedVec;
}

// 2D MATRIX
template <class T>
vector<vector<T>> transpose(const vector<vector<T>> &mat)
{
    int rows = static_cast<int>(mat.size());
    int cols = static_cast<int>(mat[0].size());
    vector<vector<T>> result(cols, vector<T>(rows));

    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            result[j][i] = mat[i][j];
        }
    }
    return result;
}

// Reorder 3D Vectors following a given order
template <typename T>
vector<vector<vector<T>>> reorder3D(const vector<vector<vector<T>>> &original)
{
    if (original.empty() || original[0].empty() || original[0][0].empty())
    {
        return {}; // Return an empty vector if the original vector is empty or not properly sized
    }

    size_t size1 = original.size();       // Original first dimension size
    size_t size2 = original[0].size();    // Original second dimension size
    size_t size3 = original[0][0].size(); // Original third dimension size

    // New vector with dimensions rearranged according to (2, 0, 1)
    vector<vector<vector<double>>> transposed(size3, vector<vector<double>>(size1, vector<double>(size2)));

    for (size_t i = 0; i < size1; ++i)
    {
        for (size_t j = 0; j < size2; ++j)
        {
            for (size_t k = 0; k < size3; ++k)
            {
                transposed[k][i][j] = original[i][j][k];
            }
        }
    }

    return transposed;
}

// Function to collect rows from a matrix based on a vector of indices
vec2 collectRows(const vec2 &matrix, const vector<int> &indices)
{
    // If no indices are provided, return empty matrix
    if (indices.empty())
    {
        return {};
    }
    vec2 result;
    for (int index : indices)
    {
        if (index < matrix.size())
        {
            result.push_back(matrix[index]);
        }
        else
        {
            // Handle the case where index is out of bounds
            throw std::out_of_range("Index out of range");
        }
    }
    return result;
}

// Function to collect rows from a Cube based on a vector of indices
vec3 collectRows(const vec3 &cube, const vector<int> &indices)
{
    // If no indices are provided, return empty matrix
    if (indices.empty())
    {
        return {};
    }
    vec3 result;
    for (int index : indices)
    {
        if (index < cube.size())
        {
            result.push_back(cube[index]);
        }
        else
        {
            // Handle the case where index is out of bounds
            throw std::out_of_range("Index out of range");
        }
    }
    return result;
}

// Element-wise exponentiation of a matrix
vec2 elementWiseExponentiation(const vec2 &matrix, double exponent)
{
    vec2 result = matrix; // Copy the original matrix to preserve its dimensions

    for (size_t i = 0; i < matrix.size(); ++i)
    { // Iterate over rows
        for (size_t j = 0; j < matrix[i].size(); ++j)
        {                                                    // Iterate over columns
            result[i][j] = std::pow(matrix[i][j], exponent); // Apply exponentiation
        }
    }

    return result;
}

void equicombsparse(int natoms, int nang1, int nang2, int nrad1, int nrad2,
                    const vector<vector<vector<vector<complex<double>>>>> &v1,
                    const vector<vector<vector<vector<complex<double>>>>> &v2,
                    const vector<double> &w3j, int llmax,
                    const vector<vector<int>> &llvec, int lam,
                    const vector<vector<complex<double>>> &c2r,
                    int featsize, int nfps, const vector<int64_t> &vfps,
                    vector<vector<vector<double>>> &p)
{
    // Initialize p with zeros
    p.assign(2 * lam + 1, vector<vector<double>>(nfps, vector<double>(natoms, 0.0)));

    // Declare variables at the beginning
    int iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2;
    double inner, normfact;

#pragma omp parallel for private(iat, n1, n2, il, imu, im1, im2, i, j, ifeat, iwig, l1, l2, mu, m1, m2, inner, normfact) default(none) shared(natoms, nang1, nang2, nrad1, nrad2, v1, v2, w3j, llmax, llvec, lam, c2r, nfps, vfps, p, featsize)
    for (iat = 0; iat < natoms; ++iat)
    {
        const int l21 = 2 * lam + 1;
        vec2 ptemp(l21, vector<double>(featsize, 0.0));
        cvec pcmplx(l21, complex<double>(0.0, 0.0));
        vec preal(l21, 0.0);
        inner = 0.0;

        ifeat = 0;
        for (n1 = 0; n1 < nrad1; ++n1)
        {
            for (n2 = 0; n2 < nrad2; ++n2)
            {
                iwig = 0;
                for (il = 0; il < llmax; ++il)
                {
                    l1 = llvec[0][il];
                    l2 = llvec[1][il];

                    fill(pcmplx.begin(), pcmplx.end(), complex<double>(0.0, 0.0));

                    for (imu = 0; imu < l21; ++imu)
                    {
                        mu = imu - lam;
                        for (im1 = 0; im1 < 2 * l1 + 1; ++im1)
                        {
                            m1 = im1 - l1;
                            m2 = m1 - mu;
                            if (abs(m2) <= l2)
                            {
                                im2 = m2 + l2;
                                pcmplx[imu] += w3j[iwig] * v1[im1][l1][n1][iat] * conj(v2[im2][l2][n2][iat]);
                                iwig++;
                            }
                        }
                    }

                    fill(preal.begin(), preal.end(), 0.0);
                    for (i = 0; i < l21; ++i)
                    {
                        for (j = 0; j < l21; ++j)
                        {
                            preal[i] += real(c2r[i][j] * pcmplx[j]);
                        }
                        inner += preal[i] * preal[i];
                        ptemp[i][ifeat] = preal[i];
                    }
                    ifeat++;
                }
            }
        }

        normfact = sqrt(inner);
        for (int n = 0; n < nfps; ++n)
        {
            for (imu = 0; imu < l21; ++imu)
            {
                p[imu][n][iat] = ptemp[imu][vfps[n]] / normfact;
            }
        }
    }
}

void equicomb(int natoms, int nang1, int nang2, int nrad1, int nrad2,
              cvec4 &v1,
              cvec4 &v2,
              vector<double> &w3j, int llmax,
              vector<vector<int>> &llvec, int lam,
              cvec2 &c2r, int featsize,
              vec3 &p)
{

    cout << "WARNING EQUICOMB IS HAS NOT BEEN TESTED, PROCEED WITH CAUTION!!!" << endl;
    // Initialize p with zeros
    p = vec3(2 * lam + 1, vec2(featsize, vector<double>(natoms, 0.0)));

    // Parallel region
#pragma omp parallel for default(none) shared(natoms, nang1, nang2, nrad1, nrad2, v1, v2, w3j, llmax, llvec, lam, c2r, featsize, p)
    for (int iat = 0; iat < natoms; ++iat)
    {
        double inner = 0.0;
        vec2 ptemp(2 * lam + 1, vec(featsize, 0.0));
        int ifeat = 0;
        for (int n1 = 0; n1 < nrad1; ++n1)
        {
            for (int n2 = 0; n2 < nrad2; ++n2)
            {
                int iwig = 0;
                for (int il = 0; il < llmax; ++il)
                {
                    int l1 = llvec[0][il];
                    int l2 = llvec[1][il];
                    cvec pcmplx(2 * lam + 1, complex<double>(0.0, 0.0));
                    for (int imu = 0; imu < 2 * lam + 1; ++imu)
                    {
                        int mu = imu - lam;
                        for (int im1 = 0; im1 < 2 * l1 + 1; ++im1)
                        {
                            int m1 = im1 - l1;
                            int m2 = m1 - mu;
                            if (std::abs(m2) <= l2)
                            {
                                int im2 = m2 + l2;
                                pcmplx[imu] += w3j[iwig] * v1[im1][l1][n1][iat] * std::conj(v2[im2][l2][n2][iat]);
                                iwig++;
                            }
                        }
                    }
                    vec preal(2 * lam + 1);
                    // Matrix multiplication c2r * pcmplx
                    for (int i = 0; i < 2 * lam + 1; ++i)
                    {
                        cdouble sum = 0.0;
                        for (int j = 0; j < 2 * lam + 1; ++j)
                        {
                            sum += c2r[i][j] * pcmplx[j];
                        }
                        preal[i] = std::real(sum);
                        inner += std::norm(sum);
                        ptemp[i][ifeat] = preal[i];
                    }
                    ifeat++;
                }
            }
        }
        const double normfact = std::sqrt(inner);
        for (ifeat = 0; ifeat < featsize; ++ifeat)
        {
            for (int imu = 0; imu < 2 * lam + 1; ++imu)
            {
                p[imu][ifeat][iat] = ptemp[imu][ifeat] / normfact;
            }
        }
    }
}

vector<cvec2> complex_to_real_transformation(vector<int> sizes)
{
    vector<cvec2> matrices{};
    for (int i = 0; i < sizes.size(); i++)
    {
        int lval = (sizes[i] - 1) / 2;
        int st = (lval & 1) ? 1 : -1;;
        ;
        cvec2 transformed_matrix(sizes[i], vector<complex<double>>(sizes[i], 0.0));
        for (int j = 0; j < lval; j++)
        {
            transformed_matrix[j][j] = complex<double>(0.0, 1.0);
            transformed_matrix[j][sizes[i] - j - 1] = complex<double>(0.0, st);
            transformed_matrix[sizes[i] - j - 1][j] = complex<double>(1.0, 0.0);
            transformed_matrix[sizes[i] - j - 1][sizes[i] - j - 1] = complex<double>(-st, 0.0);
            st = -st;
        }
        transformed_matrix[lval][lval] = sqrt(2.0);
        // Divide each element by sqrt(2.0)
        for (auto &row : transformed_matrix)
        {
            for (auto &elem : row)
            {
                elem /= sqrt(2.0);
            }
        }
        matrices.push_back(transformed_matrix);
    }
    return matrices;
}

int get_lmax_max(unordered_map<string, int> &lmax)
{
    int lmax_max = 0;
    for (auto &[key, value] : lmax)
    {
        if (value > lmax_max)
        {
            lmax_max = value;
        }
    }
    return lmax_max;
}

void set_lmax_nmax(unordered_map<string, int> &lmax, unordered_map<string, int> &nmax, array<vector<primitive>, 35> basis_set, vector<string> species)
{
    // lmax = {"C": 5, "H":2,...} with the numbers beeing the maximum angular momentum (type) for the given atom
    // nmax = {C0: 10, C1: 7, ...} with the numbers beeing the maximum number of primitives for the given atom and type

    for (auto &spe : species)
    {
        int atom_index = get_Z_from_label(spe.c_str());
        // get the last element of the basis set for the given atom
        lmax[spe] = basis_set[atom_index].back().type;
        // initialize nmax with symbol + type
        for (int i = 0; i < basis_set[atom_index].back().type + 1; i++)
        {
            nmax[spe + to_string(i)] = 0;
        }
        // count the number of primitives for the given atom and type
        for (int i = 0; i < basis_set[atom_index].size(); i++)
        {
            nmax[spe + to_string(basis_set[atom_index][i].type)] += 1;
        }
    }
}

// ALL below RASCALINE
struct RadialBasis
{
    std::string type;
    double spline_accuracy;
};

struct CutoffFunction
{
    std::string type;
    double width;
};

struct HyperParametersDensity
{
    double cutoff;
    int max_radial;
    int max_angular;
    double atomic_gaussian_width;
    double center_atom_weight;
    RadialBasis radial_basis;
    CutoffFunction cutoff_function;
};

std::string to_json(const HyperParametersDensity &params)
{
    std::ostringstream oss;
    oss << "{\n"
        << "  \"cutoff\": " << params.cutoff << ",\n"
        << "  \"max_radial\": " << params.max_radial << ",\n"
        << "  \"max_angular\": " << params.max_angular << ",\n"
        << "  \"atomic_gaussian_width\": " << params.atomic_gaussian_width << ",\n"
        << "  \"center_atom_weight\": " << params.center_atom_weight << ",\n"
        << "  \"radial_basis\": {\n"
        << "    \"Gto\": {\n"
        << "      \"spline_accuracy\": " << params.radial_basis.spline_accuracy << "\n"
        << "    }\n"
        << "  },\n"
        << "  \"cutoff_function\": {\n"
        << "    \"ShiftedCosine\": {\n"
        << "      \"width\": " << params.cutoff_function.width << "\n"
        << "    }\n"
        << "  }\n"
        << "}";
    return oss.str();
}

string gen_parameters(double rcut, int nrad, int nang, double atomic_gaussian_width = 0.3, double center_atom_weight = 1.0, double spline_accuracy = 1e-6, double cutoff_width = 0.1)
{
    std::ostringstream oss;
    HyperParametersDensity hyper_parameters_density = {
        rcut,                           // cutoff
        nrad,                           // max_radial
        nang,                           // max_angular
        atomic_gaussian_width,          // atomic_gaussian_width
        center_atom_weight,             // center_atom_weight
        {"Gto", spline_accuracy},       // radial_basis
        {"ShiftedCosine", cutoff_width} // cutoff_function
    };
    string json_string = to_json(hyper_parameters_density);
    return json_string;
}

#if defined(_WIN32) || defined(__RASCALINE__)
// RASCALINE1
metatensor::TensorMap get_feats_projs(string filepath, double rcut1, int nrad1, int nang1, double atomic_gaussian_width, std::vector<std::string> neighspe, std::vector<std::string> species)
{
    rascaline::BasicSystems system = rascaline::BasicSystems(filepath);
    // Construct the parameters for the calculator from the inputs given
    string temp_p = gen_parameters(rcut1, nrad1, nang1, atomic_gaussian_width);
    const char *parameters = temp_p.c_str();

    // size_t nspe1 = neighspe.size();
    std::vector<std::vector<int32_t>> keys_array;
    for (int l = 0; l < nang1 + 1; ++l)
    {
        for (const std::string &specen : species)
        {
            for (const std::string &speneigh : neighspe)
            {
                // Directly emplace back initializer_lists into keys_array
                keys_array.emplace_back(std::vector<int32_t>{l, 1, get_Z_from_label(specen.c_str()) + 1, get_Z_from_label(speneigh.c_str()) + 1});
            }
        }
    }

    // Assuming metatensor::Labels expects a flat sequence of integers for each label
    std::vector<int32_t> flattened_keys;
    for (const auto &subVector : keys_array)
    {
        flattened_keys.insert(flattened_keys.end(), subVector.begin(), subVector.end());
    }

    // Convert keys_array to rascaline::Labels
    std::vector<std::string> names = {"o3_lambda", "o3_sigma", "center_type", "neighbor_type"};
    metatensor::Labels keys_selection(names, flattened_keys.data(), flattened_keys.size() / names.size());

    // create the calculator with its name and parameters
    rascaline::Calculator calculator = rascaline::Calculator("spherical_expansion", parameters);

    rascaline::CalculationOptions calc_opts;
    calc_opts.selected_keys = keys_selection;
    // run the calculation
    metatensor::TensorMap descriptor = calculator.compute(system, calc_opts);

    // The descriptor is a metatensor `TensorMap`, containing multiple blocks.
    // We can transform it to a single block containing a dense representation,
    // with one sample for each atom-centered environment.
    descriptor = descriptor.keys_to_samples("center_type");
    descriptor = descriptor.keys_to_properties("neighbor_type");
    // descriptor.save("spx_pred.npy");

    return descriptor;
}

// RASCALINE2
cvec4 get_expansion_coeffs(vector<uint8_t> spx_buff, int n_atoms, int nang, int nrad, int nspe)
{
    metatensor::TensorMap spx = metatensor::TensorMap::load_buffer(spx_buff);
    vector<vector<vector<vector<complex<double>>>>> omega(nang + 1, vector<vector<vector<complex<double>>>>(n_atoms, vector<vector<complex<double>>>(2 * nang + 1, vector<complex<double>>(nspe * nrad, {0.0, 0.0}))));
    for (int l = 0; l < nang + 1; ++l)
    {
        cvec2 c2r = complex_to_real_transformation({(2 * l) + 1})[0];
        metatensor::TensorBlock spx_block = spx.block_by_id(l);
        metatensor::NDArray<double> spx_values = spx_block.values();

        // Perform the matrix multiplication and assignment
        for (int a = 0; a < n_atoms; ++a)
        {
            for (int c = 0; c < 2 * l + 1; ++c)
            {
                for (int d = 0; d < nspe * nrad; ++d)
                {
                    omega[l][a][c][d] = 0.0;
                    for (int r = 0; r < 2 * l + 1; ++r)
                    {
                        omega[l][a][c][d] += conj(c2r[r][c]) * spx_values(a, r, d);
                    }
                }
            }
        }
        c2r.clear();
    }

    vector<vector<vector<vector<complex<double>>>>> v(2 * nang + 1, vector<vector<vector<complex<double>>>>(nang + 1, vector<vector<complex<double>>>(nspe * nrad, vector<complex<double>>(n_atoms))));

    for (int a = 0; a < omega.size(); ++a)
    {
        for (int b = 0; b < omega[0].size(); ++b)
        {
            for (int c = 0; c < omega[0][0].size(); ++c)
            {
                for (int d = 0; d < omega[0][0][0].size(); ++d)
                {
                    v[c][a][d][b] = omega[a][b][c][d];
                }
            }
        }
    }
    return v;
}
#endif

// Function to filter out atoms that belong to species not available for the model selected
std::vector<std::string> filter_species(const std::vector<std::string> &atomic_symbols, const std::vector<std::string> &species)
{
    std::vector<std::string> filtered_symbols;
    std::set<std::string> excluded_species;

    // Convert species vector to a set for efficient lookup
    std::set<std::string> species_set(species.begin(), species.end());

    // Find all species that are not in the input species set
    for (const auto &symbol : atomic_symbols)
    {
        if (species_set.find(symbol) == species_set.end())
        {
            excluded_species.insert(symbol);
        }
    }

    // Filter out excluded species from atomic_symbols
    for (const auto &symbol : atomic_symbols)
    {
        if (excluded_species.find(symbol) == excluded_species.end())
        {
            filtered_symbols.push_back(symbol);
        }
    }

    return filtered_symbols;
}

struct Config
{
    string predict_filename;
    bool average;
    bool field;
    bool sparsify;
    int ncut;
    std::vector<std::string> species;
    double rcut1;
    double rcut2;
    int nang1;
    int nang2;
    int nrad1;
    int nrad2;
    double sig1;
    double sig2;
    std::vector<std::string> neighspe1;
    std::vector<std::string> neighspe2;
    double zeta;
    int Menv;
    int Ntrain;
    float trainfrac;
};

void populateConfigFromFile(const std::string &filename, Config &config)
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
                // Remove leading and trailing whitespaces from value
                value = trim(value);
                key = trim(key);

                // Populate the Config struct based on key
                if (key == "average")
                    config.average = (value == "True");
                else if (key == "field")
                    config.field = (value == "True");
                else if (key == "sparsify")
                    config.sparsify = (value == "True");
                else if (key == "ncut")
                    config.ncut = std::stoi(value);
                else if (key == "species")
                {
                    // Parse the species vector
                    size_t start = value.find("[");
                    size_t end = value.find("]");
                    std::string speciesList = value.substr(start + 1, end - start - 1);
                    // Remove quotes and spaces from individual species
                    speciesList.erase(std::remove_if(speciesList.begin(), speciesList.end(), [](char c)
                                                     { return std::isspace(c) || c == '"' || c == '\''; }),
                                      speciesList.end());
                    std::istringstream speciesStream(speciesList);
                    std::string speciesItem;
                    while (std::getline(speciesStream, speciesItem, ','))
                    {
                        config.species.push_back(speciesItem);
                    }
                }
                else if (key == "rcut1")
                    config.rcut1 = std::stod(value);
                else if (key == "rcut2")
                    config.rcut2 = std::stod(value);
                else if (key == "nang1")
                    config.nang1 = std::stoi(value);
                else if (key == "nang2")
                    config.nang2 = std::stoi(value);
                else if (key == "nrad1")
                    config.nrad1 = std::stoi(value);
                else if (key == "nrad2")
                    config.nrad2 = std::stoi(value);
                else if (key == "sig1")
                    config.sig1 = std::stod(value);
                else if (key == "sig2")
                    config.sig2 = std::stod(value);
                else if (key == "neighspe1")
                {
                    // Parse the neighspe1 vector
                    size_t start = value.find("[");
                    size_t end = value.find("]");
                    std::string neighspe1List = value.substr(start + 1, end - start - 1);
                    neighspe1List.erase(std::remove_if(neighspe1List.begin(), neighspe1List.end(), [](char c)
                                                       { return std::isspace(c) || c == '"' || c == '\''; }),
                                        neighspe1List.end());
                    std::istringstream neighspe1Stream(neighspe1List);
                    std::string neighspe1Item;
                    while (std::getline(neighspe1Stream, neighspe1Item, ','))
                    {
                        config.neighspe1.push_back(neighspe1Item);
                    }
                }
                else if (key == "neighspe2")
                {
                    // Parse the neighspe2 vector
                    size_t start = value.find("[");
                    size_t end = value.find("]");
                    std::string neighspe2List = value.substr(start + 1, end - start - 1);
                    neighspe2List.erase(std::remove_if(neighspe2List.begin(), neighspe2List.end(), [](char c)
                                                       { return std::isspace(c) || c == '"' || c == '\''; }),
                                        neighspe2List.end());
                    std::istringstream neighspe2Stream(neighspe2List);
                    std::string neighspe2Item;
                    while (std::getline(neighspe2Stream, neighspe2Item, ','))
                    {
                        config.neighspe2.push_back(neighspe2Item);
                    }
                }
                else if (key == "zeta")
                    config.zeta = std::stod(value);
                else if (key == "Menv")
                    config.Menv = std::stoi(value);
                else if (key == "Ntrain")
                    config.Ntrain = std::stoi(value);
                else if (key == "trainfrac")
                    config.trainfrac = std::stof(value);
            }
        }
    }
}

vec predict(const WFN& wavy, const string model_folder)
{
    // Read the configuration file
    cout << "Starting prediction... " << flush;
    Config config;
    string _f_path("inputs.txt");
    string _path = model_folder;
    join_path(_path, _f_path);
    populateConfigFromFile(_path, config);

    config.predict_filename = wavy.get_path();
    int nspe1 = static_cast<int>(config.neighspe1.size());
    int nspe2 = static_cast<int>(config.neighspe2.size());
    /* bool average = true;
     bool sparsify = true;
     int ncut = 1000;
     bool field = false;
     double rcut1 = 4.0;
     int nang1 = 7;
     int nang2 = 7;
     int nrad1 = 7;
     int nrad2 = 7;
     double sig1 = 0.3;

     float zeta = 2.0;
     int Ntrain = 600;
     float trainfrac = 0.8;
     std::vector<std::string> neighspe1 = { "H", "C", "N", "O", "S"};
     std::vector<std::string> species = { "H", "C", "N", "O", "S" };
     int nspe1 = neighspe1.size();
     int nspe2 = nspe2;*/

    vector<string> atomic_symbols{};
    for (int i = 0; i < wavy.atoms.size(); i++)
    {
        atomic_symbols.push_back(wavy.atoms[i].label);
    }
    //# Define system excluding atoms that belong to species not listed in SALTED input
    atomic_symbols = filter_species(atomic_symbols, config.species);
    //// Output the result
    // for (const auto& symbol : atomic_symbols) {
    //     std::cout << symbol << " ";
    // }
    // std::cout << std::endl;

    int natoms = static_cast<int>(atomic_symbols.size());
    unordered_map<string, vector<int>> atom_idx{};
    unordered_map<string, int> natom_dict{};
    for (int i = 0; i < atomic_symbols.size(); i++)
    {
        atom_idx[atomic_symbols[i]].push_back(i);
        natom_dict[atomic_symbols[i]] += 1;
    }

    unordered_map<string, int> lmax{};
    unordered_map<string, int> nmax{};
    set_lmax_nmax(lmax, nmax, QZVP_JKfit, config.species);
#if defined(_WIN32) || defined(__RASCALINE__)
    // RASCALINE (Generate descriptors)
    metatensor::TensorMap spx = get_feats_projs(config.predict_filename, config.rcut1, config.nrad1, config.nang1, config.sig1, config.neighspe1, config.species);
    vector<uint8_t> spx_buf = spx.save_buffer();
    cvec4 v1 = get_expansion_coeffs(spx_buf, natoms, config.nang1, config.nrad1, nspe1);

    metatensor::TensorMap spx_2 = get_feats_projs(config.predict_filename, config.rcut2, config.nrad2, config.nang2, config.sig2, config.neighspe2, config.species);
    vector<uint8_t> spx_buf_2 = spx.save_buffer();
    cvec4 v2 = get_expansion_coeffs(spx_buf_2, natoms, config.nang2, config.nrad2, nspe2);
    // END RASCALINE
#else
    err_not_impl_f("RASCALINE is not supported by this build", std::cout);
    cvec4 v2, v1;
#endif

    // Read Model variables
    unordered_map<string, vec2> Vmat{};
    unordered_map<string, int> Mspe{};
    unordered_map<string, vec2> power_env_sparse{};
    // Define zeta as a string with one decimal
    std::ostringstream stream;
    stream << std::fixed << std::setprecision(1) << config.zeta;
    std::string zeta_str = stream.str();
#if defined(_WIN32) || defined(__RASCALINE__)
    H5::H5File features(model_folder + "/GPR_data/FEAT_M-" + to_string(config.Menv) + ".h5", H5F_ACC_RDONLY);
    H5::H5File projectors(model_folder + "/GPR_data/projector_M" + to_string(config.Menv) + "_zeta" + zeta_str + ".h5", H5F_ACC_RDONLY);
    for (string spe : config.species)
    {
        for (int lam = 0; lam < lmax[spe] + 1; lam++)
        {
            power_env_sparse[spe + to_string(lam)] = readHDF5(features, "sparse_descriptors/" + spe + "/" + to_string(lam));
            Vmat[spe + to_string(lam)] = readHDF5(projectors, "projectors/" + spe + "/" + to_string(lam));

            if (lam == 0)
            {
                Mspe[spe] = static_cast<int>(power_env_sparse[spe + to_string(lam)].size());
            }
            if (config.zeta == 1)
            {
                vec2 Vmat_t = transpose<double>(Vmat[spe + to_string(lam)]);
                power_env_sparse[spe + to_string(lam)] = dot<double>(Vmat_t, power_env_sparse[spe + to_string(lam)]);
            }
        }
    }
    features.close();
    projectors.close();
#endif

    unordered_map<int, vector<int64_t>> vfps{};
    if (config.sparsify)
    {
        vfps = read_fps<int64_t>(model_folder + "/GPR_data/fps" + to_string(config.ncut) + "-", get_lmax_max(lmax));
    };

    int ntrain = static_cast<int>(config.Ntrain * config.trainfrac);
    vector<double> weights{};
    if (config.field)
    {
        cout << "Field" << endl;
        err_not_impl_f("Calculations using 'Field = True' are not yet supported", std::cout);
    }
    else
    {
        read_npy(model_folder + "/GPR_data/weights_N" + to_string(ntrain) + "_reg-6.npy", weights);
    }
    // END READ MODEL VARIABLES

    // Compute equivariant descriptors for each lambda value entering the SPH expansion of the electron density
    vec2 pvec_l0{};
    unordered_map<int, vec3> pvec{};
    for (int lam = 0; lam < get_lmax_max(lmax) + 1; lam++)
    {
        int llmax = 0;
        unordered_map<int, vector<int>> lvalues{};
        for (int l1 = 0; l1 < config.nang1 + 1; l1++)
        {
            for (int l2 = 0; l2 < config.nang2 + 1; l2++)
            {
                // keep only even combination to enforce inversion symmetryc
                if ((lam + l1 + l2) % 2 == 0)
                {
                    if (abs(l2 - lam) <= l1 && l1 <= (l2 + lam))
                    {
                        lvalues[llmax] = {l1, l2};
                        llmax += 1;
                    }
                }
            }
        }
        // Fill dense array from dictionary
        vector<vector<int>> llvec(llmax, vector<int>(2));
        for (int i = 0; i < llmax; i++)
        {
            llvec[i] = lvalues[i];
        }

        vector<double> wigner3j = readVectorFromFile<double>(model_folder + "/wigners/wigner_lam-" + to_string(lam) + "_lmax1-" + to_string(config.nang1) + "_lmax2-" + to_string(config.nang2) + ".dat");

        cvec2 c2r = complex_to_real_transformation({2 * lam + 1})[0];

        int featsize = nspe1 * nspe2 * config.nrad1 * config.nrad2 * llmax;
        vec3 p;
        vector<vector<int>> llvec_t = transpose<int>(llvec);
        if (config.sparsify)
        {
            int nfps = static_cast<int>(vfps[lam].size());
            equicombsparse(natoms, config.nang1, config.nang2, (nspe1 * config.nrad1), (nspe2 * config.nrad2), v1, v2, wigner3j, llmax, llvec_t, lam, c2r, featsize, nfps, vfps[lam], p);
            featsize = config.ncut;
        }
        else
        {
            equicomb(natoms, config.nang1, config.nang2, (nspe1 * config.nrad1), (nspe2 * config.nrad2), v1, v2, wigner3j, llmax, llvec_t, lam, c2r, featsize, p);
        }
        wigner3j.clear();

        vec3 p_t = reorder3D<double>(p);
        vector<double> flat_p = flatten<double>(p_t);

        if (lam == 0)
        {
            pvec_l0 = reshape<double>(flat_p, Shape2D{natoms, featsize});
        }
        else
        {
            pvec[lam] = reshape<double>(flat_p, Shape3D{natoms, 2 * lam + 1, featsize});
        }
        flat_p.clear();
    }

    unordered_map<string, vec2> psi_nm{};

    for (const string &spe : config.species)
    {
        if (atom_idx.find(spe) == atom_idx.end())
        {
            continue;
        }
        vec2 power_env_sparse_t = transpose<double>(power_env_sparse[spe + "0"]);
        vec2 kernell0_nm;
        if (config.zeta == 1)
        {
            psi_nm[spe + "0"] = dot<double>(collectRows(pvec_l0, atom_idx[spe]), power_env_sparse_t);
        }
        else
        {
            kernell0_nm = dot<double>(collectRows(pvec_l0, atom_idx[spe]), power_env_sparse_t);
            vec2 kernel_nm = elementWiseExponentiation(kernell0_nm, config.zeta);
            psi_nm[spe + "0"] = dot<double>(kernel_nm, Vmat[spe + "0"]);
        }

        for (int lam = 1; lam <= lmax[spe]; ++lam)
        {
            int featsize = static_cast<int>(pvec[lam][0][0].size());
            power_env_sparse_t = transpose<double>(power_env_sparse[spe + to_string(lam)]);
            vec3 pVec_Rows = collectRows(pvec[lam], atom_idx[spe]);
            vec2 pvec_lam = reshape<double>(flatten<double>(pVec_Rows), Shape2D{natom_dict[spe] * (2 * lam + 1), featsize});

            if (config.zeta == 1)
            {
                psi_nm[spe + to_string(lam)] = dot<double>(pvec_lam, power_env_sparse_t);
            }
            else
            {
                vec2 kernel_nm = dot<double>(pvec_lam, power_env_sparse_t);
                for (size_t i1 = 0; i1 < natom_dict[spe]; ++i1)
                {
                    for (size_t i2 = 0; i2 < Mspe[spe]; ++i2)
                    {
                        for (size_t i = 0; i < 2 * lam + 1; ++i)
                        {
                            for (size_t j = 0; j < 2 * lam + 1; ++j)
                            {
                                kernel_nm[i1 * (2 * lam + 1) + i][i2 * (2 * lam + 1) + j] *= pow(kernell0_nm[i1][i2], config.zeta - 1);
                            }
                        }
                    }
                }
                psi_nm[spe + to_string(lam)] = dot<double>(kernel_nm, Vmat[spe + to_string(lam)]);
            }
        }
    }
    
    unordered_map<string, vector<double>> C{};
    unordered_map<string, int> ispe{};
    int isize = 0;
    for (const string &spe : config.species)
    {
        if (atom_idx.find(spe) == atom_idx.end())
        {
            for (int l = 0; l <= lmax[spe]; ++l)
            {
                for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
                {
                    isize += static_cast<int>(Vmat[spe + to_string(l)].size());
                }
            }
            continue;
        }
        ispe[spe] = 0;
        for (int l = 0; l <= lmax[spe]; ++l)
        {
            for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
            {
                int Mcut = static_cast<int>(psi_nm[spe + to_string(l)][0].size());
                vector<double> weights_subset(weights.begin() + isize, weights.begin() + isize + Mcut);
                C[spe + to_string(l) + to_string(n)] = dot(psi_nm[spe + to_string(l)], weights_subset);
                isize += Mcut;
            }
        }
    }

    int Tsize = 0;
    for (int iat = 0; iat < natoms; iat++)
    {
        string spe = atomic_symbols[iat];
        for (int l = 0; l < lmax[spe] + 1; l++)
        {
            for (int n = 0; n < nmax[spe + to_string(l)]; n++)
            {
                Tsize += 2 * l + 1;
            }
        }
    }

    vector<double> Av_coeffs(Tsize, 0.0);
    unordered_map<string, vector<double>> av_coefs{};
    if (config.average)
    {
        for (string spe : config.species)
        {
            read_npy(model_folder + "/averages/averages_" + spe + ".npy", av_coefs[spe]);
        }
    }

    // fill vector of predictions
    int i = 0;
    vector<double> pred_coefs(Tsize, 0.0);
    for (int iat = 0; iat < natoms; ++iat)
    {
        string spe = atomic_symbols[iat];
        for (int l = 0; l <= lmax[spe]; ++l)
        {
            for (int n = 0; n < nmax[spe + to_string(l)]; ++n)
            {
                for (int ind = 0; ind < 2 * l + 1; ++ind)
                {
                    pred_coefs[i + ind] = C[spe + to_string(l) + to_string(n)][ispe[spe] * (2 * l + 1) + ind];
                }
                if (config.average && l == 0)
                {
                    Av_coeffs[i] = av_coefs[spe][n];
                }
                i += 2 * l + 1;
            }
        }
        ispe[spe] += 1;
    }

    if (config.average)
    {
        for (i = 0; i < Tsize; i++)
        {
            pred_coefs[i] += Av_coeffs[i];
        }
    }

    cout << "          ... done!\nNumber of predicted coefficients: " << pred_coefs.size() << endl;
    // npy::npy_data<double> coeffs;
    // coeffs.data = pred_coefs;
    // coeffs.fortran_order = false;
    // coeffs.shape = { unsigned long(pred_coefs.size()) };
    // npy::write_npy("coeffs_by_black_magic.npy", coeffs);
    return pred_coefs;
}

//Wrapper function to generate the ML density prediction
vec gen_SALTED_densities(const WFN& wave, options opt, time_point& start, time_point& end_SALTED) {
    // Run generation of tsc file
    if (opt.debug) cout << "Finished ML Density Prediction!" << endl;
    vec coefs = predict(wave, opt.SALTED_DIR);
    end_SALTED = get_time();
    if (opt.debug) {
        long long int dur = get_sec(start, end_SALTED);
        cout << "Finished ML Density Prediction!" << endl;
        if (dur < 1)
            cout << "Time for SALTED prediction: " << fixed << setprecision(0) << get_msec(start, end_SALTED) << " ms" << endl;
        else
            cout << "Time for SALTED prediction: " << fixed << setprecision(0) << dur << " s" << endl;
    }
    if (opt.debug && (opt.wfn == string("test_cysteine.xyz") || opt.wfn == string("test_sucrose.xyz"))) {
        vector<unsigned long> shape{};
        bool fortran_order;
        vec ref_coefs{};
        //Depending on the value of opt.wfn read either the cysteine or the sucrose reference coefs
        if (opt.wfn == string("test_sucrose.xyz"))
            npy::LoadArrayFromNumpy("sucrose_ref_Combined_v1.npy", shape, fortran_order, ref_coefs);
        else
            npy::LoadArrayFromNumpy("cysteine_ref_Cysteine.npy", shape, fortran_order, ref_coefs);
        // Compare coefs with the reference
        vector<double> diff_vec;
        double diff = 0.0;
        for (int i = 0; i < coefs.size(); i++)
        {
            diff += abs(coefs[i] - ref_coefs[i]);
            if (abs(coefs[i] - ref_coefs[i]) > 1e-4) {
                cout << "Difference in coef " << fixed << setprecision(3) << i << " : " << coefs[i] << " - " << ref_coefs[i] << " = " << abs(coefs[i] - ref_coefs[i]) << endl;
            }
            diff_vec.push_back((coefs[i] / ref_coefs[i]) - 1);
        }
        cout << "Difference between calculated and reference coefs: " << fixed << setprecision(3) << diff << endl;
        cout << "Maximum ((pred / ref) -1): " << fixed << setprecision(3) << *max_element(diff_vec.begin(), diff_vec.end()) << endl;
        if (diff > 0.1) {
			cout << "WARNING: The difference between the calculated and reference coefficients is too large!" << endl;
            exit(1);
        }
    }

    return coefs;
}
