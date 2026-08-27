#include "pch.h"
#include "convenience.h"
#include "cell.h"
#include "tsc_block.h"
#include "tsc_label_converter.h"
#include "test_functions.h"
#include "integrator.h"
#include "properties.h"
#include "wfn_class.h"
#include "atoms.h"
#include "basis_set.h"
#include "fchk.h"
#include "bondwise_analysis.h"

#ifdef _WIN32
#include <windows.h>
#include <commdlg.h>
#include <cderr.h>
#endif

namespace {

    bool is_valid_occ_data_path(const std::filesystem::path &base_path)
    {
        return std::filesystem::is_directory(base_path) &&
            std::filesystem::is_directory(base_path / "basis") &&
            std::filesystem::is_directory(base_path / "methods");
    }

    std::filesystem::path resolve_executable_directory(const char *argv0)
    {
        // Prefer argv0 — it comes directly from the program's own argument vector
        // and avoids reading /proc/self (which Flawfinder flags as a potential
        // information-hiding / unterminated-result risk on Linux).
        if (argv0 != nullptr)
        {
            try
            {
                auto p = std::filesystem::absolute(std::filesystem::path(argv0)).parent_path();
                if (!p.empty() && std::filesystem::exists(p))
                    return p;
            }
            catch (...)
            {
                // Filesystem operations may throw (e.g. on very unusual path strings);
                // fall through to the platform-specific fallback below.
            }
        }

        // Fall back to platform-specific APIs when argv0 is unavailable or unusable.
#ifdef _WIN32
        char module_path[MAX_PATH];
        DWORD len = GetModuleFileNameA(nullptr, module_path, MAX_PATH);
        if (len > 0 && len < MAX_PATH)
        {
            return std::filesystem::path(module_path).parent_path();
        }
#elif defined(__APPLE__)
        uint32_t size = 0;
        _NSGetExecutablePath(nullptr, &size);
        std::string exe_path(size, '\0');
        if (_NSGetExecutablePath(exe_path.data(), &size) == 0)
        {
            return std::filesystem::path(exe_path.c_str()).parent_path();
        }
#else
        char exe_path[4096];
        ssize_t len = readlink("/proc/self/exe", exe_path, sizeof(exe_path) - 1);
        if (len > 0)
        {
            exe_path[len] = '\0';
            return std::filesystem::path(exe_path).parent_path();
        }
#endif

        return {};
    }

    bool set_occ_data_path(const std::filesystem::path &path)
    {
#ifdef _WIN32
        return _putenv_s("OCC_DATA_PATH", path.string().c_str()) == 0;
#else
        return setenv("OCC_DATA_PATH", path.string().c_str(), 1) == 0;
#endif
    }

    std::filesystem::path choose_occ_data_path_from_exe_dir(const std::filesystem::path &exe_dir)
    {
        const auto parent = exe_dir.parent_path();
        const auto grandparent = parent.parent_path();

        const std::vector<std::filesystem::path> candidates = {
            exe_dir,
            exe_dir / "occ" / "share",
            exe_dir / "share",
            parent / "occ" / "share",
            parent / "share",
            grandparent / "occ" / "share"
        };

        for (const auto &candidate : candidates)
        {
            if (!candidate.empty() && is_valid_occ_data_path(candidate))
            {
                return candidate;
            }
        }
        return {};
    }

    ivec parse_rgbi_group_indices(const std::string &spec)
    {
        ivec indices;
        for (const auto &raw_part : split_string<std::string>(spec, ","))
        {
            const std::string part = trim(raw_part);
            if (part.empty())
                continue;

            const size_t range_pos = part.find('-', 1);
            if (range_pos == std::string::npos)
            {
                indices.push_back(std::stoi(part));
                continue;
            }

            const std::string first_text = trim(part.substr(0, range_pos));
            const std::string last_text = trim(part.substr(range_pos + 1));
            err_checkf(!first_text.empty() && !last_text.empty(),
                "Invalid RGBI atom index range: " + part, std::cout);

            const int first = std::stoi(first_text);
            const int last = std::stoi(last_text);
            const int step = first <= last ? 1 : -1;
            for (int atom_index = first;; atom_index += step)
            {
                indices.push_back(atom_index);
                if (atom_index == last)
                    break;
            }
        }
        return indices;
    }

    void validate_rgbi_group_set(const ivec2 &group_set)
    {
        std::map<int, int> assigned_group;
        for (int group_index = 0; group_index < static_cast<int>(group_set.size()); ++group_index)
        {
            for (int atom_index : group_set[group_index])
            {
                auto [it, inserted] = assigned_group.emplace(atom_index, group_index);
                err_checkf(inserted,
                    "Invalid RGBI groups: atom index " + std::to_string(atom_index) +
                    " is assigned to both group " + std::to_string(it->second) +
                    " and group " + std::to_string(group_index) + ".",
                    std::cout);
            }
        }
    }

    // **One definition of the geometry-aid hyperparameters, used by both the
    // single-structure and the batch flag.** Two copies of this block would be
    // the worst kind of bug available here: a descriptor of the right length
    // computed with the wrong settings is not rejected by anything downstream
    // and produces confident nonsense.
    //
    // These must match what the geometry-aid models were trained with, in
    //   geometry-aid/multi_layer_classifier/c_only_training.py :: SOAP_HP
    // and NOT the older values in geometry-aid/external_script.py.
    //
    // The feature count is the check: 11 species give 66 unique pairs, and the
    // descriptor length is 66 * (max_radial+1)^2 * (max_angular+1).
    //   old: 66 * 5^2 * 10 = 16,500
    //   new: 66 * 7^2 * 13 = 42,042   <- what every shipped model expects
    //
    // SALTED is unaffected: it builds its own FeatomicHyperParameters from
    // config.nang1 / config.nang2 in SALTED_predictor.cpp.
    // **The one field that differs between the two shipped model families.**
    //
    // geometry-aid trains two variants and saves the hyperparameters beside
    // each model, which is what settles this:
    //
    //   c_only : cutoff radius 3.5   trained on all-carbon input
    //   dirty  : cutoff radius 3.0   trained on partly *wrong* labels
    //
    // Every other field is identical -- smoothing 0.7, density width 0.2,
    // max_angular 12, max_radial 6, spline 1e-6 -- so both descriptors are
    // 42,042 long and nothing downstream needs to change shape.
    //
    // The distinction matters because only the `dirty` model can be *iterated*.
    // Dora's thesis pipeline predicts, relabels the atoms, and recomputes the
    // descriptor on the relabelled structure -- SOAP is species-aware, so
    // correcting one atom improves its neighbours' descriptors. Feeding
    // relabelled input to `c_only` is out of distribution for it, which is why
    // an earlier attempt at exactly this scored oxygen at 0.05 against 0.44 and
    // was abandoned. It was the wrong model, not the wrong idea.
    //
    // Default stays 3.5 so existing behaviour and every shipped result is
    // unchanged; `-geometry_aid_cutoff 3.0` selects the iterable variant.
    double geometry_aid_cutoff_radius = 3.5;

    SALTED_Utils::FeatomicHyperParameters geometry_aid_hyperparameters()
    {
        const std::vector<std::string> species{ "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br", "I" };

        // **Diagnostic override, never for production output.** Two guesses at
        // where the fixed ~0.7 s of a descriptor call goes have been wrong: it
        // is not the calculator constructor (0.010 s, and cached it makes no
        // difference to a batch) and it is not the 65 empty species pairs
        // (compute is flat in both atom count and element count). What is left
        // is the radial-integral splining, which depends only on
        // `spline_accuracy` -- so loosening it is the one measurement that
        // settles it. A descriptor computed at a different accuracy is **not**
        // comparable with any trained model, hence the warning: a vector of the
        // right length computed with the wrong settings is rejected by nothing
        // downstream and produces confident nonsense.
        double spline_accuracy = 1E-6;
        if (const char* override_accuracy = std::getenv("NOSPHERA2_SPLINE_ACCURACY"))
        {
            spline_accuracy = std::atof(override_accuracy);
            std::cout << "  WARNING spline_accuracy overridden to " << spline_accuracy
                      << " -- this descriptor does NOT match any trained model "
                         "and must be used for timing only" << std::endl;
        }

        // Same diagnostic, for the two parameters that set the descriptor's
        // size: 66 pairs * (max_radial+1)^2 * (max_angular+1) = 42,042 today.
        // The question these answer is whether a *retrained* model on a smaller
        // descriptor would be cheaper to evaluate, which is the only remaining
        // way to cut the fixed cost inside `compute` -- everything else has
        // been measured and ruled out. Changing either invalidates every
        // shipped model, hence the same warning.
        int max_radial = 6, max_angular = 12;
        if (const char* override_radial = std::getenv("NOSPHERA2_MAX_RADIAL"))
        {
            max_radial = std::atoi(override_radial);
            std::cout << "  WARNING max_radial overridden to " << max_radial
                      << " -- timing only" << std::endl;
        }
        if (const char* override_angular = std::getenv("NOSPHERA2_MAX_ANGULAR"))
        {
            max_angular = std::atoi(override_angular);
            std::cout << "  WARNING max_angular overridden to " << max_angular
                      << " -- timing only" << std::endl;
        }

        return SALTED_Utils::FeatomicHyperParameters{
            .cutoff_radius = geometry_aid_cutoff_radius,
            .max_radial = max_radial,
            .max_angular = max_angular,
            .atomic_gaussian_width = 0.2,
            .center_atom_weight = 1.0,
            .species = species,
            .neighspe = species,
            .radial_basis = {.type = "Gto", .spline_accuracy = spline_accuracy },
            .cutoff_function = {.type = "ShiftedCosine", .width = 0.7 }
        };
    }

    // ---- geometry-aid classifier, evaluated here instead of in Python ------
    //
    // Olex2 used to get 42,042 doubles per atom on disk and do the PCA and the
    // three dense layers itself. Measured per call: 0.006 s to write the npy,
    // 0.020 s to read it, 0.139 s to decompress the model npz, 0.030 s for the
    // projection and 0.008 s for the layers -- about 0.20 s of a 1.07 s call,
    // removable at **no cost in accuracy** because it is the same arithmetic.
    //
    // The weights come from `geometry_aid_model.bin`, produced once by
    // `make_geometry_aid_bin.py`. The `.npz` is a deflated ZIP and there is no
    // zlib here; a flat file avoids adding a zip reader to this codebase.
    struct GeometryAidModel
    {
        int n_features = 0, n_components = 0, n_layers = 0, n_classes = 0;
        bool whiten = false;
        std::vector<std::string> classes;
        vec mean;                       // n_features
        vec components;                 // n_features x n_components, transposed
        vec mean_projection;            // n_components: mean . components^T
        vec explained_variance;         // n_components, only when whiten
        std::vector<int> rows, cols;
        std::vector<vec> w, b;
    };

    template <typename T>
    void read_exact(std::istream& in, T* into, size_t count, const char* what)
    {
        in.read(reinterpret_cast<char*>(into), static_cast<std::streamsize>(count*sizeof(T)));
        err_checkf(static_cast<size_t>(in.gcount()) == count*sizeof(T),
            std::string("geometry-aid model truncated while reading ") + what, std::cout);
    }

    GeometryAidModel load_geometry_aid_model(const std::filesystem::path& path)
    {
        std::ifstream in(path, std::ios::binary);
        err_checkf(in.good(), "Cannot open the geometry-aid model: " + path.string(), std::cout);

        char magic[8] = { 0 };
        read_exact(in, magic, 8, "the magic");
        err_checkf(std::string(magic, 8) == "GEOAID01",
            "This is not a GEOAID01 file. Regenerate it with "
            "make_geometry_aid_bin.py -- a stale .bin beside a newer .npz is "
            "exactly the mismatch the magic exists to catch.", std::cout);

        GeometryAidModel m;
        int header[5] = { 0 };
        read_exact(in, header, 5, "the header");
        m.n_features = header[0];
        m.n_components = header[1];
        m.n_layers = header[2];
        m.n_classes = header[3];
        m.whiten = header[4] != 0;

        for (int c = 0; c < m.n_classes; ++c)
        {
            int length = 0;
            read_exact(in, &length, 1, "a class name length");
            std::string name(static_cast<size_t>(length), '\0');
            if (length > 0) read_exact(in, name.data(), static_cast<size_t>(length), "a class name");
            m.classes.push_back(name);
        }

        m.mean.resize(static_cast<size_t>(m.n_features));
        read_exact(in, m.mean.data(), m.mean.size(), "the PCA mean");
        m.components.resize(static_cast<size_t>(m.n_features)*m.n_components);
        read_exact(in, m.components.data(), m.components.size(), "the PCA components");
        if (m.whiten)
        {
            m.explained_variance.resize(static_cast<size_t>(m.n_components));
            read_exact(in, m.explained_variance.data(), m.explained_variance.size(), "the explained variance");
        }

        for (int l = 0; l < m.n_layers; ++l)
        {
            int shape[2] = { 0, 0 };
            read_exact(in, shape, 2, "a layer shape");
            m.rows.push_back(shape[0]);
            m.cols.push_back(shape[1]);
            vec weights(static_cast<size_t>(shape[0])*shape[1]);
            read_exact(in, weights.data(), weights.size(), "a weight matrix");
            vec bias(static_cast<size_t>(shape[1]));
            read_exact(in, bias.data(), bias.size(), "a bias vector");
            m.w.push_back(std::move(weights));
            m.b.push_back(std::move(bias));
        }

        // The constant term of the projection, computed once here so the hot
        // loop can skip the descriptor's structural zeros. See
        // `classify_descriptor` for why that is worth 8x.
        m.mean_projection.assign(static_cast<size_t>(m.n_components), 0.0);
        for (size_t f = 0; f < m.mean.size(); ++f)
        {
            const double mf = m.mean[f];
            if (mf == 0.0) continue;
            const double* comp = m.components.data() + f*static_cast<size_t>(m.n_components);
            for (size_t c = 0; c < static_cast<size_t>(m.n_components); ++c)
                m.mean_projection[c] += mf*comp[c];
        }
        return m;
    }

    const GeometryAidModel& cached_geometry_aid_model(const std::filesystem::path& path)
    {
        static std::map<std::string, GeometryAidModel> cache;
        const std::string key = path.string();
        auto found = cache.find(key);
        if (found == cache.end())
            found = cache.emplace(key, load_geometry_aid_model(path)).first;
        return found->second;
    }

    // (n_atoms, n_classes) row-major probabilities.
    //
    // Summation order is not numpy's, so the last bits will differ from the
    // Python route. That is fine and it is checked the way it matters:
    // `bench_geometry_cpp.py` requires the **argmax and the full ranking** to
    // agree on every atom, which is what the assignment actually consumes.
    vec classify_descriptor(const double* descriptor, size_t n_atoms,
        size_t n_features, const GeometryAidModel& m)
    {
        err_checkf(n_features == static_cast<size_t>(m.n_features),
            "The descriptor has " + std::to_string(n_features) + " features and "
            "the model expects " + std::to_string(m.n_features) + ". These come "
            "from different SOAP hyperparameters and the result would be "
            "meaningless rather than merely worse.", std::cout);

        const size_t k = static_cast<size_t>(m.n_components);
        vec projected(n_atoms*k, 0.0);

        // **(x - mean) . C^T  ==  x . C^T  -  mean . C^T**, and only the first
        // term touches the descriptor. That matters because the descriptor is
        // sparse in exactly the way the mean is not: an all-carbon .xyz -- what
        // Olex2 sends on the first pass -- populates one of the 66 species-pair
        // blocks, so 637 of 42,042 entries are non-zero. Centring first
        // destroys that: `row[f] - mean[f]` is non-zero wherever the *mean* is,
        // which is nearly everywhere, and the skip below never fires.
        //
        // Written the obvious way this measured 0.302 s against numpy's 0.038,
        // eight times slower than the thing it was replacing. With the constant
        // term lifted out it is a sparse row times a dense matrix.
        const double* mean_projection = m.mean_projection.data();
#pragma omp parallel for
        for (long long a = 0; a < static_cast<long long>(n_atoms); ++a)
        {
            const double* row = descriptor + static_cast<size_t>(a)*n_features;
            double* out = projected.data() + static_cast<size_t>(a)*k;
            for (size_t c = 0; c < k; ++c) out[c] = -mean_projection[c];
            for (size_t f = 0; f < n_features; ++f)
            {
                const double value = row[f];
                if (value == 0.0) continue;
                const double* comp = m.components.data() + f*k;
                for (size_t c = 0; c < k; ++c) out[c] += value*comp[c];
            }
            if (m.whiten)
                for (size_t c = 0; c < k; ++c) out[c] /= std::sqrt(m.explained_variance[c]);
        }

        vec current = std::move(projected);
        size_t width = k;
        for (int l = 0; l < m.n_layers; ++l)
        {
            const size_t out_width = static_cast<size_t>(m.cols[l]);
            vec next(n_atoms*out_width, 0.0);
            const bool last = (l == m.n_layers - 1);
#pragma omp parallel for
            for (long long a = 0; a < static_cast<long long>(n_atoms); ++a)
            {
                const double* in_row = current.data() + static_cast<size_t>(a)*width;
                double* out_row = next.data() + static_cast<size_t>(a)*out_width;
                for (size_t o = 0; o < out_width; ++o) out_row[o] = m.b[l][o];
                for (size_t i = 0; i < width; ++i)
                {
                    const double v = in_row[i];
                    if (v == 0.0) continue;
                    const double* wrow = m.w[l].data() + i*out_width;
                    for (size_t o = 0; o < out_width; ++o) out_row[o] += v*wrow[o];
                }
                if (!last)
                    for (size_t o = 0; o < out_width; ++o) out_row[o] = std::max(out_row[o], 0.0);
            }
            current = std::move(next);
            width = out_width;
        }

        // softmax, shifted by the row maximum exactly as the Python does
#pragma omp parallel for
        for (long long a = 0; a < static_cast<long long>(n_atoms); ++a)
        {
            double* row = current.data() + static_cast<size_t>(a)*width;
            double biggest = row[0];
            for (size_t o = 1; o < width; ++o) biggest = std::max(biggest, row[o]);
            double total = 0.0;
            for (size_t o = 0; o < width; ++o) { row[o] = std::exp(row[o] - biggest); total += row[o]; }
            if (total <= 0.0) total = 1.0;
            for (size_t o = 0; o < width; ++o) row[o] /= total;
        }
        return current;
    }

    void write_featomic_descriptor(const std::filesystem::path& structure,
        const std::filesystem::path& out_path,
        const SALTED_Utils::FeatomicHyperParameters& hyperparams)
    {
        const bool time_phases = std::getenv("NOSPHERA2_TIME_SOAP") != nullptr;
        auto mark = std::chrono::steady_clock::now();
        auto lap = [&mark, time_phases](const char* what) {
            if (!time_phases) return;
            const auto now = std::chrono::steady_clock::now();
            std::cout << "  SOAP_PHASE " << what << " "
                      << std::chrono::duration<double>(now - mark).count() << std::endl;
            mark = now;
        };

        featomic::SimpleSystem system = SALTED_Utils::gen_featomic_system(structure);
        lap("read_structure");
        metatensor::TensorMap descriptor = SALTED_Utils::calculate_SOAP_Powerspectrum(
            std::move(system), hyperparams);
        // Reset here or the next lap spans the whole SOAP call as well, which
        // reported the 13 MB copy below as 0.6 s when it is 15 ms.
        mark = std::chrono::steady_clock::now();

        metatensor::TensorBlock temp_block = descriptor.block_by_id(0);
        metatensor::NDArray<double> temp_values = temp_block.values();
        std::vector<size_t> sizes = temp_block.values_shape();
        vec data(sizes[0] * sizes[1]);
        std::copy(temp_values.data(), temp_values.data() + data.size(), data.data());

        npy::npy_data<double> np_descr;
        np_descr.data = data;
        np_descr.fortran_order = false;
        np_descr.shape = { static_cast<unsigned long>(sizes[0]), static_cast<unsigned long>(sizes[1]) };
        lap("copy_out");
        npy::write_npy(out_path.string(), np_descr);
        lap("write_npy");
    }

    // The same descriptor, classified here and written as (n_atoms, n_classes)
    // instead of (n_atoms, 42042). For a 40-atom structure that is 3.5 kB out
    // rather than 13.5 MB.
    void write_geometry_aid_probabilities(const std::filesystem::path& structure,
        const std::filesystem::path& out_path,
        const std::filesystem::path& model_path,
        const SALTED_Utils::FeatomicHyperParameters& hyperparams)
    {
        const bool time_phases = std::getenv("NOSPHERA2_TIME_SOAP") != nullptr;
        auto mark = std::chrono::steady_clock::now();
        auto lap = [&mark, time_phases](const char* what) {
            if (!time_phases) return;
            const auto now = std::chrono::steady_clock::now();
            std::cout << "  SOAP_PHASE " << what << " "
                      << std::chrono::duration<double>(now - mark).count() << std::endl;
            mark = now;
        };

        const GeometryAidModel& model = cached_geometry_aid_model(model_path);
        lap("load_model");

        metatensor::TensorMap descriptor = SALTED_Utils::calculate_SOAP_Powerspectrum(
            SALTED_Utils::gen_featomic_system(structure), hyperparams);
        mark = std::chrono::steady_clock::now();

        metatensor::TensorBlock temp_block = descriptor.block_by_id(0);
        metatensor::NDArray<double> temp_values = temp_block.values();
        std::vector<size_t> sizes = temp_block.values_shape();
        const vec probabilities = classify_descriptor(temp_values.data(), sizes[0], sizes[1], model);
        lap("classify");

        npy::npy_data<double> np_probs;
        np_probs.data = probabilities;
        np_probs.fortran_order = false;
        np_probs.shape = { static_cast<unsigned long>(sizes[0]),
                           static_cast<unsigned long>(model.n_classes) };
        npy::write_npy(out_path.string(), np_probs);
        lap("write_npy");
    }

}

std::string help_message =
("\n============================================================================\n"
 " NoSpherA2 -- non-spherical scattering factors and density analysis\n"
 "============================================================================\n"
 "Syntax: NoSpherA2 [options]\n"
 "Values in [brackets] are optional.  Repeatable options may be supplied more\n"
 "than once.  Input files and options may be given in any order unless noted.\n\n"
 "GETTING STARTED\n"
 "  HAR/TSC:  -cif model.cif -hkl data.hkl -wfn wavefunction.wfx -acc 2\n"
 "  IAM TSC:  -cif model.cif -hkl data.hkl -xyz model.xyz -IAM -acc 2\n\n"
 "CORE INPUT AND SCATTERING-FACTOR CALCULATION\n"
 "  -cif <file.cif>                    Structure and atom labels.\n"
 "  -wfn <file>                        Wavefunction input (.wfn, .wfx, .ffn,\n"
 "                                    .molden, .gbw, .xtb, or fchk/fch).\n"
 "  -cube_density <file.cube>          Experimental/computed electron-density\n"
 "                                    cube used for partitioned SF generation\n"
 "                                    (requires -cif and -hkl/-dmin).\n"
 "  -occ <file.toml>                   Run an OCC wavefunction calculation.\n"
 "  -xyz <file.xyz>                    Atomic positions for IAM or SALTED.\n"
 "  -hkl <file.hkl>                    Reflection list.\n"
 "  -dmin <angstrom>                   Generate reflections to this d-spacing\n"
 "                                    instead of reading -hkl.  Takes precedence\n"
 "                                    over -hkl_min_max and -hkl.\n"
 "  -hkl_min_max hmin hmax kmin kmax lmin lmax\n"
 "                                    Explicit HKL bounds; used instead of -hkl,\n"
 "                                    but only when -dmin is not given.\n"
 "  -IAM                               Use Thakkar independent-atom factors.\n"
 "  -tsc_block <n>                     Reflections per block when writing the\n"
 "                                    tsc [1000]. The table is streamed to disk\n"
 "                                    a block at a time instead of being held\n"
 "                                    whole, which is what keeps a protein\n"
 "                                    inside a small machine. 0 holds it all.\n"
 "                                    Derived from -mem when that is given\n"
 "                                    and this is not.\n"
 "  -acc <0..4>                        Numerical grid accuracy [2]; 4 is the\n"
 "                                    practical maximum.\n"
 "  -group <n ...>                     CIF disorder groups for the asymmetric\n"
 "                                    unit (prefix a number with + to invert).\n"
 "  -charge <n>  -mult <n>             Override charge/multiplicity read from\n"
 "                                    the wavefunction.\n"
 "  -method <RKS|RHF>                  Select DFT or Hartree-Fock treatment.\n"
 "  -ECP [1|2|3]                       Apply ECP correction: 1=def2, 2=xTB,\n"
 "                                    3=pTB.  Alias: -ecp.\n"
 "  -ED                                Electron-diffraction scattering factors.\n"
 "  -twin <9 values>                   Add a 3x3 twin law; may be repeated.\n"
 "  -old_tsc                           Request the legacy TSC format.\n\n"
 "PARTITIONING, BASIS, AND PERFORMANCE\n"
 "  -b <basis>  -d <directory>          Basis name/file and Tonto-style basis\n"
 "                                    directory (also needed for -fchk).\n"
 "  -Becke | -TFVC | -mbis | -embis    Select partitioning (default Hirshfeld).\n"
 "  -ri_fit [basis ...]                RI partitioning; omit a basis or use\n"
 "                                    auto_aux to generate one automatically.\n"
 "  -cpus <n>                          Maximum worker threads [all available].\n"
 "  -mem <MB>                          Memory budget for everything sliceable\n"
 "                                    [unset]. When given, the tsc block size\n"
 "                                    and the XCW I tensor window are chosen\n"
 "                                    to fit it, and anything that fits whole\n"
 "                                    is held whole, that being the fastest\n"
 "                                    arrangement. An explicit -tsc_block or\n"
 "                                    XCW i_tensor_mb overrides it.\n"
 "  -pbc <n>                           Periodic-boundary setting.\n\n"
 "TSC/TSCB TABLE UTILITIES\n"
 "  -tscb <table.tsc|table.tscb>        Convert between text .tsc and binary\n"
 "                                    .tscb, preserving the basename.\n"
 "  -tsc_labels [<table> <model.cif> [output.tsc]]\n"
 "                                    With table+cif: resolve SCATTERER_IDS\n"
 "                                    against CIF atoms and write labels .tsc\n"
 "                                    (default: <table>.labels.tsc). Without\n"
 "                                    operands: force label-based scatterers in\n"
 "                                    binary .tscb output (default behavior).\n"
 "  -merge <table ...>                 Merge .tsc/.tscb tables with checks.\n"
 "  -merge_nocheck <table ...>         Merge tables only when you know their\n"
 "                                    HKL indices are identical.\n\n"
 "MULTI-FRAGMENT WORKFLOWS\n"
 "  -mtc <wfn groups ...>              Multi-TSC calculation: alternating\n"
 "                                    wavefunction and comma-separated groups.\n"
 "  -cmtc <wfn cif groups ...>          CIF-based multi-TSC: alternating\n"
 "                                    wavefunction, CIF, and groups.\n"
 "  -mtc_mult <n ...>                  Per-input multiplicities for -mtc/-cmtc.\n"
 "  -mtc_charge <n ...>                Per-input charges (negative: n<value>).\n"
 "  -mtc_ECP <mode ...>                Per-input ECP modes.\n"
 "  -Cation <label ...>  -Anion <label ...>\n"
 "                                    Explicit ionic fragment labels.\n\n"
 "DENSITY, PROPERTIES, AND ANALYSIS\n"
 "  -rho_cube <wfn>                    Write an electron-density cube.\n"
 "  -elf  -eli  -lap  -rdg  -esp       Request ELF, ELI, Laplacian, RDG, or\n"
 "                                    electrostatic-potential properties.\n"
 "  -def  -HDEF                        Request deformation density or HDEF.\n"
 "  -MO <number|all>                   Generate a molecular-orbital property.\n"
 "  -fukui                             Fukui functions f+/f-/f0 and the dual\n"
 "                                    descriptor, in the frontier-orbital\n"
 "                                    approximation. Writes _fukui_plus,\n"
 "                                    _fukui_minus, _fukui_zero and\n"
 "                                    _dual_descriptor cubes, plus condensed\n"
 "                                    per-atom values for all five partitions.\n"
 "  -fukui_analysis <wfn>              Reactivity analysis of one wavefunction:\n"
 "                                    frontier orbitals, HOMO-LUMO gap and the\n"
 "                                    condensed Fukui functions under Hirshfeld,\n"
 "                                    Becke, TFVC, MBIS and EMBIS. No cubes, so\n"
 "                                    no grid, radius or CIF is needed.\n"
 "  -radius <angstrom>  -resolution <angstrom>\n"
 "                                    Grid settings for property calculations.\n"
 "  -hirsh <atom-index>                Hirshfeld analysis for one atom.\n"
 "  -hirshfeld_surface <wfn1> <wfn2>  Hirshfeld-surface analysis.\n"
 "  -rgbi                              Roby-Gould bond-index analysis.\n"
 "  -rgbi_no_sym                       RGBI without atomic O_h symmetrization.\n"
 "  -rgbi_basis <nao|ano>              RGBI basis: occupied NAO [nao] or ANO.\n"
 "  -rgbi-groups <range ...>           RGBI groups, e.g. 0-5,7; repeat option\n"
 "                                    for multiple group sets.\n"
 "  -promol_nci <a.xyz> <b.xyz> [rcut1 rcut2 rho_max rdg_max]\n"
 "                                    Promolecular NCI/RDG outputs. Defaults:\n"
 "                                    rcut1=0.95 and rcut2=0.75.\n"
 "  -promol_nci_single_thread          Disable NCI parallel processing.\n"
 "  -qtaim_eli <rho.cube> <eli.cube> <atoms> [background]\n"
 "  -qtaim_eli <wfn> <atoms> [resolution radius background]\n"
 "                                    Keep ELI only in QTAIM basins of 0-based,\n"
 "                                    comma-separated atom indices.\n"
 "  -laplacian_bonds <wfn>             Plot density Laplacian along bonds.\n\n"
 "CONVERSION, ML, AND SPECIALISED TOOLS\n"
 "  -gbw2wfn -wfn <file.gbw>            Convert GBW input to .wfn.\n"
 "  -convert_to_47 <wfn>               Write an NBO File47 (.47).\n"
 "  -fchk <output.fchk>                Write FCHK output (requires -b and -d).\n"
 "  -SALTED <model-dir>                Predict density with a SALTED model.\n"
 "  -SALTED_COEFS <model-dir>          Write SALTED_COEFS.npy (requires -wfn).\n"
 "  -RI_CUBE <coefficients.npy>        Write an RI density cube; use -wfn and\n"
 "                                    -ri_fit first.\n"
 "  -write_ri_coefs                    Write RI_COEFS.npy; use -wfn and\n"
 "                                    -ri_fit first.\n"
 "  -combine_mos <wfn1> <wfn2>          Combine molecular orbitals.\n"
 "  -cmos1 <MO ...>  -cmos2 <MO ...>   MO selections for -combine_mos.\n"
 "  -QCT                               Enter the legacy QCT workflow.\n\n"
 "DIAGNOSTICS AND INTERNAL TOOLS\n"
 "  -v | -v2 | -debug                  Verbose diagnostic output.\n"
 "  -profiling [tests-root]            Run the internal profiling suite\n"
 "                                    [./tests]. Alias: -profile.\n"
 "  -no-date                           Suppress date information in output.\n"
 "  -draw_orbits l,m[,resolution,radius]\n"
 "                                    Draw a spherical-harmonic orbital.\n"
 "  -eli_analysis <wfn> <resolution> <radius>\n"
 "  -ewal_sum <cube> [kmax] [accuracy] Ewald sum of a cube.\n"
 "  -atom_dens <wfn> [alpha-MOs beta-MOs]\n"
 "  -atom_dens_diff <gbw1> <gbw2>      Difference density from two GBW files.\n"
 "  -spherical_aver_fukui <wfn1> <wfn2>\n"
 "  -spherical_aver_hirsh <wfn>\n"
 "  -dipole_moments  -polarizabilities <7 wfn files>\n"
 "  -geometry_aid_cutoff <r>            SOAP cutoff for the geometry-aid descriptor:\n"
 "                                      3.5 (default) matches the c_only models,\n"
 "                                      3.0 matches the dirty models. Give it BEFORE\n"
 "                                      the descriptor flag.\n"
 "  -calc_featomic_descriptor           Write descriptor.npy (requires -wfn).\n"
 "  -calc_featomic_descriptors <list>   Same, for many structures in one run.\n"
 "                                      <list> holds one structure path per\n"
 "                                      line; each descriptor is written as\n"
 "                                      <path>.npy. The radial-basis splines\n"
 "                                      cost a fixed 0.7 s and are built once,\n"
 "                                      so a batch is far cheaper than N runs.\n"
 "  -classify_atoms <model.bin> [out]   Element probabilities per atom from the\n"
 "                                      geometry-aid model (requires -wfn).\n"
 "                                      Writes probabilities.npy, one row per\n"
 "                                      atom, instead of the 42,042-wide\n"
 "                                      descriptor. Make the .bin with\n"
 "                                      make_geometry_aid_bin.py.\n"
 "  -classify_atoms_list <list> <model.bin>\n"
 "                                      Same, for many structures in one run.\n"
 "  -SALTED_Training                    Create SALTED training data.\n"
 "  -all_charges                        Print all calculated atomic charges.\n"
 "  -atom_sfac <wfn1> <wfn2>            Compare atom-centred scattering factors.\n"
 "  -density_difference <wfn>           Use a second wavefunction for a\n"
 "                                    density-difference calculation.\n"
 "  -e_field <value>                    Apply an external electric field.\n"
 "  -fractal <name>                     Enable the legacy fractal workflow.\n"
 "  -get_g                              Enable reciprocal-space g calculation.\n"
 "  -refine [accuracy]                  Set refinement integral accuracy [0.1].\n"
 "  -rgbi_EVs                           Include RGBI eigenvectors.\n"
 "  -sfac_diffuse x y z cif wfn dmin    Calculate diffuse scattering factors.\n\n"
 "EXPERIMENTAL AND DEVELOPER COMMANDS\n"
 "  -coef <file>                        Use externally supplied SALTED\n"
 "                                    coefficients.\n"
 "  -convert_XCW <stdout> <lambda-step> Convert Tonto XCW lambda-step output.\n"
 "  -do_XCW  -calc_F  -anom_disp <file> XCW/Fcalc/anomalous-dispersion modes.\n"
 "  -XCW_settings <file>                Keywords for -do_XCW. Besides the\n"
 "                                    refinement settings, three control how\n"
 "                                    the I tensor - nr_refl blocks of\n"
 "                                    nmo(nmo+1)/2 complex doubles, quadratic\n"
 "                                    in the basis and usually the largest\n"
 "                                    thing in the process - is kept:\n"
 "                                      stream          put it on disk and\n"
 "                                                      read a window of\n"
 "                                                      reflections at a time\n"
 "                                      i_tensor_mb <n> the same, with the\n"
 "                                                      budget named in MB\n"
 "                                      safe / read     write or reuse the\n"
 "                                                      whole tensor as\n"
 "                                                      I_tensor (not\n"
 "                                                      combinable with a\n"
 "                                                      budget)\n"
 "                                    Without any of them the tensor is held\n"
 "                                    in memory, which is the fastest and what\n"
 "                                    -mem also picks whenever it fits.\n"
 "                                    Prefer a GENEROUS budget: holding 12%\n"
 "                                    of the tensor cost 8% more time than\n"
 "                                    holding all of it, while shrinking to\n"
 "                                    0.4% saved a further 7% of memory and\n"
 "                                    cost another 9% of time.\n"
 "  -partitioning_test  -NNLS_TEST      Run internal partitioning/NNLS tests.\n"
 "  -test_RI  -RI_WFN_DIFF              RI fitting diagnostics (after -wfn and\n"
 "                                    -ri_fit).\n"
 "  -wfn_cif                            Write a CIF from the wavefunction.\n"
 "  -lukas_test                          SALTED density diagnostic.\n"
 "  -calc_dens_1D [atom1 atom2 points padding]\n"
 "                                    One-dimensional density between atoms.\n"
 "  -spherical_harmonic  -spherical_atoms\n"
 "                                    Internal spherical-atom diagnostics.\n"
 "  -test                               Enable test mode.\n\n"
 "HELP\n"
 "  -h, --h, -help, --help              Show this message.\n\n"
 "EXAMPLES\n"
 "  Normal HAR\n"
 "    NoSpherA2 -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7\n"
 "  Independent-atom (Thakkar) TSC\n"
 "    NoSpherA2 -cif A.cif -hkl A.hkl -xyz A.xyz -IAM -acc 1 -cpus 7\n"
 "  Disorder / multi-TSC\n"
 "    NoSpherA2 -cif A.cif -hkl A.hkl -acc 1 -cpus 7 \\\n"
 "      -mtc 1.wfn 0,1 2.wfn 0,2 3.wfn 0,3 \\\n"
 "      -mtc_charge 0 0 0 -mtc_mult 1 1 1 -mtc_ECP 0 0 0\n"
 "  Fragment HAR\n"
 "    NoSpherA2 -cif A.cif -hkl A.hkl -acc 1 -cpus 7 \\\n"
 "      -cmtc 1.wfn 1.cif 0 2.wfn 2.cif 0 3_1.wfn 3_1.cif 0,1 \\\n"
 "      3_2.wfn 3_2.cif 0,2\n"
 "  Merge tables\n"
 "    NoSpherA2 -merge A.tsc B.tsc C.tsc\n"
 "    NoSpherA2 -merge_nocheck A.tsc B.tsc C.tsc  # only identical HKLs\n"
 "  Convert, label, and GBW conversion\n"
 "    NoSpherA2 -tscb A.tscb\n"
 "    NoSpherA2 -tsc_labels A.tscb A.cif A.labels.tsc\n"
 "    NoSpherA2 -gbw2wfn -wfn A.gbw\n"
 "  Twin law\n"
 "    NoSpherA2 -cif A.cif -hkl A.hkl -wfn A.wfx -acc 1 -cpus 7 \\\n"
 "      -twin -1 0 0 0 -1 0 0 0 -1\n");
std::string NoSpherA2_message(bool no_date)
{
    std::string t = "    _   __     _____       __              ___   ___\n";
    t.append("   / | / /___ / ___/____  / /_  ___  _____/   | |__ \\\n");
    t.append("  /  |/ / __ \\\\__ \\/ __ \\/ __ \\/ _ \\/ ___/ /| | __/ /\n");
    t.append(" / /|  / /_/ /__/ / /_/ / / / /  __/ /  / ___ |/ __/\n");
    t.append("/_/ |_/\\____/____/ .___/_/ /_/\\___/_/  /_/  |_/____/\n");
    t.append("                /_/\n");
    t.append("This software is part of the cuQCT software suite developed by Florian Kleemiss.\n");
    t.append("Please give credit and cite corresponding pieces!\nThis Software is published with BSD-2 clause license.\n");
    if (!no_date)
    {
        t.append("List of contributors of pieces of code or functionality:\n");
        t.append("      Florian Kleemiss,\n");
        t.append("      Emmanuel Hupf,\n");
        t.append("      Alessandro Genoni,\n");
        t.append("      Lukas M. Seifert,\n");
        t.append("      Daniel Bruex,\n");
        t.append("      Marti Gimferrer,\n");
        t.append("      Anker Nielsen,\n");
        t.append("      Lucas Militao,\n");
        t.append("      and many more in communications or by feedback!\n");
        t.append("NoSpherA2 uses featomic, Metatensor, and the mdspan library, as well as OCC for the calculation of wavefunctions, when required.\n");
        t.append("The used packages are published under BSD-3 clause License or explicit consent for the use in this project was given.\n");
        t.append("Please see, respectively:\n");
        t.append("   https://github.com/Luthaf/featomic\n");
        t.append("   https://github.com/lab-cosmo/metatensor\n");
        t.append("   This software utilizes Intel(c) Math Kernel Library (oneMKL), version 2025.2.0.629, for optimized mathematical computations\n");
        t.append("OCC can be found at: https://github.com/peterspackman/occ\n");
        t.append("NoSpherA2 was published at  : Kleemiss et al. Chem. Sci., 2021, 12, 1675 - 1692.\n");
        t.append("Slater IAM was published at : Kleemiss et al. J. Appl. Cryst. 2024, 57, 161 - 174.\n");
        t.append("ECP correction functions at : Kleemiss et al. J. Appl. Cryst. 2025, 58, 374 - 382.\n");
        t.append("Aux basis /RI partitioning  : Seifert et al. Z. Krist. - Cryst. Mat. 2026, 10.1515/zkri-2026-0013.\n");
        t.append("TFVC partitioning at        : Gimferrer et al. TBA.\n");
        t.append("MBIS/EMBIS partitioning at  : Nielsen et al. TBA.\n");
    }
    return t;
}

std::string build_date = ("This Executable was built on: " + std::string(__DATE__) + " " + std::string(__TIME__) + "\n");

bool ensure_occ_data_path(const char *argv0)
{
#ifdef _WIN32
    char* occ_data_path_env = nullptr;
    size_t len = 0;
    errno_t err = _dupenv_s(&occ_data_path_env, &len, "OCC_DATA_PATH");

    if (err == 0 && occ_data_path_env != nullptr)
    {
        std::filesystem::path path(occ_data_path_env);
        free(occ_data_path_env);

        if (is_valid_occ_data_path(path))
            return true;
        else
            std::cerr << "OCC DATA PATH is invalid!" << std::endl;
    }
    else if (occ_data_path_env != nullptr)
    {
        free(occ_data_path_env);
    }
#else
    const char* tmp_occ_data_path_env = std::getenv("OCC_DATA_PATH");
    if (tmp_occ_data_path_env != nullptr)
    {
        std::string occ_data_path_env(tmp_occ_data_path_env);
        if (is_valid_occ_data_path(std::filesystem::path(occ_data_path_env)))
            return true;
        else
            std::cerr << "OCC DATA PATH is invalid!" << std::endl;
    }
#endif

    std::filesystem::path exe_dir = resolve_executable_directory(argv0);
    if (exe_dir.empty())
    {
        std::cerr << "OCC_DATA_PATH not set or invalid and executable directory could not be resolved." << std::endl;
        return false;
    }

    std::filesystem::path selected_path = choose_occ_data_path_from_exe_dir(exe_dir);
    if (selected_path.empty())
    {
        std::cerr << "OCC_DATA_PATH not set or invalid. No valid OCC data directory found near executable directory: "
            << exe_dir << std::endl;
        return false;
    }

    if (!set_occ_data_path(selected_path))
    {
        std::cerr << "Failed to set OCC_DATA_PATH to resolved directory: " << selected_path << std::endl;
        return false;
    }

    std::cerr << "OCC_DATA_PATH not set or invalid. Using: " << selected_path << std::endl;
    return true;
}

bool is_similar_rel(const double &first, const double &second, const double &tolerance)
{
    double diff = abs(first - second);
    if (diff > abs((first + second + 0.01) * tolerance / 2))
        return false;
    else
        return true;
};

bool is_similar(const double &first, const double &second, const double &tolerance)
{
    double diff = abs(first - second);
    if (diff > pow(10, tolerance))
        return false;
    else
        return true;
};

bool is_similar_abs(const double &first, const double &second, const double &tolerance)
{
    double diff = abs(first - second);
    if (diff > abs(tolerance))
        return false;
    else
        return true;
};

/*
cosinus_annaeherung::cosinus_annaeherung() : mSize(0), mBase_values(nullptr), mStepwidth(1.0) {
    resize(100);
}

void cosinus_annaeherung::resize(size_t size)
{
    mSize = size;
    if (mBase_values) delete[] mBase_values;
    mBase_values = new double[mSize + 1];
#pragma omp parallel for
    for (auto i = 0; i < mSize + 1; i++)  // Fuer einen Werte mehr die Stueststellen speichern
    {
        double y = cos((MPI2 * i) / mSize);
        //std::cout << "resize: i="<<i<<" y=" << y << endl;
        mBase_values[i] = y;
    }
    mStepwidth = MPI2 / size;
}

double cosinus_annaeherung::calculate_error_at(double x) const
{
    return cos(x) - get(x);
}
*/
void copy_file(std::filesystem::path &from, std::filesystem::path &to)
{
    std::ifstream source(from.c_str(), std::ios::binary);
    std::ofstream dest(to.c_str(), std::ios::binary);

    dest << source.rdbuf();

    source.close();
    dest.close();
};

d3 vec_diff(const d3 &a, const d3 &b) {
    return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
};

d3 vec_cross(const d3 &a, const d3 &b)
{
    return {
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0] };
}

double vec_dot(const d3 &a, const d3 &b)
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};

//---------------------------Configuration files ---------------------------------------------------

std::filesystem::path get_home_path(void)
{
#ifdef _WIN32
    char *homeDrive = nullptr;
    size_t len = 0;
    std::string temp1, temp2;

    errno_t err = _dupenv_s(&homeDrive, &len, "HOMEDRIVE");
    if (err == 0 && homeDrive != nullptr) {
        temp1 = homeDrive;
        // Free the allocated memory
        free(homeDrive);
    }
    else {
        std::cerr << "Failed to retrieve the environment variable." << std::endl;
    }
    err = _dupenv_s(&homeDrive, &len, "HOMEPATH");
    if (err == 0 && homeDrive != nullptr) {
        temp2 = homeDrive;
        // Free the allocated memory
        free(homeDrive);
    }
    else {
        std::cerr << "Failed to retrieve the environment variable." << std::endl;
    }
    temp1.append(temp2);
    return temp1;
#else
    const char *home_env = getenv("HOME");
    if (home_env == nullptr) {
        std::cerr << "Warning: HOME environment variable not set." << std::endl;
        return std::filesystem::path("/tmp"); // Fallback to /tmp
    }

    std::string home = home_env;
    // Basic validation: check if it's a valid path and not empty
    if (home.empty() || home.find_first_of('\0') != std::string::npos) {
        std::cerr << "Warning: Invalid HOME environment variable." << std::endl;
        return std::filesystem::path("/tmp"); // Fallback to /tmp
    }

    return std::filesystem::path(home);
#endif
}

void write_spherical_atoms() {
#pragma omp parallel for
    for (int i = 1; i < 103; i++) {
        std::cout << "Atom " << i << std::endl;
        Thakkar atom(i);
        std::ofstream out("spherical_" + std::string(constants::atnr2letter(i)) + ".txt");
        out << std::scientific << std::setprecision(10);
        out << std::to_string(i) << " r Density" << std::endl;
        const double min_dist = 1E-7;
        const double incr = 1.025;

        // Make radial grids
        double current = 1;
        double _dist = min_dist;
        while (current > 1E-15)
        {
            current = atom.get_radial_density(_dist);
            out << _dist << " " << current << std::endl;
            _dist *= incr;
        }
    }
}

bool check_bohr(const WFN &wave, bool debug)
{
    double min_length = 300.0;
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        d3 atom1 = wave.get_atom_pos(i);
        for (int j = i + 1; j < wave.get_ncen(); j++)
        {
            d3 atom2 = wave.get_atom_pos(j);
            const double length = array_length(atom1, atom2);
            if (debug)
                std::cout << "Length for: " << i << ";" << j << ": " << length << ", min_length: " << min_length << std::endl;
            if (length < min_length)
                min_length = length;
        }
    }
    if (debug)
    {
        if (min_length < 2)
            std::cout << "Decided it's written in Angstrom" << std::endl;
        else
            std::cout << "Decided it's written in Bohr" << std::endl;
    }
    return (!(min_length < 2));
};

std::string go_get_string(std::ifstream &file, std::string search, bool rewind)
{
    if (rewind)
    {
        file.clear();
        file.seekg(0, file.beg);
    }
    std::string line;
    while (line.find(search) == std::string::npos && !file.eof() && getline(file, line))
        continue;
    if (file.eof())
        return "";
    else
        return line;
}

std::string shrink_string(std::string &input)
{
    while (input.find(" ") != -1)
    {
        input.erase(input.find(" "), 1);
    }
    while (input.find("1") != -1)
    {
        input.erase(input.find("1"), 1);
    }
    while (input.find("2") != -1)
    {
        input.erase(input.find("2"), 1);
    }
    while (input.find("3") != -1)
    {
        input.erase(input.find("3"), 1);
    }
    while (input.find("4") != -1)
    {
        input.erase(input.find("4"), 1);
    }
    while (input.find("5") != -1)
    {
        input.erase(input.find("5"), 1);
    }
    while (input.find("6") != -1)
    {
        input.erase(input.find("6"), 1);
    }
    while (input.find("7") != -1)
    {
        input.erase(input.find("7"), 1);
    }
    while (input.find("8") != -1)
    {
        input.erase(input.find("8"), 1);
    }
    while (input.find("9") != -1)
    {
        input.erase(input.find("9"), 1);
    }
    while (input.find("0") != -1)
    {
        input.erase(input.find("0"), 1);
    }
    while (input.find("(") != -1)
    {
        input.erase(input.find("("), 1);
    }
    while (input.find(")") != -1)
    {
        input.erase(input.find(")"), 1);
    }
    return input;
};

std::string shrink_string_to_atom(std::string &input, const int &atom_number)
{
    while (input.find(" ") != -1)
    {
        input.erase(input.find(" "), 1);
    }
    while (input.find("1") != -1)
    {
        input.erase(input.find("1"), 1);
    }
    while (input.find("2") != -1)
    {
        input.erase(input.find("2"), 1);
    }
    while (input.find("3") != -1)
    {
        input.erase(input.find("3"), 1);
    }
    while (input.find("4") != -1)
    {
        input.erase(input.find("4"), 1);
    }
    while (input.find("5") != -1)
    {
        input.erase(input.find("5"), 1);
    }
    while (input.find("6") != -1)
    {
        input.erase(input.find("6"), 1);
    }
    while (input.find("7") != -1)
    {
        input.erase(input.find("7"), 1);
    }
    while (input.find("8") != -1)
    {
        input.erase(input.find("8"), 1);
    }
    while (input.find("9") != -1)
    {
        input.erase(input.find("9"), 1);
    }
    while (input.find("0") != -1)
    {
        input.erase(input.find("0"), 1);
    }
    while (input.find("(") != -1)
    {
        input.erase(input.find("("), 1);
    }
    while (input.find(")") != -1)
    {
        input.erase(input.find(")"), 1);
    }
    std::string temp = constants::atnr2letter(atom_number);
    err_checkf(temp != "PROBLEM", "Problem identifying atoms!", std::cout);
    if (input.find(temp) != 1)
        return temp;
    if (temp != "PROBLEM")
        while (input.size() > temp.size())
            input.pop_back();
    return input;
};

bool read_block_from_fortran_binary(std::ifstream &file, void *Target)
{
    int size_begin = 0, size_end = 0;
    file.read(reinterpret_cast<char *>(&size_begin), sizeof(int));
    file.read(reinterpret_cast<char *>(Target), size_begin);
    file.read(reinterpret_cast<char *>(&size_end), sizeof(int));
    if (size_begin != size_end)
    {
        std::cout << "Error reading block from binary file: " << size_begin << " vs. " << size_end << std::endl;
        return false;
    }
    return true;
}
template<typename T>
bool read_block_from_fortran_binary(std::ifstream &file, std::vector<T> &Target)
{
    int size_begin = 0, size_end = 0;
    file.read(reinterpret_cast<char *>(&size_begin), sizeof(int));
    Target.resize(size_begin / sizeof(T));
    file.read(reinterpret_cast<char *>(Target.data()), size_begin);
    file.read(reinterpret_cast<char *>(&size_end), sizeof(int));
    if (size_begin != size_end)
    {
        std::cout << "Error reading block from binary file: " << size_begin << " vs. " << size_end << std::endl;
        return false;
    }
    return true;
}
template bool read_block_from_fortran_binary(std::ifstream &file, std::vector<double> &Target);
template bool read_block_from_fortran_binary(std::ifstream &file, std::vector<int> &Target);
template bool read_block_from_fortran_binary(std::ifstream &file, std::vector<float> &Target);
template bool read_block_from_fortran_binary(std::ifstream &file, std::vector<char> &Target);

primitive::primitive(int c, int t, double e, double coef) : center(c), type(t), exp(e), coefficient(coef)
{
    exp_l_plus_3_2 = pow(exp, type + 1.5);
    norm_const = pow(
        pow(2, 7 + 4 * type) * pow(exp, 3 + 2 * type) / constants::PI / pow(doublefactorial(2 * type + 1), 2),
        0.25);
    normalized_coefficient = coefficient * norm_const;
};

primitive::primitive(const SimplePrimitive &other) : center(other.center), type(other.type), exp(other.exp), coefficient(other.coefficient)
{
    exp_l_plus_3_2 = pow(other.exp, type + 1.5);
    norm_const = pow(
        pow(2, 7 + 4 * type) * pow(exp, 3 + 2 * type) / constants::PI / pow(doublefactorial(2 * type + 1), 2),
        0.25);
    normalized_coefficient = coefficient * norm_const;
};

void select_cubes(std::vector<std::vector<unsigned int>> &selection, std::vector<WFN> &wavy, unsigned int nr_of_cubes, bool wfnonly, bool debug)
{
    // asks which wfn to use, if wfnonly is set or whcih cubes up to nr of cubes to use
    // Returns values in selection[0][i] for iths selection of wavefunction and
    //  selection[1][i] for iths selection of cube
    using namespace std;
    std::cout << "Which of the following cubes to use? Need to select " << nr_of_cubes << " file";
    if (nr_of_cubes > 1)
        std::cout << "s in total." << endl;
    else
        std::cout << "." << endl;
    std::cout << endl
        << endl;
    for (int w = 0; w < wavy.size(); w++)
    {
        stringstream stream;
        std::cout << "_____________________________________________________________" << endl;
        std::cout << "WFN ";
        stream << setw(2) << w;
        std::cout << stream.str() << ") " << wavy[w].get_path().stem() << endl;
        stream.str("");
        for (int c = 0; c < wavy[w].get_cube_count(); c++)
        {
            if (c == 0)
                std::cout << "        |" << endl
                << "Cube    |" << endl;
            else
                std::cout << "        |" << endl;
            if (!wfnonly)
            {
                std::cout << setw(2) << w;
                std::cout << ".";
                std::cout << setw(2) << c;
            }
            else
                std::cout << "     ";
            std::cout << "   |_ " << wavy[w].get_cube_path(c).stem();
            if (!exists(wavy[w].get_cube_path(c)))
                std::cout << " (MEM ONLY)";
            std::cout << endl;
        }
        std::cout << "_____________________________________________________________" << endl
            << endl
            << endl;
    }
    // bool happy = false;
    unsigned int selected_cubes = 0;
    do
    {
        std::cout << "Select " << selected_cubes + 1 << ". ";
        if (wfnonly)
            std::cout << "WFN ";
        else
            std::cout << "cube ";
        std::cout << "please: ";
        string input;
        cin >> input;
        if (!wfnonly)
        {
            if (input.find('.') == string::npos)
            {
                std::cout << "no . found in input!" << endl;
                continue;
            }
        }
        else
        {
            if (input.find('.') == string::npos)
                std::cout << "Ignoring the .!" << endl;
            unsigned int nr_wave = fromString<unsigned int>(input);
            if (nr_wave < 0 || nr_wave >= wavy.size())
            {
                std::cout << "Invalid choice!" << endl;
                continue;
            }
            selected_cubes++;
            selection[0].push_back(nr_wave);
            if (selected_cubes == nr_of_cubes)
                return;
            else
                continue;
        }
        if (debug)
        {
            std::cout << "input: " << input << endl;
            std::cout << "with . found at: " << input.find('.') << endl;
            std::cout << "substr1: " << input.substr(0, input.find('.')) << endl;
            std::cout << "substr2: " << input.substr(input.find('.') + 1) << endl;
        }
        string wave(input.substr(0, input.find('.')));
        string cube(input.substr(input.find('.') + 1));
        unsigned int nr_wave = fromString<unsigned int>(wave);
        int nr_cube = fromString<int>(cube);
        if (debug)
            std::cout << "Translated: " << nr_wave << " " << nr_cube << endl;
        if (nr_wave < 0 || nr_wave >= wavy.size() || nr_cube < 0 || nr_cube >= wavy[nr_wave].get_cube_count())
        {
            std::cout << "Invalid choice!" << endl;
            continue;
        }
        selection[0][selected_cubes] = nr_wave;
        selection[1][selected_cubes] = nr_cube;
        selected_cubes++;
        if (selected_cubes == nr_of_cubes)
        {
            if (debug)
                std::cout << "Going to return!" << endl;
            return;
        }
    } while (true);
};

bool unsaved_files(std::vector<WFN> &wavy)
{
    for (int w = 0; w < wavy.size(); w++)
        for (int c = 0; c < wavy[w].get_cube_count(); c++)
            if (!exists(wavy[w].get_cube_path(c)))
                return true;
    return false;
};

void readxyzMinMax_fromWFN(
    const WFN &wavy,
    properties_options &opts,
    bool no_bohr)
{
    vec2 PosAtoms;
    PosAtoms.resize(3);
    for (int i = 0; i < 3; i++)
        PosAtoms[i].resize(wavy.get_ncen());
    bool bohrang = true;
    if (!no_bohr)
        bohrang = !check_bohr(wavy, false);

    for (int j = 0; j < wavy.get_ncen(); j++)
    {
        PosAtoms[0][j] = wavy.get_atom_coordinate(j, 0);
        PosAtoms[1][j] = wavy.get_atom_coordinate(j, 1);
        PosAtoms[2][j] = wavy.get_atom_coordinate(j, 2);
        if (!bohrang)
        {
            for (int i = 0; i < 3; i++)
                PosAtoms[i][j] = constants::ang2bohr(PosAtoms[i][j]);
        }
    }
    opts.MinMax[0] = *std::min_element(PosAtoms[0].begin(), PosAtoms[0].end());
    opts.MinMax[3] = *std::max_element(PosAtoms[0].begin(), PosAtoms[0].end());
    opts.MinMax[1] = *std::min_element(PosAtoms[1].begin(), PosAtoms[1].end());
    opts.MinMax[4] = *std::max_element(PosAtoms[1].begin(), PosAtoms[1].end());
    opts.MinMax[2] = *std::min_element(PosAtoms[2].begin(), PosAtoms[2].end());
    opts.MinMax[5] = *std::max_element(PosAtoms[2].begin(), PosAtoms[2].end());

    const double temp_rad = constants::ang2bohr(opts.radius);
    opts.MinMax[0] -= temp_rad;
    opts.MinMax[3] += temp_rad;
    opts.MinMax[1] -= temp_rad;
    opts.MinMax[4] += temp_rad;
    opts.MinMax[2] -= temp_rad;
    opts.MinMax[5] += temp_rad;

    opts.NbSteps[0] = (int)ceil(constants::bohr2ang(opts.MinMax[3] - opts.MinMax[0]) / opts.resolution);
    opts.NbSteps[1] = (int)ceil(constants::bohr2ang(opts.MinMax[4] - opts.MinMax[1]) / opts.resolution);
    opts.NbSteps[2] = (int)ceil(constants::bohr2ang(opts.MinMax[5] - opts.MinMax[2]) / opts.resolution);
}

void readxyzMinMax_fromCIF(
    std::filesystem::path cif,
    properties_options &opts,
    vec2 &cm)
{
    using namespace std;
    cell cell(cif);

    cm[0][0] = constants::ang2bohr(cell.get_a());
    cm[0][1] = constants::ang2bohr(cell.get_b() * cell.get_cg());
    cm[0][2] = constants::ang2bohr(cell.get_c() * cell.get_cb());
    cm[1][1] = constants::ang2bohr(cell.get_b() * cell.get_sg());
    cm[1][2] = constants::ang2bohr(cell.get_c() * (cell.get_ca() - cell.get_cb() * cell.get_cg()) / cell.get_sg());
    cm[2][2] = constants::ang2bohr(cell.get_V() / (cell.get_a() * cell.get_b() * cell.get_sg()));

    opts.MinMax[0] = 0.0;
    opts.MinMax[1] = 0.0;
    opts.MinMax[2] = 0.0;

    opts.MinMax[3] = (cell.get_a() + cell.get_b() * cell.get_cg() + cell.get_c() * cell.get_cb()) / 0.529177249;
    opts.MinMax[4] = (cell.get_b() * cell.get_sg() + cell.get_c() * (cell.get_ca() - cell.get_cb() * cell.get_cg()) / cell.get_sg()) / 0.529177249;
    opts.MinMax[5] = cm[2][2];

    opts.NbSteps[0] = (int)ceil(cell.get_a() / opts.resolution);
    opts.NbSteps[1] = (int)ceil(cell.get_b() / opts.resolution);
    opts.NbSteps[2] = (int)ceil(cell.get_c() / opts.resolution);

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            cm[i][j] /= opts.NbSteps[j];
}

bool generate_sph2cart_mat(vec2 &p, vec2 &d, vec2 &f, vec2 &g)
{
    //
    // From 3P: P0 P1 P2
    // To 3P : Z X Y (4 2 3, as in ORCA format)
    //
    p.resize(3);
    for (int i = 0; i < 3; i++)
    {
        p[i].resize(3, 0.0);
    }
    p[0][2] = 1.0;
    p[1][0] = 1.0;
    p[2][1] = 1.0;

    //
    // From 5D: D 0, D + 1, D - 1, D + 2, D - 2
    // To 6D : 5  6  7  8  9 10
    // XX, YY, ZZ, XY, XZ, YZ
    //
    d.resize(6);
    for (int i = 0; i < 6; i++)
    {
        d[i].resize(5, 0.0);
    }
    // XX = -0.5/SQRT(3) * D0 + 0.5 * D2
    d[0][0] = -0.5 / sqrt(3);
    d[0][3] = 0.5;
    // YY = -0.5/SQRT(3) * D0 - 0.5 * D2
    d[1][0] = -0.5 / sqrt(3);
    d[1][3] = -0.5;
    // ZZ = SQRT(1/3) * D0
    d[2][0] = sqrt(1.0 / 3.0);
    // XY = D-2
    d[3][4] = 1.0;
    // XZ = D1
    d[4][1] = 1.0;
    // YZ = D-1
    d[5][2] = 1.0;

    // From 7F: F 0, F + 1, F - 1, F + 2, F - 2, F + 3, F - 3
    // To 10F : 11   12   13   14   15   16   17   18   19  20
    // XXX, YYY, ZZZ, XXY, XXZ, YYZ, XYY, XZZ, YZZ, XYZ (AIMALL order!)
    //
    f.resize(10);
    for (int i = 0; i < 10; i++)
    {
        f[i].resize(7, 0.0);
    }
    f[0][1] = -sqrt(0.025);
    f[0][5] = -sqrt(1.0 / 24.0);

    f[1][2] = -sqrt(0.025);
    f[1][6] = sqrt(1.0 / 24.0);

    f[2][0] = sqrt(1.0 / 15.0);

    f[3][2] = -sqrt(0.025);
    f[3][6] = -sqrt(0.375);

    f[4][0] = -sqrt(0.15);
    f[4][3] = 0.5;

    f[5][0] = -sqrt(0.15);
    f[5][3] = -0.5;

    f[6][1] = -sqrt(0.025);
    f[6][5] = sqrt(0.375);

    f[7][1] = sqrt(0.4);

    f[8][2] = sqrt(0.4);

    f[9][4] = 1.0;

    g.resize(15);
    for (int i = 0; i < 15; i++)
    {
        g[i].resize(9, 0.0);
    }
    g[0][0] = 0.375 * sqrt(1.0 / 35.0);
    g[0][3] = -0.25 / sqrt(7);
    g[0][7] = -0.125;

    g[1][0] = g[0][0];
    g[1][3] = -g[0][3];
    g[1][7] = g[0][7];

    g[2][0] = sqrt(1.0 / 35.0);

    g[3][4] = -sqrt(1.0 / 28.0);
    g[3][8] = -0.5;

    g[4][1] = -0.5 / sqrt(98.0 / 63.0);
    g[4][5] = -1.0 / sqrt(8.0);

    g[5][4] = g[3][4];
    g[5][8] = -g[3][8];

    g[6][2] = g[4][1];
    g[6][6] = -g[4][5];

    g[7][1] = sqrt(2.0 / 7.0);

    g[8][2] = g[7][1];

    g[9][0] = 3.0 / sqrt(560);
    g[9][7] = 0.75;

    g[10][0] = -3.0 / sqrt(35);
    g[10][3] = 1.5 / sqrt(7);

    g[11][0] = g[10][0];
    g[11][3] = -g[10][3];

    g[12][2] = g[4][1];
    g[12][6] = -0.75 * sqrt(2);

    g[13][1] = g[4][1];
    g[13][5] = -g[12][6];

    g[14][4] = 3.0 / sqrt(7);
    return true;
}
bool generate_cart2sph_mat(vec2 &d, vec2 &f, vec2 &g, vec2 &h)
{
    //
    // From 5D: D 0, D + 1, D - 1, D + 2, D - 2
    // To 6D : 1  2  3  4  5  6
    // XX, YY, ZZ, XY, XZ, YZ
    //
    d.resize(6);
    for (int i = 0; i < 6; i++)
    {
        d[i].resize(5, 0.0);
    }
    // D0 = -0.5 * XX - 0.5 * YY + ZZ
    d[0][0] = -0.5;
    d[1][0] = -0.5;
    d[2][0] = 1.0;
    // D + 1 = XZ
    d[4][1] = 1.0;
    // D - 1 = YZ
    d[5][2] = 1.0;
    // D + 2 = SQRT(3) / 2 * (XX - YY)
    d[0][3] = sqrt(3.0) / 2.0;
    d[1][3] = -sqrt(3.0) / 2.0;
    // D - 2 = XY
    d[3][4] = 1.0;

    // From 7F: F 0, F + 1, F - 1, F + 2, F - 2, F + 3, F - 3
    // To 10F : 1   2   3   4   5   6   7   8   9  10
    // XXX, YYY, ZZZ, XYY, XXY, XXZ, XZZ, YZZ, YYZ, XYZ(Gaussian sequence, not identical to Multiwfn)
    //
    f.resize(10);
    for (int i = 0; i < 10; i++)
    {
        f[i].resize(7, 0.0);
    }
    // F 0 = -3 / (2 * sqrt5) * (XXZ + YYZ) + ZZZ
    f[2][0] = 1.0;
    f[5][0] = -1.5 / sqrt(5.0);
    f[8][0] = -1.5 / sqrt(5.0);
    // F + 1 = -sqrt(3 / 8) * XXX - sqrt(3 / 40) * XYY + sqrt(6 / 5) * XZZ
    f[0][1] = -sqrt(3.0 / 8.0);
    f[3][1] = -sqrt(3.0 / 40.0);
    f[6][1] = sqrt(6.0 / 5.0);
    // F - 1 = -sqrt(3 / 40) * XXY - sqrt(3 / 8) * YYY + sqrt(6 / 5) * YZZ
    f[1][2] = -sqrt(3.0 / 8.0);
    f[4][2] = -sqrt(3.0 / 40.0);
    f[7][2] = sqrt(6.0 / 5.0);
    // F + 2 = sqrt3 / 2 * (XXZ - YYZ)
    f[5][3] = sqrt(3.0) / 2.0;
    f[8][3] = -sqrt(3.0) / 2.0;
    // F - 2 = XYZ
    f[9][4] = 1.0;
    // F + 3 = sqrt(5 / 8) * XXX - 3 / sqrt8 * XYY
    f[0][5] = sqrt(5.0 / 8.0);
    f[3][5] = -3.0 / sqrt(8.0);
    // F - 3 = 3 / sqrt8 * XXY - sqrt(5 / 8) * YYY
    f[1][6] = -sqrt(5.0 / 8.0);
    f[4][6] = 3.0 / sqrt(8.0);

    // From 9G: G 0, G + 1, G - 1, G + 2, G - 2, G + 3, G - 3, G + 4, G - 4
    // To 15G : 1    2    3    4    5    6    7    8
    // ZZZZ, YZZZ, YYZZ, YYYZ, YYYY, XZZZ, XYZZ, XYYZ
    // 9   10   11   12   13   14   15
    // XYYY, XXZZ, XXYZ, XXYY, XXXZ, XXXY, XXXX
    //
    g.resize(15);
    for (int i = 0; i < 15; i++)
    {
        g[i].resize(9, 0.0);
    }
    // G 0 = ZZZZ + 3 / 8 * (XXXX + YYYY) - 3 * sqrt(3 / 35) * (XXZZ + YYZZ - 1 / 4 * XXYY)
    g[0][0] = 1.0;
    g[2][0] = -3.0 * sqrt(3.0 / 35.0);
    g[4][0] = 3.0 / 8.0;
    g[9][0] = -3.0 * sqrt(3.0 / 35.0);
    g[11][0] = 3.0 / 4.0 * sqrt(3.0 / 35.0);
    g[14][0] = 3.0 / 8.0;
    // G + 1 = 2 * sqrt(5 / 14) * XZZZ - 3 / 2 * sqrt(5 / 14) * XXXZ - 3 / 2 / sqrt14 * XYYZ
    g[5][1] = 2.0 * sqrt(5.0 / 14.0);
    g[7][1] = -1.5 / sqrt(14.0);
    g[12][1] = -1.5 * sqrt(5.0 / 14.0);
    // G - 1 = 2 * sqrt(5 / 14) * YZZZ - 3 / 2 * sqrt(5 / 14) * YYYZ - 3 / 2 / sqrt14 * XXYZ
    g[1][2] = 2.0 * sqrt(5.0 / 14.0);
    g[3][2] = -1.5 * sqrt(5.0 / 14.0);
    g[10][2] = -1.5 / sqrt(14.0);
    // G + 2 = 3 * sqrt(3 / 28) * (XXZZ - YYZZ) - sqrt5 / 4 * (XXXX - YYYY)
    g[2][3] = -3.0 * sqrt(3.0 / 28.0);
    g[4][3] = sqrt(5.0) / 4.0;
    g[9][3] = 3.0 * sqrt(3.0 / 28.0);
    g[14][3] = -sqrt(5.0) / 4.0;
    // G - 2 = 3 / sqrt7 * XYZZ - sqrt(5 / 28) * (XXXY + XYYY)
    g[6][4] = 3.0 / sqrt(7.0);
    g[8][4] = -sqrt(5.0 / 28.0);
    g[13][4] = -sqrt(5.0 / 28.0);
    // G + 3 = sqrt(5 / 8) * XXXZ - 3 / sqrt8 * XYYZ
    g[7][5] = -3.0 / sqrt(8.0);
    g[12][5] = sqrt(5.0 / 8.0);
    // G - 3 = -sqrt(5 / 8) * YYYZ + 3 / sqrt8 * XXYZ
    g[3][6] = -sqrt(5.0 / 8.0);
    g[10][6] = 3.0 / sqrt(8.0);
    // G + 4 = sqrt35 / 8 * (XXXX + YYYY) - 3 / 4 * sqrt3 * XXYY
    g[4][7] = sqrt(35.0) / 8.0;
    g[11][7] = -3.0 / 4.0 * sqrt(3.0);
    g[14][7] = sqrt(35.0) / 8.0;
    // G - 4 = sqrt5 / 2 * (XXXY - XYYY)
    g[8][8] = -sqrt(5.0) / 2.0;
    g[13][8] = sqrt(5.0) / 2.0;

    // From 11H: H 0, H + 1, H - 1, H + 2, H - 2, H + 3, H - 3, H + 4, H - 4, H + 5, H - 5
    // To 21H : 1     2     3     4     5     6     7     8     9    10
    // ZZZZZ YZZZZ YYZZZ YYYZZ YYYYZ YYYYY XZZZZ XYZZZ XYYZZ XYYYZ
    // 11    12    13    14    15    16    17    18    19    20    21
    // XYYYY XXZZZ XXYZZ XXYYZ XXYYY XXXZZ XXXYZ XXXYY XXXXZ XXXXY XXXXX
    //
    h.resize(21);
    for (int i = 0; i < 21; i++)
    {
        h[i].resize(11);
        std::fill(h[i].begin(), h[i].end(), 0.0);
    }
    // H 0 = ZZZZZ - 5 / sqrt21 * (XXZZZ + YYZZZ) + 5 / 8 * (XXXXZ + YYYYZ) + sqrt(15 / 7) / 4 * XXYYZ
    h[0][0] = 1.0;
    h[11][0] = -5.0 / sqrt(21.0);
    h[2][0] = -5.0 / sqrt(21.0);
    h[18][0] = 5.0 / 8.0;
    h[4][0] = 5.0 / 8.0;
    h[13][0] = sqrt(15.0 / 7.0) / 4.0;
    // H + 1 = sqrt(5 / 3) * XZZZZ - 3 * sqrt(5 / 28) * XXXZZ - 3 / sqrt28 * XYYZZ + sqrt15 / 8 * XXXXX + sqrt(5 / 3) / 8 * XYYYY + sqrt(5 / 7) / 4 * XXXYY
    h[6][1] = sqrt(5.0 / 3.0);
    h[15][1] = -3.0 * sqrt(5.0 / 28.0);
    h[8][1] = -3.0 / sqrt(28.0);
    h[20][1] = sqrt(15.0) / 8.0;
    h[10][1] = sqrt(5.0 / 3.0) / 8.0;
    h[17][1] = sqrt(5.0 / 7.0) / 4.0;
    // H - 1 = sqrt(5 / 3) * YZZZZ - 3 * sqrt(5 / 28) * YYYZZ - 3 / sqrt28 * XXYZZ + sqrt15 / 8 * YYYYY + sqrt(5 / 3) / 8 * XXXXY + sqrt(5 / 7) / 4 * XXYYY
    h[1][2] = sqrt(5.0 / 3.0);
    h[3][2] = -3.0 * sqrt(5.0 / 28.0);
    h[12][2] = -3.0 / sqrt(28.0);
    h[5][2] = sqrt(15.0) / 8.0;
    h[19][2] = sqrt(5.0 / 3.0) / 8.0;
    h[14][2] = sqrt(5.0 / 7.0) / 4.0;
    // H + 2 = sqrt5 / 2 * (XXZZZ - YYZZZ) - sqrt(35 / 3) / 4 * (XXXXZ - YYYYZ)
    h[11][3] = sqrt(5.0) / 2.0;
    h[2][3] = -sqrt(5.0) / 2.0;
    h[18][3] = -sqrt(35.0 / 3.0) / 4.0;
    h[4][3] = sqrt(35.0 / 3.0) / 4.0;
    // H - 2 = sqrt(5 / 3) * XYZZZ - sqrt(5 / 12) * (XXXYZ + XYYYZ)
    h[7][4] = sqrt(5.0 / 3.0);
    h[16][4] = -sqrt(5.0 / 12.0);
    h[9][4] = -sqrt(5.0 / 12.0);
    // H + 3 = sqrt(5 / 6) * XXXZZ - sqrt(3 / 2) * XYYZZ - sqrt(35 / 2) / 8 * (XXXXX - XYYYY) + sqrt(5 / 6) / 4 * XXXYY
    h[15][5] = sqrt(5.0 / 6.0);
    h[8][5] = -sqrt(1.5);
    h[20][5] = -sqrt(17.5) / 8.0;
    h[10][5] = sqrt(17.5) / 8.0;
    h[17][5] = sqrt(5.0 / 6.0) / 4.0;
    // H - 3 = -sqrt(5 / 6) * YYYZZ + sqrt(3 / 2) * XXYZZ - sqrt(35 / 2) / 8 * (XXXXY - YYYYY) - sqrt(5 / 6) / 4 * XXYYY
    h[3][6] = -sqrt(5.0 / 6.0);
    h[12][6] = sqrt(1.5);
    h[19][6] = -sqrt(17.5) / 8.0;
    h[5][6] = sqrt(17.5) / 8.0;
    h[14][6] = -sqrt(5.0 / 6.0) / 4.0;
    // H + 4 = sqrt35 / 8 * (XXXXZ + YYYYZ) - 3 / 4 * sqrt3 * XXYYZ
    h[18][7] = sqrt(35.0) / 8.0;
    h[4][7] = sqrt(35.0) / 8.0;
    h[13][7] = -0.75 * sqrt(3.0);
    // H - 4 = sqrt5 / 2 * (XXXYZ - XYYYZ)
    h[16][8] = sqrt(5.0) / 2.0;
    h[9][8] = -sqrt(5.0) / 2.0;
    // H + 5 = 3 / 8 * sqrt(7 / 2) * XXXXX + 5 / 8 * sqrt(7 / 2) * XYYYY - 5 / 4 * sqrt(3 / 2) * XXXYY
    h[20][9] = 3.0 / 8.0 * sqrt(3.5);
    h[10][9] = 5.0 / 8.0 * sqrt(3.5);
    h[17][9] = -1.25 * sqrt(1.5);
    // H - 5 = 3 / 8 * sqrt(7 / 2) * YYYYY + 5 / 8 * sqrt(7 / 2) * XXXXY - 5 / 4 * sqrt(3 / 2) * XXYYY
    h[5][10] = 3.0 / 8.0 * sqrt(3.5);
    h[19][10] = 5.0 / 8.0 * sqrt(3.5);
    h[14][10] = -1.25 * sqrt(1.5);
    return true;
}

bool read_fracs_ADPs_from_CIF(const std::filesystem::path& cif, WFN& wavy, cell& unit_cell, std::ofstream& log3, const bool& debug)
{
    using namespace std;
    vec2 Uij, Cijk, Dijkl;
    ifstream asym_cif_input(cif, std::ios::in);
    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    string line;
    svec labels;
    int count_fields = 0;
    int position_field[3] = { 0, 0, 0 };
    int label_field = 100;
    vec2 positions;
    positions.resize(wavy.get_ncen());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < wavy.get_ncen(); i++)
    {
        positions[i].resize(3);
    }
    bool atoms_read = false;
    while (!asym_cif_input.eof() && !atoms_read)
    {
        getline(asym_cif_input, line);
        if (line.find("loop_") != string::npos)
        {
            while (line.find("_") != string::npos)
            {
                getline(asym_cif_input, line);
                if (debug)
                    log3 << "line in loop field definition: " << line << endl;
                if (line.find("label") != string::npos)
                    label_field = count_fields;
                else if (line.find("fract_x") != string::npos)
                    position_field[0] = count_fields;
                else if (line.find("fract_y") != string::npos)
                    position_field[1] = count_fields;
                else if (line.find("fract_z") != string::npos)
                    position_field[2] = count_fields;
                else if (label_field == 100)
                {
                    if (debug)
                        log3 << "I don't think this is the atom block.. moving on!" << endl;
                    break;
                }
                count_fields++;
            }
            while (line.find("_") == string::npos && line.length() > 3)
            {
                atoms_read = true;
                stringstream s(line);
                svec fields;
                fields.resize(count_fields);
                for (int i = 0; i < count_fields; i++)
                    s >> fields[i];
                if (debug)
                    log3 << "label: " << fields[label_field] << " frac_position: " << stod(fields[position_field[0]]) << " " << stod(fields[position_field[1]]) << " " << stod(fields[position_field[2]]) << endl;
                positions[labels.size()] = unit_cell.get_coords_cartesian(stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]));
                bool found_this_one = false;
                if (debug)
                    log3 << "label: " << fields[label_field] << " cartesian position: " << positions[labels.size()][0] << " " << positions[labels.size()][1] << " " << positions[labels.size()][2] << endl;
                for (int i = 0; i < wavy.get_ncen(); i++)
                {
                    if (is_similar(positions[labels.size()][0], wavy.get_atom_coordinate(i, 0), -1) && is_similar(positions[labels.size()][1], wavy.get_atom_coordinate(i, 1), -1) && is_similar(positions[labels.size()][2], wavy.get_atom_coordinate(i, 2), -1))
                    {
                        if (debug)
                            log3 << "WFN position: " << wavy.get_atom_coordinate(i, 0) << " " << wavy.get_atom_coordinate(i, 1) << " " << wavy.get_atom_coordinate(i, 2) << endl
                            << "Found an atom: " << fields[label_field] << " Corresponding to atom charge " << wavy.get_atom_charge(i) << endl;
                        wavy.set_atom_label(i, fields[label_field]);
                        wavy.set_atom_frac_coords(i, { stod(fields[position_field[0]]), stod(fields[position_field[1]]), stod(fields[position_field[2]]) });
                        found_this_one = true;
                        break;
                    }
                }
                if (!found_this_one && debug)
                    log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
                labels.push_back(fields[label_field]);
                getline(asym_cif_input, line);
            }
        }
    }

    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    count_fields = 0;
    int ADP_field[15] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    label_field = 100;
    atoms_read = false;
    Uij.resize(wavy.get_ncen());
    while (!asym_cif_input.eof() && !atoms_read)
    {
        getline(asym_cif_input, line);
        if (line.find("loop_") != string::npos)
        {
            while (line.find("_") != string::npos)
            {
                getline(asym_cif_input, line);
                if (debug)
                    log3 << "line in loop field definition: " << line << endl;
                if (line.find("aniso_label") != string::npos)
                    label_field = count_fields;
                else if (line.find("aniso_U_11") != string::npos)
                    ADP_field[0] = count_fields;
                else if (line.find("aniso_U_22") != string::npos)
                    ADP_field[1] = count_fields;
                else if (line.find("aniso_U_33") != string::npos)
                    ADP_field[2] = count_fields;
                else if (line.find("aniso_U_12") != string::npos)
                    ADP_field[3] = count_fields;
                else if (line.find("aniso_U_13") != string::npos)
                    ADP_field[4] = count_fields;
                else if (line.find("aniso_U_23") != string::npos)
                    ADP_field[5] = count_fields;
                else if (label_field == 100)
                {
                    if (debug)
                        log3 << "I don't think this is the Uij block.. moving on!" << endl;
                    break;
                }
                count_fields++;
            }
            while (line.find("_") == string::npos && line.length() > 3)
            {
                atoms_read = true;
                stringstream s(line);
                svec fields;
                fields.resize(count_fields);
                for (int i = 0; i < count_fields; i++)
                    s >> fields[i];
                if (debug)
                    log3 << "label: " << fields[label_field] << endl;
                bool found_this_one = false;
                for (int i = 0; i < wavy.get_ncen(); i++)
                {
                    if (fields[label_field] == wavy.get_atom_label(i))
                    {
                        Uij[i].resize(6);
                        for (int j = 0; j < 6; j++)
                            Uij[i][j] = stod(fields[ADP_field[j]]);
                        found_this_one = true;
                        break;
                    }
                }
                if (!found_this_one && debug)
                    log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
                getline(asym_cif_input, line);
            }
        }
    }

    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    count_fields = 0;
    label_field = 100;
    atoms_read = false;
    Cijk.resize(wavy.get_ncen());
    while (!asym_cif_input.eof() && !atoms_read)
    {
        getline(asym_cif_input, line);
        if (line.find("loop_") != string::npos)
        {
            while (line.find("_") != string::npos)
            {
                getline(asym_cif_input, line);
                if (debug)
                    log3 << "line in loop field definition: " << line << endl;
                if (line.find("C_label") != string::npos)
                    label_field = count_fields;
                else if (line.find("C_111") != string::npos)
                    ADP_field[0] = count_fields;
                else if (line.find("C_112") != string::npos)
                    ADP_field[1] = count_fields;
                else if (line.find("C_113") != string::npos)
                    ADP_field[2] = count_fields;
                else if (line.find("C_122") != string::npos)
                    ADP_field[3] = count_fields;
                else if (line.find("C_123") != string::npos)
                    ADP_field[4] = count_fields;
                else if (line.find("C_133") != string::npos)
                    ADP_field[5] = count_fields;
                else if (line.find("C_222") != string::npos)
                    ADP_field[6] = count_fields;
                else if (line.find("C_223") != string::npos)
                    ADP_field[7] = count_fields;
                else if (line.find("C_233") != string::npos)
                    ADP_field[8] = count_fields;
                else if (line.find("C_333") != string::npos)
                    ADP_field[9] = count_fields;
                else if (label_field == 100)
                {
                    if (debug)
                        log3 << "I don't think this is the Cijk block.. moving on!" << endl;
                    break;
                }
                count_fields++;
            }
            while (line.find("_") == string::npos && line.length() > 3)
            {
                atoms_read = true;
                stringstream s(line);
                svec fields;
                fields.resize(count_fields);
                for (int i = 0; i < count_fields; i++)
                    s >> fields[i];
                if (debug)
                    log3 << "label: " << fields[label_field] << endl;
                bool found_this_one = false;
                for (int i = 0; i < wavy.get_ncen(); i++)
                {
                    if (fields[label_field] == wavy.get_atom_label(i))
                    {
                        Cijk[i].resize(10);
                        for (int j = 0; j < 6; j++)
                            Cijk[i][j] = stod(fields[ADP_field[j]]);
                        found_this_one = true;
                        break;
                    }
                }
                if (!found_this_one && debug)
                    log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
                getline(asym_cif_input, line);
            }
        }
    }

    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    count_fields = 0;
    label_field = 100;
    atoms_read = false;
    Dijkl.resize(wavy.get_ncen());
    while (!asym_cif_input.eof() && !atoms_read)
    {
        getline(asym_cif_input, line);
        if (line.find("loop_") != string::npos)
        {
            while (line.find("_") != string::npos)
            {
                getline(asym_cif_input, line);
                if (debug)
                    log3 << "line in loop field definition: " << line << endl;
                if (line.find("D_label") != string::npos)
                    label_field = count_fields;
                else if (line.find("D_1111") != string::npos)
                    ADP_field[0] = count_fields;
                else if (line.find("D_1112") != string::npos)
                    ADP_field[1] = count_fields;
                else if (line.find("D_1113") != string::npos)
                    ADP_field[2] = count_fields;
                else if (line.find("D_1122") != string::npos)
                    ADP_field[3] = count_fields;
                else if (line.find("D_1123") != string::npos)
                    ADP_field[4] = count_fields;
                else if (line.find("D_1133") != string::npos)
                    ADP_field[5] = count_fields;
                else if (line.find("D_1222") != string::npos)
                    ADP_field[6] = count_fields;
                else if (line.find("D_1223") != string::npos)
                    ADP_field[7] = count_fields;
                else if (line.find("D_1233") != string::npos)
                    ADP_field[8] = count_fields;
                else if (line.find("D_1333") != string::npos)
                    ADP_field[9] = count_fields;
                else if (line.find("D_2222") != string::npos)
                    ADP_field[10] = count_fields;
                else if (line.find("D_2223") != string::npos)
                    ADP_field[11] = count_fields;
                else if (line.find("D_2233") != string::npos)
                    ADP_field[12] = count_fields;
                else if (line.find("D_2333") != string::npos)
                    ADP_field[13] = count_fields;
                else if (line.find("D_3333") != string::npos)
                    ADP_field[14] = count_fields;
                else if (label_field == 100)
                {
                    if (debug)
                        log3 << "I don't think this is the Dijk block.. moving on!" << endl;
                    break;
                }
                count_fields++;
            }
            while (line.find("_") == string::npos && line.length() > 3)
            {
                atoms_read = true;
                stringstream s(line);
                svec fields;
                fields.resize(count_fields);
                for (int i = 0; i < count_fields; i++)
                    s >> fields[i];
                if (debug)
                    log3 << "label: " << fields[label_field] << endl;
                bool found_this_one = false;
                for (int i = 0; i < wavy.get_ncen(); i++)
                {
                    if (fields[label_field] == wavy.get_atom_label(i))
                    {
                        Dijkl[i].resize(15);
                        for (int j = 0; j < 6; j++)
                            Dijkl[i][j] = stod(fields[ADP_field[j]]);
                        found_this_one = true;
                        break;
                    }
                }
                if (!found_this_one && debug)
                    log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
                getline(asym_cif_input, line);
            }
        }
    }

    for (int i = 0; i < wavy.get_ncen(); i++)
        wavy.set_atom_ADPs(i, { Uij[i], Cijk[i], Dijkl[i] });

    return true;
};

bool read_fracs_ADPs_from_CIF(const std::filesystem::path &cif, WFN &wavy, std::ofstream &log3, const bool &debug, const bool &grown, const ivec3 &symmetry_linking_list)
{
    using namespace std;
    ifstream asym_cif_input(cif, std::ios::in);
    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    string line;
    int ncen;
    if (!grown) {
        ncen = wavy.get_ncen();
    }
    else {
        ncen = symmetry_linking_list.size();
    }
    svec labels(ncen);

    for (int i = 0; i < ncen; i++) {
        labels[i] = wavy.get_atoms()[i].get_label();
    }

    while (getline(asym_cif_input, line)) {
        if (!line.starts_with("loop_")) {
            if (debug)
                log3 << "This is not part of a loop. Moving on.";
            continue;
        }
        if (debug)
            log3 << "Found a loop!";
        getline(asym_cif_input, line);
        if (line.find("_atom_site_aniso_label") != string::npos) {
            if (debug) {
                log3 << "This loop contains anisotropic displacement parameters.";
            }
            ivec fields;
            while (line.find("_atom_site_aniso") != string::npos && line.length() > 3) {
                getline(asym_cif_input, line);
                if (line.find("U_11") != string::npos)
					fields.push_back(0);
				else if (line.find("U_22") != string::npos)
					fields.push_back(1);
				else if (line.find("U_33") != string::npos)
					fields.push_back(2);
				else if (line.find("U_12") != string::npos)
					fields.push_back(3);
				else if (line.find("U_13") != string::npos)
					fields.push_back(4);
				else if (line.find("U_23") != string::npos)
					fields.push_back(5);	
            }
            while (line.find_first_not_of(" \t\r\n") != std::string::npos) {
                std::vector<std::string> entries;
                std::istringstream iss(line);
                std::string token;
                while (entries.size() < 7 && iss >> token)
                    entries.push_back(token);
                while (entries.size() < 7 && std::getline(asym_cif_input, line)) {
                    std::istringstream nextLine(line);
                    while (entries.size() < 7 && nextLine >> token)
                        entries.push_back(token);
                }
                bool atom_found = false;
                for (int a = 0; a < ncen; a++) {
                    if (entries[0] == wavy.get_atom(a).get_label()) {
						vec2 ADPs = wavy.get_atom(a).get_ADPs();
                        if (ADPs.size() != 3)
                            ADPs.resize(3);
                        ADPs[0].resize(6);
                        for (int i = 0; i < 6; i++) {
							ADPs[0][fields[i]] = stof(entries[i + 1]);
                        }
                        wavy.set_atom_ADPs(a, ADPs);
                        if (grown) {
                            for (int b = 0; b < symmetry_linking_list[a].size(); b++) {
                                wavy.set_atom_ADPs(symmetry_linking_list[a][b][0], ADPs);
                            }
                        }
                        atom_found = true;
                        break;
                    }
                }
                if (!atom_found) throw std::runtime_error("Displacement parameters found for atom that is not recognized!");
                getline(asym_cif_input, line);
            }
        }
        else if (line.find("_atom_site_anharm_GC_C_label") != string::npos) {
			if (debug)
				log3 << "This loop contains anharmonic Gram-Charlier coefficients C.";
            ivec fields;
            while (line.find("_atom_site_anharm") != string::npos && line.length() > 3) {
                getline(asym_cif_input, line);
                if (line.find("C_111") != string::npos)
                    fields.push_back(0);
                else if (line.find("C_112") != string::npos)
                    fields.push_back(1);
                else if (line.find("C_113") != string::npos)
                    fields.push_back(2);
                else if (line.find("C_122") != string::npos)
                    fields.push_back(3);
                else if (line.find("C_123") != string::npos)
                    fields.push_back(4);
                else if (line.find("C_133") != string::npos)
                    fields.push_back(5);
                else if (line.find("C_222") != string::npos)
                    fields.push_back(6);
                else if (line.find("C_223") != string::npos)
                    fields.push_back(7);
                else if (line.find("C_233") != string::npos)
                    fields.push_back(8);
                else if (line.find("C_333") != string::npos)
                    fields.push_back(9);   
            }
            while (line.find_first_not_of(" \t\r\n") != std::string::npos) {
                std::vector<std::string> entries;
                std::istringstream iss(line);
                std::string token;
                while (entries.size() < 11 && iss >> token) {
                    const int pos = token.find('(');
                    entries.push_back(token);
                }
                while (entries.size() < 11 && std::getline(asym_cif_input, line)) {
                    std::istringstream nextLine(line);
                    while (entries.size() < 11 && nextLine >> token) {
                        entries.push_back(token);
                    }
                }
                bool atom_found = false;
                for (int a = 0; a < ncen; a++) {
                    if (entries[0] == wavy.get_atom(a).get_label()) {
                        vec2 ADPs = wavy.get_atom(a).get_ADPs();
                        if (ADPs.size() != 3)
                            ADPs.resize(3);
                        ADPs[1].resize(10);
                        for (int i = 0; i < 10; i++) {
                            ADPs[1][fields[i]] = stof(entries[i + 1]);
                        }
                        wavy.set_atom_ADPs(a, ADPs);
                        if (grown) {
                            for (int b = 0; b < symmetry_linking_list[a].size(); b++) {
                                wavy.set_atom_ADPs(symmetry_linking_list[a][b][0], ADPs);
                            }
                        }
                        atom_found = true;
                        break;
                    }
                }
                if (!atom_found) throw std::runtime_error("Displacement parameters found for atom that is not recognized!");
                getline(asym_cif_input, line);
            }
        }
        else if (line.find("_atom_site_anharm_GC_D_label") != string::npos) {
			if (debug)
				log3 << "This loop contains anharmonic Gram-Charlier coefficients D.";
            ivec fields;
            while (line.find("_atom_site_anharm") != string::npos && line.length() > 3) {
                getline(asym_cif_input, line);
                if (line.find("D_1111") != string::npos)
                    fields.push_back(0);
                else if (line.find("D_1112") != string::npos)
                    fields.push_back(1);
                else if (line.find("D_1113") != string::npos)
                    fields.push_back(2);
                else if (line.find("D_1122") != string::npos)
                    fields.push_back(3);
                else if (line.find("D_1123") != string::npos)
                    fields.push_back(4);
                else if (line.find("D_1133") != string::npos)
                    fields.push_back(5);
                else if (line.find("D_1222") != string::npos)
                    fields.push_back(6);
                else if (line.find("D_1223") != string::npos)
                    fields.push_back(7);
                else if (line.find("D_1233") != string::npos)
                    fields.push_back(8);
                else if (line.find("D_1333") != string::npos)
                    fields.push_back(9);
                else if (line.find("D_2222") != string::npos)
                    fields.push_back(10);
                else if (line.find("D_2223") != string::npos)
                    fields.push_back(11);
                else if (line.find("D_2233") != string::npos)
                    fields.push_back(12);
                else if (line.find("D_2333") != string::npos)
                    fields.push_back(13);
                else if (line.find("D_3333") != string::npos)
                    fields.push_back(14);
            }
            while (line.find_first_not_of(" \t\r\n") != std::string::npos) {
                std::vector<std::string> entries;
                std::istringstream iss(line);
                std::string token;
                while (entries.size() < 16 && iss >> token)
                    entries.push_back(token);
                while (entries.size() < 16 && std::getline(asym_cif_input, line)) {
                    std::istringstream nextLine(line);
                    while (entries.size() < 16 && nextLine >> token)
                        entries.push_back(token);
                }
                bool atom_found = false;
                for (int a = 0; a < ncen; a++) {
                    if (entries[0] == wavy.get_atom(a).get_label()) {
                        vec2 ADPs = wavy.get_atom(a).get_ADPs();
                        if (ADPs.size() != 3)
                            ADPs.resize(3);
                        ADPs[2].resize(15);
                        for (int i = 0; i < 15; i++) {
                            ADPs[2][fields[i]] = stof(entries[i + 1]);
                        }                
                        wavy.set_atom_ADPs(a, ADPs);
                        if (grown) {
                            for (int b = 0; b < symmetry_linking_list[a].size(); b++) {
                                wavy.set_atom_ADPs(symmetry_linking_list[a][b][0], ADPs);
                            }
                        }
                        atom_found = true;
                        break;
                    }
                }
                if (!atom_found) throw std::runtime_error("Displacement parameters found for atom that is not recognized!");
                getline(asym_cif_input, line);
            }
        }
        else {
            if (debug)
                log3 << "This was not the right loop. Moving on.";
            continue;
        }
    }
    return true;
    // closing function
};

vec read_U_iso_from_CIF(const std::filesystem::path &cif, WFN &wavy, cell &unit_cell, std::ofstream &log3, const bool& debug)
{
    using namespace std;
    vec U_iso;
    ifstream asym_cif_input(cif, std::ios::in);
    asym_cif_input.clear();
    asym_cif_input.seekg(0, asym_cif_input.beg);
    string line;
    svec labels;
    int count_fields = 0;
    int position_field[3] = { 0, 0, 0 };
    int label_field = 100;
    int U_iso_field = -1;          // <-- NEW: tracks column index of U_iso_or_equiv
    vec2 positions;
    positions.resize(wavy.get_ncen());

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < wavy.get_ncen(); i++)
        positions[i].resize(3);

    U_iso.resize(wavy.get_ncen(), 0.0f);  // <-- NEW: pre-fill with zeros

    bool atoms_read = false;
    while (!asym_cif_input.eof() && !atoms_read)
    {
        getline(asym_cif_input, line);
        if (line.find("loop_") != string::npos)
        {
            while (line.find("_") != string::npos)
            {
                getline(asym_cif_input, line);
                if (debug)
                    log3 << "line in loop field definition: " << line << endl;
                if (line.find("_atom_site_label") != string::npos          // be specific to avoid
                    && line.find("aniso") == string::npos)                 // matching aniso_label
                    label_field = count_fields;
                else if (line.find("fract_x") != string::npos)
                    position_field[0] = count_fields;
                else if (line.find("fract_y") != string::npos)
                    position_field[1] = count_fields;
                else if (line.find("fract_z") != string::npos)
                    position_field[2] = count_fields;
                else if (line.find("U_iso_or_equiv") != string::npos)     // <-- NEW
                    U_iso_field = count_fields;
                else if (label_field == 100)
                {
                    if (debug)
                        log3 << "I don't think this is the atom block.. moving on!" << endl;
                    break;
                }
                count_fields++;
            }
            while (line.find("_") == string::npos && line.length() > 3)
            {
                atoms_read = true;
                stringstream s(line);
                svec fields;
                fields.resize(count_fields);
                for (int i = 0; i < count_fields; i++)
                    s >> fields[i];
                if (debug)
                    log3 << "label: " << fields[label_field]
                    << " frac_position: " << stod(fields[position_field[0]])
                    << " " << stod(fields[position_field[1]])
                    << " " << stod(fields[position_field[2]]) << endl;
                positions[labels.size()] = unit_cell.get_coords_cartesian(
                    stod(fields[position_field[0]]),
                    stod(fields[position_field[1]]),
                    stod(fields[position_field[2]]));
                bool found_this_one = false;
                if (debug)
                    log3 << "label: " << fields[label_field]
                    << " cartesian position: " << positions[labels.size()][0]
                    << " " << positions[labels.size()][1]
                    << " " << positions[labels.size()][2] << endl;
                for (int i = 0; i < wavy.get_ncen(); i++)
                {
                    if (is_similar(positions[labels.size()][0], wavy.get_atom_coordinate(i, 0), -1)
                        && is_similar(positions[labels.size()][1], wavy.get_atom_coordinate(i, 1), -1)
                        && is_similar(positions[labels.size()][2], wavy.get_atom_coordinate(i, 2), -1))
                    {
                        if (debug)
                            log3 << "WFN position: "
                            << wavy.get_atom_coordinate(i, 0) << " "
                            << wavy.get_atom_coordinate(i, 1) << " "
                            << wavy.get_atom_coordinate(i, 2) << endl
                            << "Found an atom: " << fields[label_field]
                            << " Corresponding to atom charge "
                            << wavy.get_atom_charge(i) << endl;
                        wavy.set_atom_label(i, fields[label_field]);
                        wavy.set_atom_frac_coords(i, {
                            stod(fields[position_field[0]]),
                            stod(fields[position_field[1]]),
                            stod(fields[position_field[2]]) });

                        // NEW: read U_iso_or_equiv, guard against '.' / '?' placeholders
                        if (U_iso_field >= 0 && U_iso_field < (int)fields.size())
                        {
                            try { U_iso[i] = (float)stod(fields[U_iso_field]); }
                            catch (...) { U_iso[i] = 0.0f; }
                        }

                        found_this_one = true;
                        break;
                    }
                }
                if (!found_this_one && debug)
                    log3 << "I DID NOT FIND THIS ATOM IN THE CIF?! WTF?!" << endl;
                labels.push_back(fields[label_field]);
                getline(asym_cif_input, line);
            }
        }
    }
    return U_iso;
}

void swap_sort(ivec order, cvec &v)
{
    int i = 0;
    while (i < v.size() - 1)
    {
        int new_index = 0;
        for (int j = i; j < v.size(); j++)
            if (order[j] < order[i])
                new_index++;
        if (new_index > 0)
        {
            std::complex<double> temp = v[i];
            v[i] = v[i + new_index];
            v[i + new_index] = temp;
            int temp2 = order[i];
            order[i] = order[i + new_index];
            order[i + new_index] = temp2;
        }
        else
            i++;
    }
}

void swap_sort_multi(ivec order, std::vector<ivec> &v)
{
    int i = 0;
    ivec temp;
    temp.resize(v.size());
    while (i < v.size() - 1)
    {
        int new_index = 0;
#pragma omp parallel for reduction(+ : new_index)
        for (int j = i; j < v.size(); j++)
            if (order[j] < order[i])
                new_index++;
        if (new_index > 0)
        {
#pragma omp parallel for
            for (int run = 0; run < v.size(); run++)
            {
                temp[run] = v[run][i];
                v[run][i] = v[run][i + new_index];
                v[run][i + new_index] = temp[run];
            }
            int temp2 = order[i];
            order[i] = order[i + new_index];
            order[i + new_index] = temp2;
        }
        else
            i++;
    }
}

double get_lambda_1(double *a)
{
    // Returns the middle eigenvalue (lambda2) of a row-major symmetric 3x3 matrix.
    const double a00 = a[0], a01 = a[1], a02 = a[2];
    const double a11 = a[4], a12 = a[5], a22 = a[8];

    const double p1 = a01 * a01 + a02 * a02 + a12 * a12;
    if (p1 < std::numeric_limits<double>::epsilon())
    {
        double eigs[3] = { a00, a11, a22 };
        std::sort(eigs, eigs + 3);
        return eigs[1];
    }

    const double q = (a00 + a11 + a22) / 3.0;
    const double b00 = a00 - q, b11 = a11 - q, b22 = a22 - q;
    const double p2 = b00 * b00 + b11 * b11 + b22 * b22 + 2.0 * p1;
    const double p = std::sqrt(p2 / 6.0);

    const double c00 = b00 / p, c01 = a01 / p, c02 = a02 / p;
    const double c11 = b11 / p, c12 = a12 / p, c22 = b22 / p;
    const double r = 0.5 * (c00 * (c11 * c22 - c12 * c12)
                          - c01 * (c01 * c22 - c12 * c02)
                          + c02 * (c01 * c12 - c11 * c02));

    const double phi = std::acos(std::clamp(r, -1.0, 1.0)) / 3.0;
    const double eig_max = q + 2.0 * p * std::cos(phi);
    const double eig_min = q + 2.0 * p * std::cos(phi + 2.0 * constants::PI / 3.0);
    return 3.0 * q - eig_max - eig_min;
}

double get_bessel_ratio(const double nu, const double x)
{
    const double RECUR_BIG = DBL_MAX;
    const double RECUR_SMALL = DBL_MIN;
    const int maxiter = 10000;
    int n = 1;
    double Anm2 = 1.0;
    double Bnm2 = 0.0;
    double Anm1 = 0.0;
    double Bnm1 = 1.0;
    double a1 = x / (2.0 * (nu + 1.0));
    double An = Anm1 + a1 * Anm2;
    double Bn = Bnm1 + a1 * Bnm2;
    double an;
    double fn = An / Bn;
    double dn = a1;
    double s = 1.0;

    while (n < maxiter)
    {
        double old_fn;
        double del;
        n++;
        Anm2 = Anm1;
        Bnm2 = Bnm1;
        Anm1 = An;
        Bnm1 = Bn;
        an = -x * x / (4.0 * (nu + n - 1.0) * (nu + n));
        An = Anm1 + an * Anm2;
        Bn = Bnm1 + an * Bnm2;

        if (fabs(An) > RECUR_BIG || fabs(Bn) > RECUR_BIG)
        {
            An /= RECUR_BIG;
            Bn /= RECUR_BIG;
            Anm1 /= RECUR_BIG;
            Bnm1 /= RECUR_BIG;
            Anm2 /= RECUR_BIG;
        }
        else if (fabs(An) < RECUR_SMALL || fabs(Bn) < RECUR_SMALL)
        {
            An /= RECUR_SMALL;
            Bn /= RECUR_SMALL;
            Anm1 /= RECUR_SMALL;
            Bnm1 /= RECUR_SMALL;
            Anm2 /= RECUR_SMALL;
            Bnm2 /= RECUR_SMALL;
        }

        old_fn = fn;
        fn = An / Bn;
        del = old_fn / fn;

        dn = 1.0 / (2.0 * (nu + n) / x - dn);
        if (dn < 0.0)
            s = -s;

        if (fabs(del - 1.0) < 2.0 * DBL_EPSILON)
            break;
    }

    return fn;
}

double bessel_first_kind(const int l, const double x)
{
    if (l < 0 || x < 0.0)
    {
        err_not_impl_f("This is not implemented, pelase dont do this to me!", std::cout);
        return -1000;
    }
    else if (x == 0.0)
    {
        return (l > 0 ? 0.0 : 1.0);
    }
    else if (l == 0)
    {
        return sin(x) / x;
    }
    else if (l == 1)
    {
        return (sin(x) / x - cos(x)) / x;
    }
    else if (l == 2)
    {
        const double f = (3.0 / (x * x) - 1.0);
        return (f * sin(x) - 3.0 * cos(x) / x) / x;
    }
    else if (l == 3)
    {
        double x2 = x * x;
        const double f1 = (x2 - 15.0);
        const double f2 = (6. * x2 - 15.);
        return (-f2 * sin(x) + f1 * cos(x) * x) / pow(x, 4);
    }
    else if (l == 4)
    {
        double x2 = x * x;
        const double f1 = (10. * x2 - 105.0);
        const double f2 = (x2 * x2 - 45. * x2 + 105.);
        return (f2 * sin(x) + f1 * cos(x) * x) / pow(x, 5);
    }
    else if (l == 5)
    {
        double x2 = x * x;
        double x4 = x2 * x2;
        const double f1 = (-x4 + 105.0 * x2 - 945.);
        const double f2 = (15. * x4 - 420. * x2 + 945.);
        return (f2 * sin(x) + f1 * cos(x) * x) / (x4 * x2);
    }
    else if (l == 6)
    {
        double x2 = x * x;
        double x4 = x2 * x2;
        const double f1 = 21. * (x4 - 60.0 * x2 + 495.);
        const double f2 = (-x4 * x2 + 210. * x4 - 4725 * x2 + 10395.);
        return (f2 * sin(x) - f1 * cos(x) * x) / pow(x, 7);
    }
    else
    {
        double ratio = get_bessel_ratio(l + 0.5, x);
        const double smallest = DBL_MIN / DBL_EPSILON;
        double jellp1 = smallest * ratio;
        double jell = smallest;
        double jellm1;
        int ell;
        for (ell = l; ell > 0; ell--)
        {
            jellm1 = -jellp1 + (2 * ell + 1) / x * jell;
            jellp1 = jell;
            jell = jellm1;
        }

        if (fabs(jell) > fabs(jellp1))
        {
            double pre = smallest / jell;
            return bessel_first_kind(0, x) * pre;
        }
        else
        {
            double pre = smallest / jellp1;
            return bessel_first_kind(1, x) * pre;
        }
    }
}


double get_decimal_precision_from_CIF_number(std::string &given_string)
{
    int len = (int)given_string.length();
    int open_bracket = -1;
    int close_bracket = -1;
    int decimal_point = -1;
    // const char* gs = given_string.c_str();
    for (int i = 0; i < len; i++)
    {
        if (given_string[i] == '(' && open_bracket == -1)
        {
            open_bracket = i;
        }
        else if (given_string[i] == ')' && close_bracket == -1)
        {
            close_bracket = i;
        }
        else if (given_string[i] == '.' && decimal_point == -1)
        {
            decimal_point = i;
        }
    }
    double result = 0;
    int precision = 1;
    int size_of_precision = 1;
    if (open_bracket != -1 && close_bracket != -1)
    {
        size_of_precision = close_bracket - open_bracket - 1;
        std::string temp = given_string.substr(open_bracket + 1, size_of_precision);
        precision = std::stoi(temp);
    }
    int digits = 0;
    if (open_bracket != -1 && close_bracket != -1)
    {
        if (decimal_point != -1)
        {
            digits = open_bracket - decimal_point - 1;
        }
        else
        {
            digits = close_bracket - open_bracket - 1;
        }
        if (digits == 0)
        {
            return 0.001;
        }
        result = abs(precision * pow(10, -digits));
        return result;
    }
    else
        return 0.005;
};

void options::digest_options()
{
    using namespace std;
    // Lets print what was the command line, for debugging
    if (debug)
    {
        std::cout << " Recap of input:\nsize: " << arguments.size() << endl;
    }
    // This loop figures out command line options
    int argc = (int)arguments.size();
    for (int i = 0; i < arguments.size(); i++)
    {
        if (debug)
            std::cout << arguments[i] << endl;
        string temp = arguments[i];
        if (temp.find("-") > 0)
            continue;
        if (temp == "-acc")
            accuracy = stoi(arguments[i + 1]);
        if (temp == "-all_charges")
            all_charges = true;
        else if (temp == "-Anion")
        {
            int n = 1;
            string store;
            if (debug)
                std::cout << "Looking for Anions!" << endl;
            while (i + n < argc && string(arguments[i + n]).find("-") != 0)
            {
                if (i + n - 1 > arguments.size())
                    break;
                store = arguments[i + n];
                svec Z = split_string<string>(store, " ");
                for (int r = 0; r < Z.size(); r++)
                {
                    if (debug)
                        std::cout << Z[r] << endl;
                    Anions.push_back(Z[r]);
                }
                n++;
            }
        }
        else if (temp == "-laplacian_bonds")
        {
            wfn = arguments[i + 1];
            bondwise_laplacian_plots(wfn);
            exit(0);
        }
        else if (temp == "-atom_dens")
        {
            std::cout << NoSpherA2_message() << endl;
            wfn = arguments[i + 1];
            err_checkf(std::filesystem::exists(wfn), "WFN doesn't exist", std::cout);
            ivec val_MOs;
            ivec val_MOs_beta;
            if (argc >= i + 3)
            {
                val_MOs = split_string<int>(arguments[i + 2], ",");
                std::cout << "Alpha MOs to keep: ";
                for (int j = 0; j < val_MOs.size(); j++)
                    std::cout << val_MOs[j] << " ";
                std::cout << endl;
                val_MOs_beta = split_string<int>(arguments[i + 3], ",");
                std::cout << "Beta MOs to keep: ";
                for (int j = 0; j < val_MOs_beta.size(); j++)
                    std::cout << val_MOs_beta[j] << " ";
                std::cout << endl;
            }
            spherically_averaged_density(*this, val_MOs, val_MOs_beta);
            exit(0);
        }
        else if (temp == "-atom_sfac")
        {
            std::cout << NoSpherA2_message() << endl;
            wfn = arguments[i + 1];
            wfn2 = arguments[i + 2];
            err_checkf(std::filesystem::exists(wfn), "WFN doesn't exist", std::cout);
            WFN wavy(e_origin::NOT_YET_DEFINED);
            wavy.read_known_wavefunction_format(wfn, std::cout, debug);
            wavy.delete_unoccupied_MOs();

            WFN wavy2(e_origin::NOT_YET_DEFINED);
            wavy2.read_known_wavefunction_format(wfn2, std::cout, debug);
            wavy2.delete_unoccupied_MOs();

            bvec needs_grid(wavy.get_ncen(), false);
            needs_grid[0] = true;
            GridConfiguration conf;
            conf.partition_type = PartitionType::Hirshfeld;
            conf.accuracy = 4;
            GridManager grid(conf);
            vec2 d1, d2, d3, dens;
            vec2 d1_2, d2_2, d3_2, dens_2;

            cell unit_cell(10.0, 10.0, 10.0, 90, 90, 90);
            ivec asym_atom_list(1, 0);
            // Setup grids for the molecule
            auto grid2 = grid;
            grid.setup3DGridsForMolecule(wavy, asym_atom_list, needs_grid, unit_cell);
            grid.getDensityVectors(wavy, asym_atom_list, d1, d2, d3, dens);
            grid2.setup3DGridsForMolecule(wavy2, asym_atom_list, needs_grid, unit_cell);
            grid2.getDensityVectors(wavy2, asym_atom_list, d1_2, d2_2, d3_2, dens_2);

            for (int j = 0; j < d1.size(); j++)
            {
                for (int p = 0; p < d1[0].size(); p++)
                {
                    dens[j][p] -= dens_2[j][p];
                }
            }

            // Calculate partitioned charges
            PartitionResults results = grid.calculatePartitionedCharges(wavy, unit_cell);

            grid.printChargeTable({ "Fe" }, wavy, asym_atom_list, std::cout, results);

            const int points = grid.getTotalGridPoints();
            for (double k = 0.001; k < 10; k += 0.002) {
                cdouble res = calc_spherically_averaged_at_k(d1, d2, d3, dens, k);
                std::cout << "k: " << k << " sfac: " << setprecision(9) << setw(16) << scientific << res.real() << " " << setprecision(9) << setw(16) << scientific << res.imag() << endl;
            }
            exit(0);
        }
        else if (temp == "-b")
            basis_set = arguments[i + 1];
        else if (temp == "-becke" || temp == "-BECKE" || temp == "-Becke")
            partition_type = PartitionType::Becke;
        else if (temp == "-lahvatest")
        {
            //_test_lahva();
            exit(0);
        }
        else if (temp == "-Cation")
        {
            int n = 1;
            string store;
            if (debug)
                std::cout << "Looking for Cations!" << endl;
            while (i + n < argc && string(arguments[i + n]).find("-") != 0)
            {
                if (i + n - 1 > arguments.size())
                    break;
                store = arguments[i + n];
                svec Z = split_string<string>(store, " ");
                for (int r = 0; r < Z.size(); r++)
                {
                    if (debug)
                        std::cout << Z[r] << endl;
                    Cations.push_back(Z[r]);
                }
                n++;
            }
        }
        else if (temp == "-charge")
        {
            charge = stoi(arguments[i + 1]);
        }
        else if (temp == "-coef")
        {
            coef_file = arguments[i + 1];
            err_checkf(std::filesystem::exists(coef_file), "coef_file doesn't exist", std::cout);
            SALTED = true;
        }
        else if (temp == "-cif")
        {
            cif = arguments[i + 1];
            err_checkf(std::filesystem::exists(cif), "CIF doesn't exist", std::cout);
        }
        else if (temp == "-cpus")
        {
            threads = stoi(arguments[i + 1]);
            MKL_Set_Num_Threads(threads);
#ifdef _OPENMP
            omp_set_num_threads(threads);
            omp_set_dynamic(0);
#endif
        }
        else if (temp == "-cmtc")
        {
            cif_based_combined_tsc_calc = true;
            int n = 1;
            string delimiter = ",";
            groups.pop_back();
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                combined_tsc_calc_files.push_back(arguments[i + n]);
                n++;
                combined_tsc_calc_cifs.push_back(arguments[i + n]);
                n++;
                const string _temp = arguments[i + n];
                if (_temp.find(delimiter) == string::npos) {
                    if (debug)
                        std::cout << "--Delimiter not found, using ." << endl;
                    delimiter = ".";
                }
                groups.push_back(split_string<int>(_temp, delimiter));
                if (debug)
                {
                    std::cout << "--Group: " << _temp << endl << "--";
                    for (int run = 0; run < groups[groups.size() - 1].size(); run++)
                        std::cout << groups[groups.size() - 1][run] << " ";
                    std::cout << endl;
                }
                n++;
            }
        }
        else if (temp == "-combine_mos")
        {
            combine_mo.push_back(arguments[i + 1]);
            combine_mo.push_back(arguments[i + 2]);
            do_combine_mo(*this);
            exit(0);
        }
        else if (temp == "-cmos1")
        {
            int j = 1;
            while (i + j < argc && arguments[i + j].find("-") >= 1)
            {
                cmo1.push_back(stoi(arguments[i + j]));
                j++;
            }
        }
        else if (temp == "-cmos2")
        {
            int j = 1;
            while (i + j < argc && arguments[i + j].find("-") >= 1)
            {
                cmo2.push_back(stoi(arguments[i + j]));
                j++;
            }
        }
        else if (temp == "-convert_to_47") {
            err_checkf(argc >= i + 2, "Not enough arguments for -convert_to_47\nPlease provide at least stdout name!", std::cout);
            std::filesystem::path _wfn = arguments[i + 1];
            WFN wavy(e_origin::NOT_YET_DEFINED);
            wavy.read_known_wavefunction_format(_wfn, std::cout, debug);
            wavy.write_nbo(_wfn.replace_extension(".47"), debug, &std::cout);
            exit(0);
        }
        else if (temp == "-convert_XCW")
        {
            err_checkf(argc >= i + 3, "Not enough arguments for -convert_XCW\nPlease provide at least stdout name and lambda step!", std::cout);
            std::string stdo = arguments[i + 1];
            std::string step = arguments[i + 2];
            convert_tonto_XCW_lambda_steps(stdo, step, debug, *this);
            exit(0);
        }
        else if (temp == "-def" || temp == "-DEF")
            properties.def = true;
        else if (temp == "-density_difference" || temp == "-density-difference")
        {
            wfn2 = arguments[i + 1];
        }
        else if (temp == "-dmin")
            dmin = stod(arguments[i + 1]);
        else if (temp == "-d")
            basis_set_path = arguments[i + 1];
        else if (temp == "-dipole_moments")
        {
            dipole_moments(*this);
            exit(0);
        }
        // Visualize the specified orbital using spherical harmonics.
        // Call as -draw_orbits lambda,m,resolution,radius
        // Where resolution and radius are optional
        else if (temp == "-draw_orbits")
        {
            vec opts = split_string<double>(arguments[i + 1], ",");
            int l = static_cast<int>(opts[0]);
            int m = static_cast<int>(opts[1]);
            properties.resolution = 0.025;
            properties.radius = 3.5;
            if (opts.size() >= 3)
            {
                properties.resolution = opts[2];
            }
            if (opts.size() == 4)
            {
                properties.radius = opts[3];
            }

            draw_orbital(l, m, properties.resolution, properties.radius);
            exit(0);
        }
        else if (temp == "-e_field")
            efield = stod(arguments[i + 1]);
        else if (temp == "-ECP" || temp == "-ecp" || temp == "-Ecp")
        {
            ECP = true;
            if (argc >= i + 2 && string(arguments[i + 1]).find("-") != 0)
            {
                ECP_mode = stoi(arguments[i + 1]);
            }
        }
        else if (temp == "-ED")
            electron_diffraction = true;
        else if (temp == "-eli")
            properties.eli = true;
        else if (temp == "-eli_analysis") {
            err_checkf(argc >= i + 4, "Not enough arguments for -eli_analysis\nPlease provide at least wfn, resolution and radius!", std::cout);
            wfn = arguments[i + 1];
            properties.resolution = stod(arguments[i + 2]);
            properties.radius = stod(arguments[i + 3]);
            ELI_analysis(wfn, *this);
            exit(0);
        }
        else if (temp == "-qtaim_eli") {
            // Cube-files mode:  -qtaim_eli <rho.cube> <eli.cube> <atoms_csv> [<bg_value>]
            // WFN mode:         -qtaim_eli <wfn_file> <atoms_csv> [<resolution>] [<radius>] [<bg_value>]
            err_checkf(i + 2 < argc,
                "Usage:\n"
                "  -qtaim_eli <rho.cube> <eli.cube> <atoms_csv> [bg_value]\n"
                "  -qtaim_eli <wfn_file> <atoms_csv> [resolution] [radius] [bg_value]\n"
                "atoms_csv: comma-separated 0-based atom indices, e.g. 0,3,7\n"
                "bg_value:  value assigned to non-selected voxels (default 0)",
                std::cout);

            std::filesystem::path arg1 = arguments[i + 1];
            const bool cube_mode = (arg1.extension() == ".cube");

            std::filesystem::path rho_path, eli_path_arg;
            std::string atoms_csv;
            double bg_val = 0.0;

            if (cube_mode) {
                err_checkf(i + 3 < argc,
                    "Cube mode requires: -qtaim_eli <rho.cube> <eli.cube> <atoms_csv>", std::cout);
                rho_path     = arg1;
                eli_path_arg = arguments[i + 2];
                atoms_csv    = arguments[i + 3];
                if (i + 4 < argc) bg_val = stod(arguments[i + 4]);
            } else {
                rho_path  = arg1;  // wfn path
                atoms_csv = arguments[i + 2];
                if (i + 3 < argc && !std::filesystem::path(arguments[i + 3]).has_extension())
                    ; // next arg looks like atoms_csv already consumed; nothing extra
                if (i + 3 < argc) {
                    try { properties.resolution = stod(arguments[i + 3]); }
                    catch (...) { /* optional — might be bg_val */ }
                }
                if (i + 4 < argc) {
                    try { properties.radius = stod(arguments[i + 4]); }
                    catch (...) {}
                }
                if (i + 5 < argc) {
                    try { bg_val = stod(arguments[i + 5]); }
                    catch (...) {}
                }
            }

            // Parse comma-separated 0-based atom indices
            std::vector<int> indices;
            {
                std::stringstream ss(atoms_csv);
                std::string tok;
                while (std::getline(ss, tok, ',')) {
                    std::string t = trim(tok);
                    if (!t.empty()) indices.push_back(std::stoi(t));
                }
            }
            err_checkf(!indices.empty(), "No atom indices parsed from: " + atoms_csv, std::cout);

            run_QTAIM_ELI_mask(rho_path, eli_path_arg, indices, bg_val, *this, log_file);
            exit(0);
        }
        else if (temp == "-elf")
            properties.elf = true;
        else if (temp == "-embis" || temp == "-EMBIS")
        {
            partition_type = PartitionType::EMBIS;
        }
        else if (temp == "-esp")
            properties.esp = true;
        else if (temp == "-fukui" || temp == "-Fukui")
            properties.fukui = true;
        else if (temp == "-fukui_analysis")
        {
            // Only records the request; the analysis itself runs from
            // run_app_impl.
            //
            // Doing the work here instead - the pattern -dipole_moments,
            // -laplacian_bonds and friends use - would put the whole table in
            // NoSpherA2.log and nothing on the terminal, because run_app_impl
            // has already pointed std::cout at that log file by the time
            // digest_options() runs. (Verified: -dipole_moments prints nothing
            // to stdout and 2.9 kB to the log.) That is tolerable for a
            // side-effect command and useless for one whose entire output is
            // meant to be read, so this one is dispatched later, where the
            // console buffer can be restored first.
            fukui_analysis_run = true;
            // Wavefunction may be given inline (-fukui_analysis mol.gbw) or with
            // the usual -wfn. The inline form is the point of the flag, so it
            // wins, but only if the next token is not itself an option.
            if (i + 1 < argc && arguments[i + 1][0] != '-')
            {
                wfn = arguments[i + 1];
                i++;
            }
        }
        else if (temp == "-ewal_sum")
        {
            // bool read, WFN& wave, std::ostream& file,
            WFN *temp_w = new WFN(e_origin::cub);
            cube residual(arguments[i + 1], true, *temp_w, std::cout);
            if (argc >= i + 3)
            {
                int k_max = stoi(arguments[i + 2]);
                if (argc >= i + 4)
                    residual.ewald_sum(k_max, stod(arguments[i + 3]));
                else
                    residual.ewald_sum(k_max);
            }
            else
                residual.ewald_sum();
            delete (temp_w);
            exit(0);
        }
        else if (temp == "-geometry_aid_cutoff") {
            // **Selects which geometry-aid model family the descriptor is for.**
            // Must appear BEFORE the flag that computes anything, like -wfn --
            // the descriptor flags exit(0) as soon as they have run.
            //
            //   3.5  (default)  the `c_only` models, trained on all-carbon input
            //   3.0             the `dirty` models, trained on partly wrong
            //                   labels and the only ones that can be iterated
            //
            // A descriptor computed at the wrong cutoff is 42,042 long either
            // way and is rejected by nothing downstream, so the value used is
            // echoed rather than applied silently.
            // `arguments`, not `argv` -- this parser walks a vector<string>.
            err_chkf(i + 1 < arguments.size(),
                     "-geometry_aid_cutoff needs a radius in Angstrom", std::cout);
            geometry_aid_cutoff_radius = std::stod(arguments[++i]);
            std::cout << "  geometry-aid SOAP cutoff radius set to "
                      << geometry_aid_cutoff_radius << " A ("
                      << (geometry_aid_cutoff_radius > 3.25 ? "c_only" : "dirty")
                      << " model family)" << std::endl;
        }
        else if (temp == "-calc_featomic_descriptor") {
            err_chkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEFORE -calc_featomic_descriptor to specify a molecule.", std::cout);
            // Unchanged behaviour, down to the output file name: Olex2 calls
            // exactly this and reads `descriptor.npy` from the working
            // directory. The hyperparameters moved to a shared function so the
            // batch flag below cannot drift away from them.
            write_featomic_descriptor(wfn, "descriptor.npy", geometry_aid_hyperparameters());
            exit(0);
        }
        else if (temp == "-classify_atoms") {
            // -wfn <structure> -classify_atoms <model.bin> [<out.npy>]
            //
            // Writes class probabilities, one row per atom, in the class order
            // the model carries. Default output `probabilities.npy`, matching
            // the way `-calc_featomic_descriptor` defaults to descriptor.npy.
            err_chkf(!wfn.empty(), "No structure specified! Use -wfn BEFORE -classify_atoms.", std::cout);
            const std::filesystem::path model_path = arguments[i + 1];
            err_chkf(std::filesystem::exists(model_path),
                "The geometry-aid model does not exist: " + model_path.string(), std::cout);
            std::filesystem::path out_path = "probabilities.npy";
            if (i + 2 < arguments.size() && !arguments[i + 2].empty() && arguments[i + 2][0] != '-')
                out_path = arguments[i + 2];
            write_geometry_aid_probabilities(wfn, out_path, model_path,
                geometry_aid_hyperparameters());
            exit(0);
        }
        else if (temp == "-classify_atoms_list") {
            // The batch form: -classify_atoms_list <list> <model.bin>, one
            // structure path per line, each written as <path>.probs.npy.
            const std::filesystem::path list_file = arguments[i + 1];
            const std::filesystem::path model_path = arguments[i + 2];
            err_chkf(std::filesystem::exists(list_file), "The structure list does not exist: " + list_file.string(), std::cout);
            err_chkf(std::filesystem::exists(model_path), "The geometry-aid model does not exist: " + model_path.string(), std::cout);

            std::vector<std::filesystem::path> jobs;
            {
                std::ifstream list(list_file);
                std::string line;
                while (std::getline(list, line))
                {
                    const std::string entry = trim(line);
                    if (entry.empty() || entry[0] == '#') continue;
                    jobs.push_back(std::filesystem::path(entry));
                }
            }
            err_chkf(!jobs.empty(), "The structure list is empty: " + list_file.string(), std::cout);

            const SALTED_Utils::FeatomicHyperParameters hyperparams = geometry_aid_hyperparameters();
            const auto started = std::chrono::steady_clock::now();
            size_t done = 0, failed = 0;
            for (const std::filesystem::path& structure : jobs)
            {
                if (!std::filesystem::exists(structure))
                {
                    std::cout << "MISSING " << structure.string() << std::endl;
                    ++failed;
                    continue;
                }
                try
                {
                    const auto one = std::chrono::steady_clock::now();
                    std::filesystem::path out_path = structure;
                    out_path += ".probs.npy";
                    write_geometry_aid_probabilities(structure, out_path, model_path, hyperparams);
                    std::cout << "PROBABILITIES " << out_path.string() << " seconds="
                              << std::chrono::duration<double>(std::chrono::steady_clock::now() - one).count()
                              << std::endl;
                    ++done;
                }
                catch (const std::exception& e)
                {
                    std::cout << "FAILED " << structure.string() << " : " << e.what() << std::endl;
                    ++failed;
                }
            }
            const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
            std::cout << "BATCH done=" << done << " failed=" << failed
                      << " seconds=" << total
                      << " per_structure=" << (done ? total / done : 0.0) << std::endl;
            exit(failed && !done ? 1 : 0);
        }
        else if (temp == "-calc_featomic_descriptors") {
            // **Many structures in one process.** Building the featomic
            // calculator splines the radial integral for every (n, l) pair and
            // costs a fixed 0.72 s; the descriptor itself costs about 0.0009 s
            // per atom. Measured 6 August 2026: 2 atoms 0.74 s, 30 atoms 0.90 s,
            // 1000 atoms 1.57 s. One structure per process therefore spends
            // roughly ninety per cent of its time on setup that does not depend
            // on the molecule at all.
            //
            // `calculate_SOAP_Powerspectrum` keeps the calculator for the
            // lifetime of the process (see `cached_calculator`), so everything
            // after the first entry here pays only the marginal cost.
            //
            // The list file is **one structure path per line**, and the
            // descriptor is written beside it as `<path>.npy`. One field per
            // line deliberately: an output column would need quoting rules, and
            // paths with spaces in them are ordinary on Windows.
            // Blank lines and lines beginning with '#' are ignored.
            const std::filesystem::path list_file = arguments[i + 1];
            err_chkf(std::filesystem::exists(list_file), "The structure list does not exist: " + list_file.string(), std::cout);

            std::vector<std::filesystem::path> jobs;
            {
                std::ifstream list(list_file);
                std::string line;
                while (std::getline(list, line))
                {
                    const std::string entry = trim(line);
                    if (entry.empty() || entry[0] == '#')
                        continue;
                    jobs.push_back(std::filesystem::path(entry));
                }
            }
            err_chkf(!jobs.empty(), "The structure list is empty: " + list_file.string(), std::cout);

            const SALTED_Utils::FeatomicHyperParameters hyperparams = geometry_aid_hyperparameters();
            const auto started = std::chrono::steady_clock::now();
            size_t done = 0, failed = 0;
            for (const std::filesystem::path& structure : jobs)
            {
                // **One bad structure must not cost the other 999.** A batch is
                // the whole point of this flag, and a list of ten thousand
                // peaks files will contain one that cannot be parsed.
                if (!std::filesystem::exists(structure))
                {
                    std::cout << "MISSING " << structure.string() << std::endl;
                    ++failed;
                    continue;
                }
                try
                {
                    const auto one = std::chrono::steady_clock::now();
                    std::filesystem::path out_path = structure;
                    out_path += ".npy";
                    write_featomic_descriptor(structure, out_path, hyperparams);
                    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - one).count();
                    std::cout << "DESCRIPTOR " << out_path.string() << " seconds=" << seconds << std::endl;
                    ++done;
                }
                catch (const std::exception& e)
                {
                    std::cout << "FAILED " << structure.string() << " : " << e.what() << std::endl;
                    ++failed;
                }
            }
            const double total = std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count();
            std::cout << "BATCH done=" << done << " failed=" << failed
                      << " seconds=" << total
                      << " per_structure=" << (done ? total / done : 0.0) << std::endl;
            exit(failed && !done ? 1 : 0);
        }
        else if (temp == "-fchk")
            fchk = arguments[i + 1];
        else if (temp == "-fractal")
            fract = true, fract_name = arguments[i + 1];
        else if (temp == "-gbw2wfn")
            gbw2wfn = true;
        else if (temp == "-get_g")
            get_g = true;
        else if (temp == "-group")
        {
            int n = 1;
            while (i + n < argc)
            {
                const string& group_arg = arguments[i + n];
                // Olex2 can emit a bare -group followed by an empty argument.
                // Do not index or parse an empty string: doing so corrupts the
                // CRT state and later manifests as a stack-buffer-overrun.
                if (group_arg.empty() || group_arg.find("-") != string::npos)
                    break;
                int group;
                if (group_arg[0] == '+')
                    group = -stoi(group_arg);
                else
                    group = stoi(group_arg);
                groups[0].push_back(group), n++;
            }
            i += n - 1;
        }
        else if (temp == "-HDEF")
            properties.hdef = true;
        else if (temp == "-hirsh")
            properties.hirsh = true, properties.hirsh_number = stoi(arguments[i + 1]);
        else if (temp == "-hirshfeld_surface")
        {
            hirshfeld_surface = arguments[i + 1];
            hirshfeld_surface2 = arguments[i + 2];
        }
        else if (temp == "-tsc_block")
        {
            tsc_block_size = static_cast<size_t>(std::stoll(arguments[i + 1]));
            tsc_block_given = true;
        }
        else if (temp == "-hkl")
        {
            hkl = arguments[i + 1];
            err_checkf(std::filesystem::exists(hkl), "hkl doesn't exist", std::cout);
        }
        else if (temp == "-hkl_min_max")
        {
            int h_min(stoi(arguments[i + 1]));
            int h_max(stoi(arguments[i + 2]));
            int k_min(stoi(arguments[i + 3]));
            int k_max(stoi(arguments[i + 4]));
            int l_min(stoi(arguments[i + 5]));
            int l_max(stoi(arguments[i + 6]));
            hkl_min_max = { {h_min, h_max}, {k_min, k_max}, {l_min, l_max} };
        }
        else if (temp == "-IAM")
            iam_switch = true;
        else if (temp == "-lap")
            properties.lap = true;
        else if (temp == "-mbis" || temp == "-MBIS")
        {
            partition_type = PartitionType::MBIS;
        }
        else if (temp == "-mem")
        {
            mem = stod(arguments[i + 1]); // In MB
            mem_given = true;
            vec a;
            size_t vec_max_size = a.max_size();
            double doubel_max_size = static_cast<double>(vec_max_size * sizeof(double)) * 1e-6;
            if (mem > doubel_max_size)
            {
                std::cout << "Max memory set to " << mem << " MB, which is larger than the maximum allowed size of " << doubel_max_size << " MB. Setting max memory to " << 50000 << " MB." << endl;
                mem = 50000.0;
            }
        }
        else if (temp == "-method")
            method = arguments[i + 1];
        else if (temp == "-merge")
        {
            pathvec filenames;
            int n = 1;
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                filenames.push_back(arguments[i + n]);
                n++;
            }
            merge_tscs("combine", filenames, old_tsc);
            exit(0);
        }
        else if (temp == "-merge_nocheck")
        {
            pathvec filenames;
            int n = 1;
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                filenames.push_back(arguments[i + n]);
                n++;
            }
            merge_tscs_without_checks("combine", filenames, old_tsc);
            exit(0);
        }
        else if (temp == "-MO")
        {
            if (string(arguments[i + 1]) != "all")
                properties.MO_numbers.push_back(stoi(arguments[i + 1]));
            else
                properties.all_mos = true;
        }
        else if (temp == "-mtc")
        {
            combined_tsc_calc = true;
            int n = 1;
            string delimiter = ",";
            groups.pop_back();
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                combined_tsc_calc_files.push_back(arguments[i + n]);
                n++;
                const string _temp = arguments[i + n];
                if (_temp.find(delimiter) == string::npos) {
                    if (debug)
                        std::cout << "--Delimiter not found, using ." << endl;
                    delimiter = ".";
                }
                groups.push_back(split_string<int>(_temp, delimiter));
                if (debug)
                {
                    std::cout << "--Group: " << _temp << endl << "--";
                    for (int run = 0; run < groups[groups.size() - 1].size(); run++)
                        std::cout << groups[groups.size() - 1][run] << " ";
                    std::cout << endl;
                }
                n++;
            }
        }
        else if (temp == "-mtc_mult")
        {
            int n = 1;
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                combined_tsc_calc_mult.push_back(stoi(arguments[i + n]));
                n++;
            }
        }
        else if (temp == "-mtc_charge")
        {
            int n = 1;
            while (i + n < argc && string(arguments[i + n]).find("-") > 0)
            {
                if (arguments[i + n][0] == 'n')
                    combined_tsc_calc_charge.push_back(-stoi(arguments[i + n].substr(1)));
                else
                    combined_tsc_calc_charge.push_back(stoi(arguments[i + n]));
                n++;
            }
        }
        else if (temp == "-mtc_ECP")
        {
            int m = 1;
            while (i + m < argc && string(arguments[i + m]).find("-") > 0)
            {
                combined_tsc_calc_ECP.push_back(stoi(arguments[i + m]));
                m++;
            }
        }
        else if (temp == "-mult")
            mult = stoi(arguments[i + 1]);
        else if (temp == "-NNLS_TEST")
        {
            test_NNLS();
            exit(0);
        }
        else if (temp == "-no-date" || temp == "-no_date")
            no_date = true;
        else if (temp == "-pbc")
            pbc = stoi(arguments[i + 1]);
        else if (temp == "-polarizabilities")
        {
            pol_wfns = { arguments[i + 1],
                        arguments[i + 2],
                        arguments[i + 3],
                        arguments[i + 4],
                        arguments[i + 5],
                        arguments[i + 6],
                        arguments[i + 7] };
        }
        else if (temp == "-QCT" || temp == "-qct")
            qct = true;
        else if (temp == "-profiling" || temp == "-profile")
        {
            profiling = true;
            if (i + 1 < argc && arguments[i + 1].find("-") != 0)
            {
                profiling_tests_root = arguments[i + 1];
            }
        }
        else if (temp == "-promol_nci")
        {
            err_checkf(i + 2 < argc,
                "Usage: -promol_nci <frag1.xyz> <frag2.xyz> [rcut1=0.95] [rcut2=0.75] [rho_abs_max] [rdg_max]",
                std::cout);
            promol_nci = true;
            promol_nci_xyz1 = arguments[i + 1];
            promol_nci_xyz2 = arguments[i + 2];
            err_checkf(std::filesystem::exists(promol_nci_xyz1), "First XYZ file doesn't exist: " + promol_nci_xyz1.string(), std::cout);
            err_checkf(std::filesystem::exists(promol_nci_xyz2), "Second XYZ file doesn't exist: " + promol_nci_xyz2.string(), std::cout);

            double *optional_values[] = {
                &properties.promol_nci_rcut1,
                &properties.promol_nci_rcut2,
                &properties.promol_nci_rho_abs_max,
                &properties.promol_nci_rdg_max
            };
            int optional_index = 0;
            while (optional_index < 4 && i + 3 + optional_index < argc)
            {
                try
                {
                    size_t consumed = 0;
                    const std::string &candidate = arguments[i + 3 + optional_index];
                    const double value = std::stod(candidate, &consumed);
                    if (consumed != candidate.size())
                        break;
                    *optional_values[optional_index] = value;
                    optional_index++;
                }
                catch (...)
                {
                    break;
                }
            }

            i += 2 + optional_index;
        }
        else if (temp == "-promol_nci_single_thread")
            properties.promol_nci_single_threaded = true;
        else if (temp == "-radius")
            properties.radius = stod(arguments[i + 1]);
        else if (temp == "-resolution")
            properties.resolution = stod(arguments[i + 1]);
        else if (temp == "-refine")
            argc > i + 1 ? properties.integral_accuracy = stod(arguments[i + 1]) : properties.integral_accuracy = 0.1;
        else if (temp == "-rdg")
            properties.rdg = true;
        else if (temp == "-rgbi")
            rgbi = true;
        else if (temp == "-rgbi_no_sym") {
            rgbi = true;
            rgbi_no_sym = true;
        }
        else if (temp == "-rgbi_EVs") {
            rgbi_EVs = true;
        }
        else if (temp == "-rgbi_basis") {
            err_checkf(i + 1 < argc, "Not enough arguments for -rgbi_basis. Use 'nao' or 'ano'.", std::cout);
            rgbi = true;
            std::string basis = arguments[i + 1];
            std::transform(basis.begin(), basis.end(), basis.begin(),
                [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
            if (basis == "nao")
                rgbi_orbital_basis = RGBIOrbitalBasis::NAO;
            else if (basis == "ano")
                rgbi_orbital_basis = RGBIOrbitalBasis::ANO;
            else
                err_checkf(false, "Invalid -rgbi_basis value '" + basis + "'. Use 'nao' or 'ano'.", std::cout);
        }
        else if (temp == "-rgbi-groups") {
            int n = 1;
            ivec2 group_set;
            while (i + n < argc && string(arguments[i + n]).find("-") > 0) {
                group_set.push_back(parse_rgbi_group_indices(arguments[i + n]));
                n++;
            }
            if (!group_set.empty()) {
                validate_rgbi_group_set(group_set);
                rgbi_group_sets.push_back(group_set);
            }
            i += n - 1;
            rgbi = true;
        }
        else if (temp.find("-rkpts") < 1)
            read_k_pts = true;
        else if (temp == "-rho_cube")
        {
            string wfn_name = arguments[i + 1];
            std::cout << "Reading wavefunction: " << wfn_name << endl;
            WFN wavy = WFN(wfn_name);
            std::cout << "Assigning ECPs" << endl;
            if (ECP)
                wavy.set_has_ECPs(true);
            std::cout << "Starting cube calculation" << endl;
            wavy.write_rho_cube();
            exit(0);
        }
        else if (temp == "-RI_FIT" || temp == "-ri_fit")
        {
            RI_FIT = true;
            partition_type = PartitionType::RI;
            int next_basis_set = i + 1;
            // Check if next argument is a valid basis set name or a new argument starting with "-"
            while (next_basis_set < argc && arguments[next_basis_set].find("-") != 0) {
                if (arguments[next_basis_set] == "auto_aux") {
                    double beta = 2.0;
                    //Check if the next argument is a valid double
                    if (next_basis_set + 1 < argc && arguments[next_basis_set + 1].find("-") != 0) {
                        beta = std::stod(arguments[next_basis_set + 1]);
                    }
                    aux_basis.push_back(std::make_shared<BasisSet>());
                    break;
                }
                err_chkf(BasisSetLibrary::check_basis_set_exists(arguments[next_basis_set]),
                    "Basis set " + arguments[next_basis_set] + " not found in the library. Exiting.", std::cout);
                aux_basis.push_back(BasisSetLibrary::get_basis_set(arguments[next_basis_set]));
                next_basis_set++;
            }
            if (aux_basis.size() == 0) {
                cout << "No basis set specified. Falling back to automatic generation using beta = 2.0!" << endl;
                aux_basis.push_back(std::make_shared<BasisSet>());
            }
        }
        else if (temp == "-write_ri_coefs") {
            WFN wavy(wfn);
            WFN wavy_aux = generate_aux_wfn(wavy, aux_basis);
            DensityFitting::CONFIG config;
            config.analyze_quality = debug;
            //config.restrain_type = DensityFitting::RESTRAINT_TYPE::SIMPLE_AND_TIK;
            //config.charge_scheme = DensityFitting::CHARGE_SCHEME::HIRSHFELD;
            vec ri_coefs = DensityFitting::density_fit(wavy, wavy_aux, config);
            npy::npy_data<double> np_coeffs;
            np_coeffs.data = ri_coefs;
            np_coeffs.fortran_order = false;
            np_coeffs.shape = { static_cast<unsigned long>(ri_coefs.size()) };
            npy::write_npy("RI_COEFS.npy", np_coeffs);
            exit(0);
        }
        else if (temp == "-RI_CUBE" || temp == "-ri_cube")
        {
            err_chkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEVORE -SALTED_COEFS to specify a wavefunction.", std::cout);
            err_checkf(!aux_basis.empty(), "No auxiliary basis set specified! Use -RI_FIT option BEVORE -test_RI to specify an auxiliary basis set.", std::cout);
            std::string coef_file = arguments[i + 1];
            std::vector<unsigned long> shape{};
            bool fortran_order;
            vec coefs{};

            npy::LoadArrayFromNumpy(coef_file, shape, fortran_order, coefs);


            WFN wavy(wfn);
            WFN wavy_aux = generate_aux_wfn(wavy, aux_basis);

            int nr_coefs = 0;
            for (const atom &atm : wavy_aux.get_atoms()) {
                int prim = 0;
                for (int shell = 0; shell < atm.get_shellcount_size(); shell++) {
                    const int type = atm.get_basis_set_entry(prim).get_type();
                    nr_coefs += 2 * type + 1;
                    prim += atm.get_shellcount(shell);
                }
            }

            std::cout << coefs.size() << " vs. " << nr_coefs << " ceofficients" << std::endl;

            // First name of coef_file, second name of xyz file
            cube_from_coef_npy(coefs, wavy_aux);

            // std::string aux_basis = arguments[i + 1];
            //gen_CUBE_for_RI(wavy, "def2_qzvppd_rifit", this);
            //gen_CUBE_for_RI(wavy, "def2_universal_jkfit", this);
            //gen_CUBE_for_RI(wavy, "combo-basis-fit", this);
            //gen_CUBE_for_RI(wavy, "cc-pvqz-jkfit", this);

            exit(0);
        }

        else if (temp.find("-s_rho") < 1)
            properties.s_rho = true;
        else if (temp == "-SALTED" || temp == "-salted")
        {
            SALTED = true;
            salted_model_dir = arguments[i + 1];
        }
        else if (temp == "-SALTED_COEFS" || temp == "-salted_coefs")
        {
            salted_model_dir = arguments[i + 1];

            //Check that wfn is not empty
            err_chkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEVORE -SALTED_COEFS to specify a wavefunction.", std::cout);

            WFN wavy(wfn);
            SALTEDPredictor SP(wavy, *this);
            filesystem::path salted_model_path = SP.get_salted_filename();
            log_file << "Using " << salted_model_path << " for the prediction" << endl;
            if (!SP.basis_set_loaded()) {
                const string df_basis_name = SP.get_dfbasis_name();
                std::shared_ptr<BasisSet> _aux_basis = BasisSetLibrary::get_basis_set(df_basis_name);
                load_basis_into_WFN(SP.wavy, _aux_basis);
            }
            vec coefs = SP.gen_SALTED_densities();
            npy::npy_data<double> np_coeffs;
            np_coeffs.data = coefs;
            np_coeffs.fortran_order = false;
            np_coeffs.shape = { static_cast<unsigned long>(coefs.size()) };
            npy::write_npy("SALTED_COEFS.npy", np_coeffs);
        }
        else if (temp == "-SALTED_Training") {
            err_chkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEVORE -test_RI to specify a wavefunction.", std::cout);
            err_checkf(!aux_basis.empty(), "No auxiliary basis set specified! Use -RI_FIT option BEVORE -test_RI to specify an auxiliary basis set.", std::cout);

            WFN wavy(wfn);
            WFN wavy_aux = generate_aux_wfn(wavy, aux_basis);

            create_SALTED_training_data(wavy, wavy_aux);
            exit(0);
        }
        else if (temp == "-skpts")
            save_k_pts = true;
        else if (temp == "-sfac_diffuse")
        {
            sfac_diffuse[0] = fromString<double>(arguments[i + 1]);
            sfac_diffuse[1] = fromString<double>(arguments[i + 2]);
            sfac_diffuse[2] = fromString<double>(arguments[i + 3]);
            cif = arguments[i + 4];
            wfn = arguments[i + 5];
            dmin = fromString<double>(arguments[i + 6]);
            calc_sfac_diffuse(*this, std::cout);
        }
        else if (temp == "-atom_dens_diff")
        {
            filesystem::path name_wfn_1 = arguments[i + 1];
            filesystem::path name_wfn_2 = arguments[i + 2];

            subtract_dens_from_gbw(name_wfn_1, name_wfn_2, 2, 0.05);
            exit(0);
        }
        else if (temp == "-spherical_aver_fukui")
        {
            filesystem::path wfn1_name = arguments[i + 1];
            filesystem::path wfn2_name = arguments[i + 2];
            WFN *wavy1 = new WFN(wfn1_name);
            WFN *wavy2 = new WFN(wfn2_name);
            ofstream outputFile("fukui_averaged_density_wfn.dat");
            for (double r = 0.001; r < 10.0; r += 0.001)
            {
                // double dens = calc_grid_averaged_at_r_from_cube(cube_from_file, r, 360, 5800);
                double dens = calc_fukui_averaged_at_r(*wavy1, *wavy2, r, 5810, 5810);
                outputFile << r << " " << dens << "\n";
            }
            outputFile.close();
            std::cout << "Data written to output.dat" << endl;
            delete (wavy1);
            delete (wavy2);
            exit(0);
        }
        else if (temp == "-spherical_aver_hirsh")
        {
            string wfn_name = arguments[i + 1];
            std::cout << "Reading wavefunction: " << wfn_name << endl;
            WFN *wavy = new WFN(wfn_name);
            std::cout << "Assigning ECPs" << endl;
            if (ECP)
                wavy->set_has_ECPs(true);
            std::cout << "Starting spherical averaging" << endl;
            double dens;

            for (int index_atom = 0; index_atom < wavy->get_ncen(); index_atom += 1)
            {
                ofstream outputFile("hirsh_averaged_density_" + std::to_string(index_atom) + ".dat");
                for (double r = 0.001; r < 5.0; r += 0.002)
                {
                    dens = calc_hirsh_grid_averaged_at_r(*wavy, index_atom, r, 360, 5800);
                    outputFile << r << " " << dens << "\n";
                }
                outputFile.close();
            }
            std::cout << "Data written to output.dat" << endl;
            delete (wavy);
            exit(0);
        }
        else if (temp == "-spherical_harmonic")
        {
            spherical_harmonic_test();
            exit(0);
        }
        else if (temp == "-spherical_atoms") {
            write_spherical_atoms();
            exit(0);
        }
        else if (temp == "-test")
            std::cout << "Running in test mode!" << endl, test = true;
        else if (temp == "-twin")
        {
            twin_law.resize(twin_law.size() + 1);
            twin_law.back().resize(9);
            for (int twl = 0; twl < 9; twl++)
                twin_law.back()[twl] = stod(arguments[i + 1 + twl]);
            if (debug)
            {
                std::cout << "twin_law: ";
                for (int twl = 0; twl < 9; twl++)
                    std::cout << setw(7) << setprecision(2) << twin_law.back()[twl];
                std::cout << endl;
            }
            i += 9;
        }
        else if (temp == "-old_tsc")
        {
            old_tsc = true;
        }
        else if (temp == "-tfvc" || temp == "-TFVC")
        {
            partition_type = PartitionType::TFVC;
        }
        else if (temp == "-tscb")
        {
            std::filesystem::path name = arguments[i + 1];
            string cif_name = "test.cif";
            if (name.extension() == ".tscb")
            {
                tsc_block<int, cdouble> blocky = read_tsc_table(name);
                blocky.write_tsc_file(cif_name, name.replace_extension(".tsc"));
            }
            else if (name.extension() == ".tsc")
            {
                tsc_block<int, cdouble> blocky = read_tsc_table(name);
                blocky.write_tscb_file(cif_name, name.replace_extension(".tscb"));
            }
            else
                err_checkf(false, "Wrong file ending!", std::cout);
            exit(0);
        }
        else if (temp == "-tsc_labels")
        {
            const bool has_table_and_cif =
                i + 2 < argc &&
                arguments[i + 1].find('-') != 0 &&
                arguments[i + 2].find('-') != 0;

            if (has_table_and_cif)
            {
                const std::filesystem::path table = arguments.at(i + 1);
                const std::filesystem::path cif_file = arguments.at(i + 2);
                std::filesystem::path output = table;
                output.replace_extension(".labels.tsc");
                if (i + 3 < argc && arguments[i + 3].find('-') != 0)
                    output = arguments[i + 3];
                if (!convert_tsc_ids_to_labels(table, cif_file, output, std::cout))
                    exit(1);
                exit(0);
            }

            label_tsc_output = true;
        }
        else if (temp == "-test_RI")
        {
            err_chkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEVORE -test_RI to specify a wavefunction.", std::cout);
            err_checkf(!aux_basis.empty(), "No auxiliary basis set specified! Use -RI_FIT option BEVORE -test_RI to specify an auxiliary basis set.", std::cout);


            WFN wavy(wfn);
            WFN wavy_aux = generate_aux_wfn(wavy, aux_basis);
            DensityFitting::demonstrate_enhanced_density_fitting(wavy, wavy_aux);
            exit(0);

        }
        else if (temp == "-RI_WFN_DIFF") {
            err_chkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEVORE -RI_WFN_DIFF to specify a wavefunction.", std::cout);
            err_checkf(!aux_basis.empty(), "No auxiliary basis set specified! Use -RI_FIT option BEVORE -RI_WFN_DIFF to specify an auxiliary basis set.", std::cout);
            WFN wavy(wfn);
            WFN wavy_aux = generate_aux_wfn(wavy, aux_basis);
            DensityFitting::QM_RI_difference_cube(wavy, wavy_aux);
            exit(0);

            //exit(0);
        }
        else if (temp == "-wfn")
        {
            wfn = arguments[i + 1];
            err_checkf(std::filesystem::exists(wfn), "Wavefunction does not exist!", std::cout);
        }
        else if (temp == "-cube_density" || temp == "-cube")
        {
            cube_density = arguments[i + 1];
            err_checkf(std::filesystem::exists(cube_density), "Cube density file does not exist!", std::cout);
        }
        else if (temp == "-wfn_cif")
        {
            write_CIF = true;
        }
        else if (temp == "-xyz")
        {
            xyz_file = arguments[i + 1];
        }
        else if (temp == "-do_XCW") {
            do_XCW = true;
            // Optional trailing "stepsize max_value" to limit the lambda scan
            // range, e.g. for quick tests: -do_XCW 0.01 0.01
            // CURRENTLY NOT IN USE SINCE THIS IS HANDLED IN THE INPUT FILE
            //if (i + 2 < argc &&
            //    string(arguments[i + 1]).find("-") != 0 &&
            //    string(arguments[i + 2]).find("-") != 0)
            //{
            //    xcw_lambda_step = stod(arguments[i + 1]);
            //    xcw_lambda_max = stod(arguments[i + 2]);
            //}
        }
        else if (temp == "-calc_F") {
            calc_F_calc = true;
        }
        else if (temp == "-xcw_gaussian_halt") {
            xcw_gaussian_halt = true;
        }
        else if (temp == "-xcw_strong_cutoff") {
            xcw_strong_cutoff = stod(arguments[i + 1]);
        }
        else if (temp == "-XCW_settings") {
			xcw_settings_path = arguments[i + 1];
        }
        else if (temp == "-anom_disp")
        {
            anom_disp_path = arguments[i + 1];
        }
        else if (temp == "-partitioning_test")
        {
            calc_partition_densities();
        }
        else if (temp == "-occ")
        {
            occ = arguments[i + 1];
            err_checkf(std::filesystem::exists(occ), "OCC input doesn't exist!", std::cout);

        }
        else if (temp == "-lukas_test")
        {
            //Check that wfn is not empty
            err_chkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEFORE -SALTED_COEFS to specify a wavefunction.", std::cout);
            err_chkf(!salted_model_dir.empty(), "No SALTED model directory specified! Use -SALTED option BEFORE -lukas_test to specify a model directory.", std::cout);

            WFN wavy(wfn);
            SALTEDPredictor SP(wavy, *this);
            string df_basis_name = SP.get_dfbasis_name();
            filesystem::path salted_model_path = SP.get_salted_filename();
            log_file << "Using " << salted_model_path << " for the prediction" << endl;
            if (!SP.basis_set_loaded()) {
                std::shared_ptr<BasisSet> aux_basis = BasisSetLibrary::get_basis_set(df_basis_name);
                load_basis_into_WFN(SP.wavy, aux_basis);
            }
            vec coefs = SP.gen_SALTED_densities();

            cube atom_cube = calc_cube_ML(coefs, SP.wavy);
            atom_cube.write_file("DBA_total.cube");

            for (int atm_idx = 0; atm_idx < wavy.get_ncen(); atm_idx++) {
                atom_cube = calc_cube_ML(coefs, SP.wavy, atm_idx);
                atom_cube.write_file("DBA_atom_" + std::to_string(atm_idx) + ".cube");
            }
        }
        else if (temp == "-calc_dens_1D")
        {
            err_chkf(!wfn.empty(), "No wavefunction specified! Use -wfn option BEVORE -calc_dens_1D to specify a wavefunction.", std::cout);
            err_checkf(!aux_basis.empty(), "No auxiliary basis set specified! Use -RI_FIT option BEVORE -calc_dens_1D to specify an auxiliary basis set.", std::cout);
            int atom_idx_1 = 0;
            int atom_idx_2 = 1;
            int gridpoints = 1000;
            double padding = 2.0;
            if (string(arguments[i + 1]).find("-") != 0) {
                atom_idx_1 = stoi(arguments[i + 1]);
            }
            if (string(arguments[i + 2]).find("-") != 0) {
                atom_idx_2 = stoi(arguments[i + 2]);
            }
            if (string(arguments[i + 3]).find("-") != 0) {
                gridpoints = stoi(arguments[i + 3]);
            }
            if (string(arguments[i + 4]).find("-") != 0) {
                padding = stod(arguments[i + 4]);
            }

            WFN wavy(wfn);
            string out = "Calculating 1D density between atoms " + wavy.get_atom_label(atom_idx_1) + " (" + std::to_string(atom_idx_1) + ") and " + wavy.get_atom_label(atom_idx_2) + " (" + std::to_string(atom_idx_2) + ") with " + std::to_string(gridpoints) + " gridpoints and " + std::to_string(padding) + " Angstrom padding.";
            std::cout << out << std::endl;
            get1DGridData(wavy, aux_basis, atom_idx_1, atom_idx_2, gridpoints, padding);
            exit(0);
        }
    }

    // SALTED predicts a density from atom positions.  Historically its
    // calculation path read those positions through -wfn, although the CLI
    // also documents -xyz for SALTED input.  Accept the documented form and
    // preserve an explicitly supplied -wfn when both options are present.
    if (SALTED && wfn.empty() && !xyz_file.empty())
    {
        wfn = xyz_file;
        if (debug)
            log_file << "Using -xyz input as the SALTED structure: " << wfn << endl;
    }
};

namespace {
    // Captured during static initialisation, so it still refers to the console after
    // run_app() has redirected std::cout into NoSpherA2.log. A function local static
    // would only be captured on the first error, by which time the redirect has happened
    // and every error message would end up in the log file instead of on the console.
    std::streambuf *const initial_coutbuf = std::cout.rdbuf();
    std::streambuf *original_coutbuf()
    {
        return initial_coutbuf;
    }
}

void options::look_for_debug(int &argc, char **argv)
{
    // This loop figures out command line options
    for (int i = 0; i < argc; i++)
    {
        std::string temp = argv[i];
        arguments.push_back(temp);
        if (temp.find("-") > 0)
            continue;
        else if (temp == "-v" || temp == "-v2" || temp == "-debug")
            std::cout << "Turning on verbose mode!" << std::endl, debug = true,
            ProgressBar::report_counts = true;
        else if (temp == "--h" || temp == "-h" || temp == "-help" || temp == "--help")
        {
            // run_app() has already pointed std::cout at NoSpherA2.log, and the
            // exit() below unwinds nothing, so without this the entire help text
            // is written to the log file and the console stays empty.
            std::cout.rdbuf(original_coutbuf());
            std::cout << NoSpherA2_message() << help_message << build_date << std::endl;
            exit(0);
        }
    }
};

bool is_nan(const double &in)
{
    return in != in;
};
bool is_nan(const float &in)
{
    return in != in;
};
bool is_nan(const long double &in)
{
    return in != in;
};
bool is_nan(const cdouble &in)
{
    return in != in;
};

bool ends_with(const std::string &str, const std::string &suffix)
{
    if (str.length() >= suffix.length())
    {
        return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
    }
    return false;
}

cdouble hypergeometric(double a, double b, double c, cdouble x)
{
    const double TOLERANCE = 1.0e-10;
    cdouble term = a * b * x / c;
    cdouble value = 1.0 + term;
    int n = 1;

    while (std::abs(term) > TOLERANCE)
    {
        a++, b++, c++, n++;
        term *= a * b * x / c / static_cast<double>(n);
        value += term;
    }

    return value;
}

double hypergeometric(double a, double b, double c, double x)
{
    const double TOLERANCE = 1.0e-10;
    double term = a * b * x / c;
    double value = 1.0 + term;
    int n = 1;

    while (std::abs(term) > TOLERANCE)
    {
        a++, b++, c++, n++;
        term *= a * b * x / c / n;
        value += term;
    }

    return value;
}

double double_from_string_with_esd(std::string in)
{
    if (in.find('(') == std::string::npos)
        return stod(in);
    else
        return stod(in.substr(0, in.find('(')));
}

std::string trim(const std::string &s)
{
    if (s == "")
        return "";
    auto start = s.begin();
    while (start != s.end() && std::isspace(*start))
    {
        start++;
    }

    auto end = s.end();
    do
    {
        end--;
    } while (std::distance(start, end) > 0 && std::isspace(*end));

    return std::string(start, end + 1);
}

int CountWords(const char *str)
{
    if (str == NULL)
        return -1;

    bool inSpaces = true;
    int numWords = 0;

    while (*str != '\0')
    {
        if (std::isspace(*str))
        {
            inSpaces = true;
        }
        else if (inSpaces)
        {
            numWords++;
            inSpaces = false;
        }

        str++;
    }

    return numWords;
};

void print_duration(std::ostream &file, const std::string &description, const std::chrono::microseconds &duration, std::optional<std::chrono::microseconds> total_duration = std::nullopt)
{
    auto mins = std::chrono::duration_cast<std::chrono::minutes>(duration);
    auto secs = std::chrono::duration_cast<std::chrono::seconds>(duration) % 60;
    auto millisecs = std::chrono::duration_cast<std::chrono::milliseconds>(duration) % 1000;

    file << std::setw(35) << std::left << std::setfill(' ') << description << ": " << std::right
        << std::setw(2) << std::setfill('0') << mins.count() << ":"
        << std::setw(2) << std::setfill('0') << secs.count() << ":"
        << std::setw(3) << std::setfill('0') << millisecs.count();
    if (total_duration.has_value())
    {
        double percentage = (double(duration.count()) / total_duration->count()) * 100.0;
        file << "  (" << std::fixed << std::setprecision(2) << percentage << "%)";
    };
    file << std::endl;
    // Disable setfill 0 again
    file << std::setfill(' ');
}

void write_timing_to_file(std::ostream &file,
    std::vector<_time_point> time_points,
    std::vector<std::string> descriptions)
{
    using namespace std;
    // Check if either vector is empty
    if (time_points.empty() || descriptions.empty())
    {
        file << "Error: Empty vector passed to write_timing_to_file" << endl;
        return;
    }
    std::chrono::microseconds total_time = std::chrono::duration_cast<std::chrono::microseconds>(time_points.back() - time_points.front());
    file << "\n\n----------------------------- Time Breakdown! -----------------------------" << endl;
    file << "                                     mm:ss:ms" << endl;
    // Time_points.size()-1 because we are comparing each time point to the next one meaning we need to stop at the second to last element
    for (int i = 0; i < time_points.size() - 1; i++)
    {
        std::chrono::microseconds dur = std::chrono::duration_cast<std::chrono::microseconds>(time_points[i + 1] - time_points[i]);
        print_duration(file, "... for " + descriptions[i], dur, total_time);
    }
    print_duration(file, "Total Time", total_time);
    file << "---------------------------------------------------------------------------" << endl;
}

void remove_empty_elements(svec &input, const std::string &empty)
{
    for (int i = (int)input.size() - 1; i >= 0; i--)
        if (input[i] == empty || input[i] == "")
            input.erase(input.begin() + i);
}

std::chrono::high_resolution_clock::time_point get_time()
{
    // gets the current time using std chrono library
    std::chrono::high_resolution_clock::time_point time = std::chrono::high_resolution_clock::now();
    return time;
}

long long int get_musec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
    // gets the time difference in microseconds
    std::chrono::microseconds musec = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    return musec.count();
}

long long int get_msec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
    // gets the time difference in milliseconds
    std::chrono::milliseconds msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    return msec.count();
}

long long int get_sec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
    // gets the time difference in seconds
    std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);
    return sec.count();
}

const int shell2function(const int &type, const int &prim)
{
    switch (type)
    {
    case (-5):
        return -32 + prim;
    case (-4):
        return -21 + prim;
    case (-3):
        return -12 + prim;
    case (-2):
        return -5 + prim;
    case (-1):
        return 1 + prim;
    case (0):
        return 1;
    case (1):
        return 2 + prim;
    case (2):
        return 5 + prim;
    case (3):
        if (prim == 0)
            return 11;
        if (prim == 1)
            return 12;
        if (prim == 2)
            return 13;
        if (prim == 3)
            return 17;
        if (prim == 4)
            return 14;
        if (prim == 5)
            return 15;
        if (prim == 6)
            return 18;
        if (prim == 7)
            return 19;
        if (prim == 8)
            return 16;
        if (prim == 9)
            return 20;
        break;
    case (4):
        return 21 + prim;
    case (5):
        return 36 + prim;
    default:
        return 0;
    }
    return 0;
}

bool open_file_dialog(std::filesystem::path &path, bool debug, std::vector <std::string> filter, const std::string &current_path) {
#ifdef _WIN32
    char filename[1024];

    OPENFILENAMEA ofn;
    ZeroMemory(&filename, sizeof(filename));
    ZeroMemory(&ofn, sizeof(ofn));
    ofn.lStructSize = sizeof(ofn);
    ofn.lpstrFilter = "Known Formats\0*.gbw;*.fchk;*.wfx;*.wfn;*.ffn;*.molden;*.molden.input;*.xtb\0gbw Files\0*.gbw\0wfn Files\0*.wfn\0wfx Files\0*.wfx\0ffn Files\0*.ffn\0cube Files\0*.cub;*.cube\0xtb Files\0*.xtb\0Any File\0*\0";
    ofn.lpstrFile = filename;
    ofn.nMaxFile = 1024;
    ofn.lpstrTitle = "Select a File";
    ofn.Flags = OFN_DONTADDTORECENT | OFN_FILEMUSTEXIST;

    if (GetOpenFileNameA(&ofn)) {
        if (debug) std::cout << "You chose the file \"" << filename << "\"\n";
        auto p = std::filesystem::path(filename);
        if (exists(p)) {
            path = p;
            return true;
        }
    }
    else
    {
        // All this stuff below is to tell you exactly how you messed up above.
        // Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
        switch (CommDlgExtendedError())
        {
        case CDERR_DIALOGFAILURE: std::cout << "CDERR_DIALOGFAILURE\n";   break;
        case CDERR_FINDRESFAILURE: std::cout << "CDERR_FINDRESFAILURE\n";  break;
        case CDERR_INITIALIZATION: std::cout << "CDERR_INITIALIZATION\n";  break;
        case CDERR_LOADRESFAILURE: std::cout << "CDERR_LOADRESFAILURE\n";  break;
        case CDERR_LOADSTRFAILURE: std::cout << "CDERR_LOADSTRFAILURE\n";  break;
        case CDERR_LOCKRESFAILURE: std::cout << "CDERR_LOCKRESFAILURE\n";  break;
        case CDERR_MEMALLOCFAILURE: std::cout << "CDERR_MEMALLOCFAILURE\n"; break;
        case CDERR_MEMLOCKFAILURE: std::cout << "CDERR_MEMLOCKFAILURE\n";  break;
        case CDERR_NOHINSTANCE: std::cout << "CDERR_NOHINSTANCE\n";     break;
        case CDERR_NOHOOK: std::cout << "CDERR_NOHOOK\n";          break;
        case CDERR_NOTEMPLATE: std::cout << "CDERR_NOTEMPLATE\n";      break;
        case CDERR_STRUCTSIZE: std::cout << "CDERR_STRUCTSIZE\n";      break;
        case FNERR_BUFFERTOOSMALL: std::cout << "FNERR_BUFFERTOOSMALL\n";  break;
        case FNERR_INVALIDFILENAME: std::cout << "FNERR_INVALIDFILENAME\n"; break;
        case FNERR_SUBCLASSFAILURE: std::cout << "FNERR_SUBCLASSFAILURE\n"; break;
        default: std::cout << "You cancelled.\n";
        }
    }
    return false;
#else
    std::string command;
    bool use_zenity = (system("which zenity > /dev/null 2>&1") == 0);
    bool use_kdialog = false;

    if (use_zenity) {
        command = "zenity --file-selection --title=\"Select a file to load\" --filename=\"";
        command += current_path;
        command += "/\"";
        for (const auto &f : filter) {
            command += " --file-filter='";
            command += f;
            command += "'";
        }
        command += " 2> /dev/null";
    }
    else {
        use_kdialog = (system("which kdialog > /dev/null 2>&1") == 0);
        if (use_kdialog) {
            command = "kdialog --getopenfilename \"";
            command += current_path;
            command += "/\" '";
            for (const auto &f : filter) {
                command += f;
                command += " ";
            }
            command += "'";
            command += " --title \"Select a file to load\" 2> /dev/null";
        }
        else {
            std::cout << "No suitable file dialog tool found (zenity/kdialog)." << std::endl;
            std::cout << "Please enter the full path to the file: " << std::flush;
            std::string input_path;
            std::getline(std::cin, input_path);

            // Trim leading/trailing whitespace
            input_path.erase(0, input_path.find_first_not_of(" \t\n\r"));
            input_path.erase(input_path.find_last_not_of(" \t\n\r") + 1);

            if (input_path.empty()) {
                if (debug) std::cout << "No path entered." << std::endl;
                return false;
            }

            path = input_path;
            if (std::filesystem::exists(path)) {
                if (debug) std::cout << "Selected file via manual input: " << path << std::endl;
                return true;
            }
            else {
                std::cerr << "Error: File not found at path: " << path << std::endl;
                return false;
            }
        }
    }

    if (debug) {
        std::cout << "Executing command: " << command << std::endl;
    }

    FILE *f = popen(command.c_str(), "r");
    if (!f) {
        std::cerr << "Error: Failed to execute file dialog command." << std::endl;
        return false;
    }

    std::string file_str;
    char buffer[1024];
    if (fgets(buffer, sizeof(buffer), f) == NULL) {
        if (debug) std::cout << "File selection cancelled." << std::endl;
        pclose(f);
        return false;
    }
    file_str = buffer;

    int pclose_status = pclose(f);
    if (pclose_status != 0) {
        if (debug) std::cout << (use_zenity ? "Zenity" : "KDialog") << " returned non-zero status: " << pclose_status << ". User might have cancelled." << std::endl;
        // This can happen on cancel, so we check if a file was actually returned.
        if (file_str.empty()) return false;
    }

    // Clean up the path string which might have a newline
    file_str.erase(file_str.find_last_not_of(" \n\r\t") + 1);

    if (file_str.empty()) {
        if (debug) std::cout << "File selection cancelled or returned empty path." << std::endl;
        return false;
    }

    path = file_str;
    if (debug) {
        std::cout << "Selected file: " << path << std::endl;
    }

    return std::filesystem::exists(path);
#endif
};

bool save_file_dialog(std::filesystem::path &path, bool debug, const std::vector<std::string> &endings, const std::string &filename_given, const std::string &current_path) {
#ifdef _WIN32
    constexpr size_t MAX_FILENAME_SIZE = 4096;
    std::vector<char> filename_buf(MAX_FILENAME_SIZE);

    OPENFILENAMEA sfn;
    ZeroMemory(&filename_buf[0], filename_buf.size());
    ZeroMemory(&sfn, sizeof(sfn));
    sfn.lStructSize = sizeof(sfn);
    sfn.hwndOwner = NULL;  // If you have a window to center over, put its HANDLE here
    sfn.lpstrFilter = "Known Formats\0*.gbw;*.fchk;*.wfx;*.wfn;*.ffn;*.molden;*.molden.input;*.xtb\0gbw Files\0*.gbw\0wfn Files\0*.wfn\0wfx Files\0*.wfx\0ffn Files\0*.ffn\0cube Files\0*.cub;*.cube\0xtb Files\0*.xtb\0Any File\0*\0";
    sfn.lpstrFile = filename_buf.data();
    sfn.nMaxFile = static_cast<DWORD>(filename_buf.size());
    sfn.lpstrTitle = "Select a File for saving!";
    sfn.Flags = OFN_DONTADDTORECENT;
    bool end = false;
    while (!end) {
        if (GetSaveFileNameA(&sfn)) {
            std::string chosen(filename_buf.data());
            if (debug) std::cout << "You chose the file \"" << chosen << "\"\n";
            if (exists(std::filesystem::path(chosen))) {
                std::cout << chosen << " exists, do you want to overwrite it?";
                if (yesno()) {
                    path = chosen;
                    bool found = false;
                    for (int i = 0; i < endings.size(); i++) if (path.extension() == endings[i]) found = true;
                    if (found) end = true;
                }
                else return false;
            }
            else {
                path = chosen;
                bool found = false;
                for (int i = 0; i < endings.size(); i++) if (path.extension() == endings[i]) found = true;
                if (found) end = true;
            }
        }
        else
        {
            // All this stuff below is to tell you exactly how you messed up above.
            // Once you've got that fixed, you can often (not always!) reduce it to a 'user cancelled' assumption.
            switch (CommDlgExtendedError())
            {
            case CDERR_DIALOGFAILURE: std::cout << "CDERR_DIALOGFAILURE\n";   break;
            case CDERR_FINDRESFAILURE: std::cout << "CDERR_FINDRESFAILURE\n";  break;
            case CDERR_INITIALIZATION: std::cout << "CDERR_INITIALIZATION\n";  break;
            case CDERR_LOADRESFAILURE: std::cout << "CDERR_LOADRESFAILURE\n";  break;
            case CDERR_LOADSTRFAILURE: std::cout << "CDERR_LOADSTRFAILURE\n";  break;
            case CDERR_LOCKRESFAILURE: std::cout << "CDERR_LOCKRESFAILURE\n";  break;
            case CDERR_MEMALLOCFAILURE: std::cout << "CDERR_MEMALLOCFAILURE\n"; break;
            case CDERR_MEMLOCKFAILURE: std::cout << "CDERR_MEMLOCKFAILURE\n";  break;
            case CDERR_NOHINSTANCE: std::cout << "CDERR_NOHINSTANCE\n";     break;
            case CDERR_NOHOOK: std::cout << "CDERR_NOHOOK\n";          break;
            case CDERR_NOTEMPLATE: std::cout << "CDERR_NOTEMPLATE\n";      break;
            case CDERR_STRUCTSIZE: std::cout << "CDERR_STRUCTSIZE\n";      break;
            case FNERR_BUFFERTOOSMALL: std::cout << "FNERR_BUFFERTOOSMALL\n";  break;
            case FNERR_INVALIDFILENAME: std::cout << "FNERR_INVALIDFILENAME\n"; break;
            case FNERR_SUBCLASSFAILURE: std::cout << "FNERR_SUBCLASSFAILURE\n"; break;
            default: std::cout << "You cancelled.\n";
            }
            return false;
        }
    }
#else
    std::string command;
    command = "zenity --file-selection --title=\"Select where to save\" --filename=\"";
    command += current_path;
    command += filename_given;
    command += "/\" --save --confirm-overwrite 2> /dev/null";
    bool end = false;
    while (!end) {
        FILE *f = popen(command.c_str(), "r");
        if (!f) {
            std::cout << "ERROR" << std::endl;
            return false;
        }
        std::string file;
        char buf[256];
        while (fgets(buf, sizeof(buf), f)) {
            file += buf;
        }
        if (file.empty())
            return false;
        if (debug)
            std::cout << "Filename: " << file << std::endl;
        path = file;
        std::stringstream ss(path);
        std::string name = path.string();
        getline(ss, name);
        if (debug) std::cout << "Path: " << path << std::endl;
        if (pclose(f) != 0) std::cout << "Zenity returned non zero, whatever that means..." << std::endl;
        bool found = false;
        for (int i = 0; i < endings.size(); i++)
            if (path.string().find(endings[i]) != std::string::npos)
                found = true;
        if (found)
            end = true;
    }
#endif
    return true;
};

const int sht2nbas(const int &type)
{
    const int st2bas[9]{ 1, 3, 6, 10, 15, 21, 28, 36 };
    const int nst2bas[9]{ 17,15,13,11, 9, 7, 5, 4, 1 };
    if (type >= 0)
        return st2bas[type];
    else
        return nst2bas[8 + type];
};

char asciitolower(char in)
{
    if (in <= 'Z' && in >= 'A')
        return in - ('Z' - 'z');
    return in;
}

int vec_sum(const bvec &in)
{
    int sum = 0;
    for (bool val : in)
    {
        sum += val;
    }
    return sum;
}

int vec_sum(const ivec &in)
{
    int sum = 0;
    for (int val : in)
    {
        sum += val;
    }
    return sum;
}

double vec_sum(const vec &in)
{
    double sum = 0.0;
    for (double val : in)
    {
        sum += val;
    }
    return sum;
}

cdouble vec_sum(const cvec &in)
{
    cdouble res = 0.0;
    for (int i = 0; i < in.size(); i++)
        res += in[i];
    return res;
}

double vec_length(const vec &in)
{
    double sum = 0.0;
    for (double val : in)
    {
        sum += val * val;
    }
    return sqrt(sum);
}


void error_check(const bool condition, const std::source_location loc, const std::string &error_message, std::ostream &log_file)
{
    if (!condition)
    {
        log_file << "Error in " << loc.function_name() << "\n\t\tat: " << loc.file_name() << " line: " << loc.line() << "\n\t\t\t" << error_message << std::endl;
        log_file.flush();
        std::cout.rdbuf(original_coutbuf()); // reset to standard output again
        std::cout << "Error in " << loc.function_name() << " at: " << loc.file_name() << " line: " << loc.line() << " " << error_message << std::endl;
        exit(-1);
    }
};
void not_implemented(const std::source_location loc, const std::string &error_message, std::ostream &log_file)
{
    log_file << loc.function_name() << "\n\t\tat: " << loc.file_name() << " line: " << loc.line() << "\n\t\t\t" << error_message << " not yet implemented!" << std::endl;
    log_file.flush();
    std::cout.rdbuf(original_coutbuf()); // reset to standard output again
    std::cout << "Error in " << loc.function_name() << " at: " << loc.file_name() << " : " << loc.line() << " " << error_message << " not yet implemented!" << std::endl;
    exit(-1);
};

void sha::sha256_transform(uint32_t state[8], const uint8_t block[64])
{
    uint32_t w[64];
    uint32_t a, b, c, d, e, f, g, h;

    for (int i = 0; i < 16; i++)
    {
        w[i] = (block[i * 4] << 24) | (block[i * 4 + 1] << 16) |
            (block[i * 4 + 2] << 8) | (block[i * 4 + 3]);
    }

    for (int i = 16; i < 64; i++)
    {
        w[i] = SIG1(w[i - 2]) + w[i - 7] + SIG0(w[i - 15]) + w[i - 16];
    }

    a = state[0];
    b = state[1];
    c = state[2];
    d = state[3];
    e = state[4];
    f = state[5];
    g = state[6];
    h = state[7];

    for (int i = 0; i < 64; i++)
    {
        uint32_t temp1 = h + EP1(e) + CH(e, f, g) + k[i] + w[i];
        uint32_t temp2 = EP0(a) + MAJ(a, b, c);
        h = g;
        g = f;
        f = e;
        e = d + temp1;
        d = c;
        c = b;
        b = a;
        a = temp1 + temp2;
    }

    state[0] += a;
    state[1] += b;
    state[2] += c;
    state[3] += d;
    state[4] += e;
    state[5] += f;
    state[6] += g;
    state[7] += h;
}

// SHA-256 update function
void sha::sha256_update(uint32_t state[8], uint8_t buffer[64], const uint8_t *data, size_t len, uint64_t &bitlen)
{
    for (size_t i = 0; i < len; i++)
    {
        size_t buf_idx = (bitlen / 8) % 64;
        if (buf_idx >= 64) {
            std::cerr << "Buffer overflow detected in sha256_update!" << std::endl;
            return;
        }
        buffer[buf_idx] = data[i];
        bitlen += 8;
        if (bitlen % 512 == 0)
        {
            sha256_transform(state, buffer);
        }
    }
}

// SHA-256 padding and final hash computation
void sha::sha256_final(uint32_t state[8], uint8_t buffer[64], uint64_t bitlen, uint8_t hash[32])
{
    size_t i = bitlen / 8 % 64;

    if (i >= 64) {
        std::cerr << "Buffer overflow detected in sha256_final (initial)!" << std::endl;
        return;
    }
    buffer[i++] = 0x80;
    if (i > 56)
    {
        while (i < 64)
        {
            if (i >= 64) {
                std::cerr << "Buffer overflow detected in sha256_final (padding)!" << std::endl;
                return;
            }
            buffer[i++] = 0x00;
        }
        sha256_transform(state, buffer);
        i = 0;
    }

    while (i < 56)
    {
        if (i >= 64) {
            std::cerr << "Buffer overflow detected in sha256_final (final padding)!" << std::endl;
            return;
        }
        buffer[i++] = 0x00;
    }

    // memcpy to buffer + 56 is safe as buffer is 64 bytes
    bitlen = custom_bswap_64(bitlen);
    memcpy(buffer + 56, &bitlen, 8);
    sha256_transform(state, buffer);

    for (i = 0; i < 8; i++)
    {
        constexpr size_t HASH_SIZE = 32;
        const size_t off = static_cast<size_t>(i) * 4;
        if (off + 4 > HASH_SIZE) {
            std::cerr << "Buffer overflow detected in sha256_final (hash write)!" << std::endl;
            return;
        }
        // convert to big-endian value locally and copy 4 bytes
        uint32_t be = custom_bswap_32(state[i]);
        std::memcpy(hash + off, &be, sizeof(be));
    }
}

// Function to calculate SHA-256 hash
std::string sha::sha256(const std::string &input)
{
    uint32_t state[8] = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
        0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19 };

    uint8_t buffer[64] = { 0 };
    uint8_t hash[32];
    uint64_t bitlen = 0;

    sha256_update(state, buffer, reinterpret_cast<const uint8_t *>(input.c_str()), input.length(), bitlen);
    sha256_final(state, buffer, bitlen, hash);

    std::stringstream ss;
    for (int i = 0; i < 32; i++)
    {
        ss << std::hex << std::setw(2) << std::setfill('0') << (int)hash[i];
    }

    return ss.str();
}

/*
#ifdef _WIN32
// Function to convert a narrow string to a wide string
std::wstring s2ws(const std::string &s)
{
    int len;
    int slength = (int)s.length() + 1;
    len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
    std::wstring r(len, L'\0');
    MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, &r[0], len);
    return r;
}
*/

bool ProgressBar::report_counts = false;

ProgressBar::~ProgressBar()
{
    progress_ = 100.0f;
    write_progress();
    stream_ << std::endl;
    if (report_counts)
    {
        // updates is how often callers reported an item; writes is how often that
        // needed the omp critical. Every item used to take the lock, so writes
        // being far below updates is the whole point of the change.
        const unsigned long long u = update_calls_.load(), w = bar_writes_.load();
        stream_ << "[progress] " << status_text_ << ": " << u << " updates, "
                << w << " serialised writes";
        if (u > 0)
            stream_ << "  (" << (w * 100.0 / u) << "% of calls took the lock, "
                    << (u > w ? u / (w ? w : 1) : 1) << "x fewer barriers)";
        stream_ << std::endl;
    }
#ifdef _WIN32
    if (taskbarList_)
    {
        taskbarList_->SetProgressState(GetConsoleWindow(), TBPF_NOPROGRESS);
        taskbarList_->Release();
        taskbarList_ = nullptr;
    }
#endif
}

void ProgressBar::write_progress()
{
    // No need to write once progress is 100%
    if (progress_ > 100.0f)
        return;

    // Move cursor to the first position on the same line
    // Check if os is a file stream
    if (dynamic_cast<std::filebuf *>(stream_.rdbuf()))
    {
        stream_.seekp(linestart); // Is a file stream
    }
    else
    {
        stream_ << "\r" << std::flush; // Is not a file stream
    }

    // Start bar
    stream_ << "[";

    const auto completed = static_cast<size_t>(progress_ * static_cast<float>(bar_width_) / 100.0);
    for (size_t i = 0; i <= completed; i++)
    {
        stream_ << fill_;
    }

    // End bar
    if (((progress_ < 100.0f) ? progress_ : 100.0f) == 100)
    {
        stream_ << "] 100% " << std::flush;
#ifdef _WIN32
        if (taskbarList_)
        {
            taskbarList_->SetProgressValue(GetConsoleWindow(), 100, 100);
            taskbarList_->SetProgressState(GetConsoleWindow(), TBPF_NOPROGRESS);
        }
#endif
        return;
    }

    stream_ << std::flush;

#ifdef _WIN32
    // Update taskbar progress
    if (taskbarList_)
    {
        taskbarList_->SetProgressValue(GetConsoleWindow(), static_cast<ULONGLONG>(progress_), 100);
    }
#endif
}

#ifdef _WIN32
void ProgressBar::initialize_taskbar_progress()
{
    if (SUCCEEDED(CoInitialize(nullptr)))
    {
        if (SUCCEEDED(CoCreateInstance(CLSID_TaskbarList, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&taskbarList_))))
        {
            taskbarList_->HrInit();
            taskbarList_->SetProgressState(GetConsoleWindow(), TBPF_NORMAL);
        }
    }
}
#endif

void convert_tonto_XCW_lambda_steps(const std::string &str, const std::string &lambda_step, bool debug, options &opt) {
    double lambda = 0.0;
    const double ls = stod(lambda_step);
    std::string jobname, line;
    std::filesystem::path stdout_file = str;
    err_checkf(std::filesystem::exists(stdout_file), "couldn't open or find " + stdout_file.string() + ", leaving", std::cout);
    std::ifstream rf(stdout_file.string().c_str(), std::ios::in);
    rf.seekg(0);
    while (rf.good() && line.find("Name ...") == std::string::npos) {
        getline(rf, line);
    }
    jobname = split_string<std::string>(line, " ")[2];
    std::cout << "Conervting XCW wavefunctions with lambda step " + std::to_string(ls) + " and jobname: " + jobname << std::endl;

    //iterate over files in the same folder looking for .orbital_energies,restricted or .MO_energies,r and .molecular_orbitals,restricted or .MOs,r at the given lambda value
    //example filename: NiP3.molecular_orbitals,lambda=0.000000,restricted
    //make lambda value into a string with 6 decimal places
    std::stringstream ss;
    ss << std::fixed << std::setprecision(6) << lambda;
    std::string formatted_lambda = ss.str();
    std::filesystem::path energies_file = stdout_file.parent_path() / (jobname + ".orbital_energies,lambda=" + formatted_lambda + ",restricted");
    if (!std::filesystem::exists(energies_file))
        energies_file = stdout_file.parent_path() / (jobname + ".MO_energies,lambda=" + formatted_lambda + ",r");
    if (!std::filesystem::exists(energies_file))
        energies_file = stdout_file.parent_path() / (jobname + ".orbital_energies,lambda=" + formatted_lambda + ",alpha");
    if (!std::filesystem::exists(energies_file))
        energies_file = stdout_file.parent_path() / (jobname + ".MO_energies,lambda=" + formatted_lambda + ",a");

    std::filesystem::path orbitals_file = stdout_file.parent_path() / (jobname + ".molecular_orbitals,lambda=" + formatted_lambda + ",restricted");
    if (!std::filesystem::exists(orbitals_file))
        orbitals_file = stdout_file.parent_path() / (jobname + ".MOs,lambda=" + formatted_lambda + ",r");
    if (!std::filesystem::exists(orbitals_file))
        orbitals_file = stdout_file.parent_path() / (jobname + ".molecular_orbitals,lambda=" + formatted_lambda + ",alpha");
    if (!std::filesystem::exists(orbitals_file))
        orbitals_file = stdout_file.parent_path() / (jobname + ".MOs,lambda=" + formatted_lambda + ",a");

    while (std::filesystem::exists(energies_file) && std::filesystem::exists(orbitals_file)) {
        const std::filesystem::path of = orbitals_file;
        const std::filesystem::path ef = energies_file;

        err_checkf(of.string() != "", "Orbitals file name is empty?!", std::cout);
        err_checkf(ef.string() != "", "Energy file name is empty?!", std::cout);
        err_checkf(std::filesystem::exists(of), "couldn't open or find " + of.string() + ", leaving", std::cout);
        err_checkf(std::filesystem::exists(ef), "couldn't open or find " + ef.string() + ", leaving", std::cout);
        std::cout << "lambda = " + std::to_string(lambda) + "..." << std::flush;

        std::vector<WFN> wavy;
        wavy.emplace_back(e_origin::tonto);
        wavy.back().read_tonto(stdout_file, std::cout, debug, ef, of);
        std::filesystem::path basename = stdout_file.parent_path() / (jobname + "_l_" + formatted_lambda);
        wavy.back().write_wfn(basename.string() + ".wfn", debug, false);
        free_fchk(std::cout, basename.string() + ".fchk", "", wavy.back(), debug, true);

        if (opt.cif != "") {
            svec ka;
            int nr = 0;
            opt.groups[0] = { 0 };
            tsc_block<int, cdouble> result = calculate_scattering_factors<itsc_block, std::vector<WFN> &>(opt, wavy, std::cout, ka, nr);
            result.write_tscb_file(opt.cif, basename.string() + ".tscb");
        }

        lambda += ls;
        std::stringstream ss_l;
        ss_l << std::fixed << std::setprecision(6) << lambda;
        formatted_lambda = ss_l.str();
        energies_file = stdout_file.parent_path() / (jobname + ".orbital_energies,lambda=" + formatted_lambda + ",restricted");
        if (!std::filesystem::exists(energies_file))
            energies_file = stdout_file.parent_path() / (jobname + ".MO_energies,lambda=" + formatted_lambda + ",r");
        if (!std::filesystem::exists(energies_file))
            energies_file = stdout_file.parent_path() / (jobname + ".orbital_energies,lambda=" + formatted_lambda + ",alpha");
        if (!std::filesystem::exists(energies_file))
            energies_file = stdout_file.parent_path() / (jobname + ".MO_energies,lambda=" + formatted_lambda + ",a");

        orbitals_file = stdout_file.parent_path() / (jobname + ".molecular_orbitals,lambda=" + formatted_lambda + ",restricted");
        if (!std::filesystem::exists(orbitals_file))
            orbitals_file = stdout_file.parent_path() / (jobname + ".MOs,lambda=" + formatted_lambda + ",r");
        if (!std::filesystem::exists(orbitals_file))
            orbitals_file = stdout_file.parent_path() / (jobname + ".molecular_orbitals,lambda=" + formatted_lambda + ",alpha");
        if (!std::filesystem::exists(orbitals_file))
            orbitals_file = stdout_file.parent_path() / (jobname + ".MOs,lambda=" + formatted_lambda + ",a");
        std::cout << " .. done!" << std::endl;
    }
};
