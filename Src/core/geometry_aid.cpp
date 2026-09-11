#include "pch.h"
#include "geometry_aid.h"
#include "npy.h"

namespace
{
    template <typename T>
    void read_exact(std::istream& in, T* into, size_t count, const char* what)
    {
        in.read(reinterpret_cast<char*>(into), static_cast<std::streamsize>(count*sizeof(T)));
        err_checkf(static_cast<size_t>(in.gcount()) == count*sizeof(T),
            std::string("geometry-aid model truncated while reading ") + what, std::cout);
    }


    // Shared by the two batch flags: (structure, output) pairs, one unreadable
    // structure must not abort the rest. Returns the exit code.
    template <typename F>
    int run_jobs(const geometry_aid::jobvec& jobs, const char* tag, F job)
    {
        const auto started = std::chrono::steady_clock::now();
        size_t done = 0, failed = 0;
        for (const auto& [structure, out_path] : jobs)
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
                job(structure, out_path);
                std::cout << tag << " " << out_path.string() << " seconds="
                          << std::chrono::duration<double>(std::chrono::steady_clock::now() - one).count() << std::endl;
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
        return failed && !done ? 1 : 0;
    }
}

namespace geometry_aid
{
    // One definition of the geometry-aid hyperparameters, for both the single
    // structure and the batch flag. They must match what the models were
    // trained with, in geometry-aid/multi_layer_classifier/c_only_training.py
    // :: SOAP_HP, NOT the older values in geometry-aid/external_script.py.
    // A descriptor of the right length computed with the wrong settings is
    // rejected by nothing downstream. The feature count is the check: 11
    // species give 66 unique pairs and the length is
    // 66 * (max_radial+1)^2 * (max_angular+1) = 66 * 7^2 * 13 = 42,042.
    // SALTED is unaffected; it builds its own FeatomicHyperParameters from
    // config.nang1 / config.nang2 in SALTED_predictor.cpp.
    //
    // The cutoff radius is the one field that differs between the two shipped
    // model families: 3.5 for `c_only`, trained on all-carbon input, and 3.0
    // for `dirty`. Both descriptors are 42,042 long. Only `dirty` may be
    // iterated (predict, relabel, recompute); relabelled input is out of
    // distribution for `c_only`.
    SALTED_Utils::FeatomicHyperParameters hyperparameters(double cutoff_radius)
    {
        const std::vector<std::string> species{ "B", "C", "N", "O", "F", "Si", "P", "S", "Cl", "Br", "I" };

        // Diagnostic override, never for production output. A descriptor
        // computed at a different spline accuracy is not comparable with any
        // trained model, and nothing downstream rejects it -- hence the
        // warning.
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
        // Changing either invalidates every shipped model.
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
            .cutoff_radius = cutoff_radius,
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

    Model load_model(const std::filesystem::path& path)
    {
        std::ifstream in(path, std::ios::binary);
        err_checkf(in.good(), "Cannot open the geometry-aid model: " + path.string(), std::cout);

        char magic[8] = { 0 };
        read_exact(in, magic, 8, "the magic");
        err_checkf(std::string(magic, 8) == "GEOAID01",
            "This is not a GEOAID01 file. Regenerate it with "
            "make_geometry_aid_bin.py -- a stale .bin beside a newer .npz is "
            "exactly the mismatch the magic exists to catch.", std::cout);

        Model m;
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

    const Model& cached_model(const std::filesystem::path& path)
    {
        static std::map<std::string, Model> cache;
        const std::string key = path.string();
        auto found = cache.find(key);
        if (found == cache.end())
            found = cache.emplace(key, load_model(path)).first;
        return found->second;
    }

    // (n_atoms, n_classes) row-major probabilities. Summation order is not
    // numpy's, so the last bits differ from the Python route;
    // `bench_geometry_cpp.py` checks the argmax and the full ranking instead.
    vec classify_descriptor(const double* descriptor, size_t n_atoms,
        size_t n_features, const Model& m)
    {
        err_checkf(n_features == static_cast<size_t>(m.n_features),
            "The descriptor has " + std::to_string(n_features) + " features and "
            "the model expects " + std::to_string(m.n_features) + ". These come "
            "from different SOAP hyperparameters and the result would be "
            "meaningless rather than merely worse.", std::cout);

        const size_t k = static_cast<size_t>(m.n_components);
        vec projected(n_atoms*k, 0.0);

        // (x - mean) . C^T  ==  x . C^T  -  mean . C^T, and only the first term
        // touches the descriptor. Centring first destroys its sparsity: an
        // all-carbon .xyz -- what Olex2 sends on the first pass -- populates one
        // of the 66 species-pair blocks, so 637 of 42,042 entries are non-zero,
        // but `row[f] - mean[f]` is non-zero wherever the mean is and the skip
        // below never fires.
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

    void write_descriptor(const std::filesystem::path& structure,
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
    void write_probabilities(const std::filesystem::path& structure,
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

        const Model& model = cached_model(model_path);
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


    int write_descriptors(const jobvec& jobs, const SALTED_Utils::FeatomicHyperParameters& hyperparams)
    {
        return run_jobs(jobs, "DESCRIPTOR", [&hyperparams](const std::filesystem::path& structure, const std::filesystem::path& out) {
            write_descriptor(structure, out, hyperparams); });
    }

    int write_probabilities(const jobvec& jobs, const std::filesystem::path& model_path, const SALTED_Utils::FeatomicHyperParameters& hyperparams)
    {
        return run_jobs(jobs, "PROBABILITIES", [&](const std::filesystem::path& structure, const std::filesystem::path& out) {
            write_probabilities(structure, out, model_path, hyperparams); });
    }

    // One structure path per line; blank lines and lines beginning with '#' are ignored.
    pathvec read_structure_list(const std::filesystem::path& list_file)
    {
        err_checkf(std::filesystem::exists(list_file), "The structure list does not exist: " + list_file.string(), std::cout);
        pathvec structures;
        std::ifstream list(list_file);
        std::string line;
        while (getline_universal(list, line))
        {
            const std::string entry = trim(line);
            if (entry.empty() || entry[0] == '#') continue;
            structures.emplace_back(entry);
        }
        err_checkf(!structures.empty(), "The structure list is empty: " + list_file.string(), std::cout);
        return structures;
    }

    // What the flags queued, in one process; the exit code of the run.
    int run(const options& opt)
    {
        const SALTED_Utils::FeatomicHyperParameters hyperparams = hyperparameters(opt.geometry_aid_cutoff);
        jobvec descriptors, probabilities;
        if (opt.calc_featomic_descriptor)
        {
            err_checkf(!opt.wfn.empty(), "No structure given! -calc_featomic_descriptor needs -wfn.", std::cout);
            descriptors.emplace_back(opt.wfn, "descriptor.npy");
        }
        for (int i = 0; i < opt.featomic_structures.size(); i++)
            descriptors.emplace_back(opt.featomic_structures[i], opt.featomic_structures[i].string() + ".npy");
        if (!opt.classify_atoms_out.empty())
        {
            err_checkf(!opt.wfn.empty(), "No structure given! -classify_atoms needs -wfn.", std::cout);
            probabilities.emplace_back(opt.wfn, opt.classify_atoms_out);
        }
        for (int i = 0; i < opt.classify_structures.size(); i++)
            probabilities.emplace_back(opt.classify_structures[i], opt.classify_structures[i].string() + ".probs.npy");
        int code = 0;
        if (!descriptors.empty()) code |= write_descriptors(descriptors, hyperparams);
        if (!probabilities.empty()) code |= write_probabilities(probabilities, opt.geometry_aid_model, hyperparams);
        return code;
    }
}
