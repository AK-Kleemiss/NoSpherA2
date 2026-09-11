#pragma once
#include "convenience.h"
#include "SALTED_utilities.h"

// The geometry-aid ("structure solution") pipeline Olex2 drives through
// -calc_featomic_descriptor(s), -classify_atoms(_list) and -geometry_aid_cutoff:
// the SOAP power-spectrum descriptor of a structure and the classifier that
// turns it into element probabilities. The flags only queue jobs in `options`;
// `run` executes them from run_app_impl.
namespace geometry_aid
{
    typedef std::vector<std::pair<std::filesystem::path, std::filesystem::path>> jobvec;

    SALTED_Utils::FeatomicHyperParameters hyperparameters(double cutoff_radius = 3.5);

    // geometry-aid classifier: the PCA and the three dense layers Olex2 used to
    // run in Python, same arithmetic. The weights come from
    // `geometry_aid_model.bin`, produced by `make_geometry_aid_bin.py`; the
    // `.npz` it replaces is a deflated ZIP and there is no zlib here.
    struct Model
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

    Model load_model(const std::filesystem::path& path);
    const Model& cached_model(const std::filesystem::path& path);
    vec classify_descriptor(const double* descriptor, size_t n_atoms, size_t n_features, const Model& m);

    pathvec read_structure_list(const std::filesystem::path& list_file);
    void write_descriptor(const std::filesystem::path& structure, const std::filesystem::path& out_path,
        const SALTED_Utils::FeatomicHyperParameters& hyperparams);
    int write_descriptors(const jobvec& jobs, const SALTED_Utils::FeatomicHyperParameters& hyperparams);
    void write_probabilities(const std::filesystem::path& structure, const std::filesystem::path& out_path,
        const std::filesystem::path& model_path, const SALTED_Utils::FeatomicHyperParameters& hyperparams);
    int write_probabilities(const jobvec& jobs, const std::filesystem::path& model_path,
        const SALTED_Utils::FeatomicHyperParameters& hyperparams);
    int run(const options& opt);
}
