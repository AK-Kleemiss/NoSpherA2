#include "ML_predict.h"
using namespace std;

struct Config {
    string saltedname;
    string saltedpath;
    vector<string> species;
    int Menv;
    double zeta;
    double regul;
    int ncut;
    bool sparsify;
    string qmcode;
    string qmbasis;
    string dfbasis;
    bool field;
    int Ntrain;
    double trainfrac;
};

Config parse_input() {
    // Implement this function to parse the input configuration
    Config inp;
    // Dummy implementation
    return inp;
}

void get_sparse_details(unordered_map<int, vec> vfps)

int main() {
    Config inp = parse_input();

    string saltedname = inp.saltedname;
    string saltedpath = inp.saltedpath;
    vector<string> species = inp.species;
    int Menv = inp.Menv;
    double zeta = inp.zeta;
    double reg = inp.regul;
    int ncut = inp.ncut;
    bool sparsify = inp.sparsify;

    vector<int> lmax, nmax;
    if (inp.qmcode == "pyscf") {
        tie(lmax, nmax) = basiset(get_aux_basis_name(inp.qmbasis));
    }
    else {
        tie(lmax, nmax) = basiset(inp.dfbasis);
    }

    vector<int> llist;
    vector<int> nlist;
    for (const auto& spe : species) {
        llist.push_back(lmax[spe]);
        for (int l = 0; l <= lmax[spe]; ++l) {
            nlist.push_back(nmax[(spe, l)]);
        }
    }
    int lmax_max = *max_element(llist.begin(), llist.end());



    unordered_map<int, vec> vfps;
    if (sparsify) {
        for (int lam = 0; lam <= lmax_max; ++lam) {
            vector<unsigned long> shape{};
            vec sparse_data{};
            bool fortran_order;
            string filename = saltedpath + "/equirepr_" + saltedname + "/fps" + to_string(ncut) + "-" + to_string(lam) + ".npy";
            npy::LoadArrayFromNumpy(filename, shape, fortran_order, sparse_data);
            vfps[lam] = sparse_data;
        }
    }

    auto [Vmat, Mspe, power_env_sparse] = get_feats_projs(species, lmax);

    int ntrain = static_cast<int>(inp.Ntrain * inp.trainfrac);
    vec weights{};
    string weight_filename;
    if (inp.field) {
        weight_filename = saltedpath + "/regrdir_" + saltedname + "_field/M" + to_string(Menv) + "_zeta" + to_string(zeta) + "/weights_N" + to_string(ntrain) + "_reg" + to_string(static_cast<int>(log10(reg))) + ".npy";
    }
    else {
        weight_filename = saltedpath + "/regrdir_" + saltedname + "/M" + to_string(Menv) + "_zeta" + to_string(zeta) + "/weights_N" + to_string(ntrain) + "_reg" + to_string(static_cast<int>(log10(reg))) + ".npy";
    }

    unordered_map<int, vec> vfps;
    vec weights{};
    vector<unsigned long> shape{};
    bool fortran_order;
    npy::LoadArrayFromNumpy(weight_filename, shape, fortran_order, weights);

    return {lmax, nmax, lmax_max, weights, power_env_sparse, Mspe, Vmat, vfps}
}

