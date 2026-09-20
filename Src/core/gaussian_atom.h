#pragma once
#include "convenience.h"
#include "wfn_class.h"
#include "integrator.h"
#include "SALTED_utilities.h"
#include <map>
#include <optional>

using DensityFitting::Interaction_Energy;

//Atomic multipoles q_lm about the nucleus for l = 0..lmax, stored at [l*l+l+m] in the Racah normalisation
//sqrt(4pi/(2l+1)) Int rho r^l Y_lm (e bohr^l); q(0,0) is the electron count
struct Multipoles
{
    int lmax = 0;
    vec q;
    double operator()(const int l, const int m) const { return q[l * l + l + m]; }
};

//One atom whose density is a Gaussian expansion: the nucleus, its auxiliary shells and one coefficient per (shell, m).
//Everything is computed from those alone: point functions, electron count, centre multipoles, Fourier-Bessel form factors
//and the interaction energy with another Gaussian_Atom. Instances come from a Gaussian_Molecule, which slices its
//coefficients per atom.
class Gaussian_Atom
{
public:
    //by value so a temporary is moved in
    Gaussian_Atom(const atom& nucleus, vec coefficients) : at(nucleus), coefs(std::move(coefficients)), tab({ nucleus })
    {
        err_checkf(static_cast<int>(coefs.size()) == tab.n_coef, "Coefficient count does not match the auxiliary shells of " + at.get_label(), std::cout);
    }
    const atom& nucleus() const { return at; }
    //Nuclear charge less ECP electrons, what the density sees
    int Z() const { return tab.Z[0]; }
    d3 centre() const { return { tab.cx[0], tab.cy[0], tab.cz[0] }; }
    double rho(const d3& p) const { return tab(p[0], p[1], p[2], coefs.data()); }
    //Potential of this nucleus and its density
    double esp(const d3& p) const { return tab.esp(p[0], p[1], p[2], coefs.data()); }
    double values(const d3& p, d3& g, double& lap) const { return tab(p[0], p[1], p[2], coefs.data(), g[0], g[1], g[2], lap); }
    double hessian(const d3& p, d3& g, double* H) const { return tab(p[0], p[1], p[2], coefs.data(), g[0], g[1], g[2], H); }
    //Integral of the density, ECP electrons included
    double electrons() const { return calc_atomic_density({ at }, coefs)[0]; }
    //Net charge Z - electrons
    double charge() const { return at.get_charge() - electrons(); }
    //Multipoles about the nucleus, analytic: only the l shells carry the l moment
    Multipoles electrical_moments(const int lmax) const;
    //Form factor at k (bohr^-1): the Fourier-Bessel transform of the shells, i^l Y_lm(k^) (H/2)^l e^(-H^2/4a) / a^(l+3/2) per primitive
    cdouble form_factor(const double* k) const { return tab.fourier_atom(k[0], k[1], k[2], coefs.data(), 0); }
    //k_pt[3][nk] in bohr^-1, the layout of the hkl tables
    cvec scattering_factors(const vec2& k_pt) const;
    //Full interaction with B as for two molecules, see DensityFitting::interaction_energy
    Interaction_Energy interaction(const Gaussian_Atom& B, const double repulsion_K = 0.0, const int x_fun = 0) const;
    int n_aux() const { return static_cast<int>(coefs.size()); }
    const vec& coefficients() const { return coefs; }
    const aux_density_table& table() const { return tab; }
private:
    atom at;
    vec coefs;
    aux_density_table tab;
};

inline std::vector<d3> source_positions(const Gaussian_Atom& a) { return { a.centre() }; }
inline const char* source_name(const Gaussian_Atom&) { return "fitted atom density"; }

//The Gaussian_Atoms of one molecule: the RI fit of a wavefunction (-ri_fit <basis>), a SALTED prediction (-SALTED
//<model-dir>, no orbitals needed) or coefficients read from a .npy file. Keeps the auxiliary basis as a WFN, the flat
//coefficients and the flat aux_density_table the batch kernels consume, a UUID, and the interaction energies with other
//molecules by their UUID. A WFN may hold one or more of these.
class Gaussian_Molecule
{
public:
    //Coefficients from coef_file if given, else predicted with -SALTED, else the RI fit of wavy's orbitals in opt.aux_basis
    Gaussian_Molecule(const WFN& wavy, options& opt, const std::filesystem::path& coef_file = "");
    //An auxiliary basis with its coefficients already in hand; by value so temporaries are moved in
    Gaussian_Molecule(WFN aux_basis, vec coefficients);
    const std::string& uuid() const { return id; }
    const std::vector<Gaussian_Atom>& atoms() const { return ats; }
    const Gaussian_Atom& atom(const int a) const { return ats[a]; }
    int n_atoms() const { return tab.n_at; }
    int n_aux() const { return static_cast<int>(coefs.size()); }
    double rho(const d3& p) const { return tab(p[0], p[1], p[2], coefs.data()); }
    double esp(const d3& p) const { return tab.esp(p[0], p[1], p[2], coefs.data()); }
    double lap(const d3& p) const { return tab.lap(p[0], p[1], p[2], coefs.data()); }
    //Orbital-free PC07 estimate
    double eli(const d3& p) const { return tab.eli(p[0], p[1], p[2], coefs.data()); }
    //rho, its gradient and Laplacian from one shell loop
    double values(const d3& p, d3& g, double& lap) const { return tab(p[0], p[1], p[2], coefs.data(), g[0], g[1], g[2], lap); }
    //rho, its gradient and Hessian (row-major 3x3)
    double hessian(const d3& p, d3& g, double* H) const { return tab(p[0], p[1], p[2], coefs.data(), g[0], g[1], g[2], H); }
    //Electrons on every atom, ECP electrons included
    vec populations() const { return calc_atomic_density(*aux.get_atoms_ptr(), coefs); }
    double electrons() const;
    //Net charge sum Z - electrons
    double charge() const;
    //Per atom; without a partition the analytic centre moments, with one the moments of this density partitioned on a
    //Becke grid by that scheme (Hirshfeld from the promolecule, TFVC/MBIS/EMBIS from this density on the grid)
    std::vector<Multipoles> electrical_moments(const int lmax, const std::optional<DensityFitting::CHARGE_SCHEME> partition = std::nullopt) const;
    //Electrons per atom of this density partitioned by scheme
    vec populations(const DensityFitting::CHARGE_SCHEME scheme) const
    {
        vec pop;
        for (const Multipoles& M : electrical_moments(0, scheme)) pop.push_back(M.q[0]);
        return pop;
    }
    //sf[atom][k] of the listed atoms (all when empty) at k_pt[3][nk] in bohr^-1; a streaming caller passes its own bar
    cvec2 scattering_factors(const vec2& k_pt, ivec atoms = {}, ProgressBar* progress = nullptr) const;
    //Interaction with B, computed once per B.uuid() and kept
    const Interaction_Energy& interaction(const Gaussian_Molecule& B, const double repulsion_K = 0.0, const int x_fun = 0);
    const std::map<std::string, Interaction_Energy>& interactions() const { return energies; }
    //A new molecule at R x + t (bohr), the coefficients rotated with the atoms
    Gaussian_Molecule transformed(const vec2& R, const vec& t) const;
    const WFN& basis() const { return aux; }
    const vec& coefficients() const { return coefs; }
    const aux_density_table& table() const { return tab; }
    //True when the options ask for the fitted density in place of the orbitals: -ri_fit
    //<basis> on a wavefunction, or -SALTED on an xyz without orbitals
    static bool requested(const WFN& wavy, const options& opt)
    {
        return (opt.SALTED && wavy.get_nmo() == 0) || (opt.RI_FIT && !opt.aux_basis.empty() && wavy.get_nmo() > 0);
    }
private:
    void slice();
    WFN aux;
    vec coefs;
    aux_density_table tab;
    std::vector<Gaussian_Atom> ats;
    std::string id;
    std::map<std::string, Interaction_Energy> energies;
};

inline std::vector<d3> source_positions(const Gaussian_Molecule& m)
{
    std::vector<d3> c;
    for (const Gaussian_Atom& a : m.atoms()) c.push_back(a.centre());
    return c;
}
inline const char* source_name(const Gaussian_Molecule&) { return "fitted density"; }
//The batch forms of density_source.h on the fitted density's kernel: OpenMP over the shells per point, or the GPU
inline void calculate_density(const Gaussian_Molecule& m, const int n, const double* x, const double* y, const double* z, double* rho)
{
    calc_aux_density(m.table(), m.coefficients(), n, x, y, z, rho);
}
inline void calculate_density(const Gaussian_Molecule& m, const int n, const double* x, const double* y, const double* z, double* rho, double* gx, double* gy, double* gz, double* lap)
{
    calc_aux_density(m.table(), m.coefficients(), n, x, y, z, rho, gx, gy, gz, lap);
}
inline void calculate_hessian(const Gaussian_Molecule& m, const int n, const double* x, const double* y, const double* z, double* rho, double* gx, double* gy, double* gz, double* H)
{
    calc_aux_density(m.table(), m.coefficients(), n, x, y, z, rho, gx, gy, gz, nullptr, H);
}
