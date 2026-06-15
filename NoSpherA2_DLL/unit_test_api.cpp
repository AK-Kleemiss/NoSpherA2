#include "../Src/pch.h"
#include "unit_test_api.h"
#include "../Src/convenience.h"
#include "../Src/fchk.h"
#include "../Src/constants.h"
#include "../Src/wfn_class.h"
#include "../Src/tsc_block.h"
#include "../Src/atoms.h"
#include <cstring>
#include <cmath>
#include <thread>
#include <chrono>
#include <sstream>

extern "C" {

    UT_API double ut_array_length3(double x, double y, double z) {
        return std::hypot(x, y, z);
    }

    UT_API double ut_array_distance3(double x1, double y1, double z1,
                                     double x2, double y2, double z2) {
        d3 a{ x1, y1, z1 }, b{ x2, y2, z2 };
        return array_length(a, b);
    }

    UT_API void ut_vec_diff(double ax, double ay, double az,
                            double bx, double by, double bz, double out[3]) {
        d3 r = vec_diff({ ax, ay, az }, { bx, by, bz });
        out[0] = r[0]; out[1] = r[1]; out[2] = r[2];
    }

    UT_API void ut_vec_cross(double ax, double ay, double az,
                             double bx, double by, double bz, double out[3]) {
        d3 r = vec_cross({ ax, ay, az }, { bx, by, bz });
        out[0] = r[0]; out[1] = r[1]; out[2] = r[2];
    }

    UT_API double ut_vec_dot(double ax, double ay, double az,
                             double bx, double by, double bz) {
        return vec_dot({ ax, ay, az }, { bx, by, bz });
    }

    UT_API int ut_is_nan(double x) {
        // is_nan() uses "x != x" which /fp:fast optimises away; use std::isnan.
        return std::isnan(x) ? 1 : 0;
    }

    UT_API int ut_is_similar_rel(double a, double b, double tol) {
        return is_similar_rel(a, b, tol) ? 1 : 0;
    }

    UT_API int ut_is_similar_abs(double a, double b, double tol) {
        return is_similar_abs(a, b, tol) ? 1 : 0;
    }

    UT_API double ut_fast_exp_neg(double x) {
        return fast_exp_neg(x);
    }

    UT_API int ut_read_fchk_int(const char* line) {
        return read_fchk_integer(std::string(line));
    }

    UT_API double ut_read_fchk_dbl(const char* line) {
        return read_fchk_double(std::string(line));
    }

    UT_API int ut_ends_with(const char* str, const char* suffix) {
        return ends_with(std::string(str), std::string(suffix)) ? 1 : 0;
    }

    UT_API int ut_shrink_string(const char* input, char* output, int bufsize) {
        std::string s(input);
        std::string r = shrink_string(s);
        if (static_cast<int>(r.size()) >= bufsize)
            return -1;
        std::memcpy(output, r.c_str(), r.size() + 1);
        return static_cast<int>(r.size());
    }

    UT_API int ut_sha256(const char* input, char* output, int bufsize) {
        std::string r = sha::sha256(std::string(input));
        if (static_cast<int>(r.size()) >= bufsize)
            return -1;
        std::memcpy(output, r.c_str(), r.size() + 1);
        return static_cast<int>(r.size());
    }

    // --- Vector aggregates ---

    UT_API int ut_vec_sum_bool(const int* arr, int n) {
        if (n <= 0) return 0;
        bvec v(n);
        for (int i = 0; i < n; ++i) v[i] = (arr[i] != 0);
        return vec_sum(v);
    }

    UT_API int ut_vec_sum_int(const int* arr, int n) {
        if (n <= 0) return 0;
        ivec v(arr, arr + n);
        return vec_sum(v);
    }

    UT_API double ut_vec_sum_double(const double* arr, int n) {
        if (n <= 0) return 0.0;
        vec v(arr, arr + n);
        return vec_sum(v);
    }

    UT_API double ut_vec_length(const double* arr, int n) {
        if (n <= 0) return 0.0;
        vec v(arr, arr + n);
        return vec_length(v);
    }

    // --- Additional string utilities ---

    UT_API int ut_trim(const char* input, char* output, int bufsize) {
        std::string r = trim(std::string(input));
        if (static_cast<int>(r.size()) >= bufsize) return -1;
        std::memcpy(output, r.c_str(), r.size() + 1);
        return static_cast<int>(r.size());
    }

    UT_API char ut_asciitolower(char c) {
        return asciitolower(c);
    }

    UT_API double ut_double_from_esd(const char* s) {
        return double_from_string_with_esd(std::string(s));
    }

    UT_API double ut_decimal_precision_cif(const char* s) {
        std::string str(s);
        return get_decimal_precision_from_CIF_number(str);
    }

    // --- Basis-type helpers ---

    UT_API int ut_sht2nbas(int type) {
        return sht2nbas(type);
    }

    UT_API unsigned int ut_doublefactorial(int n) {
        return doublefactorial(n);
    }

    // --- Unit conversions ---

    UT_API double ut_bohr2ang(double inp)       { return constants::bohr2ang(inp); }
    UT_API double ut_ang2bohr(double inp)        { return constants::ang2bohr(inp); }
    UT_API double ut_cubic_bohr2ang(double inp)  { return constants::cubic_bohr2ang(inp); }
    UT_API double ut_cubic_ang2bohr(double inp)  { return constants::cubic_ang2bohr(inp); }

    // --- Constexpr math ---

    UT_API int    ut_const_abs(int x)              { return constants::const_abs(x); }
    UT_API double ut_constexpr_pow(double b, int e){ return constants::constexpr_pow(b, e); }
    UT_API double ut_constants_sqrt(double x)      { return constants::sqrt(x); }
    UT_API double ut_exp_approx(double x, int n)   { return constants::exp_approx(x, n); }
    UT_API double ut_log_approx(double x, int n)   { return constants::log_approx(x, n); }
    UT_API long long ut_factorial(int n)           { return static_cast<long long>(constants::ft_fun(n)); }

    // --- Orbital index mappings ---

    UT_API int ut_orca_2_pyscf(unsigned int l, int m_idx) {
        auto r = constants::orca_2_pySCF(l, m_idx);
        return r.has_value() ? static_cast<int>(r.value()) : -1;
    }

    UT_API unsigned int ut_type_2_nbo(int type) {
        return constants::type_2_nbo(type);
    }

    // --- Spherical Bessel functions ---

    UT_API double ut_bessel_j(int l, double x) {
        return bessel_first_kind(l, x);
    }

    // --- is_similar (pow-10 tolerance) ---

    UT_API int ut_is_similar_pow10(double a, double b, double tolerance) {
        return is_similar(a, b, tolerance) ? 1 : 0;
    }

    // --- shell2function ---

    UT_API int ut_shell2function(int type, int prim) {
        return shell2function(type, prim);
    }

    // --- CountWords ---

    UT_API int ut_count_words(const char* str) {
        return CountWords(str);
    }

    // --- shrink_string_to_atom ---

    UT_API int ut_shrink_string_to_atom(const char* input, int atom_number,
                                         char* output, int bufsize) {
        std::string s(input);
        std::string r = shrink_string_to_atom(s, atom_number);
        if (static_cast<int>(r.size()) >= bufsize)
            return -1;
        std::memcpy(output, r.c_str(), r.size() + 1);
        return static_cast<int>(r.size());
    }

    // --- split_string ---
    // Splits input by delim, writes up to max_out tokens into out_tokens[].
    // Returns total token count (may exceed max_out).

    UT_API int ut_split_string(const char* input, const char* delim,
                                char** out_tokens, int max_out, int tok_bufsize) {
        svec tokens = split_string<std::string>(std::string(input), std::string(delim));
        int written = std::min(static_cast<int>(tokens.size()), max_out);
        for (int i = 0; i < written; ++i) {
            strncpy_s(out_tokens[i], tok_bufsize, tokens[i].c_str(), tok_bufsize - 1);
        }
        return static_cast<int>(tokens.size());
    }

    // --- Timing helper ---

    UT_API long long ut_sleep_and_measure_us(int sleep_ms) {
        auto t0 = get_time();
        std::this_thread::sleep_for(std::chrono::milliseconds(sleep_ms));
        auto t1 = get_time();
        return get_musec(t0, t1);
    }

    // --- hypergeometric ---

    UT_API double ut_hypergeometric(double a, double b, double c, double x) {
        return hypergeometric(a, b, c, x);
    }

    // --- Atomic number → element symbol ---

    UT_API int ut_atnr2letter(int nr, char* output, int bufsize) {
        const char* sym = constants::atnr2letter(nr);
        int len = static_cast<int>(std::strlen(sym));
        if (len >= bufsize) return -1;
        std::memcpy(output, sym, len + 1);
        return len;
    }

    // --- Cartesian exponent vector from basis type ---

    UT_API void ut_type2vector(int index, int nx_out[1], int ny_out[1], int nz_out[1]) {
        int v[3];
        constants::type2vector(index, v);
        nx_out[0] = v[0];
        ny_out[0] = v[1];
        nz_out[0] = v[2];
    }

    // --- Gaussian normalization constant ---

    UT_API double ut_normgauss(int type, double exponent) {
        return constants::normgauss(type, exponent);
    }

    // --- Associated Legendre polynomial P_l^m(x) ---

    UT_API double ut_assoc_legendre(int l, int m, double x) {
        return constants::associated_legendre_polynomial(l, m, x);
    }

    // --- Cartesian → spherical coordinates ---

    UT_API void ut_cartesian_to_spherical(double x, double y, double z, double out3[3]) {
        std::vector<double> r = constants::cartesian_to_spherical(x, y, z);
        out3[0] = r[0];
        out3[1] = r[1];
        out3[2] = r[2];
    }

    // --- Median eigenvalue of 3x3 symmetric matrix ---

    UT_API double ut_get_lambda_1(const double a[9]) {
        return get_lambda_1(const_cast<double*>(a));
    }

    // -----------------------------------------------------------------------
    // WFN handle-based API
    // -----------------------------------------------------------------------

    UT_API void* ut_wfn_create() {
        return new WFN();
    }

    UT_API void ut_wfn_destroy(void* h) {
        delete static_cast<WFN*>(h);
    }

    UT_API int ut_wfn_push_atom(void* h, const char* label,
                                double x, double y, double z, int charge) {
        auto* w = static_cast<WFN*>(h);
        return w->push_back_atom(std::string(label), x, y, z, charge) ? 1 : 0;
    }

    UT_API int ut_wfn_get_ncen(void* h) {
        return static_cast<WFN*>(h)->get_ncen();
    }

    UT_API double ut_wfn_get_atom_x(void* h, int nr) {
        return static_cast<WFN*>(h)->get_atom_coordinate(nr, 0);
    }

    UT_API double ut_wfn_get_atom_y(void* h, int nr) {
        return static_cast<WFN*>(h)->get_atom_coordinate(nr, 1);
    }

    UT_API double ut_wfn_get_atom_z(void* h, int nr) {
        return static_cast<WFN*>(h)->get_atom_coordinate(nr, 2);
    }

    UT_API int ut_wfn_get_atom_charge(void* h, int nr) {
        return static_cast<WFN*>(h)->get_atom_charge(nr);
    }

    UT_API int ut_wfn_get_atom_label(void* h, int nr, char* out, int bufsize) {
        std::string s = static_cast<WFN*>(h)->get_atom_label(nr);
        if (static_cast<int>(s.size()) >= bufsize) return -1;
        std::memcpy(out, s.c_str(), s.size() + 1);
        return static_cast<int>(s.size());
    }

    UT_API double ut_wfn_atom_distance(void* h, int a, int b) {
        auto* w = static_cast<WFN*>(h);
        atom at_a = w->get_atom(a);
        atom at_b = w->get_atom(b);
        return at_a.distance_to(at_b);
    }

    UT_API int ut_wfn_add_exp(void* h, int cent, int type, double e) {
        return static_cast<WFN*>(h)->add_exp(cent, type, e) ? 1 : 0;
    }

    UT_API int ut_wfn_get_nex(void* h) {
        return static_cast<WFN*>(h)->get_nex();
    }

    UT_API int ut_wfn_push_mo(void* h, int nr, double occ, double energy) {
        return static_cast<WFN*>(h)->push_back_MO(nr, occ, energy) ? 1 : 0;
    }

    UT_API int ut_wfn_get_nmo(void* h) {
        return static_cast<WFN*>(h)->get_nmo();
    }

    UT_API double ut_wfn_get_mo_energy(void* h, int mo) {
        return static_cast<WFN*>(h)->get_MO_energy(mo);
    }

    UT_API double ut_wfn_get_mo_occ(void* h, int mo) {
        return static_cast<WFN*>(h)->get_MO_occ(mo);
    }

    UT_API int ut_wfn_set_mo_coef(void* h, int mo, int prim, double val) {
        return static_cast<WFN*>(h)->set_MO_coef(mo, prim, val) ? 1 : 0;
    }

    UT_API double ut_wfn_get_mo_coef(void* h, int mo, int prim) {
        return static_cast<WFN*>(h)->get_MO_coef(mo, prim);
    }

    UT_API int ut_wfn_add_primitive(void* h, int cent, int type, double e,
                                    const double* coefs, int n_coefs) {
        // add_primitive takes double* and uses nmo entries from it
        auto* w = static_cast<WFN*>(h);
        vec vals(coefs, coefs + n_coefs);
        return w->add_primitive(cent, type, e, vals.data()) ? 1 : 0;
    }

    UT_API int ut_wfn_delete_unoccupied_mos(void* h) {
        auto* w = static_cast<WFN*>(h);
        w->delete_unoccupied_MOs();
        return w->get_nmo();
    }

    UT_API void ut_wfn_assign_charge(void* h, int charge) {
        static_cast<WFN*>(h)->assign_charge(charge);
    }

    UT_API int ut_wfn_get_charge(void* h) {
        return static_cast<WFN*>(h)->get_charge();
    }

    UT_API void ut_wfn_assign_multi(void* h, int multi) {
        static_cast<WFN*>(h)->assign_multi(multi);
    }

    UT_API int ut_wfn_get_multi(void* h) {
        return static_cast<int>(static_cast<WFN*>(h)->get_multi());
    }

    UT_API int ut_wfn_set_method(void* h, const char* method) {
        static_cast<WFN*>(h)->set_method(std::string(method));
        return 0;
    }

    UT_API int ut_wfn_get_method(void* h, char* out, int bufsize) {
        std::string s = static_cast<WFN*>(h)->get_method();
        if (static_cast<int>(s.size()) >= bufsize) return -1;
        std::memcpy(out, s.c_str(), s.size() + 1);
        return static_cast<int>(s.size());
    }

    UT_API int ut_wfn_set_basis_set_name(void* h, const char* name) {
        static_cast<WFN*>(h)->set_basis_set_name(std::string(name));
        return 0;
    }

    UT_API int ut_wfn_get_basis_set_name(void* h, char* out, int bufsize) {
        std::string s = static_cast<WFN*>(h)->get_basis_set_name();
        if (static_cast<int>(s.size()) >= bufsize) return -1;
        std::memcpy(out, s.c_str(), s.size() + 1);
        return static_cast<int>(s.size());
    }

    UT_API int ut_wfn_get_nr_electrons(void* h) {
        return static_cast<int>(static_cast<WFN*>(h)->get_nr_electrons());
    }

    UT_API double ut_wfn_count_nr_electrons(void* h) {
        return static_cast<WFN*>(h)->count_nr_electrons();
    }

    UT_API void ut_wfn_resize_dm(void* h, int size) {
        static_cast<WFN*>(h)->resize_DM(size);
    }

    UT_API int ut_wfn_set_dm(void* h, int nr, double val) {
        return static_cast<WFN*>(h)->set_DM(nr, val) ? 1 : 0;
    }

    UT_API double ut_wfn_get_dm(void* h, int nr) {
        return static_cast<WFN*>(h)->get_DM(nr);
    }

    UT_API int ut_wfn_get_dm_size(void* h) {
        return static_cast<WFN*>(h)->get_DM_size();
    }

    // -----------------------------------------------------------------------
    // tsc_block<int, cdouble> handle-based API
    // -----------------------------------------------------------------------

    using TscBlock = tsc_block<int, cdouble>;

    UT_API void* ut_tsc_create_empty() {
        return new TscBlock();
    }

    UT_API void* ut_tsc_create(int n_scatterers, const char** labels,
                               int n_reflections,
                               const int* h_arr, const int* k_arr, const int* l_arr,
                               const double* sf_real, const double* sf_imag) {
        svec scatterers(labels, labels + n_scatterers);
        cvec2 sf(n_scatterers, cvec(n_reflections));
        for (int s = 0; s < n_scatterers; ++s)
            for (int r = 0; r < n_reflections; ++r)
                sf[s][r] = cdouble(sf_real[s * n_reflections + r],
                                   sf_imag[s * n_reflections + r]);
        std::vector<std::vector<int>> idx(3, std::vector<int>(n_reflections));
        for (int r = 0; r < n_reflections; ++r) {
            idx[0][r] = h_arr[r];
            idx[1][r] = k_arr[r];
            idx[2][r] = l_arr[r];
        }
        return new TscBlock(sf, scatterers, idx);
    }

    UT_API void ut_tsc_destroy(void* h) {
        delete static_cast<TscBlock*>(h);
    }

    UT_API int ut_tsc_is_empty(void* h) {
        return static_cast<TscBlock*>(h)->is_empty() ? 1 : 0;
    }

    UT_API int ut_tsc_scatterer_size(void* h) {
        return static_cast<int>(static_cast<TscBlock*>(h)->scatterer_size());
    }

    UT_API int ut_tsc_reflection_size(void* h) {
        return static_cast<int>(static_cast<TscBlock*>(h)->reflection_size());
    }

    UT_API int ut_tsc_get_scatterer(void* h, int nr, char* out, int bufsize) {
        std::string s = static_cast<TscBlock*>(h)->get_scatterer(nr);
        if (static_cast<int>(s.size()) >= bufsize) return -1;
        std::memcpy(out, s.c_str(), s.size() + 1);
        return static_cast<int>(s.size());
    }

    UT_API int ut_tsc_scatterers_string(void* h, char* out, int bufsize) {
        std::string s = static_cast<TscBlock*>(h)->scatterers_string();
        if (static_cast<int>(s.size()) >= bufsize) return -1;
        std::memcpy(out, s.c_str(), s.size() + 1);
        return static_cast<int>(s.size());
    }

    UT_API double ut_tsc_get_sf_real(void* h, int scatterer, int reflection) {
        cvec v = static_cast<TscBlock*>(h)->get_sf_for_scatterer(scatterer);
        return v[reflection].real();
    }

    UT_API double ut_tsc_get_sf_imag(void* h, int scatterer, int reflection) {
        cvec v = static_cast<TscBlock*>(h)->get_sf_for_scatterer(scatterer);
        return v[reflection].imag();
    }

    UT_API int ut_tsc_get_index(void* h, int dim, int refl) {
        return static_cast<TscBlock*>(h)->get_index(dim, refl);
    }

    UT_API void ut_tsc_set_ad(void* h, int val) {
        static_cast<TscBlock*>(h)->set_AD(val != 0);
    }

    UT_API int ut_tsc_get_ad(void* h) {
        return static_cast<TscBlock*>(h)->get_AD() ? 1 : 0;
    }

    UT_API int ut_tsc_append(void* h_dst, void* h_src) {
        std::ostringstream log;
        try {
            static_cast<TscBlock*>(h_dst)->append(*static_cast<TscBlock*>(h_src), log);
            return 0;
        } catch (...) {
            return -1;
        }
    }

    // -----------------------------------------------------------------------
    // primitive evaluation
    // -----------------------------------------------------------------------

    UT_API double ut_primitive_eval(int type, double alpha, double coef,
                                    double norm_const, double r) {
        primitive p;
        p.set_type(type);
        p.set_exp(alpha);
        p.set_coef(coef);
        p.set_norm_const(norm_const);
        return p.eval_gaussian(r);
    }

    UT_API double ut_primitive_eval_unnorm(int type, double alpha, double coef, double r) {
        primitive p;
        p.set_type(type);
        p.set_exp(alpha);
        p.set_coef(coef);
        return p.eval_gaussian_unnormalized(r);
    }

    // -----------------------------------------------------------------------
    // atom::distance_to (standalone)
    // -----------------------------------------------------------------------

    UT_API double ut_atom_distance(double x1, double y1, double z1,
                                   double x2, double y2, double z2) {
        atom a("A", "0", 1, x1, y1, z1, 0);
        atom b("B", "0", 2, x2, y2, z2, 0);
        return a.distance_to(b);
    }

} // extern "C"
