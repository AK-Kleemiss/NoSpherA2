#pragma once
// Lightweight export header for VS unit tests.
// Uses only primitive types so UnitTests.cpp needs no Src/ include paths.

#ifdef BUILDING_DLL
#define UT_API __declspec(dllexport)
#else
#define UT_API __declspec(dllimport)
#endif

extern "C" {
    // --- Geometry ---
    UT_API double ut_array_length3(double x, double y, double z);
    UT_API double ut_array_distance3(double x1, double y1, double z1,
                                     double x2, double y2, double z2);
    UT_API void   ut_vec_diff (double ax, double ay, double az,
                               double bx, double by, double bz, double out[3]);
    UT_API void   ut_vec_cross(double ax, double ay, double az,
                               double bx, double by, double bz, double out[3]);
    UT_API double ut_vec_dot  (double ax, double ay, double az,
                               double bx, double by, double bz);

    // --- Numeric checks ---
    UT_API int    ut_is_nan        (double x);
    UT_API int    ut_is_similar_rel(double a, double b, double tol);
    UT_API int    ut_is_similar_abs(double a, double b, double tol);
    UT_API double ut_fast_exp_neg  (double x);

    // --- FCHK line parsing ---
    // Input: a full FCHK header line; value begins at character 49.
    UT_API int    ut_read_fchk_int(const char* line);
    UT_API double ut_read_fchk_dbl(const char* line);

    // --- String utilities ---
    UT_API int    ut_ends_with    (const char* str, const char* suffix);
    // Strips digits and spaces from input; writes result to output (bufsize incl. NUL).
    // Returns result length, or -1 if buffer too small.
    UT_API int    ut_shrink_string(const char* input, char* output, int bufsize);

    // --- SHA-256 ---
    // Writes the 64-char lowercase hex digest + NUL into output (need bufsize >= 65).
    // Returns 64 on success, -1 if buffer too small.
    UT_API int    ut_sha256(const char* input, char* output, int bufsize);

    // --- Vector aggregates ---
    // Pass elements as an array + count. bools encoded as 0/1 ints.
    UT_API int    ut_vec_sum_bool  (const int* arr, int n);
    UT_API int    ut_vec_sum_int   (const int* arr, int n);
    UT_API double ut_vec_sum_double(const double* arr, int n);
    UT_API double ut_vec_length    (const double* arr, int n);

    // --- Additional string utilities ---
    UT_API int    ut_trim                   (const char* input, char* output, int bufsize);
    UT_API char   ut_asciitolower           (char c);
    UT_API double ut_double_from_esd        (const char* s);   // strips "(esd)"
    UT_API double ut_decimal_precision_cif  (const char* s);   // uncertainty from "1.234(5)"

    // --- Basis-type helpers ---
    UT_API int          ut_sht2nbas        (int type);     // shell type → Cartesian/spherical count
    UT_API unsigned int ut_doublefactorial (int n);        // n!! (double factorial)

    // --- Unit conversions (constants namespace) ---
    UT_API double ut_bohr2ang       (double inp);
    UT_API double ut_ang2bohr       (double inp);
    UT_API double ut_cubic_bohr2ang (double inp);
    UT_API double ut_cubic_ang2bohr (double inp);

    // --- Constexpr math (constants namespace) ---
    UT_API int    ut_const_abs     (int x);
    UT_API double ut_constexpr_pow (double base, int exponent);
    UT_API double ut_constants_sqrt(double x);
    UT_API double ut_exp_approx    (double x, int n);
    UT_API double ut_log_approx    (double x, int n);
    UT_API long long ut_factorial  (int n);              // n! from lookup table

    // --- Orbital index mappings (constants namespace) ---
    // orca_2_pySCF: returns index, or -1 when l is out of range.
    UT_API int          ut_orca_2_pyscf (unsigned int l, int m_idx);
    UT_API unsigned int ut_type_2_nbo   (int type);

    // --- Spherical Bessel functions ---
    UT_API double ut_bessel_j(int l, double x);

    // --- is_similar (pow-10 tolerance variant) ---
    // Returns 1 if |a-b| <= 10^tolerance, else 0.
    UT_API int    ut_is_similar_pow10(double a, double b, double tolerance);

    // --- shell2function: WFN type+prim → column index ---
    UT_API int    ut_shell2function(int type, int prim);

    // --- CountWords ---
    UT_API int    ut_count_words(const char* str);

    // --- shrink_string_to_atom: strips digits and spaces, then truncates to
    //     two characters when atom_number > 9 (matches element-symbol logic) ---
    // Returns result length, -1 if buffer too small.
    UT_API int    ut_shrink_string_to_atom(const char* input, int atom_number,
                                           char* output, int bufsize);

    // --- split_string: splits by delimiter, writes at most max_out tokens.
    //     Returns the number of tokens found (may exceed max_out). ---
    UT_API int    ut_split_string(const char* input, const char* delim,
                                  char** out_tokens, int max_out, int tok_bufsize);

    // --- Timing (monotonic; values are relative — only test sign/ordering) ---
    // Returns elapsed microseconds between two calls to ut_get_time_token().
    // Usage: int t0 = ut_get_time_token(); <work>; long long us = ut_time_diff_us(t0);
    // Since we can't pass opaque time_point objects over a C API we expose a
    // simple "sleep N ms then measure" helper that returns elapsed µs >= N*1000.
    UT_API long long ut_sleep_and_measure_us(int sleep_ms);

    // --- hypergeometric 2F1(a,b;c;x) ---
    UT_API double ut_hypergeometric(double a, double b, double c, double x);

    // --- Atomic number ↔ element symbol (constants::atnr2letter) ---
    // Returns element symbol string for atomic number 1–103; "DM" for 0; "PROBLEM" if out of range.
    // Writes into output (bufsize incl. NUL). Returns string length, -1 if buffer too small.
    UT_API int ut_atnr2letter(int nr, char* output, int bufsize);

    // --- Cartesian exponent vector from basis type (constants::type2vector) ---
    // Converts basis-function type index (1–56) to [nx, ny, nz]. Sets all to -1 if out of range.
    UT_API void ut_type2vector(int index, int nx_out[1], int ny_out[1], int nz_out[1]);

    // --- Gaussian normalization constant (constants::normgauss) ---
    UT_API double ut_normgauss(int type, double exponent);

    // --- Associated Legendre polynomial P_l^m(x) (constants::associated_legendre_polynomial) ---
    UT_API double ut_assoc_legendre(int l, int m, double x);

    // --- Cartesian → spherical coordinates (constants::cartesian_to_spherical) ---
    // Returns [r, theta, phi] written into out3[3].
    UT_API void ut_cartesian_to_spherical(double x, double y, double z, double out3[3]);

    // --- Median eigenvalue of a 3x3 symmetric matrix (get_lambda_1) ---
    // Matrix stored row-major: a[0..8]; a[1]=a[3], a[2]=a[6], a[5]=a[7].
    UT_API double ut_get_lambda_1(const double a[9]);

    // -----------------------------------------------------------------------
    // WFN handle-based API  (opaque void* owns a heap-allocated WFN object)
    // -----------------------------------------------------------------------
    UT_API void*  ut_wfn_create();
    UT_API void   ut_wfn_destroy(void* h);

    // Atoms
    UT_API int    ut_wfn_push_atom(void* h, const char* label,
                                   double x, double y, double z, int charge);
    UT_API int    ut_wfn_get_ncen(void* h);
    UT_API double ut_wfn_get_atom_x(void* h, int nr);
    UT_API double ut_wfn_get_atom_y(void* h, int nr);
    UT_API double ut_wfn_get_atom_z(void* h, int nr);
    UT_API int    ut_wfn_get_atom_charge(void* h, int nr);
    UT_API int    ut_wfn_get_atom_label(void* h, int nr, char* out, int bufsize);
    // Distance between atoms a and b (Euclidean in whatever units atoms are stored)
    UT_API double ut_wfn_atom_distance(void* h, int a, int b);

    // Primitives
    UT_API int    ut_wfn_add_exp(void* h, int cent, int type, double exp);
    UT_API int    ut_wfn_get_nex(void* h);

    // Molecular orbitals
    UT_API int    ut_wfn_push_mo(void* h, int nr, double occ, double energy);
    UT_API int    ut_wfn_get_nmo(void* h);
    UT_API double ut_wfn_get_mo_energy(void* h, int mo);
    UT_API double ut_wfn_get_mo_occ(void* h, int mo);
    UT_API int    ut_wfn_set_mo_coef(void* h, int mo, int prim, double val);
    UT_API double ut_wfn_get_mo_coef(void* h, int mo, int prim);
    // Add a full primitive (center, type, exponent) with one coefficient per MO.
    // coefs must have at least nmo elements. Returns 1 on success.
    UT_API int    ut_wfn_add_primitive(void* h, int cent, int type, double exp,
                                       const double* coefs, int n_coefs);
    // Remove virtual MOs (occ==0); returns nmo after deletion.
    UT_API int    ut_wfn_delete_unoccupied_mos(void* h);

    // Metadata
    UT_API void   ut_wfn_assign_charge(void* h, int charge);
    UT_API int    ut_wfn_get_charge(void* h);
    UT_API void   ut_wfn_assign_multi(void* h, int multi);
    UT_API int    ut_wfn_get_multi(void* h);
    UT_API int    ut_wfn_set_method(void* h, const char* method);
    UT_API int    ut_wfn_get_method(void* h, char* out, int bufsize);
    UT_API int    ut_wfn_set_basis_set_name(void* h, const char* name);
    UT_API int    ut_wfn_get_basis_set_name(void* h, char* out, int bufsize);

    // Electron counts (require atoms pushed with real Z values)
    UT_API int    ut_wfn_get_nr_electrons(void* h);
    UT_API double ut_wfn_count_nr_electrons(void* h);

    // Density matrix (upper-triangular flat storage)
    UT_API void   ut_wfn_resize_dm(void* h, int size);
    UT_API int    ut_wfn_set_dm(void* h, int nr, double val);
    UT_API double ut_wfn_get_dm(void* h, int nr);
    UT_API int    ut_wfn_get_dm_size(void* h);

    // -----------------------------------------------------------------------
    // tsc_block<int,cdouble> handle-based API
    // -----------------------------------------------------------------------
    // Create an empty block or a populated block.
    // sf_real/sf_imag are flat [n_scatterers * n_reflections] row-major arrays.
    // h_arr/k_arr/l_arr are integer arrays of length n_reflections.
    UT_API void*  ut_tsc_create_empty();
    UT_API void*  ut_tsc_create(int n_scatterers, const char** labels,
                                int n_reflections,
                                const int* h_arr, const int* k_arr, const int* l_arr,
                                const double* sf_real, const double* sf_imag);
    UT_API void   ut_tsc_destroy(void* h);

    UT_API int    ut_tsc_is_empty(void* h);
    UT_API int    ut_tsc_scatterer_size(void* h);
    UT_API int    ut_tsc_reflection_size(void* h);
    UT_API int    ut_tsc_get_scatterer(void* h, int nr, char* out, int bufsize);
    UT_API int    ut_tsc_scatterers_string(void* h, char* out, int bufsize);
    // Real and imaginary part of form factor [scatterer][reflection]
    UT_API double ut_tsc_get_sf_real(void* h, int scatterer, int reflection);
    UT_API double ut_tsc_get_sf_imag(void* h, int scatterer, int reflection);
    // Miller index: dim=0→h, 1→k, 2→l
    UT_API int    ut_tsc_get_index(void* h, int dim, int refl);
    UT_API void   ut_tsc_set_ad(void* h, int val);
    UT_API int    ut_tsc_get_ad(void* h);
    // Merge src into dst; returns 0 on success, -1 on error.
    UT_API int    ut_tsc_append(void* h_dst, void* h_src);

    // -----------------------------------------------------------------------
    // primitive evaluation (standalone — no handle needed)
    // -----------------------------------------------------------------------
    // eval_gaussian(type, exp, coef, norm_const, r) = r^type * exp(-exp*r²) * (norm*coef)
    UT_API double ut_primitive_eval(int type, double alpha, double coef,
                                    double norm_const, double r);
    // Unnormalized variant (no norm_const applied)
    UT_API double ut_primitive_eval_unnorm(int type, double alpha, double coef, double r);

    // -----------------------------------------------------------------------
    // atom::distance_to — standalone, no WFN handle needed
    // -----------------------------------------------------------------------
    UT_API double ut_atom_distance(double x1, double y1, double z1,
                                   double x2, double y2, double z2);

    // -----------------------------------------------------------------------
    // cell handle-based API
    // -----------------------------------------------------------------------
    UT_API void*  ut_cell_create(double a, double b, double c,
                                 double alpha, double beta, double gamma);
    UT_API void   ut_cell_destroy(void* h);
    UT_API double ut_cell_get_a(void* h);
    UT_API double ut_cell_get_b(void* h);
    UT_API double ut_cell_get_c(void* h);
    UT_API double ut_cell_get_angle(void* h, int i);   // 0=alpha,1=beta,2=gamma (degrees)
    UT_API double ut_cell_get_volume(void* h);
    UT_API double ut_cell_get_cm(void* h, int i, int j);
    UT_API void   ut_cell_frac_to_cart(void* h, double fx, double fy, double fz,
                                       double out3[3], int in_bohr);
    UT_API double ut_cell_d_spacing_hkl(void* h, int hh, int kk, int ll);
    UT_API double ut_cell_stl_hkl(void* h, int hh, int kk, int ll);
    UT_API int    ut_cell_get_crystal_system(void* h, char* out, int bufsize);

    // -----------------------------------------------------------------------
    // MO handle-based API
    // -----------------------------------------------------------------------
    UT_API void*  ut_mo_create(int nr, double occ, double energy);
    UT_API void   ut_mo_destroy(void* h);
    UT_API void   ut_mo_push_coef(void* h, double val);
    UT_API double ut_mo_get_coef(void* h, int nr);
    UT_API int    ut_mo_get_primitive_count(void* h);
    UT_API int    ut_mo_hdr(void* h, char* out, int bufsize);
    UT_API double ut_mo_get_occ(void* h);
    UT_API double ut_mo_get_energy(void* h);

    // -----------------------------------------------------------------------
    // lebedev_sphere quadrature
    // Returns actual number of points; writes at most max_pts into x/y/z/w.
    // -----------------------------------------------------------------------
    UT_API int    ut_lebedev_order(int order,
                                   double* x_out, double* y_out,
                                   double* z_out, double* w_out,
                                   int max_pts);

    // -----------------------------------------------------------------------
    // AtomGrid handle-based API
    // -----------------------------------------------------------------------
    UT_API void*  ut_atomgrid_create(double radial_precision,
                                     int min_angular, int max_angular,
                                     int proton_charge,
                                     double alpha_max, int max_l,
                                     const double* alpha_min, int n_alpha_min);
    UT_API void   ut_atomgrid_destroy(void* h);
    UT_API int    ut_atomgrid_get_num_points(void* h);
    UT_API int    ut_atomgrid_get_num_radial_points(void* h);
    // Radial grid (r and w arrays must each have at least get_num_radial_points() elements)
    UT_API void   ut_atomgrid_get_radial_grid(void* h, double* r_out, double* w_out);

    // -----------------------------------------------------------------------
    // WFN density evaluation (uses existing wfn handle from ut_wfn_create)
    // -----------------------------------------------------------------------
    UT_API double ut_wfn_compute_dens(void* h, double x, double y, double z);
    UT_API double ut_wfn_compute_mo(void* h, double x, double y, double z, int mo_idx);
}
