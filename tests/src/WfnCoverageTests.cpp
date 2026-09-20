#include "pch.h"
#include "core/convenience.h"
#include "core/constants.h"
#include "core/atoms.h"
#include "core/wfn_class.h"
#include "core/mo_class.h"
#include "core/cube.h"
#include <unordered_set>
#include <cstdint>

//Gaps left after WfnOpsTests / WfnReadTests / AtomSymmetryTests: the atomID hex and 64-bit views, the MO
//copy and edge headers, the hand-built basis-set shell accessors and get_norm_const, build_DM, the xdgraph
//cube writer, the unimplemented spherical MO, a synthetic open-shell pTB file and the gbw ECP branch.
namespace
{
    //scratch directory named after the test, removed by the caller on success
    static std::filesystem::path scratch_dir(const std::string& name)
    {
        const std::filesystem::path p = std::filesystem::temp_directory_path() / ("nosphera2_wfncoverage_" + name);
        std::filesystem::remove_all(p);
        std::filesystem::create_directories(p);
        return p;
    }

    //restores the working directory even when an assertion returns early
    struct cwd_guard
    {
        std::filesystem::path old = std::filesystem::current_path();
        ~cwd_guard() { std::filesystem::current_path(old); }
    };

    //one Fortran sequential record: int32 length, payload, int32 length
    static void record(std::ofstream& f, const void* data, const int bytes)
    {
        f.write(reinterpret_cast<const char*>(&bytes), sizeof(int));
        f.write(static_cast<const char*>(data), bytes);
        f.write(reinterpret_cast<const char*>(&bytes), sizeof(int));
    }
    static void record(std::ofstream& f, const int v) { record(f, &v, sizeof(int)); }
    static void record(std::ofstream& f, const double v) { record(f, &v, sizeof(double)); }
    static void record(std::ofstream& f, const vec& v) { record(f, v.data(), static_cast<int>(v.size() * sizeof(double))); }

    //the standard cartesian primitive normalisation constants for x^l with exponent a
    static double norm_s(const double a) { return std::pow(2 * a / constants::PI, 0.75); }
    static double norm_p(const double a) { return std::pow(128 * std::pow(a, 5) / constants::PI3, 0.25); }
    static double norm_d(const double a) { return std::pow(2048 * std::pow(a, 7) / (9 * constants::PI3), 0.25); }
    static double norm_f(const double a) { return std::pow(32768 * std::pow(a, 9) / (225 * constants::PI3), 0.25); }

    //He with a single s (a = 1) and a single p (b = 0.5) shell in the basis set and matching primitives in
    //tonto order; MO 0 (occ 2) is 0.6 s + 0.8 px in normalised AOs, MO 1 (occ 0) is pure py
    static WFN make_he_sp()
    {
        WFN w(e_origin::wfn);
        w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
        w.push_back_atom_basis_set(0, 1.0, 1.0, 1, 0);
        w.push_back_atom_basis_set(0, 0.5, 1.0, 2, 1);
        w.push_back_MO(1, 2.0, -0.9);
        w.push_back_MO(2, 0.0, 0.4);
        double s[2] = { 0.6 * norm_s(1.0), 0.0 };
        double px[2] = { 0.8 * norm_p(0.5), 0.0 };
        double py[2] = { 0.0, norm_p(0.5) };
        double pz[2] = { 0.0, 0.0 };
        w.add_primitive(1, 1, 1.0, s);
        w.add_primitive(1, 2, 0.5, px);
        w.add_primitive(1, 3, 0.5, py);
        w.add_primitive(1, 4, 0.5, pz);
        w.set_exp_cutoff();
        return w;
    }

    //a one-hydrogen, one-basis-function pTB file with nmomax MOs; H has no pTB core, so one alpha electron
    //and every MO after the first is beta. with_dm appends the packed density matrix record.
    static void write_ptb(const std::filesystem::path& p, const bool with_dm)
    {
        std::ofstream f(p, std::ios::binary);
        record(f, 1);
        const int infos[4] = { 1, 1, 2, 1 };
        record(f, infos, sizeof(infos));
        record(f, "H ", 2);
        record(f, 0.1);
        record(f, -0.2);
        record(f, 0.3);
        record(f, 1);
        record(f, 1);
        record(f, 1);
        record(f, 1);
        record(f, vec{ 1.0 });
        record(f, vec{ 2.0 });
        record(f, vec{ 1.0, 0.0 });
        record(f, vec{ -0.5, 0.25 });
        record(f, vec{ 0.5, 0.25 });
        if (with_dm)
            record(f, vec{ 0.9, 0.1, 0.2 });
    }

    //the 64-bit view packs frac_x | frac_y << 32 and frac_z | data << 32 | Z << 48 | reserved << 56;
    //+-16 encode to +-INT32_MAX, so the hex string is derivable by hand
    TEST(WfnCoverageAtomIdTests, HexStringRoundTripsEncodedBits)
    {
        const atomID edge(16.0, -16.0, 0.0, 0, 1);
        EXPECT_EQ(edge.to_hex_string(), "800000017fffffff0001000000000000");
        const atomID packed(0.0, 0.0, 0.0, -1, 6, 255);
        EXPECT_EQ(packed.to_hex_string(), "0000000000000000ff06ffff00000000");
        EXPECT_EQ(atomID(std::string_view(edge.to_hex_string())), edge);
        EXPECT_EQ(atomID(std::string_view(packed.to_hex_string())), packed);
        const std::array<std::uint64_t, 2> bits = packed.as_uint64();
        EXPECT_EQ(bits[0], 0u);
        EXPECT_EQ(bits[1], 0xff06ffff00000000ull);
        const atomID rebuilt(bits[0], bits[1]);
        EXPECT_EQ(rebuilt, packed);
        EXPECT_EQ(rebuilt.data(), -1);
        EXPECT_EQ(rebuilt.Z(), 6);
        EXPECT_EQ(rebuilt.reserved(), 255);
        std::ostringstream os;
        os << edge;
        EXPECT_EQ(os.str(), "atomID(frac_x: 16, frac_y: -16, frac_z: 0, Z: 1, data: 0, reserved: 0)");
        EXPECT_EQ(std::hash<atomID>{}(edge), AtomIDHash{}(edge));
        std::unordered_set<atomID> ids;
        ids.insert(edge);
        ids.insert(atomID(16.0, -16.0, 0.0, 0, 1));
        ids.insert(packed);
        EXPECT_EQ(ids.size(), 2u);
    }

    //length and character checks of the hex constructor, and the loaded-data validation of the binary one
    TEST(WfnCoverageAtomIdTests, HexAndBinaryConstructorsRejectBadInput)
    {
        EXPECT_THROW(atomID(std::string_view("0000000000000000000100000000000")), std::invalid_argument);
        EXPECT_THROW(atomID(std::string_view("000000000000000g0001000000000000")), std::invalid_argument);
        EXPECT_THROW(atomID(std::uint64_t(0), std::uint64_t(0)), std::runtime_error);
        EXPECT_THROW(atomID(std::uint64_t(0x80000000u), std::uint64_t(1) << 48), std::runtime_error);
        EXPECT_NO_THROW(atomID(std::uint64_t(0), std::uint64_t(1) << 48));
    }

    //copy construction and self-assignment keep the coefficients, insert_into_coefficients appends a range,
    //set_occ(int) stores a double, and hdr appends nothing for occ > 2 or ener == 0
    TEST(WfnCoverageMoTests, CopySelfAssignInsertAndHeaderEdges)
    {
        MO a(2, 1.0, -0.3, 1);
        a.push_back_coef(0.5);
        a.insert_into_coefficients(vec{ -0.25, 0.125 });
        ASSERT_EQ(a.get_primitive_count(), 3);
        EXPECT_DOUBLE_EQ(a.get_coefficient(2), 0.125);
        const MO b(a);
        EXPECT_EQ(b.get_op(), 1);
        EXPECT_DOUBLE_EQ(b.get_energy(), -0.3);
        EXPECT_EQ(b.get_coefficients(), (vec{ 0.5, -0.25, 0.125 }));
        MO& self = a;
        a = self;
        EXPECT_EQ(a.get_coefficients(), (vec{ 0.5, -0.25, 0.125 }));
        a.set_occ(2);
        EXPECT_DOUBLE_EQ(a.get_occ(), 2.0);
        EXPECT_EQ(MO(0, 3.0, 0.0).hdr(), std::string("MO   0") + "                  OCC NO =");
        EXPECT_EQ(MO(7, 2.0, 0.0).hdr(), std::string("MO    7") + "                  OCC NO =    2.00000000 ORB. ENERGY =");
    }

    //set_coefficient refuses a negative or past-the-end index through err_checkf
    TEST(WfnCoverageMoDeathTest, SetCoefficientOutOfRangeExits)
    {
        MO mo(1, 2.0, -1.0);
        mo.push_back_coef(0.5);
        EXPECT_EXIT(mo.set_coefficient(-1, 1.0), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
        EXPECT_EXIT(mo.set_coefficient(1, 1.0), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    }

    //shell bookkeeping on a hand-built basis: a contracted p shell (2 primitives) counts 3 * 2 cartesian
    //primitives before the next atom, the last shell ends at the basis size, out-of-range asks return -1
    TEST(WfnCoverageBasisTests, ShellAccessorsOnHandBuiltBasis)
    {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("C", 0.0, 0.0, 0.0, 6);
        w.push_back_atom("H", 2.0, 0.0, 0.0, 1);
        w.push_back_atom_basis_set(0, 3.0, 1.0, 1, 0);
        w.push_back_atom_basis_set(0, 1.0, 0.7, 2, 1);
        w.push_back_atom_basis_set(0, 0.5, 0.3, 2, 1);
        w.push_back_atom_basis_set(1, 1.2, 1.0, 1, 0);
        EXPECT_EQ(w.get_nr_basis_set_loaded(), 2);
        EXPECT_EQ(w.get_atom_shell_count(0), 2);
        EXPECT_EQ(w.get_atom_shell_count(1), 1);
        EXPECT_EQ(w.get_atom_shell_count(2), -1);
        EXPECT_EQ(w.get_atom_shell_primitives(0, 1), 2);
        EXPECT_EQ(w.get_atom_shell_primitives(0, 2), -1);
        EXPECT_EQ(w.get_shell_type(0, 0), 1);
        EXPECT_EQ(w.get_shell_type(0, 1), 2);
        EXPECT_EQ(w.get_shell_type(5, 0), -1);
        EXPECT_EQ(w.get_shell_start(0, 1), 1);
        EXPECT_EQ(w.get_shell_start(0, 9), -1);
        EXPECT_EQ(w.get_shell_end(0, 0), 0);
        EXPECT_EQ(w.get_shell_end(0, 1), 2);
        EXPECT_EQ(w.get_shell_end(1, 0), 0);
        EXPECT_EQ(w.get_shell_end(1, 1), -1);
        EXPECT_EQ(w.get_shell_start_in_primitives(0, 0), 0);
        EXPECT_EQ(w.get_shell_start_in_primitives(0, 1), 1);
        EXPECT_EQ(w.get_shell_start_in_primitives(1, 0), 7);
        EXPECT_EQ(w.get_shell_start_in_primitives(0, 9), -1);
        EXPECT_EQ(w.get_atom_primitive_count(0), 3);
        EXPECT_EQ(w.get_atom_primitive_type(0, 2), 2);
        EXPECT_EQ(w.get_atom_basis_set_exponent(0, 2), 0.5);
        EXPECT_EQ(w.get_atom_basis_set_coefficient(0, 1), 0.7);
    }

    //single-primitive s, p, d and f shells normalise to the textbook cartesian constants, with sqrt(3),
    //sqrt(5) and sqrt(15) on the mixed components, and a contracted s shell to unit self-overlap
    TEST(WfnCoverageBasisTests, NormConstMatchesCartesianConstants)
    {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
        w.push_back_atom_basis_set(0, 1.0, 1.0, 1, 0);
        w.push_back_atom_basis_set(0, 0.5, 1.0, 2, 1);
        w.push_back_atom_basis_set(0, 0.4, 1.0, 3, 2);
        w.push_back_atom_basis_set(0, 0.3, 1.0, 4, 3);
        std::ostringstream log;
        const vec n = w.get_norm_const(log, false);
        ASSERT_EQ(n.size(), 20u);
        vec expected;
        expected.push_back(norm_s(1.0));
        for (int k = 0; k < 3; k++)
            expected.push_back(norm_p(0.5));
        for (int k = 0; k < 3; k++)
            expected.push_back(norm_d(0.4));
        for (int k = 0; k < 3; k++)
            expected.push_back(std::sqrt(3.0) * norm_d(0.4));
        for (int k = 0; k < 3; k++)
            expected.push_back(norm_f(0.3));
        for (int k = 0; k < 6; k++)
            expected.push_back(std::sqrt(5.0) * norm_f(0.3));
        expected.push_back(std::sqrt(15.0) * norm_f(0.3));
        for (size_t i = 0; i < expected.size(); i++)
            EXPECT_NEAR(n[i], expected[i], 1e-12) << "entry " << i;

        WFN c(e_origin::NOT_YET_DEFINED);
        c.push_back_atom("H", 0.0, 0.0, 0.0, 1);
        c.push_back_atom_basis_set(0, 1.0, 1.0, 1, 0);
        c.push_back_atom_basis_set(0, 2.0, 0.5, 1, 0);
        const vec m = c.get_norm_const(log, true);
        ASSERT_EQ(m.size(), 2u);
        const double a[2] = { 1.0, 2.0 };
        double overlap = 0.0;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 2; j++)
                overlap += m[i] * m[j] * std::pow(constants::PI / (a[i] + a[j]), 1.5);
        EXPECT_NEAR(overlap, 1.0, 1e-12);
        EXPECT_NEAR(m[0] / m[1], norm_s(1.0) / (0.5 * norm_s(2.0)), 1e-12);
    }

    //a second atom without a basis set trips the second err_checkf of get_norm_const
    TEST(WfnCoverageBasisDeathTest, NormConstWithAnAtomLackingBasisExits)
    {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("H", 0.0, 0.0, 0.0, 1);
        w.push_back_atom("H", 1.4, 0.0, 0.0, 1);
        w.push_back_atom_basis_set(0, 1.0, 1.0, 1, 0);
        std::ostringstream log;
        EXPECT_EXIT(w.get_norm_const(log, false), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    }

    //build_DM on a wfn-origin He: check_order finds tonto p order, the basis normalises to the textbook
    //constants so the AO coefficients are (0.6, 0.8, 0, 0), and the packed DM is 2 c_u c_v over the
    //occupied MO only: 0.72, 0.96, 1.28 in the first three slots of the 4 * 5 / 2 triangle
    TEST(WfnCoverageDmTests, BuildDmFromNormalisedHeOrbitals)
    {
        WFN w = make_he_sp();
        EXPECT_EQ(w.check_order(false), 2);
        ASSERT_TRUE(w.build_DM("unused-because-basis-is-loaded", false));
        ASSERT_EQ(w.get_DM_size(), 10);
        EXPECT_NEAR(w.get_DM(0), 0.72, 1e-12);
        EXPECT_NEAR(w.get_DM(1), 0.96, 1e-12);
        EXPECT_NEAR(w.get_DM(2), 1.28, 1e-12);
        for (int i = 3; i < 10; i++)
            EXPECT_NEAR(w.get_DM(i), 0.0, 1e-12) << "slot " << i;
        EXPECT_EQ(w.get_SDM_size(), 0);
        EXPECT_EQ(w.get_DM(10), -1.0);
    }

    //the same wavefunction under an origin build_DM does not handle is refused without touching the DM
    TEST(WfnCoverageDmTests, BuildDmRefusesUnsupportedOrigin)
    {
        WFN w = make_he_sp();
        w.set_origin(e_origin::wfx);
        EXPECT_FALSE(w.build_DM("unused-because-basis-is-loaded", false));
        EXPECT_EQ(w.get_DM_size(), 0);
    }

    //without a loaded basis build_DM asks the library for one under the given path and dies when it is missing
    TEST(WfnCoverageDmDeathTest, BuildDmWithoutBasisAndMissingLibraryExits)
    {
        WFN w(e_origin::wfn);
        w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
        w.push_back_MO(1, 2.0, -0.9);
        double c = norm_s(1.0);
        w.add_primitive(1, 1, 1.0, &c);
        w.set_basis_set_name("no-such-basis");
        const std::string missing = (std::filesystem::temp_directory_path() / "nosphera2_wfncoverage_no_library").string();
        EXPECT_EXIT(w.build_DM(missing, false), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    }

    //write_cube_xdgraph: fixed header, one line per (x, z) holding the y values, no atoms when na is 0
    TEST(WfnCoverageIoTests, CubeXdgraphWritesHeaderAndYLines)
    {
        const std::filesystem::path dir = scratch_dir("xdgraph");
        const std::filesystem::path out = dir / "grid.xdgrid";
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
        cube c({ 2, 3, 4 }, 0, true);
        for (int i = 0; i < 3; i++)
        {
            c.set_origin(i, -1.0 + 0.5 * i);
            c.set_vector(i, i, 0.5);
        }
        for (int x = 0; x < 2; x++)
            for (int y = 0; y < 3; y++)
                for (int z = 0; z < 4; z++)
                    c.set_value(x, y, z, 1.0 + x + 10 * y + 100 * z);
        w.push_back_cube(c);
        w.write_cube_xdgraph(0, out, true);
        EXPECT_EQ(w.get_cube_path(0), out);
        std::ifstream in(out);
        std::vector<std::string> lines;
        std::string line;
        while (std::getline(in, line))
            lines.push_back(line);
        ASSERT_EQ(lines.size(), 20u);
        EXPECT_EQ(lines[0], "2DGRDFIL  0");
        EXPECT_EQ(lines[1], "cuQCT    FOU");
        EXPECT_EQ(lines[2], "");
        EXPECT_EQ(lines[3], "! Gridpoints, Origin, Physical Dimensions");
        EXPECT_EQ(lines[4], "             2             3             4");
        EXPECT_EQ(lines[5], "    -1.0000 -0.5000 0.0000");
        EXPECT_EQ(lines[6], "    0.500000    0.500000    0.500000");
        EXPECT_EQ(lines[7], "! Objects");
        EXPECT_EQ(lines[8], "         0");
        EXPECT_EQ(lines[9], "! Connections");
        EXPECT_EQ(lines[10], "         0");
        EXPECT_EQ(lines[11], "! Values");
        EXPECT_EQ(lines[12], "  1.0000000E+00  1.1000000E+01  2.1000000E+01");
        EXPECT_EQ(lines[13], "  1.0100000E+02  1.1100000E+02  1.2100000E+02");
        EXPECT_EQ(lines[16], "  2.0000000E+00  1.2000000E+01  2.2000000E+01");
        EXPECT_EQ(lines[19], "  3.0200000E+02  3.1200000E+02  3.2200000E+02");
        in.close();
        std::filesystem::remove_all(dir);
    }

    //the spherical MO evaluator is a stub that aborts through not_implemented
    TEST(WfnCoverageIoDeathTest, ComputeMoSphericalIsNotImplemented)
    {
        WFN w(e_origin::NOT_YET_DEFINED);
        w.push_back_atom("He", 0.0, 0.0, 0.0, 2);
        w.push_back_MO(1, 2.0, -0.9);
        double c = norm_s(1.0);
        w.add_primitive(1, 1, 1.0, &c);
        EXPECT_EXIT(w.compute_MO_spherical(d3{ 0.0, 0.0, 0.0 }, 0), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
        EXPECT_EXIT(w.write_cube_xdgraph(0, "unused.xdgrid", false), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    }

    //a synthetic pTB file: one H, one basis function, two MOs. H has one electron and no pTB core, so the
    //reader makes it a doublet with the second MO beta, primitive coefficients are momat * contr, the
    //packed density record unpacks symmetric, and debug writes the log and a wfn copy into the cwd
    TEST(WfnCoverageIoTests, PtbSyntheticOpenShellWithDensityMatrix)
    {
        const std::filesystem::path dir = scratch_dir("ptb");
        const std::filesystem::path p = dir / "wfn.xtb";
        write_ptb(p, true);
        cwd_guard guard;
        std::filesystem::current_path(dir);
        WFN w(e_origin::NOT_YET_DEFINED);
        std::ostringstream log;
        ASSERT_TRUE(w.read_ptb(p, log, true));
        EXPECT_EQ(w.get_origin(), e_origin::ptb);
        EXPECT_TRUE(w.get_isBohr());
        EXPECT_EQ(w.get_path(), p);
        ASSERT_EQ(w.get_ncen(), 1);
        EXPECT_EQ(w.get_atom_label(0), "H");
        EXPECT_EQ(w.get_atom_charge(0), 1);
        EXPECT_DOUBLE_EQ(w.get_atom_coordinate(0, 0), 0.1);
        EXPECT_DOUBLE_EQ(w.get_atom_coordinate(0, 1), -0.2);
        EXPECT_DOUBLE_EQ(w.get_atom_coordinate(0, 2), 0.3);
        EXPECT_EQ(w.get_atom_ECP_electrons(0), 0);
        EXPECT_EQ(w.get_nr_basis_set_loaded(), 1);
        EXPECT_EQ(w.get_multi(), 2u);
        ASSERT_EQ(w.get_nmo(), 2);
        ASSERT_EQ(w.get_nex(), 1);
        EXPECT_EQ(w.get_center(0), 1);
        EXPECT_EQ(w.get_type(0), 1);
        EXPECT_DOUBLE_EQ(w.get_exponent(0), 1.0);
        EXPECT_DOUBLE_EQ(w.get_MO_coef(0, 0), 1.0);
        EXPECT_DOUBLE_EQ(w.get_MO_coef(1, 0), 0.5);
        EXPECT_DOUBLE_EQ(w.get_MO_occ(0), 1.0);
        EXPECT_DOUBLE_EQ(w.get_MO_occ(1), 0.0);
        EXPECT_DOUBLE_EQ(w.get_MO_energy(0), -0.5);
        EXPECT_DOUBLE_EQ(w.get_MO_energy(1), 0.25);
        EXPECT_EQ(w.get_MO(0).get_op(), 0);
        EXPECT_EQ(w.get_MO(1).get_op(), 1);
        EXPECT_TRUE(w.get_is_unrestricted());
        EXPECT_EQ(w.get_MO_op_count(1), 1);
        EXPECT_NEAR(w.compute_dens(d3{ 0.1, -0.2, 0.3 }), 1.0, 1e-12);
        EXPECT_NEAR(w.compute_dens(d3{ 1.1, -0.2, 0.3 }), std::exp(-2.0), 1e-12);
        const dMatrix2 dm = w.get_dm();
        ASSERT_EQ(dm.extent(0), 2u);
        ASSERT_EQ(dm.extent(1), 2u);
        EXPECT_DOUBLE_EQ(dm(0, 0), 0.9);
        EXPECT_DOUBLE_EQ(dm(0, 1), 0.1);
        EXPECT_DOUBLE_EQ(dm(1, 0), 0.1);
        EXPECT_DOUBLE_EQ(dm(1, 1), 0.2);
        EXPECT_NE(log.str().find("elcount after: 1"), std::string::npos);
        EXPECT_NE(log.str().find("al/be els after:1 0"), std::string::npos);
        EXPECT_NE(log.str().find("occs: 1 0 "), std::string::npos);
        EXPECT_TRUE(std::filesystem::exists(dir / "test_convert_from_xtb.wfn"));

        //without the density record the DM stays empty and the density still comes from the MOs
        const std::filesystem::path q = dir / "nodm.xtb";
        write_ptb(q, false);
        WFN v(e_origin::NOT_YET_DEFINED);
        ASSERT_TRUE(v.read_ptb(q, log, false));
        EXPECT_EQ(v.get_dm().extent(0), 0u);
        EXPECT_NEAR(v.compute_dens(d3{ 0.1, -0.2, 0.3 }), 1.0, 1e-12);

        //a missing file is reported, not fatal
        WFN u(e_origin::NOT_YET_DEFINED);
        EXPECT_FALSE(u.read_ptb(dir / "absent.xtb", log, false));
        EXPECT_EQ(u.get_origin(), e_origin::ptb);
        std::filesystem::current_path(guard.old);
        std::filesystem::remove_all(dir);
    }

    //a pTB header announcing zero atoms, or a leading record equal to 2, is fatal
    TEST(WfnCoverageIoDeathTest, PtbBadHeaderExits)
    {
        const std::filesystem::path dir = scratch_dir("ptb_death");
        const std::filesystem::path zero = dir / "zero.xtb";
        {
            std::ofstream f(zero, std::ios::binary);
            record(f, 1);
            const int infos[4] = { 0, 1, 1, 1 };
            record(f, infos, sizeof(infos));
        }
        const std::filesystem::path two = dir / "two.xtb";
        {
            std::ofstream f(two, std::ios::binary);
            record(f, 2);
        }
        std::ostringstream log;
        WFN w(e_origin::NOT_YET_DEFINED);
        EXPECT_EXIT(w.read_ptb(zero, log, false), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
        WFN v(e_origin::NOT_YET_DEFINED);
        EXPECT_EXIT(v.read_ptb(two, log, false), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
        std::filesystem::remove_all(dir);
    }

    //the gbw ECP branch reads its section pointer from byte 32; an all-electron file of the classic layout
    //stores 0 there, which is fatal
    TEST(WfnCoverageIoDeathTest, GbwEcpPointerZeroExits)
    {
        const std::filesystem::path p = nos_test_repo_root() / "tests" / "epoxide_gbw" / "epoxide.gbw";
        ASSERT_TRUE(std::filesystem::exists(p));
        std::ostringstream log;
        WFN w(e_origin::NOT_YET_DEFINED);
        EXPECT_EXIT(w.read_gbw(p, log, false, true), ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    }

    //HgH2 with the def2 ECP on Hg: 60 core electrons should arrive on the mercury atom
    //suspected defect: Src/core/wfn_io.cpp:1442 the ECP pointer is read from byte 32 whatever the magic number says; a magic -1 file keeps its geometry pointer there, so the atom block is parsed as ECP entries and err_checkf(Z > 0) exits, while every magic 40 fixture (Rb.gbw with an ECP included) holds 0 there and exits on the pointer check
    TEST(WfnCoverageIoTests, DISABLED_GbwEcpBlockOnCurrentFormat)
    {
        const std::filesystem::path p = nos_test_repo_root() / "tests" / "ELI_heavy" / "hgh2_ecp.gbw";
        ASSERT_TRUE(std::filesystem::exists(p));
        std::ostringstream log;
        WFN w(e_origin::NOT_YET_DEFINED);
        ASSERT_TRUE(w.read_gbw(p, log, false, true));
        ASSERT_EQ(w.get_ncen(), 3);
        EXPECT_TRUE(w.get_has_ECPs());
        EXPECT_EQ(w.get_atom_ECP_electrons(0), 60);
        EXPECT_EQ(w.get_nr_ECP_electrons(), 60u);
    }
}
