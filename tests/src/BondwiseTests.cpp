#include "pch.h"
#include "core/convenience.h"
#include "core/wfn_class.h"
#include "core/cube.h"
#include "core/atoms.h"
#include "core/constants.h"
#include "core/bondwise_analysis.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
    //captures std::cout for the lifetime of the object; the analysis code reports through cout
    struct CoutCapture
    {
        std::ostringstream buffer;
        std::streambuf* old;
        CoutCapture() : old(std::cout.rdbuf(buffer.rdbuf())) {}
        ~CoutCapture() { std::cout.rdbuf(old); }
        std::string str() const { return buffer.str(); }
    };

    std::filesystem::path bondwise_temp(const std::string& name)
    {
        return std::filesystem::temp_directory_path() / ("bondwise_tests_" + name);
    }

    dMatrix2 diagonal_matrix(const std::vector<double>& diag)
    {
        const int n = static_cast<int>(diag.size());
        dMatrix2 m(n, n);
        std::fill(m.container().begin(), m.container().end(), 0.0);
        for (int i = 0; i < n; i++)
            m(i, i) = diag[i];
        return m;
    }

    double trace(const dMatrix2& m)
    {
        double t = 0.0;
        for (int i = 0; i < static_cast<int>(m.extent(0)); i++)
            t += m(i, i);
        return t;
    }

    int count_occurrences(const std::string& text, const std::string& key)
    {
        int n = 0;
        for (size_t pos = text.find(key); pos != std::string::npos; pos = text.find(key, pos + key.size()))
            n++;
        return n;
    }

    //the number printed right after key, e.g. "Population of atom 0: 9.42047"
    double value_after(const std::string& text, const std::string& key)
    {
        const size_t pos = text.find(key);
        if (pos == std::string::npos)
            return std::numeric_limits<double>::quiet_NaN();
        std::istringstream iss(text.substr(pos + key.size()));
        double v = std::numeric_limits<double>::quiet_NaN();
        iss >> v;
        return v;
    }

    //the numeric columns of the first table row that contains key, e.g. "N - Li   9.420 ..."
    std::vector<double> row_numbers_after(const std::string& text, const std::string& key)
    {
        std::vector<double> out;
        const size_t pos = text.find(key);
        if (pos == std::string::npos)
            return out;
        const size_t end = text.find('\n', pos);
        std::istringstream iss(text.substr(pos + key.size(), end - pos - key.size()));
        double v;
        while (iss >> v)
            out.push_back(v);
        return out;
    }

    //three atoms in a right angle, no orbitals; enough for every validation branch of do_bonds
    WFN three_atom_wfn()
    {
        WFN wavy(e_origin::NOT_YET_DEFINED);
        wavy.push_back_atom("C", 0.0, 0.0, 0.0, 6);
        wavy.push_back_atom("N", 0.0, 0.0, 2.0, 7);
        wavy.push_back_atom("O", 1.0, 0.0, 2.0, 8);
        return wavy;
    }

    void write_text(const std::filesystem::path& p, const std::string& text)
    {
        std::ofstream f(p);
        f << text;
    }

    //runs autobonds on the given input text with a 3-atom wfn and returns {return value, captured stdout}
    std::pair<int, std::string> run_autobonds(const std::string& name, const std::string& text)
    {
        WFN wavy = three_atom_wfn();
        const auto input = bondwise_temp(name + ".inp");
        write_text(input, text);
        int rc;
        std::string out;
        {
            CoutCapture cap;
            rc = autobonds(false, wavy, input, true);
            out = cap.str();
        }
        std::filesystem::remove(input);
        return { rc, out };
    }

    const char* const autobonds_header = " 1    1    1    1\n";

    //two atoms 3 bohr apart on a 15 x 7 x 7 grid of 0.5 bohr: rho is two equal gaussians, eli encodes the
    //voxel index as 100 x + 10 y + z + 1 so the kept voxels can be decoded after the mask
    constexpr int grid_nx = 15, grid_ny = 7, grid_nz = 7;
    constexpr int atom_a_ix = 4, atom_b_ix = 10, atom_iy = 3, atom_iz = 3;
    constexpr int mid_ix = (atom_a_ix + atom_b_ix) / 2; //the x column exactly between the two atoms

    void make_pair_grid(WFN& parent, cube& rho, cube& eli)
    {
        parent.push_back_atom("C", 0.0, 0.0, 0.0, 6);
        parent.push_back_atom("N", 3.0, 0.0, 0.0, 7);
        rho = cube({ grid_nx, grid_ny, grid_nz }, 2, true);
        eli = cube({ grid_nx, grid_ny, grid_nz }, 2, true);
        rho.give_parent_wfn(parent);
        eli.give_parent_wfn(parent);
        const double origin[3] = { -2.0, -1.5, -1.5 };
        for (int i = 0; i < 3; i++) {
            rho.set_origin(i, origin[i]);
            eli.set_origin(i, origin[i]);
            for (int j = 0; j < 3; j++) {
                rho.set_vector(i, j, i == j ? 0.5 : 0.0);
                eli.set_vector(i, j, i == j ? 0.5 : 0.0);
            }
        }
        rho.calc_dv();
        eli.calc_dv();
        for (int x = 0; x < grid_nx; x++)
            for (int y = 0; y < grid_ny; y++)
                for (int z = 0; z < grid_nz; z++) {
                    const auto p = rho.get_pos(x, y, z);
                    const double ra2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
                    const double rb2 = (p[0] - 3.0) * (p[0] - 3.0) + p[1] * p[1] + p[2] * p[2];
                    rho.set_value(x, y, z, std::exp(-ra2 / 4.0) + std::exp(-rb2 / 4.0));
                    eli.set_value(x, y, z, 100.0 * x + 10.0 * y + z + 1.0);
                }
    }

    double encoded_eli(int x, int y, int z)
    {
        return 100.0 * x + 10.0 * y + z + 1.0;
    }

    //reads a masked cube back and lists the decoded (x,y,z) indices of every kept voxel
    struct MaskedCube
    {
        std::array<int, 3> size;
        std::array<double, 3> origin;
        std::vector<std::array<int, 3>> kept;
        int background_voxels = 0;
        std::string comment1, comment2;
    };

    MaskedCube read_masked(const std::filesystem::path& p, double background, std::ostream& log)
    {
        WFN dummy;
        cube c(p, true, dummy, log);
        MaskedCube m;
        m.comment1 = c.get_comment1();
        m.comment2 = c.get_comment2();
        for (int i = 0; i < 3; i++) {
            m.size[i] = c.get_size(i);
            m.origin[i] = c.get_origin(i);
        }
        for (int x = 0; x < m.size[0]; x++)
            for (int y = 0; y < m.size[1]; y++)
                for (int z = 0; z < m.size[2]; z++) {
                    const double v = c.get_value(x, y, z);
                    if (std::abs(v - background) < 1e-9) {
                        m.background_voxels++;
                        continue;
                    }
                    const int code = static_cast<int>(std::lround(v)) - 1;
                    m.kept.push_back({ code / 100, (code / 10) % 10, code % 10 });
                }
        return m;
    }
}

//a cartesian p shell is one irrep of O_h, so any diagonal p block averages to (mean) x identity
TEST(BondwiseSymmetrizeTests, CartesianPShellAveragesToScalar)
{
    dMatrix2 m = diagonal_matrix({ 1.0, 2.0, 3.0 });
    symmetrize_atomic_matrix_oh(m, { 1 }, false);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_NEAR(m(i, j), i == j ? 2.0 : 0.0, 1e-12) << i << " " << j;
}

//cartesian d splits into {xx,yy,zz} and {xy,xz,yz}; each set averages separately and the trace survives
TEST(BondwiseSymmetrizeTests, CartesianDShellSplitsIntoTwoSets)
{
    dMatrix2 m = diagonal_matrix({ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 });
    const double before = trace(m);
    symmetrize_atomic_matrix_oh(m, { 2 }, false);
    for (int i = 0; i < 3; i++)
        EXPECT_NEAR(m(i, i), 2.0, 1e-12) << i;
    for (int i = 3; i < 6; i++)
        EXPECT_NEAR(m(i, i), 5.0, 1e-12) << i;
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 6; j++)
            if (i != j)
                EXPECT_NEAR(m(i, j), 0.0, 1e-12) << i << " " << j;
    EXPECT_NEAR(trace(m), before, 1e-12);
}

//one xx-yy coupling is spread by the axis permutations over the three pairs of {xx,yy,zz}: 1/3 each
TEST(BondwiseSymmetrizeTests, CartesianDOffDiagonalSpreadsOverAxisPairs)
{
    dMatrix2 m = diagonal_matrix({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
    m(0, 1) = 1.0;
    m(1, 0) = 1.0;
    symmetrize_atomic_matrix_oh(m, { 2 }, false);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_NEAR(m(i, j), i == j ? 0.0 : 1.0 / 3.0, 1e-12) << i << " " << j;
    for (int i = 3; i < 6; i++)
        for (int j = 0; j < 6; j++) {
            EXPECT_NEAR(m(i, j), 0.0, 1e-12) << i << " " << j;
            EXPECT_NEAR(m(j, i), 0.0, 1e-12) << j << " " << i;
        }
}

//an s-p coupling changes sign under inversion, so the cross block between different l averages to zero
TEST(BondwiseSymmetrizeTests, CartesianCrossShellBlocksVanish)
{
    dMatrix2 m = diagonal_matrix({ 4.0, 1.0, 1.0, 1.0 });
    for (int i = 1; i < 4; i++) {
        m(0, i) = 0.7;
        m(i, 0) = 0.7;
    }
    symmetrize_atomic_matrix_oh(m, { 0, 1 }, false);
    EXPECT_NEAR(m(0, 0), 4.0, 1e-12);
    for (int i = 1; i < 4; i++) {
        EXPECT_NEAR(m(0, i), 0.0, 1e-12) << i;
        EXPECT_NEAR(m(i, 0), 0.0, 1e-12) << i;
        EXPECT_NEAR(m(i, i), 1.0, 1e-12) << i;
    }
}

//group averaging is a projector: a second pass changes nothing, the result is symmetric and the trace is kept
//(f and g shells, which also checks the type_vector exponent lookup for l = 3 and 4)
TEST(BondwiseSymmetrizeTests, CartesianFgAveragingIsIdempotentAndTracePreserving)
{
    const int n = 10 + 15;
    dMatrix2 m(n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            m(i, j) = std::cos(0.37 * i + 0.11 * j) + std::cos(0.37 * j + 0.11 * i);
    const double before = trace(m);
    symmetrize_atomic_matrix_oh(m, { 3, 4 }, false);
    dMatrix2 twice = m;
    symmetrize_atomic_matrix_oh(twice, { 3, 4 }, false);
    EXPECT_NEAR(trace(m), before, 1e-10);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            EXPECT_NEAR(m(i, j), m(j, i), 1e-10) << i << " " << j;
            EXPECT_NEAR(twice(i, j), m(i, j), 1e-10) << i << " " << j;
        }
}

//a diagonal f or g matrix averages within each O_h orbit of components: the orbit of a component is the
//set with the same multiset of exponents (xxx,yyy,zzz | xxy,xxz,... | xyz for f), and a diagonal input
//stays diagonal because every operation maps one component onto one component
TEST(BondwiseSymmetrizeTests, CartesianFgDiagonalAveragesOverExponentOrbits)
{
    for (const int l : { 3, 4 }) {
        const int n = (l + 1) * (l + 2) / 2;
        const int first_type = l * (l + 1) * (l + 2) / 6 + 1;
        std::vector<std::array<int, 3>> orbit(n);
        std::vector<double> d(n);
        for (int i = 0; i < n; i++) {
            int e[3];
            constants::type2vector(first_type + i, e);
            ASSERT_EQ(e[0] + e[1] + e[2], l) << l << " " << i;
            std::sort(e, e + 3);
            orbit[i] = { e[0], e[1], e[2] };
            d[i] = i + 1.0;
        }
        dMatrix2 m = diagonal_matrix(d);
        symmetrize_atomic_matrix_oh(m, { l }, false);
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            int members = 0;
            for (int j = 0; j < n; j++)
                if (orbit[j] == orbit[i]) {
                    sum += d[j];
                    members++;
                }
            EXPECT_NEAR(m(i, i), sum / members, 1e-12) << l << " " << i;
            for (int j = 0; j < n; j++)
                if (i != j)
                    EXPECT_NEAR(m(i, j), 0.0, 1e-12) << l << " " << i << " " << j;
        }
    }
}

//spherical p (libcint order py, pz, px) is one irrep too: diagonal averages to the mean
TEST(BondwiseSymmetrizeTests, SphericalPShellAveragesToScalar)
{
    dMatrix2 m = diagonal_matrix({ 1.0, 2.0, 3.0 });
    symmetrize_atomic_matrix_oh(m, { 1 }, true);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            EXPECT_NEAR(m(i, j), i == j ? 2.0 : 0.0, 1e-10) << i << " " << j;
}

//spherical d in libcint order (dxy, dyz, dz2, dxz, dx2-y2): t2g at 0,1,3 -> 7/3, eg at 2,4 -> 4
TEST(BondwiseSymmetrizeTests, SphericalDShellSplitsIntoT2gAndEg)
{
    dMatrix2 m = diagonal_matrix({ 1.0, 2.0, 3.0, 4.0, 5.0 });
    symmetrize_atomic_matrix_oh(m, { 2 }, true);
    const double expected[5] = { 7.0 / 3.0, 7.0 / 3.0, 4.0, 7.0 / 3.0, 4.0 };
    for (int i = 0; i < 5; i++)
        for (int j = 0; j < 5; j++)
            EXPECT_NEAR(m(i, j), i == j ? expected[i] : 0.0, 1e-10) << i << " " << j;
    EXPECT_NEAR(trace(m), 15.0, 1e-10);
}

//shell offsets: an s + p + d spherical atom (9 functions) keeps s, averages p and d blocks in place
//and zeroes every cross block
TEST(BondwiseSymmetrizeTests, SphericalMultiShellOffsets)
{
    std::vector<double> d(9);
    for (int i = 0; i < 9; i++)
        d[i] = i + 1.0;
    dMatrix2 m = diagonal_matrix(d);
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++)
            if (i != j)
                m(i, j) = 0.25;
    symmetrize_atomic_matrix_oh(m, { 0, 1, 2 }, true);
    EXPECT_NEAR(m(0, 0), 1.0, 1e-10);
    for (int i = 1; i < 4; i++)
        EXPECT_NEAR(m(i, i), 3.0, 1e-10) << i;
    const double d_expect[5] = { (5.0 + 6.0 + 8.0) / 3.0, (5.0 + 6.0 + 8.0) / 3.0, 8.0, (5.0 + 6.0 + 8.0) / 3.0, 8.0 };
    for (int i = 0; i < 5; i++)
        EXPECT_NEAR(m(4 + i, 4 + i), d_expect[i], 1e-10) << i;
    for (int i = 0; i < 9; i++)
        for (int j = 0; j < 9; j++)
            if (i != j)
                EXPECT_NEAR(m(i, j), 0.0, 1e-10) << i << " " << j;
    EXPECT_NEAR(trace(m), 45.0, 1e-10);
}

//two s shells are invariant under every operation: the full 2x2 matrix, off-diagonal included, is untouched
TEST(BondwiseSymmetrizeTests, SOnlyMatrixIsUnchanged)
{
    for (const bool spherical : { false, true }) {
        dMatrix2 m = diagonal_matrix({ 1.5, -0.5 });
        m(0, 1) = 0.3;
        m(1, 0) = 0.3;
        symmetrize_atomic_matrix_oh(m, { 0, 0 }, spherical);
        EXPECT_NEAR(m(0, 0), 1.5, 1e-12);
        EXPECT_NEAR(m(1, 1), -0.5, 1e-12);
        EXPECT_NEAR(m(0, 1), 0.3, 1e-12);
        EXPECT_NEAR(m(1, 0), 0.3, 1e-12);
    }
}

//an empty shell list on an empty matrix is a no-op rather than an error
TEST(BondwiseSymmetrizeTests, EmptyMatrixIsNoOp)
{
    dMatrix2 m(0, 0);
    symmetrize_atomic_matrix_oh(m, {}, false);
    EXPECT_EQ(m.extent(0), 0u);
    symmetrize_atomic_matrix_oh(m, {}, true);
    EXPECT_EQ(m.extent(0), 0u);
}

//selecting the first atom keeps exactly its basin. The density is separable, f(x) g(y,z), with f sampled
//at x = -2 + 0.5 i: f(5) = 1.14902 > f(6) = 1.14668 > f(7) = 1.13957 < f(8) = f(6) (the merge persistence
//of 0.8 % is above the 0.5 % threshold), so the grid maxima are columns 5 and 9 and an ascent from a
//column <= 6 can never cross column 7. Columns 0..6 (7 x 49 = 343 voxels) are therefore all kept, nothing
//beyond column 7 is, and only the exactly-tied column 7 is left to the tie-break. The written cube keeps
//the comments and the full-grid origin because column 0 is kept.
TEST(BondwiseQtaimMaskTests, SelectedAtomKeepsOnlyItsBasin)
{
    WFN parent(e_origin::NOT_YET_DEFINED);
    cube rho, eli;
    make_pair_grid(parent, rho, eli);
    const auto out = bondwise_temp("mask_first.cube");
    std::ostringstream log;
    QTAIM_ELI_mask(rho, eli, parent, parent.get_atoms(), { 0 }, -1.0, out, false, log);
    ASSERT_TRUE(std::filesystem::exists(out));
    const MaskedCube m = read_masked(out, -1.0, log);
    std::filesystem::remove(out);

    EXPECT_EQ(m.comment1, "QTAIM-masked ELI");
    EXPECT_EQ(m.comment2, "Selected atoms: 0");
    EXPECT_EQ(m.size[1], grid_ny);
    EXPECT_EQ(m.size[2], grid_nz);
    EXPECT_GE(m.size[0], mid_ix);
    EXPECT_LE(m.size[0], mid_ix + 1);
    EXPECT_NEAR(m.origin[0], -2.0, 1e-6);
    EXPECT_NEAR(m.origin[1], -1.5, 1e-6);
    EXPECT_NEAR(m.origin[2], -1.5, 1e-6);
    int left_of_mid = 0;
    bool own_voxel = false;
    for (const auto& k : m.kept) {
        EXPECT_LE(k[0], mid_ix);
        if (k[0] < mid_ix)
            left_of_mid++;
        if (k[0] == atom_a_ix && k[1] == atom_iy && k[2] == atom_iz)
            own_voxel = true;
    }
    EXPECT_EQ(left_of_mid, mid_ix * grid_ny * grid_nz);
    EXPECT_EQ(static_cast<int>(m.kept.size()) + m.background_voxels, m.size[0] * m.size[1] * m.size[2]);
    EXPECT_TRUE(own_voxel);
}

//the second atom selected: the mirror image, columns 8..14 all kept, nothing below column 7, and the
//shrunk origin moves to the first kept x column. Together with the first atom's basin the two selections
//partition the grid: every voxel is kept by exactly one of them
TEST(BondwiseQtaimMaskTests, SelectedSecondAtomShiftsOriginAndPartitionsGrid)
{
    WFN parent(e_origin::NOT_YET_DEFINED);
    cube rho, eli;
    make_pair_grid(parent, rho, eli);
    const auto out = bondwise_temp("mask_second.cube");
    std::ostringstream log;
    QTAIM_ELI_mask(rho, eli, parent, parent.get_atoms(), { 1 }, 0.0, out, false, log);
    const MaskedCube m = read_masked(out, 0.0, log);
    std::filesystem::remove(out);

    EXPECT_EQ(m.comment2, "Selected atoms: 1");
    int xmin = grid_nx, right_of_mid = 0;
    bool own_voxel = false;
    for (const auto& k : m.kept) {
        EXPECT_GE(k[0], mid_ix);
        if (k[0] > mid_ix)
            right_of_mid++;
        xmin = std::min(xmin, k[0]);
        if (k[0] == atom_b_ix && k[1] == atom_iy && k[2] == atom_iz)
            own_voxel = true;
    }
    EXPECT_TRUE(own_voxel);
    EXPECT_EQ(right_of_mid, (grid_nx - mid_ix - 1) * grid_ny * grid_nz);
    EXPECT_GE(xmin, mid_ix);
    EXPECT_LE(xmin, mid_ix + 1);
    EXPECT_NEAR(m.origin[0], -2.0 + 0.5 * xmin, 1e-6);
    EXPECT_EQ(m.size[0], grid_nx - xmin);
    EXPECT_EQ(static_cast<int>(m.kept.size()) + m.background_voxels, m.size[0] * m.size[1] * m.size[2]);

    const auto out_first = bondwise_temp("mask_first_again.cube");
    QTAIM_ELI_mask(rho, eli, parent, parent.get_atoms(), { 0 }, 0.0, out_first, false, log);
    const MaskedCube first = read_masked(out_first, 0.0, log);
    std::filesystem::remove(out_first);
    std::vector<int> owners(grid_nx * grid_ny * grid_nz, 0);
    for (const auto& k : first.kept)
        owners[(k[0] * grid_ny + k[1]) * grid_nz + k[2]]++;
    for (const auto& k : m.kept)
        owners[(k[0] * grid_ny + k[1]) * grid_nz + k[2]]++;
    EXPECT_EQ(std::count(owners.begin(), owners.end(), 1), static_cast<std::ptrdiff_t>(owners.size()));
}

//no basin matches an empty selection: the output is a single background voxel and a warning is logged
TEST(BondwiseQtaimMaskTests, NoSelectionGivesSingleBackgroundVoxel)
{
    WFN parent(e_origin::NOT_YET_DEFINED);
    cube rho, eli;
    make_pair_grid(parent, rho, eli);
    const auto out = bondwise_temp("mask_none.cube");
    std::ostringstream log;
    QTAIM_ELI_mask(rho, eli, parent, parent.get_atoms(), {}, 7.5, out, false, log);
    EXPECT_NE(log.str().find("no QTAIM basins matched"), std::string::npos);
    const MaskedCube m = read_masked(out, 7.5, log);
    std::filesystem::remove(out);
    EXPECT_EQ(m.size[0], 1);
    EXPECT_EQ(m.size[1], 1);
    EXPECT_EQ(m.size[2], 1);
    EXPECT_EQ(m.background_voxels, 1);
    EXPECT_TRUE(m.kept.empty());
    EXPECT_EQ(m.comment2, "Selected atoms: ");
}

//both atoms selected: every voxel belongs to a kept basin, so the output is the full eli grid unchanged
TEST(BondwiseQtaimMaskTests, AllAtomsSelectedKeepsWholeGrid)
{
    WFN parent(e_origin::NOT_YET_DEFINED);
    cube rho, eli;
    make_pair_grid(parent, rho, eli);
    const auto out = bondwise_temp("mask_all.cube");
    std::ostringstream log;
    QTAIM_ELI_mask(rho, eli, parent, parent.get_atoms(), { 0, 1 }, -1.0, out, true, log);
    EXPECT_NE(log.str().find("selected."), std::string::npos);
    const MaskedCube m = read_masked(out, -1.0, log);
    std::filesystem::remove(out);
    EXPECT_EQ(m.size[0], grid_nx);
    EXPECT_EQ(m.size[1], grid_ny);
    EXPECT_EQ(m.size[2], grid_nz);
    EXPECT_EQ(m.background_voxels, 0);
    ASSERT_EQ(static_cast<int>(m.kept.size()), grid_nx * grid_ny * grid_nz);
    EXPECT_EQ(m.comment2, "Selected atoms: 0,1");
    //decoded indices come back in storage order, which proves each voxel kept its own eli value
    int idx = 0;
    for (int x = 0; x < grid_nx; x++)
        for (int y = 0; y < grid_ny; y++)
            for (int z = 0; z < grid_nz; z++, idx++) {
                EXPECT_EQ(m.kept[idx][0], x);
                EXPECT_EQ(m.kept[idx][1], y);
                EXPECT_EQ(m.kept[idx][2], z);
            }
}

//cube-files mode: rho and eli read from disk, atoms taken from the cube header, output placed next to the
//eli file as <stem>_qtaim_masked.cube with the selected basin's voxels
TEST(BondwiseQtaimMaskTests, RunCubeFilesModeWritesNextToEli)
{
    WFN parent(e_origin::NOT_YET_DEFINED);
    cube rho, eli;
    make_pair_grid(parent, rho, eli);
    const auto rho_path = bondwise_temp("run_rho.cube");
    const auto eli_path = bondwise_temp("run_eli.cube");
    rho.set_path(rho_path);
    eli.set_path(eli_path);
    ASSERT_TRUE(rho.write_file(true));
    ASSERT_TRUE(eli.write_file(true));
    const auto out = bondwise_temp("run_eli_qtaim_masked.cube");

    options opt;
    std::ostringstream log;
    run_QTAIM_ELI_mask(rho_path, eli_path, { 1 }, -1.0, opt, log);
    std::filesystem::remove(rho_path);
    std::filesystem::remove(eli_path);
    ASSERT_TRUE(std::filesystem::exists(out));
    EXPECT_NE(log.str().find("Reading density cube"), std::string::npos);
    const MaskedCube m = read_masked(out, -1.0, log);
    std::filesystem::remove(out);

    EXPECT_EQ(m.comment2, "Selected atoms: 1");
    //same partition as SelectedSecondAtomShiftsOriginAndPartitionsGrid: every column right of the midplane
    int right_of_mid = 0;
    bool own_voxel = false;
    for (const auto& k : m.kept) {
        EXPECT_GE(k[0], mid_ix);
        if (k[0] > mid_ix)
            right_of_mid++;
        if (k[0] == atom_b_ix && k[1] == atom_iy && k[2] == atom_iz)
            own_voxel = true;
    }
    EXPECT_EQ(right_of_mid, (grid_nx - mid_ix - 1) * grid_ny * grid_nz);
    EXPECT_TRUE(own_voxel);
}

//the eli values written by the mask survive the cube text format exactly (integers up to 1500 in 6 digits),
//which the decoding in the other tests relies on
TEST(BondwiseQtaimMaskTests, EncodedEliValuesRoundTripThroughCubeFile)
{
    WFN parent(e_origin::NOT_YET_DEFINED);
    cube rho, eli;
    make_pair_grid(parent, rho, eli);
    const auto p = bondwise_temp("roundtrip_eli.cube");
    eli.set_path(p);
    ASSERT_TRUE(eli.write_file(true));
    WFN dummy;
    std::ostringstream log;
    cube back(p, true, dummy, log);
    std::filesystem::remove(p);
    ASSERT_EQ(back.get_size(0), grid_nx);
    EXPECT_EQ(dummy.get_ncen(), 2);
    EXPECT_NEAR(dummy.get_atom_coordinate(1, 0), 3.0, 1e-6);
    for (int x = 0; x < grid_nx; x += 7)
        for (int y = 0; y < grid_ny; y += 3)
            for (int z = 0; z < grid_nz; z += 2)
                EXPECT_DOUBLE_EQ(back.get_value(x, y, z), encoded_eli(x, y, z));
}

//a missing input file makes autobonds write input.example into the working directory and return 0
TEST(BondwiseAutobondsTests, MissingInputWritesExample)
{
    const auto dir = bondwise_temp("example_dir");
    std::filesystem::create_directories(dir);
    const auto old_cwd = std::filesystem::current_path();
    std::filesystem::current_path(dir);
    WFN wavy = three_atom_wfn();
    int rc;
    std::string out;
    {
        CoutCapture cap;
        rc = autobonds(false, wavy, dir / "does_not_exist.inp", true);
        out = cap.str();
    }
    std::filesystem::current_path(old_cwd);
    EXPECT_EQ(rc, 0);
    EXPECT_NE(out.find("input.example"), std::string::npos);
    const auto example = dir / "input.example";
    ASSERT_TRUE(std::filesystem::exists(example));
    std::ifstream f(example);
    std::string first, line;
    std::getline(f, first);
    int lines = 1;
    while (std::getline(f, line))
        lines++;
    f.close();
    EXPECT_EQ(first.rfind("!COMMENT", 0), 0u);
    EXPECT_GT(lines, 20);
    std::filesystem::remove_all(dir);
}

//a header shorter than ten characters is rejected before parsing
TEST(BondwiseAutobondsTests, ShortHeaderReturnsZero)
{
    const auto r = run_autobonds("short", "1 1 1 1\n 2 0 0 5 5 5 1 2 2 2 1 2 3\n");
    EXPECT_EQ(r.first, 0);
    EXPECT_EQ(r.second.find("Finished all calculations"), std::string::npos);
}

//a long header that does not hold four integers is rejected too
TEST(BondwiseAutobondsTests, UnparsableHeaderReturnsZero)
{
    const auto r = run_autobonds("unparsable", "rho rdg eli lap switches\n 2 0 0 5 5 5 1 2 2 2 1 2 3\n");
    EXPECT_EQ(r.first, 0);
}

//atom out of range, repeated atoms and mode_sel outside 1..4 are each counted as one failed run
TEST(BondwiseAutobondsTests, InvalidAtomSelectionsCountAsFailures)
{
    const std::string text = std::string(autobonds_header) +
        " 2 0 0 5 5 5 1 2 2 2 1 2 4\n" //atom 4 of 3
        " 2 0 0 5 5 5 1 2 2 2 1 1 3\n" //atom1 == atom2
        " 0 0 0 5 5 5 1 2 2 2 1 2 3\n" //mode_sel 0
        " 5 0 0 5 5 5 1 2 2 2 1 2 3\n" //mode_sel 5
        " 2 0 0 5 5 5 1 2 2 2 0 2 3\n"; //atom 0
    const auto r = run_autobonds("invalid_atoms", text);
    EXPECT_EQ(r.first, 1);
    EXPECT_EQ(count_occurrences(r.second, "Invalid selections of atoms or mode_sel"), 5);
    EXPECT_EQ(count_occurrences(r.second, "problem somewhere during the calculations"), 5);
    EXPECT_NE(r.second.find("0 out of 5 were successful"), std::string::npos);
}

//a negative resolution entry is refused after the geometry was set up
TEST(BondwiseAutobondsTests, NegativeResolutionRejected)
{
    const auto r = run_autobonds("neg_res", std::string(autobonds_header) + " 2 0 0 5 -5 5 1 2 2 2 1 2 3\n");
    EXPECT_EQ(r.first, 1);
    EXPECT_NE(r.second.find("Wrong input in res!"), std::string::npos);
    EXPECT_NE(r.second.find("0 out of 1 were successful"), std::string::npos);
}

//a negative box scaling is refused
TEST(BondwiseAutobondsTests, NegativeBoxRejected)
{
    const auto r = run_autobonds("neg_box", std::string(autobonds_header) + " 2 0 0 5 5 5 1 2 -2 2 1 2 3\n");
    EXPECT_EQ(r.first, 1);
    EXPECT_NE(r.second.find("Wrong input for box scaling!"), std::string::npos);
}

//a box of 50 or more bond lengths is refused with the "be realistic" message
TEST(BondwiseAutobondsTests, OversizedBoxRejected)
{
    const auto r = run_autobonds("big_box", std::string(autobonds_header) + " 3 0 0 5 5 5 0 50 2 2 1 2 3\n");
    EXPECT_EQ(r.first, 1);
    EXPECT_NE(r.second.find("Come on, be realistic!"), std::string::npos);
    EXPECT_NE(r.second.find("0 out of 1 were successful"), std::string::npos);
}

//comment lines before the header and short lines between entries are skipped, so only one run is counted
TEST(BondwiseAutobondsTests, CommentAndShortLinesAreSkipped)
{
    const std::string text = "!comment one\n!comment two\n" + std::string(autobonds_header) +
        "\n 1 2 3\n 4 0 0 5 5 5 1 2 2 2 1 2 4\n\n";
    const auto r = run_autobonds("comments", text);
    EXPECT_EQ(r.first, 1);
    EXPECT_NE(r.second.find("0 out of 1 were successful"), std::string::npos);
    EXPECT_EQ(count_occurrences(r.second, "Invalid selections"), 1);
}

//every mode 1..4 reaches the resolution/box validation, i.e. the orientation vectors are built for each
TEST(BondwiseAutobondsTests, EveryModeReachesBoxValidation)
{
    std::string text = autobonds_header;
    for (int mode = 1; mode <= 4; mode++)
        text += " " + std::to_string(mode) + " 1 1 0.5 0.5 0.5 0 60 2 2 1 2 3\n";
    const auto r = run_autobonds("modes", text);
    EXPECT_EQ(r.first, 1);
    EXPECT_EQ(count_occurrences(r.second, "Come on, be realistic!"), 4);
    EXPECT_EQ(count_occurrences(r.second, "Invalid selections"), 0);
    EXPECT_NE(r.second.find("0 out of 4 were successful"), std::string::npos);
}

//a valid bond line writes the rho cube with the wavefunction's atoms in its header
TEST(BondwiseAutobondsTests, ValidBondWritesRhoCube)
{
    WFN wavy = three_atom_wfn();
    wavy.push_back_MO(1, 2.0, -0.5);
    std::vector<primitive> prims;
    prims.emplace_back(1, 1, 1.0, 1.0);
    wavy.push_back_spherical_shell(0, 0, vec2{ vec{ 1.0 } }, prims, 0, 1);
    wavy.set_exp_cutoff();
    const auto stem = bondwise_temp("valid_bond");
    wavy.set_path(stem);
    const auto input = bondwise_temp("valid_bond.inp");
    write_text(input, " 1    0    0    0\n 2 0 1 5 5 5 1 2 2 2 1 2 3\n"); //mres=1: 5 grid points per axis
    int rc;
    std::string out;
    {
        CoutCapture cap;
        rc = autobonds(false, wavy, input, true);
        out = cap.str();
    }
    std::filesystem::remove(input);
    EXPECT_EQ(rc, 1);
    EXPECT_NE(out.find("1 out of 1 were successful"), std::string::npos);
    const auto rho_file = std::filesystem::path(stem.generic_string() + "_C1_N2_O3_1_rho.cube");
    ASSERT_TRUE(std::filesystem::exists(rho_file));
    WFN dummy;
    std::ostringstream log;
    cube back(rho_file, true, dummy, log);
    std::filesystem::remove(rho_file);
    EXPECT_EQ(dummy.get_ncen(), 3);
    EXPECT_EQ(back.get_size(0), 5);
    EXPECT_EQ(back.get_size(1), 5);
    EXPECT_EQ(back.get_size(2), 5);
    EXPECT_GT(back.sum(), 0.0);
}

namespace
{
    std::filesystem::path nh3li_fixture()
    {
        const auto root = nos_test_repo_root();
        if (root.empty())
            return {};
        const auto p = root / "tests" / "RGBI_groups" / "nh3li.gbw";
        return std::filesystem::exists(p) ? p : std::filesystem::path{};
    }

    //loads nh3li.gbw and runs the Roby analysis with captured stdout; empty string when the fixture is absent
    std::string roby_output(const ivec3& groups, bool symmetrize, bool ano, bool EVs, bool theta)
    {
        const auto p = nh3li_fixture();
        if (p.empty())
            return {};
        CoutCapture cap;
        WFN wavy(p);
        Roby_information roby(wavy, groups, symmetrize, ano, EVs, theta);
        return cap.str();
    }
}

//the NAO populations and the printed total reproduce tests/RGBI/nh3li_nao.good (the atom rows sum to 16.84,
//the total counts 12.92 electrons, so the two are checked separately)
TEST(BondwiseRobyTests, NaoPopulationsMatchGolden)
{
    const std::string out = roby_output({}, true, false, false, false);
    if (out.empty())
        GTEST_SKIP() << "tests/RGBI_groups/nh3li.gbw not found";
    const double golden[5] = { 9.42047, 1.4373, 1.4353, 1.43756, 3.1097 };
    for (int i = 0; i < 5; i++)
        EXPECT_NEAR(value_after(out, "Population of atom " + std::to_string(i) + ": "), golden[i], 2e-3) << i;
    EXPECT_NEAR(value_after(out, "Total Population: "), 12.9218, 2e-3);
}

//the N-Li and N-H rows match the golden table, and the printed Tot. and Pyth. columns follow from Cov. and Ion.
TEST(BondwiseRobyTests, NaoBondTableMatchesGolden)
{
    const std::string out = roby_output({}, true, false, false, false);
    if (out.empty())
        GTEST_SKIP() << "tests/RGBI_groups/nh3li.gbw not found";
    const std::vector<double> li = row_numbers_after(out, "N - Li");
    ASSERT_EQ(li.size(), 9u);
    const double golden_li[9] = { 9.420, 3.110, 12.393, 0.137, 0.184, 0.421, 0.459, 15.972, 26.173 };
    for (int i = 0; i < 9; i++)
        EXPECT_NEAR(li[i], golden_li[i], 3e-3) << i;
    EXPECT_NEAR(li[6], std::sqrt(li[4] * li[4] + li[5] * li[5]), 2e-3);
    EXPECT_NEAR(li[7], 100.0 * li[4] * li[4] / (li[6] * li[6]), 0.3);
    EXPECT_NEAR(li[8], 200.0 * std::asin(li[4] / li[6]) / constants::PI, 0.3);

    const std::vector<double> h = row_numbers_after(out, "N -  H");
    ASSERT_EQ(h.size(), 9u);
    const double golden_h[9] = { 9.420, 1.437, 9.615, 1.243, 0.905, 0.296, 0.952, 90.322, 79.861 };
    for (int i = 0; i < 9; i++)
        EXPECT_NEAR(h[i], golden_h[i], 3e-3) << i;
    EXPECT_EQ(count_occurrences(out, "N -  H"), 3);
}

//theta_info prints one theta-subspace table per bonded pair; each row's Total is sqrt(Cov^2 + Ion^2) and the
//rows summed with the bond loop's rules (pair != i, 0.573 deg < theta < 89.427 deg, i.e. 1e-2 rad) reproduce
//the Cov. and Ion. of the RGBI table, which is computed on an independent path from the same populations
TEST(BondwiseRobyTests, ThetaInfoReportsEveryBond)
{
    const std::string out = roby_output({}, true, false, false, true);
    if (out.empty())
        GTEST_SKIP() << "tests/RGBI_groups/nh3li.gbw not found";
    EXPECT_NE(out.find("RGBI theta-subspace reports enabled."), std::string::npos);
    const std::string section_key = "Roby-Gould theta subspaces for";
    EXPECT_EQ(count_occurrences(out, section_key), 4);
    EXPECT_EQ(count_occurrences(out, " Pair    theta/degrees"), 4);
    //the sections come in bond-loop order (0,1), (0,2), (0,3), (0,4); the table rows are looked up by index
    //because the RGBI table is sorted after the loop
    const char* const table_keys[4] = { "   0 -   1    N -  H", "   0 -   2    N -  H", "   0 -   3    N -  H",
                                        "   0 -   4    N - Li" };
    const double cutoff = 1e-2 * constants::INV_PI_180;
    size_t pos = out.find(section_key);
    int rows = 0;
    for (int b = 0; b < 4 && pos != std::string::npos; b++) {
        const size_t next = out.find(section_key, pos + section_key.size());
        std::istringstream lines(out.substr(pos, next == std::string::npos ? std::string::npos : next - pos));
        std::string line;
        double cov_sum = 0.0, ion_sum = 0.0;
        int summed = 0;
        while (std::getline(lines, line)) {
            //rows look like "   1,  2          12.345    0.500    0.100     0.200    0.300    0.100     0.100     0.224"
            const size_t comma = line.find(',');
            if (line.size() < 90 || comma == std::string::npos || comma > 6)
                continue;
            std::istringstream iss(line);
            int i, pair;
            char c;
            double theta, cp, cm, cov, ip, im, ion, total;
            if (!(iss >> i >> c >> pair >> theta >> cp >> cm >> cov >> ip >> im >> ion >> total))
                continue;
            rows++;
            EXPECT_GE(theta, -1e-3) << line;
            EXPECT_LE(theta, 90.0 + 1e-3) << line;
            EXPECT_NEAR(total, std::sqrt(cov * cov + ion * ion), 2e-3) << line;
            if (pair != i) {
                EXPECT_NEAR(cov, 0.5 * (cp - cm), 1.5e-3) << line;
                EXPECT_NEAR(ion, 0.5 * (ip - im), 1.5e-3) << line;
                if (theta > cutoff && theta < 90.0 - cutoff) {
                    cov_sum += cov;
                    ion_sum += ion;
                    summed++;
                }
            }
        }
        const std::vector<double> row = row_numbers_after(out, table_keys[b]);
        ASSERT_EQ(row.size(), 9u) << table_keys[b];
        //the printed rows carry three decimals, so half a unit in the last place per summed row
        EXPECT_NEAR(row[4], cov_sum, 5e-4 * summed + 2e-3) << table_keys[b];
        EXPECT_NEAR(row[5], ion_sum, 5e-4 * summed + 2e-3) << table_keys[b];
        pos = next;
    }
    EXPECT_GT(rows, 4);
    //the populations are unaffected by the extra report
    EXPECT_NEAR(value_after(out, "Population of atom 0: "), 9.42047, 2e-3);
}

//EVs=true prints the unsorted projected-density eigenvalues per atom without changing the numbers
TEST(BondwiseRobyTests, EigenvaluePrintsLeavePopulationsUnchanged)
{
    const std::string out = roby_output({}, true, false, true, false);
    if (out.empty())
        GTEST_SKIP() << "tests/RGBI_groups/nh3li.gbw not found";
    EXPECT_GE(count_occurrences(out, "Eigenvalues of projected density P (unsorted):"), 5);
    EXPECT_NE(out.find("theta_I after Ionic"), std::string::npos);
    EXPECT_NEAR(value_after(out, "Population of atom 0: "), 9.42047, 2e-3);
    EXPECT_NEAR(value_after(out, "Population of atom 4: "), 3.1097, 2e-3);
}

//singleton groups {N} and {Li} take the partial-basis projection path and must give the atom populations
//and the N-Li pair population of the plain table; the group Tot. is sqrt(Cov^2 + Ion^2)
TEST(BondwiseRobyTests, SingletonGroupsReproduceAtomPopulations)
{
    const std::string out = roby_output({ { { 0 }, { 4 } } }, true, false, false, false);
    if (out.empty())
        GTEST_SKIP() << "tests/RGBI_groups/nh3li.gbw not found";
    EXPECT_NE(out.find("G0: atoms 0(N)"), std::string::npos);
    EXPECT_NE(out.find("G1: atoms 4(Li)"), std::string::npos);
    EXPECT_NEAR(value_after(out, "Population of G0 (N)"), 9.42047, 2e-3);
    EXPECT_NEAR(value_after(out, "Population of G1 (Li)"), 3.1097, 2e-3);
    const std::vector<double> g = row_numbers_after(out, "G0 - G1");
    ASSERT_EQ(g.size(), 9u);
    EXPECT_NEAR(g[0], 9.420, 3e-3);
    EXPECT_NEAR(g[1], 3.110, 3e-3);
    EXPECT_NEAR(g[2], 12.393, 5e-3);
    //s = n_G1 + n_G2 - n_G1G2 is the golden 0.137 of the N-Li row
    EXPECT_NEAR(g[3], 0.137, 3e-3);
    EXPECT_NEAR(g[3], g[0] + g[1] - g[2], 3e-3);
    EXPECT_NEAR(g[6], std::sqrt(g[4] * g[4] + g[5] * g[5]), 2e-3);
    EXPECT_GT(g[6], 0.0);
}

//NH3 vs Li as groups: G0 is the projection onto the union of the four atomic subspaces, so it is bounded
//by the N population from below and by the electron count (a projector expectation of the density) from
//above; the atom sum (13.73) is not a valid expectation because the atomic subspaces overlap. The pair
//spans the full basis, so n_G1G2 is the total population and s = n_G1 + n_G2 - n_G1G2 holds exactly
TEST(BondwiseRobyTests, Nh3LiGroupPairSpansFullBasis)
{
    const std::string out = roby_output({ { { 0, 1, 2, 3 }, { 4 } } }, true, false, false, false);
    if (out.empty())
        GTEST_SKIP() << "tests/RGBI_groups/nh3li.gbw not found";
    const double g0 = value_after(out, "Population of G0 (NH3)");
    EXPECT_GT(g0, 9.42047);
    EXPECT_LE(g0, 13.0 + 1e-3);
    EXPECT_NEAR(value_after(out, "Population of G1 (Li)"), 3.1097, 2e-3);
    const std::vector<double> g = row_numbers_after(out, "G0 - G1");
    ASSERT_EQ(g.size(), 9u);
    EXPECT_NEAR(g[0], g0, 3e-3);
    EXPECT_NEAR(g[1], 3.110, 3e-3);
    EXPECT_NEAR(g[2], value_after(out, "Total Population: "), 5e-3);
    EXPECT_NEAR(g[3], g[0] + g[1] - g[2], 3e-3);
    EXPECT_NEAR(g[6], std::sqrt(g[4] * g[4] + g[5] * g[5]), 2e-3);
}

//switching symmetrization off changes the NAOs but not the electron count: every population is a projector
//expectation of the 13-electron density, so the total can never exceed 13 and must stay close to it, and
//the same three N-H bonds are found. No golden: tests/RGBI_groups/NoSpherA2_no_sym.log is not tracked
TEST(BondwiseRobyTests, NoSymmetrizeKeepsTotalPopulation)
{
    const std::string out = roby_output({}, false, false, false, false);
    if (out.empty())
        GTEST_SKIP() << "tests/RGBI_groups/nh3li.gbw not found";
    const double total = value_after(out, "Total Population: ");
    EXPECT_GT(total, 12.8);
    EXPECT_LE(total, 13.0 + 1e-3);
    const double n = value_after(out, "Population of atom 0: ");
    const double li = value_after(out, "Population of atom 4: ");
    EXPECT_NEAR(n, 9.42047, 0.3);
    EXPECT_NEAR(li, 3.1097, 0.3);
    EXPECT_LT(n + li, total + 1e-3); //the three H atoms carry the rest
    EXPECT_EQ(count_occurrences(out, "N -  H"), 3);
}

//the ANO basis runs an atomic SCF per atom; on nh3li every atom converges (no fallback) and the populations
//reproduce tests/RGBI/nh3li_ano.good. The five atomic SCFs are the slowest part of this file; if the run
//exceeds the ~2 s budget this is the test to move to the integration suite
TEST(BondwiseRobyTests, AnoBasisMatchesGoldenWithoutFallback)
{
    const std::string out = roby_output({}, true, true, false, false);
    if (out.empty())
        GTEST_SKIP() << "tests/RGBI_groups/nh3li.gbw not found";
    EXPECT_NE(out.find("Calculating ANOs for all atoms"), std::string::npos);
    EXPECT_NE(out.find("ANO fallback summary: no atom-level fallbacks were needed."), std::string::npos);
    const double golden[5] = { 9.37397, 1.44853, 1.44718, 1.44872, 3.10977 };
    for (int i = 0; i < 5; i++)
        EXPECT_NEAR(value_after(out, "Population of atom " + std::to_string(i) + ": "), golden[i], 5e-3) << i;
    EXPECT_NEAR(value_after(out, "Total Population: "), 12.8475, 5e-3);
    const std::vector<double> h = row_numbers_after(out, "N -  H");
    ASSERT_EQ(h.size(), 9u);
    EXPECT_NEAR(h[4], 0.885, 5e-3);
    EXPECT_NEAR(h[5], 0.295, 5e-3);
}
