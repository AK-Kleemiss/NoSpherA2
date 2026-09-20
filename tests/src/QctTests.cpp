#include "pch.h"

#include "core/convenience.h"
#include "core/wfn_class.h"

int QCT(options& opt, std::vector<WFN>& wavy);

namespace
{
    //feeds the menu a script through std::cin and captures what it prints; both streams go back when the test ends
    struct MenuRun
    {
        std::istringstream script;
        std::ostringstream out;
        std::streambuf* old_in;
        std::streambuf* old_out;
        explicit MenuRun(const std::string& keys)
            : script(keys), old_in(std::cin.rdbuf(script.rdbuf())), old_out(std::cout.rdbuf(out.rdbuf())) {}
        ~MenuRun()
        {
            std::cin.rdbuf(old_in);
            std::cin.clear(); //the script ends in EOF, the failbit would otherwise stick to the real cin
            std::cout.rdbuf(old_out);
        }
        std::string str() const { return out.str(); }
    };

    //a scratch directory named after the test; the cwd goes back and the directory is removed when the test passed
    struct Scratch
    {
        std::filesystem::path dir;
        std::filesystem::path old_cwd;
        explicit Scratch(const std::string& name) : old_cwd(std::filesystem::current_path())
        {
            dir = std::filesystem::temp_directory_path() / ("NoSpherA2_Qct_" + name);
            std::filesystem::remove_all(dir);
            std::filesystem::create_directories(dir);
            std::filesystem::current_path(dir);
        }
        ~Scratch()
        {
            std::error_code ec;
            std::filesystem::current_path(old_cwd, ec);
            if (!::testing::Test::HasFailure())
                std::filesystem::remove_all(dir, ec);
        }
    };

    //H2 with the bonding and antibonding s combination both doubly occupied (as in BondwiseCoverageTests)
    std::filesystem::path write_h2(const std::filesystem::path& dir)
    {
        WFN wavy(e_origin::NOT_YET_DEFINED);
        wavy.push_back_atom("H", -1.0, 0.0, 0.0, 1);
        wavy.push_back_atom("H", 1.0, 0.0, 0.0, 1);
        wavy.push_back_MO(1, 2.0, -0.5);
        wavy.push_back_MO(2, 2.0, -0.4);
        std::vector<primitive> prims;
        prims.emplace_back(1, 1, 1.0, 1.0);
        prims.emplace_back(2, 1, 1.0, 1.0);
        for (int mo = 0; mo < 2; mo++)
        {
            wavy.push_back_spherical_shell(mo, 0, vec2{ vec{ 1.0 } }, prims, 0, 1);
            wavy.push_back_spherical_shell(mo, 0, vec2{ vec{ mo == 0 ? 1.0 : -1.0 } }, prims, 1, 1);
        }
        const auto path = dir / "h2.wfn";
        EXPECT_TRUE(wavy.write_wfn(path, false, false));
        return path;
    }

    int run(const std::string& keys, options& opt, std::vector<WFN>& wavy, std::string& printed)
    {
        MenuRun m(keys);
        const int code = QCT(opt, wavy);
        printed = m.str();
        return code;
    }
}

//a closed stdin must end the menu, not spin on the default branch
TEST(QctTests, ClosedInputLeavesTheMenu)
{
    options opt;
    std::vector<WFN> wavy;
    std::string printed;
    EXPECT_EQ(run("", opt, wavy, printed), 0);
    EXPECT_NE(printed.find("NoSpherA2 -- QCT interactive menu"), std::string::npos);
    EXPECT_NE(printed.find("No wavefunction loaded"), std::string::npos);
    EXPECT_NE(printed.find("Input closed"), std::string::npos);
}

TEST(QctTests, QuitKeyAndUnknownKey)
{
    options opt;
    std::vector<WFN> wavy;
    std::string printed;
    EXPECT_EQ(run("zzz\nq\n", opt, wavy, printed), 0);
    EXPECT_NE(printed.find("Sorry, I did not get that"), std::string::npos);
    EXPECT_NE(printed.find("Bye!"), std::string::npos);
    EXPECT_EQ(printed.find("Input closed"), std::string::npos);
}

TEST(QctTests, SettingsKeys)
{
    options opt;
    std::vector<WFN> wavy;
    std::string printed;
    EXPECT_EQ(run("L\n4\nL\n0\nL\nabc\nT\n1.5\n0.25\nE\nD\nQ\n", opt, wavy, printed), 0);
    EXPECT_EQ(opt.threads, 4);
    EXPECT_NE(printed.find("Number of threads set to 4"), std::string::npos);
    EXPECT_NE(printed.find("Invalid value, keeping 4"), std::string::npos);
    EXPECT_NE(printed.find("'abc' is not a valid number"), std::string::npos);
    EXPECT_DOUBLE_EQ(opt.properties.radius, 1.5);
    EXPECT_DOUBLE_EQ(opt.properties.resolution, 0.25);
    EXPECT_NE(printed.find("EXPERT MODE!"), std::string::npos);
    EXPECT_TRUE(opt.debug);
    EXPECT_NE(printed.find("expert on  debug on"), std::string::npos);
}

//every wavefunction key refuses politely while nothing is loaded, and bad paths never reach the reader
TEST(QctTests, KeysWithoutWavefunctionAndBadPaths)
{
    Scratch s("BadPaths");
    std::ofstream(s.dir / "notes.txt") << "hello\n";
    options opt;
    std::vector<WFN> wavy;
    std::string printed;
    EXPECT_EQ(run("S\nM\nO\nU\nF\nR\nmissing.wfn\nR\nnotes.txt\nR\n\nC\n1\nQ\n", opt, wavy, printed), 0);
    EXPECT_TRUE(wavy.empty());
    size_t pos = 0, count = 0;
    while ((pos = printed.find("First you need to read a wavefunction!", pos)) != std::string::npos) { count++; pos++; }
    EXPECT_EQ(count, 5u);
    EXPECT_NE(printed.find("No such file: missing.wfn"), std::string::npos);
    EXPECT_NE(printed.find("Unknown extension '.txt'"), std::string::npos);
    EXPECT_NE(printed.find("No cubes loaded"), std::string::npos);
}

//read a .wfn, check its units, convert it to wfx / xyz / occupied wfn, compute a rho cube, read the cube back and integrate it
TEST(QctTests, ReadConvertPropertyCubeAndIntegrate)
{
    Scratch s("ReadConvert");
    const auto wfn = write_h2(s.dir);
    options opt;
    opt.properties.radius = 1.5;
    opt.properties.resolution = 0.3;
    std::vector<WFN> wavy;
    std::string printed;
    const std::string keys =
        "R\nh2.wfn\n"
        "U\n"
        "S\n3\n\n"          //.wfx, default name
        "S\n6\nh2_out.xyz\n"
        "S\n2\nh2_occ.wfn\n"
        "P\nrho\n"
        "R\nh2_rho.cube\n"
        "C\n1\n1\n"         //integrate cube 1
        "C\n99\n"
        "Q\n";
    EXPECT_EQ(run(keys, opt, wavy, printed), 0);
    ASSERT_EQ(wavy.size(), 1u);
    EXPECT_EQ(wavy[0].get_ncen(), 2);
    EXPECT_EQ(wavy[0].get_cube_count(), 1);
    EXPECT_NE(printed.find("Read h2.wfn: 2 atoms, 2 MOs"), std::string::npos);
    EXPECT_NE(printed.find("Active [0] h2.wfn  (wfn)"), std::string::npos);
    EXPECT_NE(printed.find("Appears to be in bohr!"), std::string::npos);
    EXPECT_TRUE(std::filesystem::exists(s.dir / "h2.wfx"));
    EXPECT_TRUE(std::filesystem::exists(s.dir / "h2_out.xyz"));
    EXPECT_TRUE(std::filesystem::exists(s.dir / "h2_occ.wfn"));
    EXPECT_TRUE(std::filesystem::exists(s.dir / "h2_rho.cube"));
    EXPECT_NE(printed.find("Attached h2_rho.cube as cube 0"), std::string::npos);
    const auto at = printed.find("Integrated value: ");
    ASSERT_NE(at, std::string::npos);
    EXPECT_GT(std::stod(printed.substr(at + 18)), 0.0);
    EXPECT_NE(printed.find("cubes 1"), std::string::npos);
    EXPECT_EQ(opt.wfn, std::filesystem::path()); //the property run borrows opt.wfn and gives it back
    //the converted files read back
    WFN wfx(s.dir / "h2.wfx", false);
    EXPECT_EQ(wfx.get_nmo(), 2);
}

//modify marks the wavefunction, quitting then asks; A and X switch and close
TEST(QctTests, ModifyActivateCloseAndUnsavedPrompt)
{
    Scratch s("Modify");
    write_h2(s.dir);
    options opt;
    std::vector<WFN> wavy;
    std::string printed;
    const std::string keys =
        "R\nh2.wfn\nR\nh2.wfn\n"
        "A\n0\n"
        "M\n5\nHe\n0\n0\n5\n2\n"  //add an atom to wavefunction 0
        "M\n3\n99\n"              //no such centre
        "M\n1\n"
        "X\nn\n"                  //modified: refuse to close
        "A\n1\nX\n"               //unmodified copy closes without a question
        "Q\nn\n"                  //unsaved -> stay
        "Q\ny\n";
    EXPECT_EQ(run(keys, opt, wavy, printed), 0);
    ASSERT_EQ(wavy.size(), 1u);
    EXPECT_EQ(wavy[0].get_ncen(), 3);
    EXPECT_TRUE(wavy[0].get_modified());
    EXPECT_NE(printed.find("Added He"), std::string::npos);
    EXPECT_NE(printed.find("No such centre."), std::string::npos);
    EXPECT_NE(printed.find("MODIFIED"), std::string::npos);
    EXPECT_NE(printed.find("Unsaved changes - close anyway (y/n)"), std::string::npos);
    EXPECT_NE(printed.find("Closed."), std::string::npos);
    EXPECT_NE(printed.find("There are unsaved wavefunctions - quit anyway (y/n)"), std::string::npos);
    EXPECT_NE(printed.find("Bye!"), std::string::npos);
}

//the bonding submenu: a bond-plane rho cube around three atoms lands as a cube of the active wavefunction,
//an invalid atom choice is reported and stays in the menu
TEST(QctTests, BondingMenuBondPlaneCube)
{
    Scratch s("Bonding");
    {
        //C at the origin, N two bohr up z, O one bohr off it, one s Gaussian on C (as in BondwiseCoverageTests)
        WFN wavy(e_origin::NOT_YET_DEFINED);
        wavy.push_back_atom("C", 0.0, 0.0, 0.0, 6);
        wavy.push_back_atom("N", 0.0, 0.0, 2.0, 7);
        wavy.push_back_atom("O", 1.0, 0.0, 2.0, 8);
        wavy.push_back_MO(1, 2.0, -0.5);
        std::vector<primitive> prims;
        prims.emplace_back(1, 1, 1.0, 1.0);
        wavy.push_back_spherical_shell(0, 0, vec2{ vec{ 1.0 } }, prims, 0, 1);
        ASSERT_TRUE(wavy.write_wfn(s.dir / "cno.wfn", false, false));
    }
    options opt;
    opt.properties.resolution = 0.5;
    std::vector<WFN> wavy;
    std::string printed;
    const std::string keys =
        "R\ncno.wfn\n"
        "N\n0\n"
        "N\n1\n2\n1\n1\n1\nrho\n"   //atom 1 three times: refused by do_bonds
        "N\n1\n2\n1\n2\n3\nrho\n"   //bond 1-2, plane through 3
        "N\n42\n"
        "Q\n";
    EXPECT_EQ(run(keys, opt, wavy, printed), 0);
    ASSERT_EQ(wavy.size(), 1u);
    EXPECT_NE(printed.find("BONDING ANALYSIS"), std::string::npos);
    EXPECT_NE(printed.find("Bond plane calculation failed"), std::string::npos);
    EXPECT_EQ(wavy[0].get_cube_count(), 1);
    EXPECT_TRUE(std::filesystem::exists(wavy[0].get_cube_path(0)));
    EXPECT_NE(printed.find("attached to the active wavefunction"), std::string::npos);
    EXPECT_NE(printed.find("Sorry, I did not get that."), std::string::npos);
}
