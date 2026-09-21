
#include "pch.h"

#include "core/convenience.h"
#include "core/constants.h"
#include "core/fchk.h"
#include "core/AtomGrid.h"
#include "core/SALTED_utilities.h"
#include "core/scattering_factors.h"
#include "core/nos_math.h"
#include "core/GridManager.h"
#include "core/atoms.h"
#include "core/tsc_block.h"
#include "core/cell.h"
#include "core/wfn_class.h"
#include "core/properties.h"
#include "core/integrator.h"
#include "core/integration_params.h"
#include "core/basis_set.h"
#include "core/geometry_aid.h"
#include "core/crystal_energies.h"
#include "core/NoSpherA2.h"
#include "core/isosurface.h"
#include "core/npy.h"
#include "core/libCintMain.h"
#include "core/throughput.h"
#include "core/i_tensor_stream.h"
#include "core/tsc_label_converter.h"
#include "core/sphere_lebedev_rule.h"
#include <occ/qm/hf.h>
#include <occ/qm/scf.h>
#include <spdlog/spdlog.h>
#undef I
#ifdef NOSPHERA2_USE_GPU
#include "core/blas_gpu.h"
#include "core/aux_density_gpu.h"
#endif

static constexpr double PI_VAL = 3.14159265358979323846;


namespace NoSpherA2UnitTests
{
    // -----------------------------------------------------------------------
    // Atom Tests
    // -----------------------------------------------------------------------
    class AtomTest : public ::testing::Test {
    protected:
        static atom make_atom()
        {
            return atom{
                "C1",
                {},
                6,
                1.25,
                -2.5,
                3.75,
                0
            };
        }

        static atom make_atom_fractional(
            const int charge,
            const double x,
            const double y,
            const double z,
            const int group_nr = 0)
        {
            atom result{
                "Test",
                {},
                charge,
                0.0,
                0.0,
                0.0,
                charge
            };

            result.set_frac_coords(d3{ x, y, z });
            result.set_group_nr(group_nr);
            return result;
        }

        static atomID binary_roundtrip(const atomID& original)
        {
            std::stringstream buffer(
                std::ios::in |
                std::ios::out |
                std::ios::binary
            );

            original.write_atom_id(buffer);
            buffer.seekg(0);

            return atomID(buffer);
        }
    };
    TEST_F(AtomTest, ID_WriteAndRead)
    {
        atom value =
            make_atom_fractional(6, 0.1, 0.2, 0.3, 0);

        const atomID originalId = value.get_ID();

        std::stringstream buffer(
            std::ios::in |
            std::ios::out |
            std::ios::binary
        );

        originalId.write_atom_id(buffer);

        EXPECT_EQ(buffer.str().size(), sizeof(atomID));
        EXPECT_EQ(buffer.str().size(), 16);

        buffer.seekg(0);
        const atomID readId(buffer);

        EXPECT_EQ(originalId, readId);
    }

    TEST_F(AtomTest, ID_WriteAndReadPreservesNegativeCoordinates)
    {
        atom value =
            make_atom_fractional(6, -1.25, 2.5, -3.75, -4);

        const atomID originalId = value.get_ID();
        const atomID readId = binary_roundtrip(originalId);

        EXPECT_EQ(originalId, readId);
    }

    TEST_F(AtomTest, ID_IsDeterministic)
    {
        atom first =
            make_atom_fractional(8, -0.125, 1.75, 12.5, 3);

        atom second =
            make_atom_fractional(8, -0.125, 1.75, 12.5, 3);

        EXPECT_EQ(first.get_ID(), second.get_ID());
    }

    TEST_F(AtomTest, ID_DifferentCoordinatesProduceDifferentIDs)
    {
        atom first =
            make_atom_fractional(6, 0.123456, 0.2, 0.3);

        atom second =
            make_atom_fractional(6, 0.123457, 0.2, 0.3);

        /*
         * The difference is 1e-6, which is comfortably larger than the
         * approximately 7.45e-9 resolution of the signed 32-bit encoding
         * over the range [-16, 16].
         */
        EXPECT_NE(first.get_ID(), second.get_ID());
    }

    TEST_F(AtomTest, ID_CoordinateSignAffectsID)
    {
        atom positive =
            make_atom_fractional(6, 1.25, 2.5, 3.75);

        atom negative =
            make_atom_fractional(6, -1.25, 2.5, 3.75);

        EXPECT_NE(positive.get_ID(), negative.get_ID());
    }

    TEST_F(AtomTest, ID_DifferentAtomicNumbersProduceDifferentIDs)
    {
        atom carbon =
            make_atom_fractional(6, 0.1, 0.2, 0.3);

        atom oxygen =
            make_atom_fractional(8, 0.1, 0.2, 0.3);

        EXPECT_NE(carbon.get_ID(), oxygen.get_ID());
    }

    TEST_F(AtomTest, ID_IsAvailableWhenNoGroupWasAssigned)
    {
        /*
         * group_nr feeds the int16_t data field of atomID, which throws on
         * anything that field cannot hold. An atom only gets a group when the
         * CIF reader matches it, so the default has to be usable on its own;
         * every other test here sets one and so never exercises it.
         */
        atom value{ "Test", {}, 6, 0.0, 0.0, 0.0, 6 };
        value.set_frac_coords(d3{ 0.1, 0.2, 0.3 });

        atomID id;
        EXPECT_NO_THROW({ id = value.get_ID(); });
        EXPECT_EQ(id, atomID(0.1, 0.2, 0.3, 0, 6));
    }

    TEST_F(AtomTest, ID_DifferentGroupsProduceDifferentIDs)
    {
        atom first =
            make_atom_fractional(6, 0.1, 0.2, 0.3, -1);

        atom second =
            make_atom_fractional(6, 0.1, 0.2, 0.3, 1);

        EXPECT_NE(first.get_ID(), second.get_ID());
    }

    TEST_F(AtomTest, ID_IsRebuiltWhenCIFPartChanges)
    {
        atom value = make_atom_fractional(6, 0.1, 0.2, 0.3, 1);
        const atomID part_one_id = value.get_ID();

        // CIF matching can assign PART after an ID has already been requested.
        value.set_group_nr(2);

        // set_group_nr() updates the cached value, rather than leaving it empty
        // for a later get_ID() call to reconstruct.
        EXPECT_EQ(value.get_ID(), atomID(0.1, 0.2, 0.3, 2, 6));
        EXPECT_NE(value.get_ID(), part_one_id);
    }

    TEST_F(AtomTest, ID_SupportsCoordinateRangeBoundaries)
    {
        EXPECT_NO_THROW({
            const atomID minimum(-16.0, -16.0, -16.0, 0, 6);
            const atomID restored = binary_roundtrip(minimum);
            EXPECT_EQ(minimum, restored);
            });

        EXPECT_NO_THROW({
            const atomID maximum(16.0, 16.0, 16.0, 0, 6);
            const atomID restored = binary_roundtrip(maximum);
            EXPECT_EQ(maximum, restored);
            });
    }

    TEST_F(AtomTest, ID_RejectsCoordinatesOutsideSupportedRange)
    {
        EXPECT_THROW(
            (atomID{ 16.000001, 0.0, 0.0, 0, 6 }),
            std::out_of_range
        );

        EXPECT_THROW(
            (atomID{ -16.000001, 0.0, 0.0, 0, 6 }),
            std::out_of_range
        );
    }

    TEST_F(AtomTest, ID_RejectsNonFiniteCoordinates)
    {
        const double infinity =
            std::numeric_limits<double>::infinity();

        const double nan =
            std::numeric_limits<double>::quiet_NaN();

        EXPECT_THROW(
            (atomID{ infinity, 0.0, 0.0, 0, 6 }),
            std::invalid_argument
        );

        EXPECT_THROW(
            (atomID{ nan, 0.0, 0.0, 0, 6 }),
            std::invalid_argument
        );
    }

    TEST_F(AtomTest, ID_RejectsInvalidAtomicNumber)
    {
        EXPECT_THROW(
            (atomID{ 0.1, 0.2, 0.3, 0, 0 }),
            std::out_of_range
        );

        EXPECT_THROW(
            (atomID{ 0.1, 0.2, 0.3, 0, 256 }),
            std::out_of_range
        );
    }

    TEST_F(AtomTest, ID_RejectsDataOutsideInt16Range)
    {
        EXPECT_THROW(
            (atomID{
                0.1,
                0.2,
                0.3,
                static_cast<int>(
                    std::numeric_limits<std::int16_t>::max()
                ) + 1,
                6
                }),
            std::out_of_range
        );

        EXPECT_THROW(
            (atomID{
                0.1,
                0.2,
                0.3,
                static_cast<int>(
                    std::numeric_limits<std::int16_t>::min()
                ) - 1,
                6
                }),
            std::out_of_range
        );
    }

    TEST_F(AtomTest, ID_DefaultConstructedObjectIsNotInitialized)
    {
        const atomID id;

        EXPECT_FALSE(id.is_initialized());
    }

    TEST_F(AtomTest, ID_CannotWriteUninitializedObject)
    {
        const atomID id;

        std::stringstream buffer(
            std::ios::in |
            std::ios::out |
            std::ios::binary
        );

        EXPECT_THROW(
            id.write_atom_id(buffer),
            std::runtime_error
        );
    }

    TEST_F(AtomTest, ID_RejectsTruncatedBinaryInput)
    {
        /*
         * A valid atomID requires 16 bytes, but this stream contains only 8.
         */
        const std::string incompleteData(8, '\0');

        std::istringstream input(
            incompleteData,
            std::ios::in | std::ios::binary
        );

        EXPECT_THROW(
            (atomID{ input }),
            std::runtime_error
        );
    }
    // -----------------------------------------------------------------------
    // Non-trivial atom behavior
    // -----------------------------------------------------------------------

    TEST_F(AtomTest, DistanceToOtherAtomIsEuclideanDistance)
    {
        const atom first{
            "A",
            {},
            1,
            1.0,
            2.0,
            3.0,
            0
        };

        const atom second{
            "B",
            {},
            1,
            4.0,
            6.0,
            3.0,
            0
        };

        EXPECT_NEAR(first.distance_to(second), 5.0, 1e-12);
        EXPECT_NEAR(second.distance_to(first), 5.0, 1e-12);
    }

    TEST_F(AtomTest, BasisSetSupportsAddingModifyingAndErasingEntries)
    {
        atom value = make_atom();

        ASSERT_TRUE(value.push_back_basis_set(10.0, 0.1, 1, 0));
        ASSERT_TRUE(value.push_back_basis_set(20.0, 0.2, 2, 1));
        ASSERT_TRUE(value.push_back_basis_set(30.0, 0.3, 3, 2));

        value.set_basis_set_exponent(1, 25.0);
        value.set_basis_set_coefficient(1, 0.25);

        EXPECT_DOUBLE_EQ(value.get_basis_set_exponent(1), 25.0);
        EXPECT_DOUBLE_EQ(value.get_basis_set_coefficient(1), 0.25);

        value.erase_basis_set(0);

        ASSERT_EQ(value.get_basis_set_size(), 2u);
        EXPECT_DOUBLE_EQ(value.get_basis_set_exponent(0), 25.0);
        EXPECT_DOUBLE_EQ(value.get_basis_set_exponent(1), 30.0);
    }

    TEST_F(AtomTest, IndexedShellCountSetterExpandsAndZeroInitializesVector)
    {
        atom value = make_atom();

        value.set_shellcount(3u, 9u);

        ASSERT_EQ(value.get_shellcount_size(), 4u);
        EXPECT_EQ(value.get_shellcount(0u), 0u);
        EXPECT_EQ(value.get_shellcount(1u), 0u);
        EXPECT_EQ(value.get_shellcount(2u), 0u);
        EXPECT_EQ(value.get_shellcount(3u), 9u);
    }

    TEST_F(AtomTest, AssignmentPerformsDeepCopy)
    {
        atom source{
            "O1",
            {},
            8,
            1.0,
            2.0,
            3.0,
            -2,
            2
        };

        source.set_frac_coords(d3{ 0.1, 0.2, 0.3 });
        source.set_shellcount(std::vector<unsigned int>{2u, 3u});
        ASSERT_TRUE(source.push_back_basis_set(25.0, 0.75, 2, 1));

        atom destination;
        destination = source;

        destination.set_label("Changed");
        destination.set_basis_set_exponent(0, 999.0);
        destination.set_shellcount(0u, 99u);

        EXPECT_EQ(source.get_label(), "O1");
        EXPECT_DOUBLE_EQ(source.get_basis_set_exponent(0), 25.0);
        EXPECT_EQ(source.get_shellcount(0u), 2u);

        EXPECT_EQ(destination.get_label(), "Changed");
        EXPECT_DOUBLE_EQ(destination.get_basis_set_exponent(0), 999.0);
        EXPECT_EQ(destination.get_shellcount(0u), 99u);
    }

    TEST_F(AtomTest, EqualityDetectsMeaningfulDifference)
    {
        atom first = make_atom();
        atom second = make_atom();

        EXPECT_TRUE(first == second);

        ASSERT_TRUE(second.push_back_basis_set(10.0, 0.5, 1, 0));

        EXPECT_FALSE(first == second);
    }

    // ParseSymopTests — cell::parse_symop, the CIF symmetry operation reader
    struct parsed_symop
    {
        int rot[3][3]{};
        double trans[3]{};
    };

    static parsed_symop parse(const std::string& operation)
    {
        parsed_symop result;
        std::ostringstream sink;
        cell::parse_symop(operation, "test.cif", result.rot, result.trans, sink);
        return result;
    }

    static void expect_rot(const parsed_symop& op, const int expected[3][3])
    {
        for (int comp = 0; comp < 3; comp++)
            for (int axis = 0; axis < 3; axis++)
                EXPECT_EQ(op.rot[comp][axis], expected[comp][axis])
                    << "component " << comp << ", axis " << axis;
    }

    TEST(ParseSymopTest, Identity)
    {
        const parsed_symop op = parse("x,y,z");
        const int expected[3][3] = { {1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
        expect_rot(op, expected);
        EXPECT_DOUBLE_EQ(op.trans[0], 0.0);
        EXPECT_DOUBLE_EQ(op.trans[1], 0.0);
        EXPECT_DOUBLE_EQ(op.trans[2], 0.0);
    }

    // The operations of P 2_1 2_1 2_1, which used to abort the process because the
    // translation was read with stof("x+1") instead of being split off the axis term.
    TEST(ParseSymopTest, ScrewAxisWithTrailingFraction)
    {
        const parsed_symop op = parse("x+1/2,-y+1/2,-z");
        const int expected[3][3] = { {1, 0, 0}, {0, -1, 0}, {0, 0, -1} };
        expect_rot(op, expected);
        EXPECT_DOUBLE_EQ(op.trans[0], 0.5);
        EXPECT_DOUBLE_EQ(op.trans[1], 0.5);
        EXPECT_DOUBLE_EQ(op.trans[2], 0.0);
    }

    TEST(ParseSymopTest, TranslationBeforeAxis)
    {
        const parsed_symop op = parse("1/2+X,1/2-Y,-Z");
        const int expected[3][3] = { {1, 0, 0}, {0, -1, 0}, {0, 0, -1} };
        expect_rot(op, expected);
        EXPECT_DOUBLE_EQ(op.trans[0], 0.5);
        EXPECT_DOUBLE_EQ(op.trans[1], 0.5);
        EXPECT_DOUBLE_EQ(op.trans[2], 0.0);
    }

    TEST(ParseSymopTest, DecimalTranslationsAndWhitespace)
    {
        const parsed_symop op = parse(" 0.5 - x , y , 0.25 + z ");
        const int expected[3][3] = { {-1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
        expect_rot(op, expected);
        EXPECT_DOUBLE_EQ(op.trans[0], 0.5);
        EXPECT_DOUBLE_EQ(op.trans[1], 0.0);
        EXPECT_DOUBLE_EQ(op.trans[2], 0.25);
    }

    // Rhombohedral obverse setting: mixed axes in one component and thirds
    TEST(ParseSymopTest, MixedAxesAndThirds)
    {
        const parsed_symop op = parse("-y+2/3,x-y+1/3,z+1/3");
        const int expected[3][3] = { {0, -1, 0}, {1, -1, 0}, {0, 0, 1} };
        expect_rot(op, expected);
        EXPECT_NEAR(op.trans[0], 2.0 / 3.0, 1e-12);
        EXPECT_NEAR(op.trans[1], 1.0 / 3.0, 1e-12);
        EXPECT_NEAR(op.trans[2], 1.0 / 3.0, 1e-12);
    }

    TEST(ParseSymopTest, Inversion)
    {
        const parsed_symop op = parse("-x,-y,-z");
        const int expected[3][3] = { {-1, 0, 0}, {0, -1, 0}, {0, 0, -1} };
        expect_rot(op, expected);
        EXPECT_DOUBLE_EQ(op.trans[0], 0.0);
        EXPECT_DOUBLE_EQ(op.trans[1], 0.0);
        EXPECT_DOUBLE_EQ(op.trans[2], 0.0);
    }


    // A malformed operation has to leave through error_check's exit(-1), not abort the
    // process with a fail-fast. The message itself cannot be matched here because
    // error_check reports on stdout while death tests only see stderr.
    TEST(ParseSymopDeathTest, MalformedOperationExitsCleanly)
    {
        parsed_symop result;
        EXPECT_EXIT(cell::parse_symop("x,y", "test.cif", result.rot, result.trans, std::cout),
            ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
        EXPECT_EXIT(cell::parse_symop("x+1/0,y,z", "test.cif", result.rot, result.trans, std::cout),
            ::testing::ExitedWithCode(ERROR_CHECK_EXIT_CODE), ".*");
    }

    // ------------------------------------------------------------------
    // atom: basis-set entries, ADP orders, IDs and coordinates

    TEST(AtomStateTests, BasisSetAdpsAndCoordinates)
    {
        atom a("C1", atomID(), 1, 0.1, 0.2, 0.3, 6);
        EXPECT_FALSE(a.get_basis_set_loaded());
        EXPECT_TRUE(a.push_back_basis_set(1.5, 0.7, 1, 0));
        EXPECT_TRUE(a.push_back_basis_set(0.5, 0.3, 1, 0));
        EXPECT_TRUE(a.push_back_basis_set(0.2, 1.0, 2, 1));
        EXPECT_FALSE(a.push_back_basis_set(0.2, 1.0, -1, 2)) << "negative type";
        EXPECT_TRUE(a.get_basis_set_loaded());
        EXPECT_EQ(a.get_basis_set_size(), 3u);
        EXPECT_EQ(a.get_shellcount(0), 2u);
        EXPECT_EQ(a.get_shellcount(1), 1u);
        EXPECT_EQ(a.get_basis_set_exponent(1), 0.5);
        EXPECT_EQ(a.get_basis_set_type(2), 2);
        a.print_values_long();
        EXPECT_FALSE(a.is_anharm());
        double uiso = 0.05;
        a.assign_ADPs(uiso);
        EXPECT_FALSE(a.is_anharm());
        vec bad(5, 0.0), second(6, 0.01), third(10, 0.001), fourth(15, 0.0001), none;
        a.assign_ADPs(bad);
        a.assign_ADPs(second);
        EXPECT_FALSE(a.is_anharm());
        a.assign_ADPs(second, none, none);
        EXPECT_FALSE(a.is_anharm());
        a.assign_ADPs(bad, third, fourth);
        a.assign_ADPs(second, bad, fourth);
        a.assign_ADPs(second, third, bad);
        EXPECT_FALSE(a.is_anharm());
        a.assign_ADPs(second, third, fourth);
        EXPECT_TRUE(a.is_anharm());
        // the ID is derived from the fractional position, group and charge on first use
        a.set_frac_coords({ 0.25, 0.5, 0.75 });
        EXPECT_EQ(a.get_frac_coordinate(1), 0.5);
        EXPECT_EQ(a.get_frac_coordinate(3), 0.0);
        const atomID id = a.get_ID();
        EXPECT_TRUE(id.is_initialized());
        EXPECT_EQ(id, atomID(0.25, 0.5, 0.75, 0, 6));
        EXPECT_EQ(a.get_ID(), id);
        a.set_ID(atomID(0.1, 0.1, 0.1, 0, 6));
        EXPECT_FALSE(a.get_ID() == id);
        a.set_coordinate(2, 9.0);
        a.set_coordinate(5, 1.0);
        EXPECT_EQ(a.get_coordinate(2), 9.0);
        EXPECT_EQ(a.get_coordinate(5), 0.0);
        EXPECT_EQ(a.get_pos()[0], 0.1);
        const atom copy(a);
        EXPECT_EQ(copy.get_basis_set_size(), 3u);
        EXPECT_TRUE(copy.is_anharm());
        atom assigned;
        assigned = a;
        EXPECT_EQ(assigned.get_label(), "C1");
        EXPECT_EQ(basis_set_entry().get_exponent(), 0.0);
    }

} // namespace NoSpherA2UnitTests
