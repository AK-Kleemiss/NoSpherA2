#pragma once
#include "convenience.h"

//-----------------Definition of atoms and basis sets--------------------

class basis_set_entry {
private:
    double coefficient;
    double exponent;
    unsigned int type; //1=S; 2=P; 3=D; 4=F; 5=G
    unsigned int shell;
    primitive p;
public:
    // Default assignment operator
    basis_set_entry& operator=(const basis_set_entry& rhs) = default;
    basis_set_entry();
    basis_set_entry(double g_coefficient, double g_exponent, unsigned int g_type, unsigned int g_shell);
    // Getter functions

    const double& get_coefficient() const { return coefficient; };
    const double& get_exponent() const { return exponent; };
    const unsigned int& get_type() const { return type; };
    const unsigned int& get_shell() const { return shell; };
    const primitive& get_primitive() const { return p; };
    // Setter functions

    void set_coefficient(const double& value) { coefficient = value; };
    void set_exponent(const double& value) { exponent = value; };
    void set_type(const unsigned int& value) { type = value; };
    void set_shell(const unsigned int& value) { shell = value; };
    void set_primitive(const primitive& value) { p = value; };

    // Equality operator
    bool operator==(const basis_set_entry& other) const {
        return coefficient == other.get_coefficient() &&
            exponent == other.get_exponent() &&
            type == other.get_type() &&
            shell == other.get_shell() &&
            p == other.get_primitive();
    }
};

namespace atomID_utils {

    inline constexpr double max_fractional_coordinate = 16.0;

    inline constexpr double integer_max =
        static_cast<double>(std::numeric_limits<std::int32_t>::max());

    inline constexpr double encode_scale =
        integer_max / max_fractional_coordinate;

    inline constexpr double decode_scale =
        max_fractional_coordinate / integer_max;

} // namespace atomID_utils
class atomID {
private:
    int32_t frac_x_int;
    int32_t frac_y_int;
    int32_t frac_z_int;
    int16_t data_;
    uint8_t Z_;
    uint8_t reserved_;

    static std::int32_t encode_coordinate(double coordinate)
    {
        if (!std::isfinite(coordinate)) {
            throw std::invalid_argument("Coordinate must be finite");
        }

        if (coordinate < -atomID_utils::max_fractional_coordinate ||
            coordinate > atomID_utils::max_fractional_coordinate) {
            throw std::out_of_range(
                "Coordinate must be between -16 and 16"
            );
        }

        return static_cast<std::int32_t>(
            std::llround(coordinate * atomID_utils::encode_scale)
         );
    }

    static double decode_coordinate(std::int32_t encoded)
    {
        return static_cast<double>(encoded) * atomID_utils::decode_scale;
    }

    void validate_loaded_data() const
    {
        if (Z_ == 0) {
            throw std::runtime_error(
                "Invalid atom ID: atomic number is zero");
        }

        /*
         * The encoder maps -16 to -INT32_MAX, not INT32_MIN.
         * Therefore, INT32_MIN cannot be produced by a valid encoder.
         */
        if (frac_x_int == std::numeric_limits<std::int32_t>::min() ||
            frac_y_int == std::numeric_limits<std::int32_t>::min() ||
            frac_z_int == std::numeric_limits<std::int32_t>::min()) {
            throw std::runtime_error(
                "Invalid atom ID: encoded coordinate is out of range");
        }

        /*
         * If reserved_ is currently required to be zero, enable this:
         *
         * if (reserved_ != 0) {
         *     throw std::runtime_error(
         *         "Unsupported atom ID format version");
         * }
         */
    }
public:
    // Default constructor initializes to an uninitialized state
    atomID() : frac_x_int(0), frac_y_int(0), frac_z_int(0), data_(0), Z_(0), reserved_(0) {}

    //Construct from fractional coordinates, data, and atomic number
    atomID(const double frac_x, const double frac_y, const double frac_z, const int data, const int Z, const int reserved = 0)
    {
        if (Z < 1 || Z > 255) {
            throw std::out_of_range("Z must be between 1 and 255");
        }
        if (data < std::numeric_limits<std::int16_t>::min() ||
            data > std::numeric_limits<std::int16_t>::max()) {
            throw std::out_of_range("data does not fit into int16_t");
        }
        if (reserved < 0 || reserved > 255) {
            throw std::out_of_range(
                "reserved must be between 0 and 255");
        }

        frac_x_int = encode_coordinate(frac_x);
        frac_y_int = encode_coordinate(frac_y);
        frac_z_int = encode_coordinate(frac_z);
        data_ = static_cast<int16_t>(data);
        Z_ = static_cast<uint8_t>(Z);
        reserved_ = static_cast<uint8_t>(reserved);
    }

    //Construct from two 64-bit integers (binary representation)
    atomID(const uint64_t first, const uint64_t second) {
        *this = from_uint64(first, second);
        validate_loaded_data();
    }

    //Construct from a hexadecimal string representation
    atomID(const std::string_view hex) {
        if (hex.size() != 32) {
            throw std::invalid_argument(
                "atomID hexadecimal string must contain exactly 32 characters");
        }

        std::uint64_t first{};
        std::uint64_t second{};

        const auto first_part = hex.substr(0, 16);
        const auto second_part = hex.substr(16, 16);

        const auto first_result = std::from_chars(
            first_part.data(),
            first_part.data() + first_part.size(),
            first,
            16);

        const auto second_result = std::from_chars(
            second_part.data(),
            second_part.data() + second_part.size(),
            second,
            16);

        if (first_result.ec != std::errc{} ||
            first_result.ptr != first_part.data() + first_part.size() ||
            second_result.ec != std::errc{} ||
            second_result.ptr != second_part.data() + second_part.size()) {
            throw std::invalid_argument(
                "atomID string contains invalid hexadecimal characters");
        }

        *this = from_uint64(first, second);
    }

    //Construct from a binary input stream
    atomID(std::istream& input) {
        input.read(reinterpret_cast<char*>(this), sizeof(atomID));
        if (!input) {
            throw std::runtime_error("Failed to read atom ID from stream");
        }
        validate_loaded_data();
    }

    std::array<std::uint64_t, 2> as_uint64() const noexcept
    {
        // First 64-bit integer:
        // bits  0-31: frac_x_int
        // bits 32-63: frac_y_int
        const std::uint64_t first =
            static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(frac_x_int))
            |
            (static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(frac_y_int)) << 32);

        // Second 64-bit integer:
        // bits  0-31: frac_z_int
        // bits 32-47: data_
        // bits 48-55: Z_
        // bits 56-63: reserved_
        const std::uint64_t second =
            static_cast<std::uint64_t>(
                static_cast<std::uint32_t>(frac_z_int))
            |
            (static_cast<std::uint64_t>(
                static_cast<std::uint16_t>(data_)) << 32)
            |
            (static_cast<std::uint64_t>(Z_) << 48)
            |
            (static_cast<std::uint64_t>(reserved_) << 56);

        return { first, second };
    }

    static atomID from_uint64(
        const std::uint64_t first,
        const std::uint64_t second) noexcept
    {
        atomID result;

        const auto frac_x_bits =
            static_cast<std::uint32_t>(first);

        const auto frac_y_bits =
            static_cast<std::uint32_t>(first >> 32);

        const auto frac_z_bits =
            static_cast<std::uint32_t>(second);

        const auto data_bits =
            static_cast<std::uint16_t>(second >> 32);

        result.frac_x_int = std::bit_cast<std::int32_t>(frac_x_bits);
        result.frac_y_int = std::bit_cast<std::int32_t>(frac_y_bits);
        result.frac_z_int = std::bit_cast<std::int32_t>(frac_z_bits);
        result.data_ = std::bit_cast<std::int16_t>(data_bits);

        result.Z_ =
            static_cast<std::uint8_t>(second >> 48);

        result.reserved_ =
            static_cast<std::uint8_t>(second >> 56);

        return result;
    }


    [[nodiscard]]
    std::string to_hex_string() const
    {
        const auto [first, second] = as_uint64();

        std::ostringstream stream;
        stream << std::hex
            << std::setfill('0')
            << std::setw(16) << first
            << std::setw(16) << second;

        return stream.str();
    }

    // Equality operator
    bool operator==(const atomID& other) const {
        return frac_x_int == other.frac_x_int &&
            frac_y_int == other.frac_y_int &&
            frac_z_int == other.frac_z_int &&
            data_ == other.data_ &&
            Z_ == other.Z_ &&
            reserved_ == other.reserved_;
    }

    //<< operator for printing atomID
    friend std::ostream& operator<<(std::ostream& os, const atomID& id) {
        os << "atomID(";
        os << "frac_x: " << id.frac_x() << ", ";
        os << "frac_y: " << id.frac_y() << ", ";
        os << "frac_z: " << id.frac_z() << ", ";
        os << "Z: " << static_cast<int>(id.Z()) << ", ";
        os << "data: " << id.data() << ", ";
        os << "reserved: " << static_cast<int>(id.reserved());
        os << ")";
        return os;
    }

    const bool is_initialized() const {
        return Z_ != 0;
    }

    void write_atom_id(std::ostream& os) const {
        if (!is_initialized()) {
            throw std::runtime_error("Atom ID is not initialized");
        }

        os.write(
            reinterpret_cast<const char*>(this),
            sizeof(atomID)
        );

        if (!os) {
            throw std::runtime_error("Failed to write atom ID to stream");
        }

    }

    double frac_x() const noexcept
    {
        return decode_coordinate(frac_x_int);
    }

    double frac_y() const noexcept
    {
        return decode_coordinate(frac_y_int);
    }

    double frac_z() const noexcept
    {
        return decode_coordinate(frac_z_int);
    }

    std::int16_t data() const noexcept
    {
        return data_;
    }

    std::uint8_t Z() const noexcept
    {
        return Z_;
    }

    std::uint8_t reserved() const noexcept
    {
        return reserved_;
    }
};

//We have to define a custom hash function for atomID because it is a user-defined type and the standard library does not provide a hash function for it.
struct AtomIDHash {
    [[nodiscard]]
    std::size_t operator()(const atomID& id) const noexcept
    {
        const auto [first, second] = id.as_uint64();

        const std::size_t h1 = std::hash<std::uint64_t>{}(first);
        const std::size_t h2 = std::hash<std::uint64_t>{}(second);

        return h1 ^ (
            h2
            + static_cast<std::size_t>(0x9e3779b97f4a7c15ULL)
            + (h1 << 6)
            + (h1 >> 2));
    }
};
namespace std {

    template<>
    struct hash<atomID> {
        [[nodiscard]]
        std::size_t operator()(const atomID& id) const noexcept
        {
            return AtomIDHash{}(id);
        }
    };

}

//Test that atomID is exactly 16 bytes, trivially copyable, and standard layout.
//To make sure every compiler will have the same binary representation of atomID.
static_assert(
    sizeof(atomID) == 16,
    "atomID binary representation must be exactly 16 bytes");

static_assert(
    std::is_trivially_copyable_v<atomID>,
    "atomID must remain trivially copyable for raw binary I/O");

static_assert(
    std::is_standard_layout_v<atomID>,
    "atomID must remain a standard-layout class");


class atom {
private:
    std::string label;
    atomID ID;
    int nr, charge, ECP_electrons;
    double x, y, z;
    d3 frac_coords;
    std::vector<basis_set_entry> basis_set;
    int basis_set_id;
    std::vector<unsigned int> shellcount;
    //The Order is:
    //[0] = second order (U11, U22, U33, U12, U13, U23)
    //[1] = third order  (C111, C112, C113, C122, C123, C133, C222, C223, C233, C333)
    //[2] = fourth order (D1111, D1112, D1113, D1122, D1123, D1133, D1222, D1223, D1233, D1333, D2222, D2223, D2233, D2333, D3333)
    vec2 ADPs;
    bool is_asym = false;
    /*
     * Feeds the data field of atomID, which is an int16_t and rejects anything
     * outside its range, so the default has to be a value that field can hold.
     * Zero is what the rest of the code already means by "no group": it is what
     * the CIF reader uses when the file carries no group column.
     */
    int group_nr = 0;
public:
    atom();
    atom(const std::string& l, const atomID& id, const int& n, const double& c1, const double& c2, const double& c3, const int& ch);
    atom(const std::string& l, const atomID& id, const int& n, const double& c1, const double& c2, const double& c3, const int& ch, const int& ECP_els);
    atom& operator=(const atom& rhs);
    void print_values() const;
    bool push_back_basis_set(const double & exponent, const double &coefficient, const int &type, const int &shell);
    void print_values_long() const;
    bool get_basis_set_loaded() const;
    bool is_anharm() const;
    void assign_ADPs(vec& second, vec& third, vec& fourth);
    void assign_ADPs(vec& second);
    void assign_ADPs(double& Uiso);
    void set_ID(const atomID& id);
    atomID get_ID();

    const basis_set_entry& get_basis_set_entry(const int& _nr) const { return basis_set[_nr]; };
    // Non-const overload to obtain a mutable reference
    basis_set_entry& get_basis_set_entry(const int& _nr) { return basis_set[_nr]; }
    double get_coordinate(const unsigned int& axis) const;
    double get_frac_coordinate(const unsigned int& axis) const;
    int get_charge() const { return charge; };
    int get_ECP_electrons() const { return ECP_electrons; };
    void set_charge(const int& ch) { charge = ch; };
    void set_ECP_electrons(const int& ECP_els) { ECP_electrons = ECP_els; };
    void set_coordinate(const unsigned int& axis, const double& value);
    const d3 get_pos() const { return { x,y,z }; };
    void set_frac_coords(const d3& frac);
    std::string get_label() const { return label; };
    int get_nr() const { return nr; };
    void set_label(const std::string& l) { label = l; };
    void set_nr(const int& n) { nr = n; };
    void clear_basis_set() { basis_set.clear(); };
    void set_basis_set(std::vector<basis_set_entry> new_basis) { basis_set = new_basis; };
    unsigned int get_basis_set_size() const { return (unsigned int)basis_set.size(); };
    double get_basis_set_exponent(const int& _nr) const { return basis_set[_nr].get_exponent(); };
    double get_basis_set_coefficient(const int& _nr) const { return basis_set[_nr].get_coefficient(); };
    void set_basis_set_exponent(const int& _nr, const double& value) { basis_set[_nr].set_exponent(value); };
    void set_basis_set_coefficient(const int& _nr, const double& value) { basis_set[_nr].set_coefficient(value); };
    void set_basis_set_type(const int& _nr, const int& value) { basis_set[_nr].set_type(value); };
    void set_basis_set_shell(const int& _nr, const int& value) { basis_set[_nr].set_shell(value); };
    std::vector<basis_set_entry> get_basis_set() const { return basis_set; };
    const basis_set_entry* get_basis_set_ptr() const { return &(basis_set[0]); };
    void erase_basis_set(const unsigned int& _nr) { basis_set.erase(basis_set.begin() + _nr); };
    int get_basis_set_type(const int& _nr) const { return basis_set[_nr].get_type(); };
    int get_basis_set_shell(const int& _nr) const { return basis_set[_nr].get_shell(); };
    void set_basis_set_id(const int& id) { basis_set_id = id; };
    int get_basis_set_id() const { return basis_set_id; };
    void set_shellcount(const std::vector<unsigned int> sc) { shellcount.resize(sc.size()); shellcount = sc; };
    std::vector<unsigned int> get_shellcount() const { return shellcount; };
    unsigned int get_shellcount_size() const { return (unsigned int)shellcount.size(); };
    // Make indexed setter safe
    void set_shellcount(const unsigned int& _nr, const unsigned int& value) {
        if (_nr >= shellcount.size()) shellcount.resize(_nr + 1u, 0u);
        shellcount[_nr] = value;
    }
    unsigned int get_shellcount(const unsigned int& _nr) const { return shellcount[_nr]; };
    void push_back_shell(const unsigned int& type) { shellcount.push_back(type); };
    void clear_shellcount() { shellcount.clear(); };
    void set_ADPs(const vec2& adps) { ADPs = adps; };
    vec2 get_ADPs() const { return ADPs; };
    bool get_is_asym() const { return is_asym; };
    void set_is_asym(const bool& val) { is_asym = val; };
    // Fractional coordinates and disorder group form the CIF-derived
    // SCATTERER_ID. Rebuild it when either value changes so a subsequent MTC
    // comparison cannot reuse an ID from another PART.
    void set_group_nr(const int group_nr) {
        this->group_nr = group_nr;
        if (charge >= 1 && charge <= 255)
            ID = atomID(frac_coords[0], frac_coords[1], frac_coords[2], this->group_nr, charge);
    };
    double distance_to(const atom& other) const;

    bool operator==(const atom& other) const;
};

atomID get_atom_ID(const int Z, const d3& frac_coords, const int dat = 0);

struct scatterer_id_masks_d5 {
    const static uint64_t
        z_mask = 0x00000000000000FF,
        a_mask = 0x0000000000FFFF00,
        b_mask = 0x000000FFFF000000,
        c_mask = 0x00FFFF0000000000,
        a_sig = 0x0100000000000000,
        b_sig = 0x0200000000000000,
        c_sig = 0x0400000000000000,
        d_mask = 0xF800000000000000,
        mask_m = 0x000000000000FFFF; // max crd value
    const static int a_shift = 16; 
    //The groups in refinement are defined positive and negative, so we need to shift the group number by 16 to make it positive and fit.
    // E.g. group 0 becomes 16, group 1 becomes 17, group -1 becomes 15, etc.
    const static int group_shift = 16;
};
