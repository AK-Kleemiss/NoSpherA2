#pragma once
#include "pch.h"
#include "convenience.h"
#include "constants.h"
#include "wfn_class.h"

struct BasisSetMetadata {
    std::string_view name;
    const SimplePrimitive* primitives;       // Pointer to this basis set's primitives
    uint32_t primitiveCount = 0;
    const std::array<int, 118> elementCounts; // Size 118
    const std::array<int, 118> elementOffsets; // Size 118
};

extern const BasisSetMetadata basis_sets[];
extern const std::size_t basis_set_count;

class BasisSet {
public:
    //Constructor loading data from a predefined basis set into the class
    BasisSet(const BasisSetMetadata& data) : _name(data.name), _primitives(const_cast<SimplePrimitive*>(data.primitives)), _primitiveCount(data.primitiveCount), _elementCounts(data.elementCounts), _elementOffsets(data.elementOffsets) {};

    //Empty consturctor for simple initialization
    BasisSet() {
        _elementCounts.fill(0);
        _elementOffsets.fill(0);
    }

    //Constructor for the generation of a basis set from a given WFN
    BasisSet(const WFN& wavy) : BasisSet() {
        gen_auto_aux(wavy);
    }

    std::shared_ptr<std::array<std::vector<primitive>, 118>> get_data();
    const std::span<const SimplePrimitive> operator[](const int& element) const;
    void operator+=(const BasisSet& other);

    void set_name(const std::string& name) { _name = name; }
    void set_primitive_count(uint32_t count) { _primitiveCount = count; }
    uint32_t get_primitive_count() const { return _primitiveCount; }
    void set_element_ranges(const std::array<int, 118>& counts) { _elementCounts = counts; }
    void set_count_for_element(int element, int count) {
        _elementCounts[element] = count;
        _elementOffsets[element] = _primitiveCount;
    }

    void add_owned_primitive(const SimplePrimitive& primitive) {
        _ownedPrimitives.emplace_back(primitive);
        _primitives = _ownedPrimitives.data();
        _primitiveCount++;
    }

    const size_t get_owned_primitive_count() const {
        return _ownedPrimitives.size();
    }

    const SimplePrimitive* get_primitives() const {
        return _primitives;
    }

    void gen_auto_aux(const WFN& orbital_wfn);
    void gen_auto_aux_for_element(const atom& atm);

    bool has_element(const int& element) const {
        return _elementCounts[element - 1] != 0;
    }

    occ::qm::AOBasis to_AOBasis(const std::vector<occ::core::Atom>& atoms) const;
private:
    std::vector<SimplePrimitive> _ownedPrimitives;
    SimplePrimitive* _primitives = nullptr;
    std::string _name;
    uint32_t _primitiveCount = 0;
    std::array<int, 118> _elementCounts; // Size 118
    std::array<int, 118> _elementOffsets; // Size 118
    std::array<std::vector<primitive>, 118> _convertedData;
};


namespace BasisSetLibrary {
    // Access basis sets
    std::shared_ptr<BasisSet> get_basis_set(std::string basis_name);
    bool check_basis_set_exists(std::string basis_name);

    bool read_basis_set_vanilla(const std::filesystem::path& basis_set_path, WFN& wave, const bool& debug);

    bool read_basis_set_missing(const std::filesystem::path& basis_set_path, WFN& wave, bool debug);
}

int load_basis_into_WFN(WFN& wavy, std::shared_ptr<BasisSet> b, bool decontract = true, bool complete = false);
WFN generate_aux_wfn(const WFN& orbital_wfn, std::vector<std::shared_ptr<BasisSet>>& aux_basis);

