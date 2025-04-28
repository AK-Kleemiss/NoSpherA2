#pragma once
#include <array>
#include <vector>
#include "pch.h"
#include "convenience.h"
#include "constants.h"
#include "wfn_class.h"


struct ElementRange {
    uint32_t start = 0;  // Start index in this basis set’s primitive array
    uint32_t count = 0;  // Number of primitives for the element
};

struct BasisSetMetadata {
    std::string_view name;
    const SimplePrimitive* primitives;       // Pointer to this basis set's primitives
    uint32_t primitiveCount = 0;
    const std::array<ElementRange,118> elementRanges; // Size 118
};

extern const BasisSetMetadata aux_basis_sets[];
extern const std::size_t aux_basis_set_count;

class BasisSet {
public:
    BasisSet(const BasisSetMetadata& data) : _name(data.name), _primitives(const_cast<SimplePrimitive*>(data.primitives)), _primitiveCount(data.primitiveCount), _elementRanges(data.elementRanges) {};
    BasisSet() = default;

    const std::array<std::vector<primitive>, 118>& get_data();
    const std::span<const SimplePrimitive> operator[](int element) const;

    void set_name(const std::string& name) { _name = name; }
    void set_primitive_count(uint32_t count) { _primitiveCount = count; }
    void set_element_ranges(const std::array<ElementRange, 118>& ranges) { _elementRanges = ranges; }
    void set_element_range_for_element(int element, uint32_t start, uint32_t count) {
        _elementRanges[element].start = start;
        _elementRanges[element].count = count;
    }
    void add_owned_primitive(const SimplePrimitive& primitive) {
        _ownedPrimitives.push_back(primitive);
        _primitives = _ownedPrimitives.data();
    }

    size_t get_owned_primitive_count() const {
        return _ownedPrimitives.size();
    }

    SimplePrimitive* get_primitives() {
        return _primitives;
    }
    
    
private:
    std::vector<SimplePrimitive> _ownedPrimitives;
    SimplePrimitive* _primitives = nullptr;
    std::string _name;
    uint32_t _primitiveCount = 0;
    std::array<ElementRange, 118> _elementRanges;
    std::array<std::vector<primitive>, 118> _converted_data;
};


class BasisSetLibrary {
public:
    // Constructor
    BasisSetLibrary() = default;
    // Access basis set
    BasisSet get_basis_set(std::string basis_name);
    bool check_basis_set_exists(std::string basis_name);
    BasisSet gen_aux(const WFN& orbital_wfn, double& n);
private:

    BasisSet auxiliary_basis;
};
