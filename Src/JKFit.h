#pragma once
#include <array>
#include <vector>
#include "pch.h"
#include "convenience.h"
#include "constants.h"
#include "wfn_class.h"


//struct ElementRange {
//    uint32_t start = 0;  // Start index in this basis sets primitive array
//    uint32_t count = 0;  // Number of primitives for the element
//};

struct BasisSetMetadata {
    std::string_view name;
    const SimplePrimitive* primitives;       // Pointer to this basis set's primitives
    uint32_t primitiveCount = 0;
    const std::array<int, 118> elementCounts; // Size 118
    const std::array<int, 118> elementOffsets; // Size 118
};

extern const BasisSetMetadata aux_basis_sets[];
extern const std::size_t aux_basis_set_count;

class BasisSet {
public:
    //Constructor loading data from a predefined basis set into the class
    BasisSet(const BasisSetMetadata& data) : _name(data.name), _primitives(const_cast<SimplePrimitive*>(data.primitives)), _primitiveCount(data.primitiveCount), _elementCounts(data.elementCounts), _elementOffsets(data.elementOffsets) {};
    
    //Empty consturctor for simple initialization
    BasisSet() {
        _elementCounts.fill(0);
        _elementOffsets.fill(0);
    }

    //Constructor for the generation of a basis set from a given WFN and beta value
    BasisSet(const WFN& wavy, const double beta) : BasisSet() {
        _beta = beta;
        gen_aux(wavy);
    }

    //Constructor for the initialization of a auto_aux basis set
    //By calling the gen_aux function on this object, the basis set will be generated
    BasisSet(const double beta) : BasisSet()  {
        _beta = beta;
    }


    std::shared_ptr<std::array<std::vector<primitive>, 118>> get_data();
    const std::span<const SimplePrimitive> operator[](const int&element) const ;

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
    
    void gen_aux(const WFN& orbital_wfn);
    
private:
    std::vector<SimplePrimitive> _ownedPrimitives;
    SimplePrimitive* _primitives = nullptr;
    std::string _name;
    uint32_t _primitiveCount = 0;
    std::array<int, 118> _elementCounts; // Size 118
    std::array<int, 118> _elementOffsets; // Size 118
    std::array<std::vector<primitive>, 118> _convertedData;
    double _beta = 2.0;
};


class BasisSetLibrary {
public:
    // Constructor
    BasisSetLibrary() = default;
    // Access basis set
    std::shared_ptr<BasisSet> get_basis_set(std::string basis_name);
    bool check_basis_set_exists(std::string basis_name);
private:

    BasisSet auxiliary_basis;
};
