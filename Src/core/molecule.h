#pragma once
#ifndef __MOLECULE_H__
#define __MOLECULE_H__

#include <iostream>
#include <vector>
#include <string>

#include "atoms.h"

// UNFINISHED!

class molecule
{
private:
    std::vector<int> atoms;
    std::vector<std::vector<int>> bonds;
    std::vector<int> plane_definition;

public:
    molecule();
    molecule(std::vector<int> &atoms, std::vector<std::vector<int>> &bonds, std::vector<int> &plane);

    std::vector<int> get_bond(int nr) { return bonds[nr]; };
    bool read_molecule(std::string path);
    bool input_molecule();
};

#endif
