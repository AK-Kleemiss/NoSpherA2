/*
 * molecule.cpp
 *
 *  Created on: May 6, 2019
 *      Author: florian
 */
#include "pch.h"

#include "convenience.h"
#include "molecule.h"
using namespace std;

// UNFINISHED WORKFLOW AND FILE!

bool molecule::read_molecule(std::string path)
{
    if (!std::filesystem::exists(path))
    {
        cout << "ERROR: Could not open file for reading!" << endl;
        return false;
    }
    ifstream file(path.c_str());
    string line;
    if (!file.good())
    {
        cout << "ERROR: Could not open file for reading!" << endl;
        return false;
    }
    file.seekg(0);
    getline_universal(file, line);
    int count = 0;
    while (!file.eof())
    {
        stringstream streamy(line);
        string junk;
        bonds.push_back(vector<int>());
        bonds[bonds.size() - 1].resize(3);
        streamy >> junk >> bonds[count][0] >> bonds[count][1] >> bonds[count][2];
        count++;
        getline_universal(file, line);
    }
    return true;
}

bool molecule::input_molecule()
{
    return true;
}
