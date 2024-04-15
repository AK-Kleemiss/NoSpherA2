/*
 * molecule.cpp
 *
 *  Created on: May 6, 2019
 *      Author: florian
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#ifdef _WIN32
#include <io.h>
#else
#include <unistd.h>
#include <sys/wait.h>
#endif
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <vector>

#include "convenience.h"
#include "molecule.h"
using namespace std;

// UNFINISHED WORKFLOW AND FILE!

bool molecule::read_molecule(std::string path)
{
    if (!exists(path))
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
    getline(file, line);
    size_t length;
    char tempchar[200];
    int count = 0;
    stringstream streamy;
    while (!file.eof())
    {
        length = line.copy(tempchar, line.size(), 0);
        tempchar[length] = '\0';
        string junk;
        streamy << tempchar;
        bonds.push_back(vector<int>());
        bonds[bonds.size() - 1].resize(3);
        streamy >> junk >> bonds[count][0] >> bonds[count][1] >> bonds[count][2];
    }
    return true;
}

bool molecule::input_molecule()
{
    return true;
}
