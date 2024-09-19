#include "basis_set.h"
#include "convenience.h"
struct basis_set_entry;
class WFN;

//-------------Reading basis sets and determining the density matrix from the wfn coefficients -----------------------

// read the basis set and give the info to the atoms of the wavefunction
/**
 * Reads a basis set from a specified file path and populates the given WFN object with the basis set information.
 *
 * @param basis_set_path The file path of the basis set.
 * @param wave The WFN object to populate with the basis set information.
 * @param debug A boolean flag indicating whether to enable debug mode.
 * @return Returns true if the basis set is successfully read and populated, false otherwise.
 */
bool read_basis_set(const std::string &basis_set_path, WFN &wave, bool debug)
{
    using namespace std;
    string basis_set_name;
    string temp;
    bool end = false;
    bool manual = false;
    while (!end)
    {
        if (wave.get_basis_set_name().length() < 3)
        {
            cout << "Which basis set do you want to load?" << endl;
            cin >> basis_set_name;
        }
        else if (!manual)
            basis_set_name = wave.get_basis_set_name();
        // assemble basis set name and look if file exists
        temp = basis_set_path;
        temp.append(basis_set_name);
        if (exists(temp))
        {
            if (debug)
                cout << "basis set is valid, continueing..." << endl;
            end = true;
        }
        else
        {
            cout << "sorry, could not find this basis set in the basis set directory specified in the programs.config file!" << endl;
            return false;
            // manual = true;
            // cout << "What is the name of the basis set in the directory: ";
            // cin >> basis_set_name;
        }
    }
    if (debug)
        cout << "File of basis set to load: " << temp << endl;
    ifstream ifile(temp.c_str(), ios::in);
    //  Looking for all the types of atoms we need to find
    svec elements_list;
    bool found = false;
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        if (debug)
            cout << "i: " << i << endl;
        for (int j = 0; j < elements_list.size(); j++)
        {
            if (elements_list[j].compare(wave.get_atom_label(i)) == 0)
                found = true;
            if (debug)
                cout << "   j: " << j << " Atom label: " << wave.get_atom_label(i) << endl;
        }
        if (!found)
        {
            elements_list.push_back(wave.get_atom_label(i));
            if (debug)
                cout << "Added an atom which was not there yet!" << wave.get_atom_label(i) << endl;
        }
        found = false;
    }
    if (debug)
    {
        cout << "Number of elements in elements_list: " << elements_list.size() << endl;
        cout << "This is the elements list:" << endl;
        for (int l = 0; l < elements_list.size(); l++)
            cout << l << ": " << elements_list[l] << endl;
    }
    // int found_counter = 0;
    for (int i = 0; i < elements_list.size(); i++)
    {
        if (debug)
            cout << "before: " << elements_list[i] << " " << i << endl;
        while (elements_list[i].find(" ") != -1)
            elements_list[i].erase(elements_list[i].find(" "), 1);
        elements_list[i].append(":");
        if (debug)
        {
            cout << "after: " << elements_list[i] << " " << i << endl;
        }
        // scan the tonto style basis set file for the entries we are looking or:
        string line;
        getline(ifile, line);
        int file_type = 0;
        // check if we support that type of basis set
        while (line.find("keys=") == string::npos && !ifile.eof())
            getline(ifile, line);
        if (debug)
        {
            cout << "Line after looking for keys=: " << line << endl;
        }
        if (line.find("keys=") < line.size() && debug)
            cout << "Found keys=!" << endl;
        if (line.find("turbomole") < line.size())
        {
            file_type = 1;
            if (debug)
                cout << "This file is written in turbomole type!" << endl;
        }
        else if (line.find("gamess-us") < line.size())
        {
            file_type = 2;
            if (debug)
                cout << "This file is written in gamess-us type!" << endl;
        }
        else if (line.find("gaussian") < line.size())
        {
            file_type = 3;
            if (debug)
                cout << "This file is written in gaussian type!" << endl;
        }
        else if (line.find("CRYSTAL") < line.size())
        {
            file_type = 1;
            wave.set_d_f_switch(true);
            if (debug)
                cout << "This file is written in CRYSTAL type!" << endl;
        }
        else
        {
            cout << "This type of file is not supported, please provide another basis set!" << endl;
            return false;
        }
        if (ifile.eof())
        {
            cout << "Please provide a basis set in the turbomole, gaussian or gamess-us format compatible with tonto."
                 << "Look at the example files \"examble.basis\" and \"examble2.basis\" in the wfn_cpp folder if you want to see how it has to look like" << endl;
            return false;
        }
        while (!(line.find(elements_list[i]) < line.size()) && !ifile.eof())
        {
            getline(ifile, line);
            if (debug)
                cout << "line while search for " << elements_list[i] << " :" << line << endl;
        }
        if (debug && line.find(elements_list[i]) != -1)
        {
            cout << "I found an entry i know from the element list!" << endl;
            cout << "The line is: " << line << endl;
        }
        if (ifile.eof())
        {
            cout << "Could not find the atom you were looking for in the basis set file... " << endl;
            return false;
        }
        unsigned int shell = 0;
        if (line.find("{") == -1)
        {
            getline(ifile, line);
            if (debug)
            {
                cout << "I read an additional line!" << endl;
                // Enter();
            }
        }
        while (line.find("}") == string::npos && !ifile.eof())
        {
            getline(ifile, line);
            stringstream stream;
            stream << line;
            if (line.find("}") != string::npos)
                break;
            int count = 0;
            //int nr_exp = 0;
            double pi = 3.14159265358979;
            char c_temp = '?';
            double temp_vals[2]{0, 0};
            int dum = 0;
            if (file_type == 1)
            {
                stream >> count >> c_temp;
                if (debug)
                    cout << "count: " << count << " type: " << c_temp << endl;
            }
            else if (file_type == 2)
            {
                stream >> c_temp >> count;
                if (debug)
                    cout << "count: " << count << " type: " << c_temp << endl;
            }
            else if (file_type == 3)
            {
                stream >> c_temp >> count;
                if (debug)
                    cout << "count: " << count << " type: " << c_temp << endl;
            }
            for (int j = 0; j < count; j++)
            {
                getline(ifile, line);
                if (debug)
                {
                    cout << "read the " << j << ". line: " << line << endl;
                }
                stringstream stream2;
                stream2 << line;
                if (file_type == 1)
                    stream2 >> temp_vals[0] >> temp_vals[1];
                else if (file_type == 2)
                    stream2 >> dum >> temp_vals[0] >> temp_vals[1];
                else if (file_type == 3)
                    stream2 >> temp_vals[0] >> temp_vals[1];
                switch (c_temp)
                {
                case 's':
                case 'S':
                    temp_vals[1] *= pow((2 * temp[0] / pi), 0.75);
                    break;
                case 'p':
                case 'P':
                    temp_vals[1] *= pow((128 * pow(temp[0], 5) / pow(pi, 3)), 0.25);
                    break;
                case 'd':
                case 'D':
                    temp_vals[1] *= (pow(2048 * pow(temp[0], 7) / (9 * pow(pi, 3)), 0.25));
                    break;
                case 'f':
                case 'F':
                    temp_vals[1] *= pow((32768 * pow(temp[0], 9) / (225 * pow(pi, 3))), 0.25);
                    break;
                default:
                    cout << "Sorry, orbital types higher than f-type are not yet supported!" << endl;
                    return false;
                } // end switch of types
                // this is where i started copying
                for (int h = 0; h < wave.get_ncen(); h++)
                {
                    string temp_label;
                    temp_label = wave.get_atom_label(h);
                    while (temp_label.find(" ") != -1)
                        temp_label.erase(temp_label.find(" "), 1);
                    temp_label.append(":");
                    if (elements_list[i].find(temp_label) != -1)
                    {
                        if (debug)
                        {
                            cout << "It's a match!" << endl;
                            cout << "element_label: " << elements_list[i] << " temp_label: " << temp_label << endl;
                        }
                        switch (c_temp)
                        {
                        case 's':
                        case 'S':
                            if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 1, shell))
                            {
                                cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                                cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                                     << " and exp: " << temp[0] << " and type S" << endl;
                            break;
                        case 'p':
                        case 'P':
                            if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 2, shell))
                            {
                                cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                                cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                                     << " and exp: " << temp[0] << " and type P" << endl;
                            break;
                        case 'd':
                        case 'D':
                            if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 3, shell))
                            {
                                cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                                cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                                     << " and exp: " << temp[0] << " and type D" << endl;
                            break;
                        case 'f':
                        case 'F':
                            if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 4, shell))
                            {
                                cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                                cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                                     << " and exp: " << temp[0] << " and type F" << endl;
                            break;
                        default:
                            cout << "Sorry, orbital types higher than f-type are not yet supported!" << endl;
                            return false;
                        } // end switch of types
                    }     // end if(find atom_label + :
                }         // end for h = ncen
                //nr_exp++;
                if (debug)
                    cout << "recapitulation[" << j << "]... type: " << c_temp << " coef: " << temp[0] << " exp: " << temp[1] << endl;
                if (dum > count)
                {
                    cout << "this should not happen, lets stop before i do something silly!" << endl;
                    return false;
                }
            }
            shell++;
        } // end while line != }
        if (debug)
            cout << "I found }: " << line << endl;
        ifile.seekg(0);
    } // end for element_list.size()
    if (debug)
    {
        cout << "FINISHED WITH READING BASIS SET!" << endl;
    }
    ifile.close();
    return true;
};

/**
 * Reads a basis set from a file and sets it in the given WFN object.
 * Does not perform basis set changes or checks.
 *
 * @param basis_set_path The path to the directory containing the basis set files.
 * @param wave The WFN object to set the basis set in.
 * @param debug A boolean flag indicating whether to enable debug mode.
 * @param manual A boolean flag indicating whether to manually input the basis set name.
 * @return Returns true if the basis set is successfully read and set, false otherwise.
 */
bool read_basis_set_vanilla(const std::string& basis_set_path, WFN& wave, const bool& debug, bool manual)
{
    using namespace std;
    string basis_set_name;
    string temp_name;
    bool end = false;
    while (!end)
    {
        if (wave.get_basis_set_name().length() < 3 || manual)
        {
            cout << "Which basis set do you want to load?" << endl;
            cin >> basis_set_name;
        }
        else if (!manual)
            basis_set_name = wave.get_basis_set_name();
        // assemble basis set name and look if file exists
        temp_name = basis_set_path;
        temp_name.append(basis_set_name);
        if (exists(temp_name))
        {
            if (debug)
                cout << "basis set is valid, continuing..." << endl;
            end = true;
        }
        else
        {
            cout << "sorry, could not find this basis set in the basis set directory specified in the programs.config file!" << endl;
            return false;
            // cout << "What is the name of the basis set in the directory: ";
            // cin >> basis_set_name;
        }
    }
    wave.set_basis_set_name(basis_set_name);
    if (debug)
        cout << "File of basis set to load: " << temp_name << endl;
    ifstream ifile(temp_name.c_str(), ios::in);
    //  Looking for all the types of atoms we need to find
    svec elements_list;
    bool found = false;
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        if (debug)
            cout << "i: " << i << endl;
        for (int j = 0; j < elements_list.size(); j++)
        {
            if (elements_list[j].compare(wave.get_atom_label(i)) == 0)
                found = true;
            if (debug)
                cout << "   j: " << j << " Atom label: " << wave.get_atom_label(i) << endl;
        }
        if (!found)
        {
            elements_list.push_back(wave.get_atom_label(i));
            if (debug)
                cout << "Added an atom which was not there yet! " << wave.get_atom_label(i) << "!" << endl;
        }
        found = false;
    }
    if (debug)
    {
        cout << "Number of elements in elements_list: " << elements_list.size() << endl;
        cout << "This is the elements list:" << endl;
        for (int l = 0; l < elements_list.size(); l++)
            cout << l << ": " << elements_list[l] << "," << endl;
    }
    // int found_counter = 0;
    for (int i = 0; i < elements_list.size(); i++)
    {
        if (debug)
            cout << "before: " << elements_list[i] << " " << i << endl;
        while (elements_list[i].find(" ") != -1)
        {
            elements_list[i].erase(elements_list[i].find(" "), 1);
        } // Cut out spaces
        elements_list[i].append(":");
        if (debug)
        {
            cout << "after: " << elements_list[i] << " " << i << endl;
        }
        // scan the tonto style basis set file for the entries we are looking or:
        string line;
        getline(ifile, line);
        int file_type = 0;
        // check if we support that type of basis set
        while (line.find("keys=") == -1 && !ifile.eof())
            getline(ifile, line);
        if (debug)
        {
            cout << "Line after looking for keys=: " << line << endl;
        }
        if (line.find("keys=") < line.size() && debug)
            cout << "Found keys=!" << endl;
        if (line.find("turbomole") < line.size())
        {
            file_type = 1;
            if (debug)
                cout << "This file is written in turbomole type!" << endl;
        }
        else if (line.find("gamess-us") < line.size())
        {
            file_type = 2;
            if (debug)
                cout << "This file is written in gamess-us type!" << endl;
        }
        else if (line.find("gaussian") < line.size())
        {
            file_type = 3;
            if (debug)
                cout << "This file is written in gaussian type!" << endl;
        }
        else if (line.find("CRYSTAL") < line.size())
        {
            file_type = 1;
            wave.set_d_f_switch(true);
            if (debug)
                cout << "This file is written in CRYSTAL type!" << endl;
        }
        else
        {
            cout << "This type of file is not supported, please provide another basis set!" << endl;
            return false;
        }
        if (ifile.eof())
        {
            cout << "Please provide a basis set in the turbomole, gaussian or gamess-us format compatible with tonto."
                << "Look at the example files \"examble.basis\" and \"examble2.basis\" in the wfn_cpp folder if you want to see how it has to look like" << endl;
            return false;
        }
        while (!(line.find(elements_list[i]) < line.size()) && !ifile.eof())
            getline(ifile, line);
        if (debug)
            cout << "line while search for " << elements_list[i] << " :" << line << endl;
        if (debug && line.find(elements_list[i]) != -1)
        {
            cout << "I found an entry i know from the element list!" << endl;
            cout << "The line is: " << line << endl;
        }
        if (ifile.eof())
        {
            cout << "Could not find the atom you were looking for in the basis set file... " << endl;
            return false;
        }
        unsigned int shell = 0;
        if (line.find("{") == -1)
        {
            getline(ifile, line);
            if (debug)
            {
                cout << "I read an additional line!" << endl;
            }
        }
        while (line.find("}") == string::npos && !ifile.eof())
        {
            getline(ifile, line);
            stringstream stream;
            stream << line;
            if (line.find("}") != string::npos)
                break;
            int count = 0;
            //int nr_exp = 0;
            char c_temp = '?';
            double temp_vals[2]{ 0, 0 };
            int dum = 0;
            if (file_type == 1)
            {
                stream >> count >> c_temp;
                if (debug)
                    cout << "count: " << count << " type: " << c_temp << endl;
            }
            else if (file_type == 2 || file_type == 3)
            {
                stream >> c_temp >> count;
                if (debug)
                    cout << "count: " << count << " type: " << c_temp << endl;
            }
            for (int j = 0; j < count; j++)
            {
                getline(ifile, line);
                if (debug)
                {
                    cout << "read the " << j << ". line: " << line << endl;
                }
                stringstream stream2;
                stream2 << line;
                if (file_type == 1)
                    stream2 >> temp_vals[0] >> temp_vals[1];
                else if (file_type == 2)
                    stream2 >> dum >> temp_vals[0] >> temp_vals[1];
                else if (file_type == 3)
                    stream2 >> temp_vals[0] >> temp_vals[1];
                // this is where i started copying
                for (int h = 0; h < wave.get_ncen(); h++)
                {
                    string temp_label;
                    temp_label = wave.get_atom_label(h);
                    while (temp_label.find(" ") != -1)
                        temp_label.erase(temp_label.find(" "), 1);
                    temp_label.append(":");
                    if (elements_list[i].find(temp_label) != -1)
                    {
                        int type = 0;
                        switch (c_temp)
                        {
                        case 's':
                        case 'S':
                            type = 1;
                            break;
                        case 'p':
                        case 'P':
                            type = 2;
                            break;
                        case 'd':
                        case 'D':
                            type = 3;
                            break;
                        case 'f':
                        case 'F':
                            type = 4;
                            break;
                        default:
                            cout << "Sorry, orbital types higher than f-type are not yet supported!" << endl;
                            return false;
                        } // end switch of types
                        if (!wave.push_back_atom_basis_set(h, temp_vals[0], temp_vals[1], type, shell))
                        {
                            cout << "ERROR while pushing back atoms basis set" << endl;
                        }
                        if (debug)
                            cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp_vals[1]
                            << " and exp: " << temp_vals[0] << " and type " << type << endl;
                    } // end if(find atom_label + :
                }     // end for h = ncen
                //nr_exp++;
                if (debug)
                    cout << "recapitulation[" << j << "]... type: " << c_temp << " coef: " << temp_vals[0] << " exp: " << temp_vals[1] << endl;
                if (dum > count)
                {
                    cout << "this should not happen, lets stop before i do something silly!" << endl;
                    return false;
                }
            }
            shell++;
        } // end while line != }
        if (debug)
            cout << "I found }: " << line << endl;
        ifile.seekg(0);
    } // end for element_list.size()
    if (debug)
    {
        cout << "FINISHED WITH READING BASIS SET!" << endl;
    }
    ifile.close();
    return true;
};

/**
 * Reads a basis set from a specified file path and updates the basis set information in the given WFN object.
 * Only applies it to atoms wihtout basis set info.
 *
 * @param basis_set_path The file path of the basis set.
 * @param wave The WFN object to update with the basis set information.
 * @param debug Flag indicating whether to enable debug mode.
 * @return Returns true if the basis set was successfully read and updated, false otherwise.
 */
bool read_basis_set_missing(const std::string &basis_set_path, WFN &wave, bool debug)
{
    using namespace std;
    string basis_set_name;
    string temp;
    bool end = false;
    bool manual = false;
    while (!end)
    {
        if (wave.get_basis_set_name().length() < 3)
        {
            cout << "Which basis set do you want to load?" << endl;
            cin >> basis_set_name;
        }
        else if (!manual)
            basis_set_name = wave.get_basis_set_name();
        // assemble basis set name and look if file exists
        temp = basis_set_path;
        temp.append(basis_set_name);
        if (exists(temp))
        {
            if (debug)
                cout << "basis set is valid, continueing..." << endl;
            end = true;
        }
        else
        {
            cout << "sorry, could not find this basis set in the basis set directory specified in the programs.config file!" << endl;
            return false;
            // manual = true;
            // cout << "What is the name of the basis set in the directory: ";
            // cin >> basis_set_name;
        }
    }
    if (debug)
        cout << "File of basis set to load: " << temp << endl;
    ifstream ifile(temp.c_str(), ios::in);
    //  Looking for all the types of atoms we need to find
    svec elements_list;
    bool found = false;
    for (int i = 0; i < wave.get_ncen(); i++)
    {
        if (debug)
            cout << "i: " << i << endl;
        for (int j = 0; j < elements_list.size(); j++)
        {
            if (elements_list[j].find(wave.get_atom_label(i)) != -1)
                found = true;
            if (debug)
                cout << "   j: " << j << endl;
        }
        if (!found)
        {
            elements_list.push_back(wave.get_atom_label(i));
            if (debug)
                cout << "Added an atom which was not there yet! " << wave.get_atom_label(i) << endl;
        }
        found = false;
    }
    if (debug)
    {
        cout << "Number of elements in elements_list: " << elements_list.size() << endl;
        cout << "This is the elements list:" << endl;
        for (int l = 0; l < elements_list.size(); l++)
            cout << l << ": " << elements_list[l] << endl;
    }
    // int found_counter = 0;
    for (int i = 0; i < elements_list.size(); i++)
    {
        if (debug)
            cout << "before: " << elements_list[i] << " " << i << endl;
        if (elements_list[i].find(" "))
            elements_list[i].erase(elements_list[i].find(" "), 1);
        elements_list[i].append(":");
        if (debug)
        {
            cout << "after: " << elements_list[i] << " " << i << endl;
        }
        // scan the tonto style basis set file for the entries we are looking or:
        string line;
        getline(ifile, line);
        int file_type = 0;
        // check if we support that type of basis set
        while (line.find("keys=") == -1 && !ifile.eof())
        {
            if (debug)
            {
                cout << "line.size of first line: " << line.size() << "line.find(\"keys=\"): " << line.find("keys=") << endl;
            }
            getline(ifile, line);
        }
        if (debug)
        {
            cout << "Line after looking for keys=: " << line << endl;
        }
        if (line.find("keys=") < line.size() && debug)
            cout << "Found keys=!" << endl;
        if (line.find("turbomole") < line.size())
        {
            file_type = 1;
            if (debug)
                cout << "This file is written in turbomole type!" << endl;
        }
        else if (line.find("gamess-us") < line.size())
        {
            file_type = 2;
            if (debug)
                cout << "This file is written in gamess-us type!" << endl;
        }
        else if (line.find("gaussian") < line.size())
        {
            file_type = 3;
            if (debug)
                cout << "This file is written in gaussian type!" << endl;
        }
        else
        {
            cout << "This type of file is not supported, please provide another basis set!" << endl;
            return false;
        }
        if (ifile.eof())
        {
            cout << "Please provide a basis set in the turbomole, gaussian or gamess-us format compatible with tonto."
                 << "Look at the example files \"examble.basis\" and \"examble2.basis\" in the wfn_cpp folder if you want to see how it has to look like" << endl;
            return false;
        }
        while (!(line.find(elements_list[i]) < line.size()) && !ifile.eof())
        {
            getline(ifile, line);
            if (debug)
                cout << "line while search for " << elements_list[i] << " :" << line << endl;
        }
        if (debug && line.find(elements_list[i]) != -1)
        {
            cout << "I found an entry i know from the element list!" << endl;
            cout << "The line is: " << line << endl;
        }
        if (ifile.eof())
        {
            cout << "Could not find the atom you were looking for in the basis set file... " << endl;
            return false;
        }
        unsigned int shell = 0;
        if (line.find("{") == -1)
        {
            getline(ifile, line);
            if (debug)
                cout << "I read an additional line!" << endl;
        }
        while (line.find("}") == -1 && !ifile.eof())
        {
            getline(ifile, line);
            stringstream stream;
            stream << line;
            int count = 0;
            //int nr_exp = 0;
            char c_temp = '?';
            double temp_num[2]{0, 0};
            int dum = 0;
            if (file_type == 1)
            {
                stream >> count >> c_temp;
                if (debug)
                    cout << "count: " << count << " type: " << c_temp << endl;
            }
            else if (file_type == 2)
            {
                stream >> c_temp >> count;
                if (debug)
                    cout << "count: " << count << " type: " << c_temp << endl;
            }
            else if (file_type == 3)
            {
                stream >> c_temp >> count;
                if (debug)
                    cout << "count: " << count << " type: " << c_temp << endl;
            }
            for (int j = 0; j < count; j++)
            {
                getline(ifile, line);
                if (debug)
                {
                    cout << "read the " << j << ". line: " << line << endl;
                }
                stringstream stream2;
                stream2 << line;
                if (file_type == 1)
                    stream2 >> temp_num[0] >> temp_num[1];
                else if (file_type == 2)
                    stream2 >> dum >> temp_num[0] >> temp_num[1];
                else if (file_type == 3)
                    stream2 >> temp_num[0] >> temp_num[1];
                // this is where i started copying
                for (int h = 0; h < wave.get_ncen(); h++)
                {
                    // skip atoms taht already have a basis set!
                    if (wave.get_atom_basis_set_loaded(h))
                        continue;
                    string temp_label;
                    temp_label = wave.get_atom_label(h);
                    if (temp_label.find(" "))
                        temp_label.erase(temp_label.find(" "), 1);
                    temp_label.append(":");
                    if (elements_list[i].find(temp_label) != -1)
                    {
                        if (debug)
                        {
                            cout << "It's a match!" << endl;
                            cout << "element_label: " << elements_list[i] << " temp_label: " << temp_label << endl;
                        }
                        switch (c_temp)
                        {
                        case 's':
                        case 'S':
                            if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 1, shell))
                            {
                                cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                                cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                                     << " and exp: " << temp[0] << " and type S" << endl;
                            break;
                        case 'p':
                        case 'P':
                            if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 2, shell))
                            {
                                cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                                cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                                     << " and exp: " << temp[0] << " and type P" << endl;
                            break;
                        case 'd':
                        case 'D':
                            if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 3, shell))
                            {
                                cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                                cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                                     << " and exp: " << temp[0] << " and type D" << endl;
                            break;
                        case 'f':
                        case 'F':
                            if (!wave.push_back_atom_basis_set(h, temp[0], temp[1], 4, shell))
                            {
                                cout << "ERROR while pushing back atoms basis set" << endl;
                            }
                            if (debug)
                                cout << "Pushing back on atom: " << h + 1 << " with coef: " << temp[1]
                                     << " and exp: " << temp[0] << " and type F" << endl;
                            break;
                        default:
                            cout << "Sorry, orbital types higher than f-type are not yet supported!" << endl;
                            return false;
                        } // end switch of types
                    }     // end if(find atom_label + :
                }         // end for h = ncen
                //nr_exp++;
                if (debug)
                    cout << "recapitulation[" << j << "]... type: " << c_temp << " coef: " << temp[0] << " exp: " << temp[1] << endl;
                if (dum > count)
                {
                    cout << "this should not happen, lets stop before i do something silly!" << endl;
                    return false;
                }
            }
            shell++;
        } // end while line != }
        if (debug)
            cout << "I found }!" << endl;
        ifile.seekg(0);
    } // end for element_list.size()
    if (debug)
    {
        cout << "FINISHED WITH READING BASIS SET!" << endl;
    }
    ifile.close();
    return true;
};

#include "wfn_class.h"
#include "atoms.h"
