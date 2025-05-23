#include "integration_params.h"
#include "JKFit.h"

Int_Params::Int_Params()
{
    wfn_origin = 0;
    ncen = 0;
}

Int_Params::Int_Params(const WFN &wavy)
{
    // Constructor, copying the atoms and basis set from the WFN object
    atoms = wavy.get_atoms();
    wfn_origin = wavy.get_origin();
    ncen = wavy.get_ncen();
    calc_integration_parameters();
}

vec Int_Params::normalize_gto(vec coef, const vec exp, const int l)
{
    // GTO norm Ref: H. B. Schlegel and M. J. Frisch, Int. J. Quant.  Chem., 54(1995), 83-87.
    for (int i = 0; i < coef.size(); i++)
    {
        coef[i] *= 1.0 / std::sqrt(gaussian_int(l * 2 + 2, 2 * exp[i]));
    }

    // Normalize contracted GTO
    // #ee = numpy.empty((nprim, nprim))
    //     #for i in range(nprim) :
    //     #    for j in range(i + 1) :
    //     #        ee[i, j] = ee[j, i] = gaussian_int(angl * 2 + 2, es[i] + es[j])
    //     #s1 = 1 / numpy.sqrt(numpy.einsum('pi,pq,qi->i', cs, ee, cs))
    //     return numpy.einsum('pi,i->pi', cs, s1)
    vec2 ee(coef.size(), vec(coef.size(), 0.0));
    for (int i = 0; i < coef.size(); i++)
    {
        for (int j = 0; j < i + 1; j++)
        {
            ee[i][j] = ee[j][i] = gaussian_int(l * 2 + 2, exp[i] + exp[j]);
        }
    }
    // Do the einsum by hand
    vec s1(coef.size(), 0.0);
    for (int i = 0; i < coef.size(); i++)
    {
        for (int j = 0; j < coef.size(); j++)
        {
            for (int k = 0; k < coef.size(); k++)
            {
                s1[i] += coef[k] * ee[k][j] * coef[j];
            }
        }
    }
    for (int i = 0; i < s1.size(); i++)
    {
        s1[i] = 1.0 / std::sqrt(s1[i]);
    }
    for (int i = 0; i < coef.size(); i++)
    {
        coef[i] *= s1[i];
    }
    return coef;
}

void Int_Params::collect_basis_data()
{
    for (int atom_idx = 0; atom_idx < ncen; atom_idx++)
    {
        // Collect number of basis functions per atom
        nbas += atoms[atom_idx].get_shellcount_size();
        // Skip if the basis set for this element has already been collected
        if (basis_sets.find(atoms[atom_idx].get_charge()) != basis_sets.end())  continue;


        std::vector<basis_set_entry> basis = atoms[atom_idx].get_basis_set();
        // Populate coefficients and exponents vectors
        vec coefficients, exponents;
        for (int shell = 0; shell < basis.size(); shell++)
        {
            coefficients.push_back(basis[shell].get_coefficient());
            exponents.push_back(basis[shell].get_exponent());
        }
        // Normalize the GTOs depending on the context
        if (wfn_origin == 9 || wfn_origin == 6)
        {
            for (int i = 0; i < coefficients.size(); i++)
            {
                int l = basis[i].get_type() - 1;
                coefficients[i] *= std::sqrt(constants::PI * 4 / constants::double_ft[2*l+1]); // Conversion factor from GBW to libcint  ... something something, spherical harmonics...
            }
        }
        else if (wfn_origin == 0)
        {
            int coef_idx = 0;
            for (unsigned int shell = 0; shell < atoms[atom_idx].get_shellcount_size(); shell++)
            {
                vec shell_coefs(coefficients.begin() + coef_idx, coefficients.begin() + atoms[atom_idx].get_shellcount(shell) + coef_idx);
                vec shell_exp(exponents.begin() + coef_idx, exponents.begin() + atoms[atom_idx].get_shellcount(shell) + coef_idx);

                shell_coefs = normalize_gto(shell_coefs, shell_exp, basis[shell].get_type());

                // Place the new coefs at the correct place in the coefficients vector
                std::copy(shell_coefs.begin(), shell_coefs.end(), coefficients.begin() + coef_idx);
                coef_idx += atoms[atom_idx].get_shellcount(shell);
            }
        }
        else
        {
            std::cout << "WFN Origin not recognized, thread carefully! No normalisation was performed!" << std::endl;
        }

        int max_l = 1;
        for (int func = 0; func < basis.size(); func++)
        {
            int new_l = 0;
            if (wfn_origin == 0)      new_l = basis[func].get_type();
            else if (wfn_origin == 9 || wfn_origin == 6) new_l = basis[func].get_type() -1;
            else {
                std::cout << "THIS WFN ORIGIN IS UNTESTED, THREAD CAREFULLY!!!!!" << std::endl;
                new_l = basis[func].get_type() -1;
            }

            if (new_l > max_l) max_l = new_l;
        }

        vec coefficients_new(coefficients.size(), 0.0), exponents_new(exponents.size(), 0.0);
        ivec shellcount_new, shelltype;
        size_t pos_in_new_coeffs = 0;
        for (unsigned int l = 0; l <= max_l; l++) {
            int n_funcs = 0;
            for (unsigned int shell_idx = 0; shell_idx < atoms[atom_idx].get_shellcount_size(); shell_idx++) {
                int curr_funcs = (int)atoms[atom_idx].get_shellcount()[shell_idx];

                if (((basis[n_funcs].get_type()-1 != l) && (wfn_origin == 9 || wfn_origin == 6))  || ((basis[n_funcs].get_type() != l) && (wfn_origin == 0))) { //Sort functions regarding the angular momentum
                    if (wfn_origin != 0 && wfn_origin != 9 && wfn_origin != 6) std::cout << "THIS WFN ORIGIN IS UNTESTED, THREAD CAREFULLY!!!!!" << std::endl;
                    n_funcs += curr_funcs;
                    continue;
                }

                std::copy(coefficients.begin() + n_funcs, coefficients.begin() + n_funcs + curr_funcs, coefficients_new.begin() + pos_in_new_coeffs);
                std::copy(exponents.begin() + n_funcs, exponents.begin() + n_funcs + curr_funcs, exponents_new.begin() + pos_in_new_coeffs);
                shellcount_new.push_back(curr_funcs);
                shelltype.push_back(l);
                pos_in_new_coeffs += curr_funcs;
                n_funcs += curr_funcs;


            }
        }
        //Populate basis_sets dictionary  (Element:[coefficients, exponents, starting point in _env vector, shellcount])
        basis_sets.insert({ atoms[atom_idx].get_charge(), LibCintBasis{coefficients_new, exponents_new, 0, shellcount_new, shelltype}});
    }
}

void Int_Params::populate_atm()
{
    // atom: (Z, ptr, nuclear_model=1, ptr+3, 0,0)  nuclear_model seems to be 1 for all atoms  USURE!!! Nested arrays per atom
    _atm.resize(ncen * 6, 0);
    int ptr = 20;
    for (int atom_idx = 0; atom_idx < ncen; atom_idx++)
    {
        // Assign_atoms
        _atm[atom_idx * 6 + 0] = atoms[atom_idx].get_charge(); // Z
        _atm[atom_idx * 6 + 1] = ptr;                          // ptr
        _atm[atom_idx * 6 + 2] = 1;                            // nuclear_model
        ptr += 3;
        _atm[atom_idx * 6 + 3] = ptr; // ptr+3
        ptr += 1;
    }
}

void Int_Params::populate_env()
{
    // env: (leading 20*0, (x,y,z,zeta=0)*ncen, (coefficients_l=0, exponentsl=0, ...)*ncen )  flat array for everything
    _env.resize(20 + ncen * 4, 0);
    int env_offset = 20; // Do not know if this is needed
    for (int atom_idx = 0; atom_idx < ncen; atom_idx++)
    {
        _env[env_offset + atom_idx * 4] = atoms[atom_idx].get_coordinate(0);
        _env[env_offset + atom_idx * 4 + 1] = atoms[atom_idx].get_coordinate(1);
        _env[env_offset + atom_idx * 4 + 2] = atoms[atom_idx].get_coordinate(2);
        _env[env_offset + atom_idx * 4 + 3] = 0; // zeta
    }
    int bas_ptr = _env.size();
    for (int charge = 1; charge < 118; charge++) {
        //Check if charge is key in basis_sets
        if (basis_sets.find(charge) == basis_sets.end()) {
            continue;
        }

        basis_sets[charge].env_idx = (double)bas_ptr;

        LibCintBasis basis_data = basis_sets[charge];
        vec coefficients = basis_data.coefficients;
        vec exponents = basis_data.exponents;

        int func_count = 0;
        for (int shell = 0; shell < basis_sets[charge].shellcount.size(); shell++) {
            int n_funcs = basis_sets[charge].shellcount[shell];
            for (int func = 0; func < n_funcs; func++) {
                _env.push_back(exponents[func + func_count]);
            }
            for (int func = 0; func < n_funcs; func++)
            {
                _env.push_back(coefficients[func + func_count]);
            }
            func_count += n_funcs;
        }
        bas_ptr += func_count*2;
    }
}

void Int_Params::populate_bas()
{
    // basis: (atom_id, l, nprim, ncentr, kappa=0, ptr, ptr+nprim, 0)  atom_id has to be consistent with order in other arrays  |  Nested arrays for each contracted basis function and then stacked for all atoms
    _bas.resize(8 * nbas, 0);
    int index = 0;
    for (int atom_idx = 0; atom_idx < ncen; atom_idx++)
    {
        LibCintBasis basis = basis_sets[atoms[atom_idx].get_charge()];
        int bas_ptr = basis.env_idx;
        for (int shell = 0; shell < basis.shellcount.size(); shell++)
        {
            _bas[index * 8 + 0] = atom_idx; // atom_id
            _bas[index * 8 + 1] = (int)basis.shelltypes[shell];  // l s=0, p=1, d=2 ...
            _bas[index * 8 + 2] = (int)basis.shellcount[shell];                  // nprim
            _bas[index * 8 + 3] = 1;                                     // ncentr    Not sure
            _bas[index * 8 + 5] = bas_ptr;                               // Pointer to the exponents of a shell
            bas_ptr += (int)basis.shellcount[shell];                     // Add the number of contracted functions per shell to gain the pointer to the coefficient
            _bas[index * 8 + 6] = bas_ptr;                               // Pointer to the coefficients of a shell
            //_bas 8 is set to 0
            bas_ptr += _bas[index * 8 + 2] * _bas[atom_idx * 8 + 3];

            nao += (_bas[index * 8 + 1] * 2 + 1) * _bas[index * 8 + 3]; // 2l+1 + n_centr
            index++;
        }
    }
}

void Int_Params::calc_integration_parameters()
{
    collect_basis_data();
    populate_atm();
    populate_env();
    populate_bas();
}

Int_Params Int_Params::operator+(const Int_Params &other)
{
    // Combine two Int_Params objects
    // Create a new Int_Params object and resize the arrays to fit the new size
    Int_Params combined;
    combined._atm.resize(_atm.size() + other._atm.size(), 0);
    combined._bas.resize(_bas.size() + other._bas.size(), 0);
    combined._env.resize(_env.size() + other._env.size(), 0);
    combined.ncen = ncen + other.ncen;
    combined.nbas = nbas + other.nbas;

    // Copy all the data from the first object into the respective containers in the combined object
    std::copy(_atm.begin(), _atm.end(), combined._atm.begin());
    std::copy(_bas.begin(), _bas.end(), combined._bas.begin());
    std::copy(_env.begin(), _env.end(), combined._env.begin());
    // Also copy the _env data from the second object to the end of the _env data in the combined object
    std::copy(other._env.begin(), other._env.end(), combined._env.begin() + _env.size());

    // Update the pointers in the second object to point to the correct place in the combined _env vector
    unsigned int off = _env.size();
    int natm_off = _atm.size() / 6;
    ivec atm2 = other._atm;
    ivec bas2 = other._bas;
    for (int a = 0; a < natm_off; a++)
    {
        atm2[a * 6 + 1] += off;
        atm2[a * 6 + 3] += off;
    }
    for (int b = 0; b < other.nbas; b++)
    {
        bas2[b * 8 + 0] += natm_off;
        bas2[b * 8 + 5] += off;
        bas2[b * 8 + 6] += off;
    }
    // Copy the data from the second object into the combined objects
    std::copy(atm2.begin(), atm2.end(), combined._atm.begin() + _atm.size());
    std::copy(bas2.begin(), bas2.end(), combined._bas.begin() + _bas.size());
    combined.nao = nao + other.nao;

    return combined;
}

void Int_Params::print_data(std::string name) {
    std::cout << "Printing data for " << name << std::endl;
    std::ofstream file(name + ".txt");
    file << "ATM:" << std::endl;
    for (int a = 0; a < _atm.size()/6; a++) {
        for (int i = 0; i < 6; i++) {
            file << _atm[a * 6 + i] << " ";
        }
        file << std::endl;
    }
    file << "\n\nBAS:" << std::endl;
    for (int b = 0; b < _bas.size() / 8; b++) {
        for (int i = 0; i < 8; i++) {
            file << _bas[b * 8 + i] << " ";
        }
        file << std::endl;
    }
    file << "\n\nENV:" << std::endl;
    for (int e = 0; e < _env.size(); e++) {
        file << _env[e] << " ";
    }
    file.close();
}


ivec make_loc(ivec& bas, int nbas) {
    ivec dims(nbas, 0);
    // Calculate (2*l + 1) * nctr for spherical harmonics
    for (size_t i = 0; i < nbas; i++)
    {
        dims[i] = (2 * bas(ANG_OF, i) + 1) * bas(NCTR_OF, i);
    }

    // Create the ao_loc array
    ivec ao_loc(nbas + 1, 0);

    // Compute the cumulative sum
    std::partial_sum(dims.begin(), dims.end(), ao_loc.begin() + 1);

    return ao_loc;
}

double CINTcommon_fac_sp(int l)
{
    switch (l)
    {
    case 0:
        return 0.282094791773878143;
    case 1:
        return 0.488602511902919921;
    default:
        return 1;
    }
}
