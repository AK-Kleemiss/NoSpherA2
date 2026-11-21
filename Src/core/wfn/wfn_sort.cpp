#include "pch.h"
#include "wfn_class.h"
#include "convenience.h"
#include "mo_class.h"
#include "cube.h"
#include "constants.h"
#include "fchk.h"
#include "basis_set.h"
#include "nos_math.h"
#include "libCintMain.h"
#include "integrator.h"
#include "cell.h"

int WFN::check_order(const bool &debug) const
{
	for (int i = 0; i < ncen; i++)
	{
		if (!get_atom_basis_set_loaded(i))
		{
			std::cout << "Sorry, consistency check only works if basis set is loaded for all atoms!\n"
				<< "Failing atom: " << i << " " << get_atom_label(i) << "\n";
			return -1;
		}
	}
	//---------------------------check if the wfn is in the right order----------------------
	int order = 0;   // 1=gaussian (P=2222 3333 4444) 2=tonto (234 234 234 234) 3=ORCA (423 423 423 423)
	int f_order = 0; // 1=gaussian (F=11 12 13 17 14 15 18 19 16 20) 2=tonto=ORCA 3=ORCA (11 12 13 14 15 17 16 18 19 20) 4=natural (11 12 13 14 15 16 17 18 19 20)
	int primcounter = 0;
	bool order_found = false;
	for (int a = 0; a < get_ncen(); a++)
	{
		for (int s = 0; s < get_atom_shell_count(a); s++)
		{
			int type = get_shell_type(a, s);
			switch (type)
			{
			case 1: // S Orbital
			{
				for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
				{
					if (types[primcounter] != 1)
					{
						order = -1;
						if (debug)
						{
							std::cout << "This should not happen, the order of your file is not ok for S-types! Checked #:" << primcounter << "\n";
						}
					}
					else
						primcounter++;
				}
				break;
			}
			case 2: // P Orbital
			{
				if (order_found)
				{
					if (order == 1)
					{
						for (int r = 0; r < get_atom_shell_primitives(a, s); r++)
						{
							if (types[primcounter] != 2 || types[primcounter + get_atom_shell_primitives(a, s)] != 3 || types[primcounter + 2 * get_atom_shell_primitives(a, s)] != 4)
							{
								if (debug)
								{
									std::cout << "The found order does not match all type entries! primcounter: " << primcounter << "\n";
								}
							}
							primcounter++;
						}
						primcounter += 2 * get_atom_shell_primitives(a, s);
					}
					else if (order == 2)
					{
						for (int r = 0; r < get_atom_shell_primitives(a, s); r++)
						{
							if (types[primcounter] != 2 || types[primcounter + 1] != 3 || types[primcounter + 2] != 4)
							{
								if (debug)
								{
									std::cout << "The found order does not match all type entries! primcounter: " << primcounter << "\n";
								}
								return -1;
							}
							primcounter += 3;
						}
					}
					else if (order == 3)
					{
						for (int r = 0; r < get_atom_shell_primitives(a, s); r++)
						{
							if (types[primcounter] != 4 || types[primcounter + 1] != 2 || types[primcounter + 2] != 3)
							{
								if (debug)
								{
									std::cout << "The found order does not match all type entries! primcounter: " << primcounter << "\n";
								}
								return -1;
							}
							primcounter += 3;
						}
					}
				}
				else
				{
					if (types[primcounter] == 2)
					{
						if (debug && a == 0)
						{
							std::cout << "Seems to be either tonto or gaussian file...\n";
						}
						if (types[primcounter + 1] == 3)
						{
							order = 2;
							order_found = true;
							if (debug)
							{
								std::cout << "This wfn file is in tonto order!\n";
							}
						}
						else if (types[primcounter + 1] == 2 && get_atom_shell_primitives(a, s) > 1)
						{
							order = 1;
							order_found = true;
							if (debug)
							{
								std::cout << "This wfn file is in gaussian order!\n";
							}
						}
						else
						{
							order = 1;
							order_found = true;
							if (debug)
							{
								std::cout << "Either some error or this shell only has 1 p-primitive and "
									<< "i didn't find any order yet... assuming gaussian\n";
							}
						}
					}
					else if (types[primcounter] == 4)
					{
						if (debug && a == 0)
						{
							std::cout << "Seems as if this file was ORCA ordered...\n";
						}
						order = 3;
						if (types[primcounter + 1] == 2)
						{
							if (debug)
							{
								std::cout << "Yep, that's right! making it permanent now!\n";
							}
							order_found = true;
						}
					}
					else
					{
						std::cout << "I can't recognize this order of the .wfn file...\n";
						return -1;
					}
					s--;
				}
				break;
			}
			case 3: // D Orbital
			{
				if (order_found)
				{
					switch (order)
					{
					case 1:
					{
						for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
						{
							if (types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] != 5 || types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] != 6 || types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] != 7 || types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] != 8 || types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] != 9 || types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] != 10)
							{
								order = -1;
								if (debug)
								{
									std::cout << "The checked types are 6 from #" << primcounter << " and are:\n"
										<< types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] << "\n";
								}
								std::cout << "Something seems to be wrong in the order of your D-Types...\n";
							}
							else
								primcounter += 6;
						}
					}
					break;
					case 2:
					case 3:
					{
						for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
						{
							if (types[get_shell_start_in_primitives(a, s) + 0 + 6 * i] != 5 || types[get_shell_start_in_primitives(a, s) + 1 + 6 * i] != 6 || types[get_shell_start_in_primitives(a, s) + 2 + 6 * i] != 7 || types[get_shell_start_in_primitives(a, s) + 3 + 6 * i] != 8 || types[get_shell_start_in_primitives(a, s) + 4 + 6 * i] != 9 || types[get_shell_start_in_primitives(a, s) + 5 + 6 * i] != 10)
							{
								order = -1;
								if (debug)
								{
									std::cout << "The checked types are 6 from #" << primcounter << " and are:\n"
										<< types[get_shell_start_in_primitives(a, s) + 0 + 6 * i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 1 + 6 * i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 2 + 6 * i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 3 + 6 * i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 4 + 6 * i] << " "
										<< types[get_shell_start_in_primitives(a, s) + 5 + 6 * i] << "\n";
								}
								std::cout << "Something seems to be wrong in the order of your D-Types...\n";
							}
							else
								primcounter += 6;
						}
					}
					break;
					}
				}
				else
				{
					std::cout << "That's highly suspicious, no order for P-type but D-type functions? Let me stop before i do something stupid!\n";
					return -1;
				}
			}
			break;
			case 4:
				for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
				{
					if (types[get_shell_start_in_primitives(a, s) + 0 * get_atom_shell_primitives(a, s) + i] != 11 || types[get_shell_start_in_primitives(a, s) + 1 * get_atom_shell_primitives(a, s) + i] != 12 || types[get_shell_start_in_primitives(a, s) + 2 * get_atom_shell_primitives(a, s) + i] != 13 || types[get_shell_start_in_primitives(a, s) + 3 * get_atom_shell_primitives(a, s) + i] != 17 || types[get_shell_start_in_primitives(a, s) + 4 * get_atom_shell_primitives(a, s) + i] != 14 || types[get_shell_start_in_primitives(a, s) + 5 * get_atom_shell_primitives(a, s) + i] != 15 || types[get_shell_start_in_primitives(a, s) + 6 * get_atom_shell_primitives(a, s) + i] != 18 || types[get_shell_start_in_primitives(a, s) + 7 * get_atom_shell_primitives(a, s) + i] != 19 || types[get_shell_start_in_primitives(a, s) + 8 * get_atom_shell_primitives(a, s) + i] != 16 || types[get_shell_start_in_primitives(a, s) + 9 * get_atom_shell_primitives(a, s) + i] != 20)
					{
						if (types[primcounter] == 11 || types[primcounter + 1] == 12 || types[primcounter + 2] == 13 || types[primcounter + 3] == 14 || types[primcounter + 4] == 15 || types[primcounter + 5] == 16 || types[primcounter + 6] == 17 || types[primcounter + 7] == 18 || types[primcounter + 8] == 19 || types[primcounter + 9] == 20)
						{
							f_order = 4;
							if (debug)
							{
								std::cout << "The checked types are 10 from #" << primcounter << " and are:\n"
									<< types[primcounter] << " " << types[primcounter + 1] << " " << types[primcounter + 2] << " "
									<< types[primcounter + 3] << " " << types[primcounter + 4] << " " << types[primcounter + 5] << " "
									<< types[primcounter + 6] << " " << types[primcounter + 7] << " " << types[primcounter + 8] << " "
									<< types[primcounter + 9] << "\n"
									<< "Appears to be already okay...\n";
							}
							primcounter += 10;
						}
						else if (types[primcounter] != 11 || types[primcounter + 1] != 12 || types[primcounter + 2] != 13 || types[primcounter + 3] != 14 || types[primcounter + 4] != 15 || types[primcounter + 5] != 17 || types[primcounter + 6] != 16 || types[primcounter + 7] != 18 || types[primcounter + 8] != 19 || types[primcounter + 9] != 20)
						{
							order = -1;
							if (debug)
							{
								std::cout << "The checked types are 10 from #" << primcounter << " and are:\n"
									<< types[primcounter] << " " << types[primcounter + 1] << " " << types[primcounter + 2] << " "
									<< types[primcounter + 3] << " " << types[primcounter + 4] << " " << types[primcounter + 5] << " "
									<< types[primcounter + 6] << " " << types[primcounter + 7] << " " << types[primcounter + 8] << " "
									<< types[primcounter + 9] << "\n"
									<< "Something seems to be wrong in the order of your F-Types...\n";
							}
						}
						else
						{
							f_order = 3;
							primcounter += 10;
						}
					}
					else
					{
						f_order = 1;
						primcounter += 10;
					}
				}
				break;
			default:
				std::cout << "ERROR in type assignement!!\n";
				return -1;
				break;
			}
		}
	}
	/*if(debug){
	   std::cout << "Going to return " << f_order*10+order << endl;
		Enter();
	}*/
	return (f_order * 10 + order);
};

bool WFN::sort_wfn(const int &g_order, const bool &debug)
{
	set_modified();
	int primcounter = 0;
	int f_order = 0;
	int order = g_order;
	// Sorry for this way of forwarding the order, i think 2 switches would have been more nicely
	while (order >= 10)
	{
		f_order++;
		order -= 10;
	}
	switch (order)
	{
	case 1:
		for (int a = 0; a < get_ncen(); a++)
			for (int s = 0; s < get_atom_shell_count(a); s++)
				switch (get_shell_type(a, s))
				{
				case 1:
					primcounter += get_atom_shell_primitives(a, s);
					break;
				case 2:
				{
					if (get_atom_shell_primitives(a, s) > 1)
					{
						ivec temp_centers;
						temp_centers.resize(3 * get_atom_shell_primitives(a, s));
						ivec temp_types;
						temp_types.resize(3 * get_atom_shell_primitives(a, s));
						vec temp_exponents;
						temp_exponents.resize(3 * get_atom_shell_primitives(a, s));
						vec2 temp_MO_coefficients;
						temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 3);
						for (int i = 0; i < get_atom_shell_primitives(a, s) * 3; i++)
							temp_MO_coefficients[i].resize(nmo);
						for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
							for (int c = 0; c < 3; c++)
							{
								temp_centers[3 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_types[3 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_exponents[3 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
								for (int m = 0; m < nmo; m++)
									temp_MO_coefficients[3 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
							}
						for (int i = 0; i < 3 * get_atom_shell_primitives(a, s); i++)
						{
							centers[primcounter + i] = temp_centers[i];
							types[primcounter + i] = temp_types[i];
							exponents[primcounter + i] = temp_exponents[i];
							for (int m = 0; m < nmo; m++)
								err_checkf(set_MO_coef(m, primcounter + i, temp_MO_coefficients[i][m]), "Error while assigning new MO coefficient!", std::cout);
						}
					}
					primcounter += get_atom_shell_primitives(a, s) * 3;
					break;
				}
				case 3:
				{
					if (get_atom_shell_primitives(a, s) > 1)
					{
						ivec temp_centers;
						temp_centers.resize(6 * get_atom_shell_primitives(a, s));
						ivec temp_types;
						temp_types.resize(6 * get_atom_shell_primitives(a, s));
						vec temp_exponents;
						temp_exponents.resize(get_atom_shell_primitives(a, s) * 6);
						vec2 temp_MO_coefficients;
						temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 6);
						for (int i = 0; i < get_atom_shell_primitives(a, s) * 6; i++)
							temp_MO_coefficients[i].resize(nmo);
						for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
							for (int c = 0; c < 6; c++)
							{
								temp_centers[6 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_types[6 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_exponents[6 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
								for (int m = 0; m < nmo; m++)
									temp_MO_coefficients[6 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
							}
						for (int i = 0; i < 6 * get_atom_shell_primitives(a, s); i++)
						{
							centers[primcounter + i] = temp_centers[i];
							types[primcounter + i] = temp_types[i];
							exponents[primcounter + i] = temp_exponents[i];
							for (int m = 0; m < nmo; m++)
								if (!set_MO_coef(m, primcounter + i, temp_MO_coefficients[i][m]))
								{
									std::cout << "Error while assigning new MO coefficient!\n";
									return false;
								}
						}
					}
					primcounter += get_atom_shell_primitives(a, s) * 6;
					break;
				}
				case 4:
				{
					if (get_atom_shell_primitives(a, s) > 1)
					{
						ivec temp_centers;
						temp_centers.resize(10 * get_atom_shell_primitives(a, s));
						ivec temp_types;
						temp_types.resize(10 * get_atom_shell_primitives(a, s));
						vec temp_exponents;
						temp_exponents.resize(get_atom_shell_primitives(a, s) * 10);
						vec2 temp_MO_coefficients;
						temp_MO_coefficients.resize(get_atom_shell_primitives(a, s) * 10);
						for (int i = 0; i < get_atom_shell_primitives(a, s) * 10; i++)
							temp_MO_coefficients[i].resize(nmo);
						for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
							for (int c = 0; c < 10; c++)
							{
								temp_centers[10 * i + c] = centers[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_types[10 * i + c] = types[primcounter + i + c * get_atom_shell_primitives(a, s)];
								temp_exponents[10 * i + c] = exponents[primcounter + i + c * get_atom_shell_primitives(a, s)];
								for (int m = 0; m < nmo; m++)
									temp_MO_coefficients[10 * i + c][m] = MOs[m].get_coefficient(primcounter + i + c * get_atom_shell_primitives(a, s));
							}
						for (int i = 0; i < 10 * get_atom_shell_primitives(a, s); i++)
						{
							centers[primcounter + i] = temp_centers[i];
							types[primcounter + i] = temp_types[i];
							exponents[primcounter + i] = temp_exponents[i];
							for (int m = 0; m < nmo; m++)
								if (!set_MO_coef(m, primcounter + i, temp_MO_coefficients[i][m]))
								{
									std::cout << "Error while assigning new MO coefficient!\n";
									return false;
								}
						}
					}
					primcounter += get_atom_shell_primitives(a, s) * 10;
					break;
				}
				}
		break;
	case 2:
		if (debug)
		{
			std::cout << "Nothing to be done here, i like tonto type..." << std::endl;
		}
		break;
	case 3:
		for (int a = 0; a < get_ncen(); a++)
		{
			for (int s = 0; s < get_atom_shell_count(a); s++)
			{
				switch (get_shell_type(a, s))
				{
				case 1:
					primcounter += get_atom_shell_primitives(a, s);
					break;
				case 2:
				{
					int temp_center;
					int temp_type;
					double temp_exponent;
					vec temp_MO_coefficients;
					temp_MO_coefficients.resize(nmo);
					for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
					{
						temp_center = centers[primcounter];
						temp_type = types[primcounter];
						temp_exponent = exponents[primcounter];
						for (int m = 0; m < nmo; m++)
							temp_MO_coefficients[m] = MOs[m].get_coefficient(primcounter);
						for (int j = 0; j < 2; j++)
						{
							centers[primcounter + j] = centers[primcounter + 1 + j];
							types[primcounter + j] = types[primcounter + 1 + j];
							exponents[primcounter + j] = exponents[primcounter + 1 + j];
							for (int m = 0; m < nmo; m++)
							{
								err_checkf(set_MO_coef(m, primcounter + j, MOs[m].get_coefficient(primcounter + 1 + j)), "Error while assigning new MO coefficient!", std::cout);
							}
						}
						centers[primcounter + 2] = temp_center;
						types[primcounter + 2] = temp_type;
						exponents[primcounter + 2] = temp_exponent;
						for (int m = 0; m < nmo; m++)
							if (!set_MO_coef(m, primcounter + 2, temp_MO_coefficients[m]))
							{
								std::cout << "Error while assigning new MO coefficient!\n";
								return false;
							}
						primcounter += 3;
					}
					break;
				}
				case 3:
					primcounter += get_atom_shell_primitives(a, s) * 6;
					break;
				case 4:
					primcounter += get_atom_shell_primitives(a, s) * 10;
					break;
				}
			}
		}
		break;
	default:
		std::cout << "order type: " << f_order << " " << order << " not supported!" << std::endl;
		return false;
		break;
	}
	primcounter = 0;
	switch (f_order)
	{
	case 1:
		for (int a = 0; a < get_ncen(); a++)
			for (int s = 0; s < get_atom_shell_count(a); s++)
				switch (get_shell_type(a, s))
				{
				case 1:
					primcounter += get_atom_shell_primitives(a, s);
					break;
				case 2:
					primcounter += get_atom_shell_primitives(a, s) * 3;
					break;
				case 3:
					primcounter += get_atom_shell_primitives(a, s) * 6;
					break;
				case 4:
				{
					ivec temp_center;
					temp_center.resize(10);
					ivec temp_type;
					temp_type.resize(10);
					vec temp_exponent;
					temp_exponent.resize(10);
					vec2 temp_MO_coefficients;
					temp_MO_coefficients.resize(nmo);
					for (int m = 0; m < nmo; m++)
						temp_MO_coefficients[m].resize(10);
					for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
					{
						for (int j = 0; j < 10; j++)
						{
							temp_center[j] = centers[get_shell_start_in_primitives(a, s) + 10 * i + j];
							temp_type[j] = types[get_shell_start_in_primitives(a, s) + 10 * i + j];
							temp_exponent[j] = exponents[get_shell_start_in_primitives(a, s) + 10 * i + j];
							for (int m = 0; m < nmo; m++)
								temp_MO_coefficients[m][j] = MOs[m].get_coefficient(get_shell_start_in_primitives(a, s) + 10 * i + j);
						}
						// mask[j] is where gaussian's j-th f function (11 12 13 17 14 15 18 19 16 20) lands in natural order
						ivec mask{ 0, 1, 2, 6, 3, 4, 7, 8, 5, 9 };
						for (int j = 0; j < 10; j++)
						{
							centers[get_shell_start_in_primitives(a, s) + 10 * i + mask[j]] = temp_center[j];
							types[get_shell_start_in_primitives(a, s) + 10 * i + mask[j]] = temp_type[j];
							exponents[get_shell_start_in_primitives(a, s) + 10 * i + mask[j]] = temp_exponent[j];
							for (int m = 0; m < nmo; m++)
								set_MO_coef(m, get_shell_start_in_primitives(a, s) + 10 * i + mask[j], temp_MO_coefficients[m][j]);
						}
					}
					primcounter += 10;
				}
				break;
				}
		break;
	case 2:
	case 3:
		for (int a = 0; a < get_ncen(); a++)
			for (int s = 0; s < get_atom_shell_count(a); s++)
				switch (get_shell_type(a, s))
				{
				case 1:
					primcounter += get_atom_shell_primitives(a, s);
					break;
				case 2:
					primcounter += get_atom_shell_primitives(a, s) * 3;
					break;
				case 3:
					primcounter += get_atom_shell_primitives(a, s) * 6;
					break;
				case 4:
				{
					ivec temp_center;
					temp_center.resize(10);
					ivec temp_type;
					temp_type.resize(10);
					vec temp_exponent;
					temp_exponent.resize(10);
					vec2 temp_MO_coefficients;
					temp_MO_coefficients.resize(nmo);
					for (int m = 0; m < nmo; m++)
						temp_MO_coefficients[m].resize(10);
					for (int i = 0; i < get_atom_shell_primitives(a, s); i++)
					{
						for (int j = 0; j < 10; j++)
						{
							temp_center[j] = centers[get_shell_start_in_primitives(a, s) + 10 * i + j];
							temp_type[j] = types[get_shell_start_in_primitives(a, s) + 10 * i + j];
							temp_exponent[j] = exponents[get_shell_start_in_primitives(a, s) + 10 * i + j];
							for (int m = 0; m < nmo; m++)
								temp_MO_coefficients[m][j] = MOs[m].get_coefficient(get_shell_start_in_primitives(a, s) + 10 * i + j);
						}
						ivec mask{ 0, 1, 2, 3, 4, 6, 5, 7, 8, 9 };
						for (int j = 0; j < 10; j++)
						{
							centers[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_center[mask[j]];
							types[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_type[mask[j]];
							exponents[get_shell_start_in_primitives(a, s) + 10 * i + j] = temp_exponent[mask[j]];
							for (int m = 0; m < nmo; m++)
								set_MO_coef(m, get_shell_start_in_primitives(a, s) + 10 * i + mask[j], temp_MO_coefficients[m][j]);
						}
					}
					primcounter += 10;
				}
				break;
				}
		break;
	case 4:
		if (debug)
		{
			std::cout << "This is fine, i like them well ordered!" << std::endl;
		}
		break;
	case 0:
		if (debug)
		{
			std::cout << "There seems to be no f-functions!" << std::endl;
		}
		break;
	default:
		std::cout << "f-order type: " << f_order << " not supported!" << std::endl;
		return false;
		break;
	}
	return true;
};
