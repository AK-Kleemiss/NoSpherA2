#include "pch.h"
#include "convenience.h"
#include "wfn_class.h"

bool is_similar_rel(const double &first, const double &second, const double &tolerance)
{
	double diff = abs(first - second);
	if (diff > abs((first + second + 0.01) * tolerance / 2))
		return false;
	else
		return true;
};

bool is_similar(const double &first, const double &second, const double &tolerance)
{
	double diff = abs(first - second);
	if (diff > pow(10, tolerance))
		return false;
	else
		return true;
};

bool is_similar_abs(const double &first, const double &second, const double &tolerance)
{
	double diff = abs(first - second);
	if (diff > abs(tolerance))
		return false;
	else
		return true;
};

/*
cosinus_annaeherung::cosinus_annaeherung() : mSize(0), mBase_values(nullptr), mStepwidth(1.0) {
	resize(100);
}

void cosinus_annaeherung::resize(size_t size)
{
	mSize = size;
	if (mBase_values) delete[] mBase_values;
	mBase_values = new double[mSize + 1];
#pragma omp parallel for
	for (auto i = 0; i < mSize + 1; i++)  // Fuer einen Werte mehr die Stueststellen speichern
	{
		double y = cos((MPI2 * i) / mSize);
		//std::cout << "resize: i="<<i<<" y=" << y << endl;
		mBase_values[i] = y;
	}
	mStepwidth = MPI2 / size;
}

double cosinus_annaeherung::calculate_error_at(double x) const
{
	return cos(x) - get(x);
}
*/
void copy_file(std::filesystem::path &from, std::filesystem::path &to)
{
	std::ifstream source(from.c_str(), std::ios::binary);
	std::ofstream dest(to.c_str(), std::ios::binary);

	dest << source.rdbuf();

	source.close();
	dest.close();
};

d3 vec_diff(const d3 &a, const d3 &b) {
	return { a[0] - b[0], a[1] - b[1], a[2] - b[2] };
};

d3 vec_cross(const d3 &a, const d3 &b)
{
	return {
		a[1] * b[2] - a[2] * b[1],
		a[2] * b[0] - a[0] * b[2],
		a[0] * b[1] - a[1] * b[0] };
}

double vec_dot(const d3 &a, const d3 &b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
};

//---------------------------Configuration files ---------------------------------------------------

bool check_bohr(const WFN &wave, bool debug)
{
	double min_length = 300.0;
	for (int i = 0; i < wave.get_ncen(); i++)
	{
		d3 atom1 = wave.get_atom_pos(i);
		for (int j = i + 1; j < wave.get_ncen(); j++)
		{
			d3 atom2 = wave.get_atom_pos(j);
			const double length = array_length(atom1, atom2);
			if (debug)
				std::cout << "Length for: " << i << ";" << j << ": " << length << ", min_length: " << min_length << std::endl;
			if (length < min_length)
				min_length = length;
		}
	}
	if (debug)
	{
		if (min_length < 2)
			std::cout << "Decided it's written in Angstrom" << std::endl;
		else
			std::cout << "Decided it's written in Bohr" << std::endl;
	}
	return (!(min_length < 2));
};

std::string go_get_string(std::ifstream &file, std::string search, bool rewind)
{
	if (rewind)
	{
		file.clear();
		file.seekg(0, file.beg);
	}
	std::string line;
	while (line.find(search) == std::string::npos && !file.eof() && getline_universal(file, line))
		continue;
	if (file.eof())
		return "";
	else
		return line;
}

std::string shrink_string(std::string &input)
{
	while (input.find(" ") != -1)
	{
		input.erase(input.find(" "), 1);
	}
	while (input.find("1") != -1)
	{
		input.erase(input.find("1"), 1);
	}
	while (input.find("2") != -1)
	{
		input.erase(input.find("2"), 1);
	}
	while (input.find("3") != -1)
	{
		input.erase(input.find("3"), 1);
	}
	while (input.find("4") != -1)
	{
		input.erase(input.find("4"), 1);
	}
	while (input.find("5") != -1)
	{
		input.erase(input.find("5"), 1);
	}
	while (input.find("6") != -1)
	{
		input.erase(input.find("6"), 1);
	}
	while (input.find("7") != -1)
	{
		input.erase(input.find("7"), 1);
	}
	while (input.find("8") != -1)
	{
		input.erase(input.find("8"), 1);
	}
	while (input.find("9") != -1)
	{
		input.erase(input.find("9"), 1);
	}
	while (input.find("0") != -1)
	{
		input.erase(input.find("0"), 1);
	}
	while (input.find("(") != -1)
	{
		input.erase(input.find("("), 1);
	}
	while (input.find(")") != -1)
	{
		input.erase(input.find(")"), 1);
	}
	return input;
};

std::string shrink_string_to_atom(std::string &input, const int &atom_number)
{
	while (input.find(" ") != -1)
	{
		input.erase(input.find(" "), 1);
	}
	while (input.find("1") != -1)
	{
		input.erase(input.find("1"), 1);
	}
	while (input.find("2") != -1)
	{
		input.erase(input.find("2"), 1);
	}
	while (input.find("3") != -1)
	{
		input.erase(input.find("3"), 1);
	}
	while (input.find("4") != -1)
	{
		input.erase(input.find("4"), 1);
	}
	while (input.find("5") != -1)
	{
		input.erase(input.find("5"), 1);
	}
	while (input.find("6") != -1)
	{
		input.erase(input.find("6"), 1);
	}
	while (input.find("7") != -1)
	{
		input.erase(input.find("7"), 1);
	}
	while (input.find("8") != -1)
	{
		input.erase(input.find("8"), 1);
	}
	while (input.find("9") != -1)
	{
		input.erase(input.find("9"), 1);
	}
	while (input.find("0") != -1)
	{
		input.erase(input.find("0"), 1);
	}
	while (input.find("(") != -1)
	{
		input.erase(input.find("("), 1);
	}
	while (input.find(")") != -1)
	{
		input.erase(input.find(")"), 1);
	}
	std::string temp = constants::atnr2letter(atom_number);
	err_checkf(temp != "PROBLEM", "Problem identifying atoms!", std::cout);
	if (input.find(temp) != 0)
		return temp;
	if (temp != "PROBLEM")
		while (input.size() > temp.size())
			input.pop_back();
	return input;
};

bool read_block_from_fortran_binary(std::ifstream &file, void *Target, const size_t capacity)
{
	int size_begin = 0, size_end = 0;
	file.read(reinterpret_cast<char *>(&size_begin), sizeof(int));
	//A damaged record length used to be written straight past the end of Target
	if (!file.good() || size_begin < 0 || static_cast<size_t>(size_begin) > capacity)
	{
		std::cout << "Error reading block from binary file: record of " << size_begin << " bytes, expected at most " << capacity << std::endl;
		return false;
	}
	file.read(reinterpret_cast<char *>(Target), size_begin);
	file.read(reinterpret_cast<char *>(&size_end), sizeof(int));
	if (!file.good() || size_begin != size_end)
	{
		std::cout << "Error reading block from binary file: " << size_begin << " vs. " << size_end << std::endl;
		return false;
	}
	return true;
}
template<typename T>
bool read_block_from_fortran_binary(std::ifstream &file, std::vector<T> &Target)
{
	int size_begin = 0, size_end = 0;
	file.read(reinterpret_cast<char *>(&size_begin), sizeof(int));
	if (!file.good() || size_begin < 0)
	{
		std::cout << "Error reading block from binary file: record of " << size_begin << " bytes" << std::endl;
		return false;
	}
	Target.resize(size_begin / sizeof(T));
	file.read(reinterpret_cast<char *>(Target.data()), size_begin);
	file.read(reinterpret_cast<char *>(&size_end), sizeof(int));
	if (!file.good() || size_begin != size_end)
	{
		std::cout << "Error reading block from binary file: " << size_begin << " vs. " << size_end << std::endl;
		return false;
	}
	return true;
}
template bool read_block_from_fortran_binary(std::ifstream &file, vec &Target);
template bool read_block_from_fortran_binary(std::ifstream &file, std::vector<int> &Target);
template bool read_block_from_fortran_binary(std::ifstream &file, std::vector<float> &Target);
template bool read_block_from_fortran_binary(std::ifstream &file, std::vector<char> &Target);

primitive::primitive(int c, int t, double e, double coef) : center(c), type(t), exp(e), coefficient(coef)
{
	exp_l_plus_3_2 = pow(exp, type + 1.5);
	norm_const = pow(
		pow(2, 7 + 4 * type) * pow(exp, 3 + 2 * type) / constants::PI / pow(doublefactorial(2 * type + 1), 2),
		0.25);
	normalized_coefficient = coefficient * norm_const;
};

primitive::primitive(const SimplePrimitive &other) : primitive(other.center, other.type, other.exp, other.coefficient) {};

void select_cubes(std::vector<std::vector<unsigned int>> &selection, std::vector<WFN> &wavy, unsigned int nr_of_cubes, bool wfnonly, bool debug)
{
	// asks which wfn to use, if wfnonly is set or whcih cubes up to nr of cubes to use
	// Returns values in selection[0][i] for iths selection of wavefunction and
	//  selection[1][i] for iths selection of cube
	using namespace std;
	std::cout << "Which of the following cubes to use? Need to select " << nr_of_cubes << " file";
	if (nr_of_cubes > 1)
		std::cout << "s in total." << endl;
	else
		std::cout << "." << endl;
	std::cout << endl
		<< endl;
	for (int w = 0; w < wavy.size(); w++)
	{
		stringstream stream;
		std::cout << "_____________________________________________________________" << endl;
		std::cout << "WFN ";
		stream << setw(2) << w;
		std::cout << stream.str() << ") " << wavy[w].get_path().stem() << endl;
		stream.str("");
		for (int c = 0; c < wavy[w].get_cube_count(); c++)
		{
			if (c == 0)
				std::cout << "        |" << endl
				<< "Cube    |" << endl;
			else
				std::cout << "        |" << endl;
			if (!wfnonly)
			{
				std::cout << setw(2) << w;
				std::cout << ".";
				std::cout << setw(2) << c;
			}
			else
				std::cout << "     ";
			std::cout << "   |_ " << wavy[w].get_cube_path(c).stem();
			if (!exists(wavy[w].get_cube_path(c)))
				std::cout << " (MEM ONLY)";
			std::cout << endl;
		}
		std::cout << "_____________________________________________________________" << endl
			<< endl
			<< endl;
	}
	// bool happy = false;
	unsigned int selected_cubes = 0;
	do
	{
		std::cout << "Select " << selected_cubes + 1 << ". ";
		if (wfnonly)
			std::cout << "WFN ";
		else
			std::cout << "cube ";
		std::cout << "please: ";
		string input;
		cin >> input;
		if (!wfnonly)
		{
			if (input.find('.') == string::npos)
			{
				std::cout << "no . found in input!" << endl;
				continue;
			}
		}
		else
		{
			if (input.find('.') == string::npos)
				std::cout << "Ignoring the .!" << endl;
			unsigned int nr_wave = fromString<unsigned int>(input);
			if (nr_wave < 0 || nr_wave >= wavy.size())
			{
				std::cout << "Invalid choice!" << endl;
				continue;
			}
			selected_cubes++;
			selection[0].push_back(nr_wave);
			if (selected_cubes == nr_of_cubes)
				return;
			else
				continue;
		}
		if (debug)
		{
			std::cout << "input: " << input << endl;
			std::cout << "with . found at: " << input.find('.') << endl;
			std::cout << "substr1: " << input.substr(0, input.find('.')) << endl;
			std::cout << "substr2: " << input.substr(input.find('.') + 1) << endl;
		}
		string wave(input.substr(0, input.find('.')));
		string cube(input.substr(input.find('.') + 1));
		unsigned int nr_wave = fromString<unsigned int>(wave);
		int nr_cube = fromString<int>(cube);
		if (debug)
			std::cout << "Translated: " << nr_wave << " " << nr_cube << endl;
		if (nr_wave < 0 || nr_wave >= wavy.size() || nr_cube < 0 || nr_cube >= wavy[nr_wave].get_cube_count())
		{
			std::cout << "Invalid choice!" << endl;
			continue;
		}
		selection[0][selected_cubes] = nr_wave;
		selection[1][selected_cubes] = nr_cube;
		selected_cubes++;
		if (selected_cubes == nr_of_cubes)
		{
			if (debug)
				std::cout << "Going to return!" << endl;
			return;
		}
	} while (true);
};

bool unsaved_files(std::vector<WFN> &wavy)
{
	for (int w = 0; w < wavy.size(); w++)
		for (int c = 0; c < wavy[w].get_cube_count(); c++)
			if (!exists(wavy[w].get_cube_path(c)))
				return true;
	return false;
};

void swap_sort(ivec order, cvec &v)
{
	int i = 0;
	while (i < v.size() - 1)
	{
		int new_index = 0;
		for (int j = i; j < v.size(); j++)
			if (order[j] < order[i])
				new_index++;
		if (new_index > 0)
		{
			std::complex<double> temp = v[i];
			v[i] = v[i + new_index];
			v[i + new_index] = temp;
			int temp2 = order[i];
			order[i] = order[i + new_index];
			order[i + new_index] = temp2;
		}
		else
			i++;
	}
}

void swap_sort_multi(ivec order, std::vector<ivec> &v)
{
	int i = 0;
	ivec temp;
	temp.resize(v.size());
	while (i < order.size() - 1)
	{
		int new_index = 0;
#pragma omp parallel for reduction(+ : new_index)
		for (int j = i; j < order.size(); j++)
			if (order[j] < order[i])
				new_index++;
		if (new_index > 0)
		{
#pragma omp parallel for
			for (int run = 0; run < v.size(); run++)
			{
				temp[run] = v[run][i];
				v[run][i] = v[run][i + new_index];
				v[run][i + new_index] = temp[run];
			}
			int temp2 = order[i];
			order[i] = order[i + new_index];
			order[i + new_index] = temp2;
		}
		else
			i++;
	}
}

double get_lambda_1(double *a)
{
	// Returns the middle eigenvalue (lambda2) of a row-major symmetric 3x3 matrix.
	const double a00 = a[0], a01 = a[1], a02 = a[2];
	const double a11 = a[4], a12 = a[5], a22 = a[8];

	const double p1 = a01 * a01 + a02 * a02 + a12 * a12;
	if (p1 < std::numeric_limits<double>::epsilon())
	{
		double eigs[3] = { a00, a11, a22 };
		std::sort(eigs, eigs + 3);
		return eigs[1];
	}

	const double q = (a00 + a11 + a22) / 3.0;
	const double b00 = a00 - q, b11 = a11 - q, b22 = a22 - q;
	const double p2 = b00 * b00 + b11 * b11 + b22 * b22 + 2.0 * p1;
	const double p = std::sqrt(p2 / 6.0);

	const double c00 = b00 / p, c01 = a01 / p, c02 = a02 / p;
	const double c11 = b11 / p, c12 = a12 / p, c22 = b22 / p;
	double r = 0.5 * (c00 * (c11 * c22 - c12 * c12)
						  - c01 * (c01 * c22 - c12 * c02)
						  + c02 * (c01 * c12 - c11 * c02));
	//r = +-1 is a double eigenvalue, where acos turns an ulp of rounding in r into sqrt(eps) in lambda2
	//(MSVC ARM64 fuses the products above and lands 1 ulp off for {4,1,1}); snap rounding-level misses.
	if (std::abs(r) > 1.0 - 16 * std::numeric_limits<double>::epsilon())
		r = std::copysign(1.0, r);

	const double phi = std::acos(std::clamp(r, -1.0, 1.0)) / 3.0;
	const double eig_max = q + 2.0 * p * std::cos(phi);
	const double eig_min = q + 2.0 * p * std::cos(phi + 2.0 * constants::PI / 3.0);
	return 3.0 * q - eig_max - eig_min;
}

double get_decimal_precision_from_CIF_number(std::string &given_string)
{
	int len = (int)given_string.length();
	int open_bracket = -1;
	int close_bracket = -1;
	int decimal_point = -1;
	// const char* gs = given_string.c_str();
	for (int i = 0; i < len; i++)
	{
		if (given_string[i] == '(' && open_bracket == -1)
		{
			open_bracket = i;
		}
		else if (given_string[i] == ')' && close_bracket == -1)
		{
			close_bracket = i;
		}
		else if (given_string[i] == '.' && decimal_point == -1)
		{
			decimal_point = i;
		}
	}
	double result = 0;
	int precision = 1;
	int size_of_precision = 1;
	if (open_bracket != -1 && close_bracket != -1)
	{
		size_of_precision = close_bracket - open_bracket - 1;
		std::string temp = given_string.substr(open_bracket + 1, size_of_precision);
		precision = std::stoi(temp);
	}
	int digits = 0;
	if (open_bracket != -1 && close_bracket != -1)
	{
		if (decimal_point != -1)
		{
			digits = open_bracket - decimal_point - 1;
			if (digits == 0)
			{
				return 0.001;
			}
		}
		result = abs(precision * pow(10, -digits));
		return result;
	}
	else
		return 0.005;
};

namespace {
	// Captured during static initialisation, so it still refers to the console after
	// run_app() has redirected std::cout into NoSpherA2.log. A function local static
	// would only be captured on the first error, by which time the redirect has happened
	// and every error message would end up in the log file instead of on the console.
	std::streambuf *const initial_coutbuf = std::cout.rdbuf();
	std::streambuf *original_coutbuf()
	{
		return initial_coutbuf;
	}
}

void options::look_for_debug(int &argc, char **argv)
{
	//-v sets it here, before digest_options; per process, so it starts off each run
	ProgressBar::report_counts = false;
	// This loop figures out command line options
	for (int i = 0; i < argc; i++)
	{
		std::string temp = argv[i];
		arguments.push_back(temp);
		if (temp.find("-") > 0)
			continue;
		else if (temp == "-v" || temp == "-v2" || temp == "-debug")
			std::cout << "Turning on verbose mode!" << std::endl, debug = true,
			ProgressBar::report_counts = true;
		else if (temp == "--h" || temp == "-h" || temp == "-help" || temp == "--help")
		{
			// run_app() has already pointed std::cout at NoSpherA2.log, and the
			// exit() below unwinds nothing, so without this the entire help text
			// is written to the log file and the console stays empty.
			std::cout.rdbuf(original_coutbuf());
			std::cout << NoSpherA2_message() << help_message << build_date << std::endl;
			exit(0);
		}
	}
};

bool ends_with(const std::string &str, const std::string &suffix)
{
	if (str.length() >= suffix.length())
	{
		return str.compare(str.length() - suffix.length(), suffix.length(), suffix) == 0;
	}
	return false;
}

double double_from_string_with_esd(std::string in)
{
	if (in.find('(') == std::string::npos)
		return stod(in);
	else
		return stod(in.substr(0, in.find('(')));
}

std::string trim(const std::string &s)
{
	if (s == "")
		return "";
	auto start = s.begin();
	while (start != s.end() && std::isspace(*start))
	{
		start++;
	}

	auto end = s.end();
	do
	{
		end--;
	} while (std::distance(start, end) > 0 && std::isspace(*end));

	return std::string(start, end + 1);
}

int CountWords(const char *str)
{
	if (str == NULL)
		return -1;

	bool inSpaces = true;
	int numWords = 0;

	while (*str != '\0')
	{
		if (std::isspace(*str))
		{
			inSpaces = true;
		}
		else if (inSpaces)
		{
			numWords++;
			inSpaces = false;
		}

		str++;
	}

	return numWords;
};

void print_duration(std::ostream &file, const std::string &description, const std::chrono::microseconds &duration, std::optional<std::chrono::microseconds> total_duration = std::nullopt)
{
	auto mins = std::chrono::duration_cast<std::chrono::minutes>(duration);
	auto secs = std::chrono::duration_cast<std::chrono::seconds>(duration) % 60;
	auto millisecs = std::chrono::duration_cast<std::chrono::milliseconds>(duration) % 1000;

	file << std::setw(35) << std::left << std::setfill(' ') << description << ": " << std::right
		<< std::setw(2) << std::setfill('0') << mins.count() << ":"
		<< std::setw(2) << std::setfill('0') << secs.count() << ":"
		<< std::setw(3) << std::setfill('0') << millisecs.count();
	if (total_duration.has_value())
	{
		double percentage = (double(duration.count()) / total_duration->count()) * 100.0;
		file << "  (" << std::fixed << std::setprecision(2) << percentage << "%)";
	};
	file << std::endl;
	// Disable setfill 0 again
	file << std::setfill(' ');
}

void write_timing_to_file(std::ostream &file,
	std::vector<_time_point> time_points,
	std::vector<std::string> descriptions)
{
	using namespace std;
	// Check if either vector is empty
	if (time_points.empty() || descriptions.empty())
	{
		file << "Error: Empty vector passed to write_timing_to_file" << endl;
		return;
	}
	std::chrono::microseconds total_time = std::chrono::duration_cast<std::chrono::microseconds>(time_points.back() - time_points.front());
	file << "\n\n----------------------------- Time Breakdown! -----------------------------" << endl;
	file << "                                     mm:ss:ms" << endl;
	// Time_points.size()-1 because we are comparing each time point to the next one meaning we need to stop at the second to last element
	for (int i = 0; i < time_points.size() - 1; i++)
	{
		std::chrono::microseconds dur = std::chrono::duration_cast<std::chrono::microseconds>(time_points[i + 1] - time_points[i]);
		print_duration(file, "... for " + descriptions[i], dur, total_time);
	}
	print_duration(file, "Total Time", total_time);
	file << "---------------------------------------------------------------------------" << endl;
}

void remove_empty_elements(svec &input, const std::string &empty)
{
	for (int i = (int)input.size() - 1; i >= 0; i--)
		if (input[i] == empty || input[i] == "")
			input.erase(input.begin() + i);
}

std::chrono::high_resolution_clock::time_point get_time()
{
	// gets the current time using std chrono library
	std::chrono::high_resolution_clock::time_point time = std::chrono::high_resolution_clock::now();
	return time;
}

long long int get_musec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
	// gets the time difference in microseconds
	std::chrono::microseconds musec = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
	return musec.count();
}

long long int get_msec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
	// gets the time difference in milliseconds
	std::chrono::milliseconds msec = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	return msec.count();
}

long long int get_sec(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end)
{
	// gets the time difference in seconds
	std::chrono::seconds sec = std::chrono::duration_cast<std::chrono::seconds>(end - start);
	return sec.count();
}

const int shell2function(const int &type, const int &prim)
{
	switch (type)
	{
	case (-5):
		return -32 + prim;
	case (-4):
		return -21 + prim;
	case (-3):
		return -12 + prim;
	case (-2):
		return -5 + prim;
	case (-1):
		return 1 + prim;
	case (0):
		return 1;
	case (1):
		return 2 + prim;
	case (2):
		return 5 + prim;
	case (3):
		if (prim == 0)
			return 11;
		if (prim == 1)
			return 12;
		if (prim == 2)
			return 13;
		if (prim == 3)
			return 17;
		if (prim == 4)
			return 14;
		if (prim == 5)
			return 15;
		if (prim == 6)
			return 18;
		if (prim == 7)
			return 19;
		if (prim == 8)
			return 16;
		if (prim == 9)
			return 20;
		break;
	case (4):
		return 21 + prim;
	case (5):
		return 36 + prim;
	default:
		return 0;
	}
	return 0;
}

const int sht2nbas(const int &type)
{
	const int st2bas[9]{ 1, 3, 6, 10, 15, 21, 28, 36 };
	const int nst2bas[9]{ 17,15,13,11, 9, 7, 5, 4, 1 };
	if (type >= 0)
		return st2bas[type];
	else
		return nst2bas[8 + type];
};

char asciitolower(char in)
{
	if (in <= 'Z' && in >= 'A')
		return in - ('Z' - 'z');
	return in;
}

int vec_sum(const bvec &in)
{
	int sum = 0;
	for (bool val : in)
	{
		sum += val;
	}
	return sum;
}

int vec_sum(const ivec &in)
{
	int sum = 0;
	for (int val : in)
	{
		sum += val;
	}
	return sum;
}

double vec_sum(const vec &in)
{
	double sum = 0.0;
	for (double val : in)
	{
		sum += val;
	}
	return sum;
}

cdouble vec_sum(const cvec &in)
{
	cdouble res = 0.0;
	for (int i = 0; i < in.size(); i++)
		res += in[i];
	return res;
}

double vec_length(const vec &in)
{
	double sum = 0.0;
	for (double val : in)
	{
		sum += val * val;
	}
	return sqrt(sum);
}


void error_check(const bool condition, const std::source_location loc, const std::string &error_message, std::ostream &log_file)
{
	if (!condition)
	{
		log_file << "Error in " << loc.function_name() << "\n\t\tat: " << loc.file_name() << " line: " << loc.line() << "\n\t\t\t" << error_message << std::endl;
		std::cout.rdbuf(original_coutbuf()); // reset to standard output again
		std::cout << "Error in " << loc.function_name() << " at: " << loc.file_name() << " line: " << loc.line() << " " << error_message << std::endl;
		exit(-1);
	}
};
void not_implemented(const std::source_location loc, const std::string &error_message, std::ostream &log_file)
{
	log_file << loc.function_name() << "\n\t\tat: " << loc.file_name() << " line: " << loc.line() << "\n\t\t\t" << error_message << " not yet implemented!" << std::endl;
	std::cout.rdbuf(original_coutbuf()); // reset to standard output again
	std::cout << "Error in " << loc.function_name() << " at: " << loc.file_name() << " : " << loc.line() << " " << error_message << " not yet implemented!" << std::endl;
	exit(-1);
};

/*
#ifdef _WIN32
// Function to convert a narrow string to a wide string
std::wstring s2ws(const std::string &s)
{
	int len;
	int slength = (int)s.length() + 1;
	len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0);
	std::wstring r(len, L'\0');
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, &r[0], len);
	return r;
}
*/
