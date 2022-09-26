#include "spherical_density.h"
#include "convenience.h"

const int Thakkar::first_ex() {
	if (atomic_number == 1) return 0;
	else if (atomic_number > 113) return 200000000;
	int ex = 0;
	for (int i = 0; i < atomic_number - 1; i++)
		ex += nex[i];
	return ex;
};

const int Thakkar::previous_element_coef() {
	if (atomic_number <= 2) return 0;
	int counter = 0;
	for (int temp = atomic_number - 2; temp >= 0; temp--) {
		for (int m = 0; m < 7; m++)
			if (occ[temp * 19 + 0 + m] != 0)
				counter += ns[temp];
		for (int m = 0; m < 6; m++)
			if (occ[temp * 19 + 7 + m] != 0)
				counter += np[temp];
		for (int m = 0; m < 4; m++)
			if (occ[temp * 19 + 13 + m] != 0)
				counter += nd[temp];
		for (int m = 0; m < 2; m++)
			if (occ[temp * 19 + 17 + m] != 0)
				counter += nf[temp];
	}
	return counter;
};

const double Thakkar::get_radial_density(double dist){
	//Speedup things for H
	if (atomic_number == 1)
		return pow(exp(-dist) * 2.0, 2) / (4 * PI);

	double Rho = 0.0;
	int nr_ex = first_ex();
	if (nr_ex == 200000000)
		return -20;
	int nr_coef = previous_element_coef();

	double Orb[19] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	const int offset = (atomic_number - 1) * 19;

	double exponent;
	for (int ex = 0; ex < ns[atomic_number - 1]; ex++) {
		for (int m = 0; m < 7; m++) {
			if (occ[offset + m] == 0) continue;
			exponent = -z[nr_ex] * dist;
			if(exponent > -46.5) { // Corresponds to at least 1E-20
				if (n[nr_ex] == 1)
					Orb[m] += c[nr_coef] * exp(exponent);
				else
					Orb[m] += c[nr_coef] * pow(dist, n[nr_ex] - 1) * exp(exponent);
			}
			nr_coef++;
		}
		nr_ex++;
	}
	for (int ex = 0; ex < np[atomic_number - 1]; ex++) {
		for (int m = 0; m < 6; m++) {
			if (occ[offset + m + 7] == 0) continue;
			exponent = -z[nr_ex] * dist;
			if (exponent > -46.5) {
				if (n[nr_ex] == 1)
					Orb[m + 7] += c[nr_coef] * exp(exponent);
				else
					Orb[m + 7] += c[nr_coef] * pow(dist, n[nr_ex] - 1) * exp(exponent);
			}
			nr_coef++;
		}
		nr_ex++;
	}
	for (int ex = 0; ex < nd[atomic_number - 1]; ex++) {
		for (int m = 0; m < 4; m++) {
			if (occ[offset + m + 13] == 0) continue;
			exponent = -z[nr_ex] * dist;
			if (exponent > -46.5) {
				if (n[nr_ex] == 1)
					Orb[m + 13] += c[nr_coef] * exp(exponent);
				else
					Orb[m + 13] += c[nr_coef] * pow(dist, n[nr_ex] - 1) * exp(exponent);
			}
			nr_coef++;
		}
		nr_ex++;
	}
	for (int ex = 0; ex < nf[atomic_number - 1]; ex++) {
		for (int m = 0; m < 2; m++) {
			if (occ[offset + m + 17] == 0) continue;
			exponent = -z[nr_ex] * dist;
			if (ex > -46.5) {
				if (n[nr_ex] == 1)
					Orb[m + 17] += c[nr_coef] * exp(exponent);
				else
					Orb[m + 17] += c[nr_coef] * pow(dist, n[nr_ex] - 1) * exp(exponent);
			}
			nr_coef++;
		}
		nr_ex++;
	}

	for (int m = 0; m < 19; m++) {
		if (Orb[m] == 0 || occ[offset + m] == 0)
			continue;
		Rho += occ[offset + m] * pow(Orb[m], 2);
	}
	return Rho / (4 * PI);
};

double cosinus_integral(const int N, const double z, const double k);

double sinus_integral(const int N, const double z, const double k) {
	//Calculates the integral 0 - inf r ^ N e ^ -zr sin(kr) dr through recursion using the general integral int[0-inf] e^ax sin(bx) dx = e^ax/(a^2+b^2) * (a sin(bx) - b cos(bx)) and partial integration
	if (N == 0)
		return k / (z * z + k * k);
	else
		return N / (z * z + k * k) * (z * sinus_integral(N - 1, z, k) + k * cosinus_integral(N - 1, z, k));
};

double cosinus_integral(const int N, const double z, const double k) {
	//Calculates the integral 0 - inf r ^ N e ^ -zr cos(kr) dr through recursion using the general integral int[0-inf] e^ax cos(bx) dx = e^ax/(a^2+b^2) * (a cos(bx) - b sin(bx)) and partial integration
	if (N == 0)
		return z / (z * z + k * k);
	else
		return N / (z * z + k * k) * (z * cosinus_integral(N - 1, z, k) - k * sinus_integral(N - 1, z, k));
};

const double Thakkar::get_form_factor(const double k_vector, std::ofstream& file, bool debug) {
	double result(0.0);
	using namespace std;
	const double local_k = cubic_ang2bohr(k_vector); // since the radial exponents are in a.u.

	const int l_ns = ns[atomic_number - 1];
	const int l_np = np[atomic_number - 1];
	const int l_nd = nd[atomic_number - 1];
	const int l_nf = nf[atomic_number - 1];
	
	int nr_ex = first_ex();
	if (nr_ex == 200000000)
		return -20;
	int nr_coef = previous_element_coef();
	if (atomic_number == 1) nr_coef = 0;//return 16.0/pow(4.0+pow(k_vector*2.13*PI,2),2);
	const int offset = (atomic_number - 1) * 19;
	double temp;
	int i_j_distance = 0;
	for (int m = 0; m < 7; m++)
		if (occ[offset + m] != 0)
			i_j_distance++;
	if (debug) file << "s type i-j-dist: " << i_j_distance << endl;
	for (int m = 0; m < 7; m++) {
		if (occ[offset + m] == 0) continue;
		if (debug) file << endl << "m=" << m << " occ= " << occ[offset + m] << endl;
		for (int i = 0; i < l_ns; i++) {
			if (debug) file << "i=" << i << "\n  c = " << c[nr_coef + m + i * i_j_distance] << " n = " << n[nr_ex + i] << " z = " << z[nr_ex + i] << endl;
			for (int j = 0; j < l_ns-i; j++) {
				temp = occ[offset + m] 
					* c[nr_coef + m + i * i_j_distance] * c[nr_coef + m + (i + j) * i_j_distance]
					* sinus_integral(n[nr_ex + i] + n[nr_ex + i + j] - 1, (z[nr_ex + i] + z[nr_ex + i + j]), local_k);
				if (j != 0)
					result += 2 * temp;
				else
					result += temp;
				if (debug) file << " result=" << result / (k_vector * 2* PI) << endl;
			}
		}
	}
	if (debug) file << "old coef:" << nr_coef << " nr_coef += " << l_ns * i_j_distance;
	nr_coef += i_j_distance * l_ns;
	if (debug) file << " new coef:" << nr_coef << endl;
	nr_ex += l_ns;
	i_j_distance = 0;
	for (int m = 7; m < 13; m++)
		if (occ[offset + m] != 0)
			i_j_distance++;
	if (debug) file << "p type i-j-dist: " << i_j_distance << " np=" << l_np << endl;
	for (int m = 7; m < 13; m++) {
		if (occ[offset + m] == 0) continue;
		if (debug) file << endl << "m=" << m << " occ= " << occ[offset + m]  << endl;
		for (int i = 0; i < l_np; i++) {
			if (debug) file << "i=" << i << "\n  c = " << c[nr_coef + m-7 + i * i_j_distance] << " n = " << n[nr_ex + i] << " z = " << z[nr_ex + i] << endl;
			for (int j = 0; j < l_np - i; j++) {
				if (debug) file << "    j=" << j << "\n      c = " << c[nr_coef + m-7 + (i + j) * i_j_distance] << " n = " << n[nr_ex + i + j] << " z = " << z[nr_ex + i + j] << endl;
				temp = occ[offset + m] 
					* c[nr_coef + m-7 + i*i_j_distance] * c[nr_coef + m-7 + (i + j) * i_j_distance] 
					* sinus_integral(n[nr_ex + i] + n[nr_ex + i + j] - 1, z[nr_ex + i] + z[nr_ex + i + j], local_k);
				if (j != 0)
					result += 2 * temp;
				else
					result += temp;
			}
		}
	}
	if (debug) file << "old coef:" << nr_coef << " nr_coef += " << l_np * i_j_distance;
	nr_coef += i_j_distance * l_np;
	if (debug) file << " new coef:" << nr_coef << endl;
	nr_ex += l_np;
	i_j_distance = 0;
	for (int m = 13; m < 17; m++)
		if (occ[offset + m] != 0)
			i_j_distance++;
	if (debug) file << "d type i-j-dist: " << i_j_distance << " nd=" << l_nd << endl;
	for (int m = 13; m < 17; m++) {
		if (occ[offset + m] == 0) continue;
		if (debug) file << endl << "m=" << m << " occ= " << occ[offset + m] << endl;
		for (int i = 0; i < l_nd; i++) {
			if (debug) file << "i=" << i << "\n  c = " << c[nr_coef + m - 13 + i * i_j_distance] << " n = " << n[nr_ex + i] << " z = " << z[nr_ex + i] << endl;
			for (int j = 0; j < l_nd - i; j++) {
				if (debug) file << "    j=" << j << "\n      c = " << c[nr_coef + m - 13 + (i + j) * i_j_distance] << " n = " << n[nr_ex + i + j] << " z = " << z[nr_ex + i + j] << endl;
				temp = occ[offset + m] 
					* c[nr_coef + m - 13 + i * i_j_distance] * c[nr_coef + m - 13 + (i + j) * i_j_distance] 
					* sinus_integral(n[nr_ex + i] + n[nr_ex + i + j] - 1, z[nr_ex + i] + z[nr_ex + i + j], local_k);
				if (j != 0)
					result += 2 * temp;
				else
					result += temp;
			}
		}
	}
	if (debug) file << "old coef:" << nr_coef << " nr_coef += " << l_nd * i_j_distance;
	nr_coef += i_j_distance * l_nd;
	if (debug) file << " new coef:" << nr_coef << endl;
	nr_ex += l_nd;
	i_j_distance = 0;
	for (int m = 17; m < 19; m++)
		if (occ[offset + m] != 0)
			i_j_distance++;
	if (debug) file << "d type i-j-dist: " << i_j_distance << " nd=" << l_nd << endl;
	for (int m = 17; m < 19; m++) {
		if (occ[offset + m] == 0) continue;
		if (debug) file << endl << "m=" << m << " occ= " << occ[offset + m] << endl;
		for (int i = 0; i < l_nf; i++) {
			if (debug) file << "i=" << i << "\n  c = " << c[nr_coef + m - 17 + i * i_j_distance] << " n = " << n[nr_ex + i] << " z = " << z[nr_ex + i] << endl;
			for (int j = 0; j < l_nf - i; j++) {
				if (debug) file << "    j=" << j << "\n      c = " << c[nr_coef + m - 17 + (i + j) * i_j_distance] << " n = " << n[nr_ex + i + j] << " z = " << z[nr_ex + i + j] << endl;
				temp = occ[offset + m] 
					* c[nr_coef + m - 17 + i * i_j_distance] * c[nr_coef + m - 17 + (i + j) * i_j_distance] 
					* sinus_integral(n[nr_ex + i] + n[nr_ex + i + j] - 1, z[nr_ex + i] + z[nr_ex + i + j], local_k);
				if (j != 0)
					result += 2 * temp;
				else
					result += temp;
			}
		}
	}

	return result / local_k;
};

const double Thakkar::get_core_form_factor(const double &k_vector, const int &core_els, std::ofstream& file, bool debug) {
	double result(0.0);
	using namespace std;
	const double local_k = cubic_ang2bohr(k_vector); // since the radial exponents are in a.u.

	const int l_ns = ns[atomic_number - 1];
	const int l_np = np[atomic_number - 1];
	const int l_nd = nd[atomic_number - 1];
	const int l_nf = nf[atomic_number - 1];

	int nr_ex = first_ex();
	if (nr_ex == 200000000)
		return -20;
	int nr_coef = previous_element_coef();
	if (atomic_number == 1) nr_coef = 0;//return 16.0/pow(4.0+pow(k_vector*2.13*PI,2),2);
	const int offset = (atomic_number - 1) * 19;
	double temp;
	int i_j_distance = 0;
	for (int m = 0; m < 7; m++)
		if (occ[offset + m] != 0)
			i_j_distance++;
	int max_s=0, max_p=0, max_d=0, max_f=0;
	if (core_els == 28) {
		max_s = 3; max_p = 2; max_d = 1; max_f = 0;
	}
	else if (core_els == 46) {
		max_s = 4; max_p = 3; max_d = 2; max_f = 0;
	}
	else if (core_els == 60) {
		max_s = 4; max_p = 3; max_d = 2; max_f = 1;
	}
	for (int m = 0; m < max_s; m++) {
		if (occ[offset + m] == 0) continue;
		for (int i = 0; i < l_ns; i++) {
			for (int j = 0; j < l_ns - i; j++) {
				temp = occ[offset + m]
					* c[nr_coef + m + i * i_j_distance] * c[nr_coef + m + (i + j) * i_j_distance]
					* sinus_integral(n[nr_ex + i] + n[nr_ex + i + j] - 1, (z[nr_ex + i] + z[nr_ex + i + j]), local_k);
				if (j != 0)
					result += 2 * temp;
				else
					result += temp;
			}
		}
	}
	nr_coef += i_j_distance * l_ns;
	nr_ex += l_ns;
	i_j_distance = 0;
	for (int m = 7; m < 13; m++)
		if (occ[offset + m] != 0)
			i_j_distance++;
	for (int m = 7; m < 7+max_p; m++) {
		if (occ[offset + m] == 0) continue;
		for (int i = 0; i < l_np; i++) {
			for (int j = 0; j < l_np - i; j++) {
				temp = occ[offset + m]
					* c[nr_coef + m - 7 + i * i_j_distance] * c[nr_coef + m - 7 + (i + j) * i_j_distance]
					* sinus_integral(n[nr_ex + i] + n[nr_ex + i + j] - 1, z[nr_ex + i] + z[nr_ex + i + j], local_k);
				if (j != 0)
					result += 2 * temp;
				else
					result += temp;
			}
		}
	}
	nr_coef += i_j_distance * l_np;
	nr_ex += l_np;
	i_j_distance = 0;
	for (int m = 13; m < 17; m++)
		if (occ[offset + m] != 0)
			i_j_distance++;
	for (int m = 13; m < 13+max_d; m++) {
		if (occ[offset + m] == 0) continue;
		for (int i = 0; i < l_nd; i++) {
			for (int j = 0; j < l_nd - i; j++) {
				temp = occ[offset + m]
					* c[nr_coef + m - 13 + i * i_j_distance] * c[nr_coef + m - 13 + (i + j) * i_j_distance]
					* sinus_integral(n[nr_ex + i] + n[nr_ex + i + j] - 1, z[nr_ex + i] + z[nr_ex + i + j], local_k);
				if (j != 0)
					result += 2 * temp;
				else
					result += temp;
			}
		}
	}
	nr_coef += i_j_distance * l_nd;
	nr_ex += l_nd;
	i_j_distance = 0;
	for (int m = 17; m < 19; m++)
		if (occ[offset + m] != 0)
			i_j_distance++;
	for (int m = 17; m < 17+max_f; m++) {
		if (occ[offset + m] == 0) continue;
		for (int i = 0; i < l_nf; i++) {
			for (int j = 0; j < l_nf - i; j++) {
				temp = occ[offset + m]
					* c[nr_coef + m - 17 + i * i_j_distance] * c[nr_coef + m - 17 + (i + j) * i_j_distance]
					* sinus_integral(n[nr_ex + i] + n[nr_ex + i + j] - 1, z[nr_ex + i] + z[nr_ex + i + j], local_k);
				if (j != 0)
					result += 2 * temp;
				else
					result += temp;
			}
		}
	}

	return result / local_k;
};

Thakkar_Anion::Thakkar_Anion(int g_atom_number) {
	if (g_atom_number != 1 && g_atom_number != 6 && g_atom_number != 8 && g_atom_number != 15) err_not_impl_f("Only selected anions are currently defined!", std::cout);
	atomic_number = g_atom_number;
	nex = &(Anion_nex[0]);
	ns = &(Anion_ns[0]);
	np = &(Anion_np[0]);
	nd = &(Anion_nd[0]);
	nf = &(Thakkar_nf[0]);
	occ = &(Anion_occ[0]);
	n = &(Anion_n[0]);
	z = &(Anion_z[0]);
	c = &(Anion_c[0]);
	charge = -1;
};

Thakkar_Cation::Thakkar_Cation(int g_atom_number) {
	if (g_atom_number < 3 || g_atom_number > 55) err_not_impl_f("Atoms with Z < 3 or bigger than 55 are not defined!", std::cout);
	atomic_number = g_atom_number;
	nex = &(Cation_nex[0]);
	ns = &(Cation_ns[0]);
	np = &(Cation_np[0]);
	nd = &(Cation_nd[0]);
	nf = &(Thakkar_nf[0]);
	occ = &(Cation_occ[0]);
	n = &(Cation_n[0]);
	z = &(Cation_z[0]);
	c = &(Cation_c[0]);
	charge = +1;
};