#include "wfn_class.h"
#include "convenience.h"
#include "spherical_density.h"
#include "cell.h"
#include "cube.h"

using namespace std;

void computeValues(
	double* PosGrid,		// [3] vector with current position on te grid
	WFN &wavy,
	double* Rho,			// Value of Electron Density
	double* normGrad,		// Gradiant Vector
	double* Hess,			// Hessian Matrix, later used to determine lambda2
	double* Elf,			// Value of the ELF
	double* Eli,			// Value of the ELI
	double* Lap				// Value for the Laplacian
)
{
	double* phi = (double*)malloc(sizeof(double) * 10 * wavy.get_nmo());
	if (phi == NULL) {
		free(phi);
		return;
	}
	for (int i = 0; i < 10 * wavy.get_nmo(); i++) phi[i] = 0.0;

	double chi[10];
	double d[3];
	int iat;
	int l[3];
	double ex;
	double xl[3][3];

	for (int j = 0; j < wavy.get_nex(); j++) {
		iat = wavy.get_center(j) - 1;

		type2vector(wavy.get_type(j), l);
		d[0] = PosGrid[0] - wavy.atoms[iat].x;
		d[1] = PosGrid[1] - wavy.atoms[iat].y;
		d[2] = PosGrid[2] - wavy.atoms[iat].z;
		double temp = -wavy.get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++) {
			if (l[k] == 0) {
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1) {
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2) {
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else {
				return;
			}
		}
		double ex2 = 2 * wavy.get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		double temp_ex = pow(ex2, 2);
		chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
		chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
		chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;
		chi[7] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[2][0] * ex;
		chi[8] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[1][0] * ex;
		chi[9] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * ex;
		for (int mo = 0; mo < wavy.get_nmo(); mo++)
			for (int i = 0; i < 10; i++)
				//				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
				phi[mo * 10 + i] += wavy.get_MO_coef(mo,j,false) * chi[i];      //build MO values at this point
	}

	double Grad[3];
	double tau;

	Rho[0] = 0; tau = 0; Elf[0] = 0;
	Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;
	Hess[0] = 0; Hess[1] = 0; Hess[2] = 0;
	Hess[8] = 0; Hess[4] = 0; Hess[5] = 0;

	for (int mo = 0; mo < wavy.get_nmo(); mo++) {
		Rho[0] += wavy.get_MO_occ(mo) * pow(phi[mo * 10], 2);
		Grad[0] += wavy.get_MO_occ(mo) * 2 * phi[mo * 10] * phi[mo * 10 + 1];
		Grad[1] += wavy.get_MO_occ(mo) * 2 * phi[mo * 10] * phi[mo * 10 + 2];
		Grad[2] += wavy.get_MO_occ(mo) * 2 * phi[mo * 10] * phi[mo * 10 + 3];
		Hess[0] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 10] * phi[mo * 10 + 4] + pow(phi[mo * 10 + 1], 2));
		Hess[4] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 10] * phi[mo * 10 + 5] + pow(phi[mo * 10 + 2], 2));
		Hess[8] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 10] * phi[mo * 10 + 6] + pow(phi[mo * 10 + 3], 2));
		Hess[1] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 10] * phi[mo * 10 + 7] + phi[mo * 10 + 1] * phi[mo * 10 + 2]);
		Hess[2] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 10] * phi[mo * 10 + 8] + phi[mo * 10 + 1] * phi[mo * 10 + 3]);
		Hess[5] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 10] * phi[mo * 10 + 9] + phi[mo * 10 + 2] * phi[mo * 10 + 3]);
		tau += wavy.get_MO_occ(mo) * (pow(phi[mo * 10 + 1], 2) + pow(phi[mo * 10 + 2], 2) + pow(phi[mo * 10 + 3], 2));
	}

	free(phi);

	Hess[3] = Hess[1];
	Hess[6] = Hess[2];
	Hess[7] = Hess[5];
	if (Rho[0] > 0) {
		normGrad[0] = alpha * sqrt(Grad[0] * Grad[0] + Grad[1] * Grad[1] + Grad[2] * Grad[2]) / pow(Rho[0], c_43);
		Elf[0] = 1 / (1 + pow(ctelf * pow(Rho[0], c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho[0]), 2));
		Eli[0] = Rho[0] * pow(12 / (Rho[0] * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), c_38);
	}
	Lap[0] = Hess[0] + Hess[4] + Hess[8];

};

void computeELIELF(
	double* PosGrid,		// [3] vector with current position on te grid
	WFN& wavy,
	double* Elf,			// Value of the ELF
	double* Eli				// Value of the ELI
)
{
	double* phi = (double*)malloc(sizeof(double) * 4 * wavy.get_nmo());
	if (phi == NULL) {
		free(phi);
		return;
	}
	for (int i = 0; i < 4 * wavy.get_nmo(); i++) phi[i] = 0.0;

	double chi[4];
	double d[3];
	int iat;
	int l[3];
	double ex;
	double xl[3][3];

	for (int j = 0; j < wavy.get_nex(); j++) {
		iat = wavy.get_center(j) - 1;

		type2vector(wavy.get_type(j), l);
		d[0] = PosGrid[0] - wavy.atoms[iat].x;
		d[1] = PosGrid[1] - wavy.atoms[iat].y;
		d[2] = PosGrid[2] - wavy.atoms[iat].z;
		double temp = -wavy.get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++) {
			if (l[k] == 0) {
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1) {
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2) {
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else {
				return;
			}
		}
		double ex2 = 2 * wavy.get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		for (int mo = 0; mo < wavy.get_nmo(); mo++)
			for (int i = 0; i < 4; i++)
				//				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
				phi[mo * 4 + i] += wavy.get_MO_coef(mo, j, false) * chi[i];      //build MO values at this point
	}

	double Grad[3];
	double tau;
	double Rho;

	Rho = 0; tau = 0;
	Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;

	for (int mo = 0; mo < wavy.get_nmo(); mo++) {
		Rho += wavy.get_MO_occ(mo) * pow(phi[mo * 4], 2);
		Grad[0] += wavy.get_MO_occ(mo) * 2 * phi[mo * 4] * phi[mo * 4 + 1];
		Grad[1] += wavy.get_MO_occ(mo) * 2 * phi[mo * 4] * phi[mo * 4 + 2];
		Grad[2] += wavy.get_MO_occ(mo) * 2 * phi[mo * 4] * phi[mo * 4 + 3];
		tau += wavy.get_MO_occ(mo) * (pow(phi[mo * 4 + 1], 2) + pow(phi[mo * 4 + 2], 2) + pow(phi[mo * 4 + 3], 2));
	}

	free(phi);

	Elf[0] = 1 / (1 + pow(ctelf * pow(Rho, c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
	Eli[0] = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), c_38);

};

void computeELI(
	double* PosGrid,		// [3] vector with current position on te grid
	WFN& wavy,
	double* Eli				// Value of the ELI
)
{
	double* phi = (double*)malloc(sizeof(double) * 4 * wavy.get_nmo());
	if (phi == NULL) {
		free(phi);
		return;
	}
	for (int i = 0; i < 4 * wavy.get_nmo(); i++) phi[i] = 0.0;

	double chi[4];
	double d[3];
	int iat;
	int l[3];
	double ex;
	double xl[3][3];

	for (int j = 0; j < wavy.get_nex(); j++) {
		iat = wavy.get_center(j) - 1;

		type2vector(wavy.get_type(j), l);
		d[0] = PosGrid[0] - wavy.atoms[iat].x;
		d[1] = PosGrid[1] - wavy.atoms[iat].y;
		d[2] = PosGrid[2] - wavy.atoms[iat].z;
		double temp = -wavy.get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++) {
			if (l[k] == 0) {
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1) {
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2) {
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else {
				return;
			}
		}
		double ex2 = 2 * wavy.get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		for (int mo = 0; mo < wavy.get_nmo(); mo++)
			for (int i = 0; i < 4; i++)
				//				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
				phi[mo * 4 + i] += wavy.get_MO_coef(mo, j, false) * chi[i];      //build MO values at this point
	}

	double Grad[3];
	double tau;
	double Rho;

	Rho = 0; tau = 0;
	Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;

	for (int mo = 0; mo < wavy.get_nmo(); mo++) {
		Rho += wavy.get_MO_occ(mo) * pow(phi[mo * 4], 2);
		Grad[0] += wavy.get_MO_occ(mo) * 2 * phi[mo * 4] * phi[mo * 4 + 1];
		Grad[1] += wavy.get_MO_occ(mo) * 2 * phi[mo * 4] * phi[mo * 4 + 2];
		Grad[2] += wavy.get_MO_occ(mo) * 2 * phi[mo * 4] * phi[mo * 4 + 3];
		tau += wavy.get_MO_occ(mo) * (pow(phi[mo * 4 + 1], 2) + pow(phi[mo * 4 + 2], 2) + pow(phi[mo * 4 + 3], 2));
	}

	free(phi);

	Eli[0] = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), c_38);

};

void computeRho(
	double* PosGrid,		// [3] vector with current position on te grid
	WFN& wavy,
	double* Rho				// Value of Rho
)
{
	double* phi = (double*)malloc(sizeof(double) * wavy.get_nmo());
	if (phi == NULL) {
		free(phi);
		return;
	}
	for (int i = 0; i < wavy.get_nmo(); i++) phi[i] = 0.0;

	double d[3];
	int iat;
	int l[3];
	double ex;
	double xl = 1.0;

	for (int j = 0; j < wavy.get_nex(); j++) {
		iat = wavy.get_center(j) - 1;

		type2vector(wavy.get_type(j), l);
		d[0] = PosGrid[0] - wavy.atoms[iat].x;
		d[1] = PosGrid[1] - wavy.atoms[iat].y;
		d[2] = PosGrid[2] - wavy.atoms[iat].z;
		double temp = -wavy.get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
			continue;
		ex = exp(temp);
		xl = 1.0;
		for (int k = 0; k < 3; k++) {
			if (l[k] == 0)		continue;
			else if (l[k] == 1) xl *= d[k];
			else if (l[k] == 2)	xl *= d[k] * d[k];
			else if (l[k] == 3)	xl *= pow(d[k],3);
			else if (l[k] == 4)	xl *= pow(d[k],4);
		}
		for (int mo = 0; mo < wavy.get_nmo(); mo++)
			phi[mo] += wavy.get_MO_coef(mo, j, false) * xl * ex;      //build MO values at this point
	}

	double rho=0;
	for (int mo = 0; mo < wavy.get_nmo(); mo++)
		rho += wavy.get_MO_occ(mo) * pow(phi[mo], 2);
	free(phi);
	Rho[0] = rho;

};

void computeELF(
	double* PosGrid,		// [3] vector with current position on te grid
	WFN& wavy,
	double* Elf				// Value of the ELF
)
{
	double* phi = (double*)malloc(sizeof(double) * 4 * wavy.get_nmo());
	if (phi == NULL) {
		free(phi);
		return;
	}
	for (int i = 0; i < 4 * wavy.get_nmo(); i++) phi[i] = 0.0;

	double chi[4];
	double d[3];
	int iat;
	int l[3];
	double ex;
	double xl[3][3];

	for (int j = 0; j < wavy.get_nex(); j++) {
		iat = wavy.get_center(j) - 1;

		type2vector(wavy.get_type(j), l);
		d[0] = PosGrid[0] - wavy.atoms[iat].x;
		d[1] = PosGrid[1] - wavy.atoms[iat].y;
		d[2] = PosGrid[2] - wavy.atoms[iat].z;
		double temp = -wavy.get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++) {
			if (l[k] == 0) {
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1) {
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2) {
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else {
				return;
			}
		}
		double ex2 = 2 * wavy.get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		for (int mo = 0; mo < wavy.get_nmo(); mo++)
			for (int i = 0; i < 4; i++)
				//				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
				phi[mo * 4 + i] += wavy.get_MO_coef(mo, j, false) * chi[i];      //build MO values at this point
	}

	double Grad[3];
	double tau;
	double Rho;

	Rho = 0; tau = 0; Elf[0] = 0;
	Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;

	for (int mo = 0; mo < wavy.get_nmo(); mo++) {
		Rho += wavy.get_MO_occ(mo) * pow(phi[mo * 4], 2);
		Grad[0] += wavy.get_MO_occ(mo) * 2 * phi[mo * 4] * phi[mo * 4 + 1];
		Grad[1] += wavy.get_MO_occ(mo) * 2 * phi[mo * 4] * phi[mo * 4 + 2];
		Grad[2] += wavy.get_MO_occ(mo) * 2 * phi[mo * 4] * phi[mo * 4 + 3];
		tau += wavy.get_MO_occ(mo) * (pow(phi[mo * 4 + 1], 2) + pow(phi[mo * 4 + 2], 2) + pow(phi[mo * 4 + 3], 2));
	}

	free(phi);

	Elf[0] = 1 / (1 + pow(ctelf * pow(Rho, c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));

};

void computeLapELIELF(
	double* PosGrid,		// [3] vector with current position on te grid
	WFN& wavy,
	double* Elf,			// Value of the ELF
	double* Eli,			// Value of the ELI
	double* Lap				// Value for the Laplacian
)
{
	double* phi = (double*)malloc(sizeof(double) * 7 * wavy.get_nmo());
	if (phi == NULL) {
		free(phi);
		return;
	}
	for (int i = 0; i < 7 * wavy.get_nmo(); i++) phi[i] = 0.0;

	double chi[7];
	double d[3];
	int iat;
	int l[3];
	double ex;
	double xl[3][3];

	for (int j = 0; j < wavy.get_nex(); j++) {
		iat = wavy.get_center(j) - 1;

		type2vector(wavy.get_type(j), l);
		d[0] = PosGrid[0] - wavy.atoms[iat].x;
		d[1] = PosGrid[1] - wavy.atoms[iat].y;
		d[2] = PosGrid[2] - wavy.atoms[iat].z;
		double temp = -wavy.get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++) {
			if (l[k] == 0) {
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1) {
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2) {
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else {
				return;
			}
		}
		double ex2 = 2 * wavy.get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		double temp_ex = pow(ex2, 2);
		chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
		chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
		chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;
		for (int mo = 0; mo < wavy.get_nmo(); mo++)
			for (int i = 0; i < 7; i++)
				//				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
				phi[mo * 7 + i] += wavy.get_MO_coef(mo, j, false) * chi[i];      //build MO values at this point
	}

	double Grad[3];
	double Hess[3];
	double tau;
	double Rho;

	Rho = 0; tau = 0; Elf[0] = 0;
	Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;
	Hess[0] = 0; Hess[1] = 0; Hess[2] = 0;

	for (int mo = 0; mo < wavy.get_nmo(); mo++) {
		Rho += wavy.get_MO_occ(mo) * pow(phi[mo * 7], 2);
		Grad[0] += wavy.get_MO_occ(mo) * 2 * phi[mo * 7] * phi[mo * 7 + 1];
		Grad[1] += wavy.get_MO_occ(mo) * 2 * phi[mo * 7] * phi[mo * 7 + 2];
		Grad[2] += wavy.get_MO_occ(mo) * 2 * phi[mo * 7] * phi[mo * 7 + 3];
		Hess[0] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 7] * phi[mo * 7 + 4] + pow(phi[mo * 7 + 1], 2));
		Hess[1] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 7] * phi[mo * 7 + 5] + pow(phi[mo * 7 + 2], 2));
		Hess[2] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 7] * phi[mo * 7 + 6] + pow(phi[mo * 7 + 3], 2));
		tau += wavy.get_MO_occ(mo) * (pow(phi[mo * 7 + 1], 2) + pow(phi[mo * 7 + 2], 2) + pow(phi[mo * 7 + 3], 2));
	}

	free(phi);

	Elf[0] = 1 / (1 + pow(ctelf * pow(Rho, c_m53) * (tau * 0.5 - 0.125 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2)) / Rho), 2));
	Eli[0] = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), c_38);
	Lap[0] = Hess[0] + Hess[1] + Hess[2];

};

void computeLapELI(
	double* PosGrid,		// [3] vector with current position on te grid
	WFN& wavy,
	double* Eli,			// Value of the ELI
	double* Lap				// Value for the Laplacian
)
{
	double* phi = (double*)malloc(sizeof(double) * 7 * wavy.get_nmo());
	if (phi == NULL) {
		free(phi);
		return;
	}
	for (int i = 0; i < 7 * wavy.get_nmo(); i++) phi[i] = 0.0;

	double chi[7];
	double d[3];
	int iat;
	int l[3];
	double ex;
	double xl[3][3];

	for (int j = 0; j < wavy.get_nex(); j++) {
		iat = wavy.get_center(j) - 1;

		type2vector(wavy.get_type(j), l);
		d[0] = PosGrid[0] - wavy.atoms[iat].x;
		d[1] = PosGrid[1] - wavy.atoms[iat].y;
		d[2] = PosGrid[2] - wavy.atoms[iat].z;
		double temp = -wavy.get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < -34.5388) //corresponds to cutoff of ex < 1E-15
			continue;
		ex = exp(temp);
		for (int k = 0; k < 3; k++) {
			if (l[k] == 0) {
				xl[k][0] = 1.0;
				xl[k][1] = 0.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 1) {
				xl[k][0] = d[k];
				xl[k][1] = 1.0;
				xl[k][2] = 0.0;
			}
			else if (l[k] == 2) {
				xl[k][0] = d[k] * d[k];
				xl[k][1] = 2 * d[k];
				xl[k][2] = 2;
			}
			else if (l[k] == 3) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d[k];
				xl[k][1] = 3 * d2;
				xl[k][2] = 6 * d[k];
			}
			else if (l[k] == 4) {
				double d2 = d[k] * d[k];
				xl[k][0] = d2 * d2;
				xl[k][1] = 4 * d2 * d[k];
				xl[k][2] = 12 * d2;
			}
			else {
				return;
			}
		}
		double ex2 = 2 * wavy.get_exponent(j);
		chi[0] = xl[0][0] * xl[1][0] * xl[2][0] * ex;
		chi[1] = (xl[0][1] - ex2 * pow(d[0], l[0] + 1)) * xl[1][0] * xl[2][0] * ex;
		chi[2] = (xl[1][1] - ex2 * pow(d[1], l[1] + 1)) * xl[0][0] * xl[2][0] * ex;
		chi[3] = (xl[2][1] - ex2 * pow(d[2], l[2] + 1)) * xl[0][0] * xl[1][0] * ex;
		double temp_ex = pow(ex2, 2);
		chi[4] = (xl[0][2] - ex2 * (2 * l[0] + 1) * xl[0][0] + temp_ex * pow(d[0], l[0] + 2)) * xl[1][0] * xl[2][0] * ex;
		chi[5] = (xl[1][2] - ex2 * (2 * l[1] + 1) * xl[1][0] + temp_ex * pow(d[1], l[1] + 2)) * xl[2][0] * xl[0][0] * ex;
		chi[6] = (xl[2][2] - ex2 * (2 * l[2] + 1) * xl[2][0] + temp_ex * pow(d[2], l[2] + 2)) * xl[0][0] * xl[1][0] * ex;
		for (int mo = 0; mo < wavy.get_nmo(); mo++)
			for (int i = 0; i < 7; i++)
				//				if( abs(chi[i]) * coefficients[nprim*maxc[j]+j] > pow(10.0,-10) )
				phi[mo * 7 + i] += wavy.get_MO_coef(mo, j, false) * chi[i];      //build MO values at this point
	}

	double Grad[3];
	double Hess[3];
	double tau;
	double Rho;

	Rho = 0; tau = 0;
	Grad[0] = 0; Grad[1] = 0; Grad[2] = 0;
	Hess[0] = 0; Hess[1] = 0; Hess[2] = 0;

	for (int mo = 0; mo < wavy.get_nmo(); mo++) {
		Rho += wavy.get_MO_occ(mo) * pow(phi[mo * 7], 2);
		Grad[0] += wavy.get_MO_occ(mo) * 2 * phi[mo * 7] * phi[mo * 7 + 1];
		Grad[1] += wavy.get_MO_occ(mo) * 2 * phi[mo * 7] * phi[mo * 7 + 2];
		Grad[2] += wavy.get_MO_occ(mo) * 2 * phi[mo * 7] * phi[mo * 7 + 3];
		Hess[0] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 7] * phi[mo * 7 + 4] + pow(phi[mo * 7 + 1], 2));
		Hess[1] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 7] * phi[mo * 7 + 5] + pow(phi[mo * 7 + 2], 2));
		Hess[2] += wavy.get_MO_occ(mo) * 2 * (phi[mo * 7] * phi[mo * 7 + 6] + pow(phi[mo * 7 + 3], 2));
		tau += wavy.get_MO_occ(mo) * (pow(phi[mo * 7 + 1], 2) + pow(phi[mo * 7 + 2], 2) + pow(phi[mo * 7 + 3], 2));
	}

	free(phi);

	Eli[0] = Rho * pow(12 / (Rho * tau - 0.25 * (pow(Grad[0], 2) + pow(Grad[1], 2) + pow(Grad[2], 2))), c_38);
	Lap[0] = Hess[0] + Hess[1] + Hess[2];

};

void computeMO(
	double* PosGrid,		// [3] vector with current position on te grid
	WFN& wavy,
	double* MO,			// Value of Electron Density
	int mo
)
{
	MO[0] = 0.0;

	double d[3];
	int iat;
	int l[3];
	double ex;
	double xl;

	for (int j = 0; j < wavy.get_nex(); j++) {
		iat = wavy.get_center(j) - 1;

		type2vector(wavy.get_type(j), l);
		d[0] = PosGrid[0] - wavy.atoms[iat].x;
		d[1] = PosGrid[1] - wavy.atoms[iat].y;
		d[2] = PosGrid[2] - wavy.atoms[iat].z;
		double temp = -wavy.get_exponent(j) * (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
		if (temp < -46.0517) //corresponds to cutoff of ex < 1E-20
			continue;
		ex = exp(temp);
		xl = 1.0;
		for (int k = 0; k < 3; k++) {
			if (l[k] == 0) continue;
			else if (l[k] == 1)	xl *= d[k];
			else if (l[k] == 2)	xl *= d[k] * d[k];
			else if (l[k] == 3)	xl *= d[k] * d[k] * d[k];
			else if (l[k] == 4)	xl *= d[k] * d[k] * d[k] * d[k];
		}
		MO[0] += wavy.get_MO_coef(mo, j, false) * xl * ex;      //build MO values at this point
	}
};

double Integrate(int &m, double i, double &expn) {
	int x;
	if (i <= 10) {
		if (expn == 0.0)
			return 0.0;
		double a = m + 0.5;
		double term = 1 / a;
		double partsum = term;
		for (x = 2; x < 50; x++) {
			a++;
			term *= (i / a);
			partsum += term;
			if (term / partsum < 1E-8)
				return 0.5 * partsum * expn;
		}
	}
	else {
		double a = m;
		double b = a + 0.5;
		a -= 0.5;
		double id = 1 / i;
		double approx = 0.88622692 * sqrt(id) * pow(id, m);
		for (x = 0; x < m; x++) {
			b--;
			approx *= b;
		}
		double mult = 0.5 * expn * id;
		if (mult == 0)
			return approx;
		double prop = mult / approx;
		double term = 1;
		double partsum = 1;
		for (x = 1; x < i + m; x++) {
			term *= a * id;
			partsum += term;
			if (abs(term * prop / partsum) < 1E-8)
				return approx - mult * partsum;
			a--;
		}
	}
	return -1;
};

//precomputed factors for [j][l][m][i] values
int pre[9][5][5][9];
void fill_pre() {
	for (int j = 0; j < 9; j++)
		for (int l = 0; l < 5; l++)
			for (int m = 0; m < 5; m++) {
				int imax = min(j, l);
				int imin = max(0, j - m);
				for (int i = imin; i <= imax; i++)
					pre[j][l][m][i] = ft[j] * ft[l] / ft[l - i] / ft[i] * ft[m] / ft[m - j + i] / ft[j - i];
			}
}
double fj(int &j, int &l, int &m, double &aa, double &bb) {
	double temp = 0.0;
	double temp2;
	int a, b;
	for (int i = max(0, j - m); i <= min(j, l); i++) {
		//pre = factorial[l] / factorial[l - i] / factorial[i] * factorial[m] / factorial[m - j + i] / factorial[j - i];
		temp2 = pre[j][l][m][i];
		a = l - i;
		b = m + i - j;
		if (a != 0)	temp2 *= pow(aa, a);
		if (b != 0)	temp2 *= pow(bb, b);
		temp += temp2;
	}
	return temp;
};

int Afac_pre[9][5][9];
void fill_Afac_pre() {
	for (int l = 0; l < 9; l++)
		for (int r = 0; r <= l / 2; r++)
			for (int s = 0; s <= (l - 2 * r) / 2; s++)
				Afac_pre[l][r][s] = ft[r] * ft[s] * ft[l - 2 * r - 2 * s];
}

double Afac(int &l, int &r, int &i, double &PC, double &gamma, double &fjtmp) {
	double temp = fjtmp * pow(0.25 / gamma, r + i) / Afac_pre[l][r][i];
	int num = l - 2 * r - 2 * i;
	if (num != 0)
		temp *= pow(PC, num);
	if (i % 2 == 1)
		return -temp;
	else
		return temp;
}

double computeESP(double* PosGrid, vector<vector<double> > &d2, WFN& wavy)
{
	double ESP = 0;
	double P[3];
	double Pi[3];
	double Pj[3];
	double PC[3];
	double Fn[11];
	double Al[54];
	double Am[54];
	double An[54];
	int maplrsl[54];
	int maplrsm[54];
	int maplrsn[54];
	int l_i[3];
	int l_j[3];
	int iat, jat, MaxFn;
	double ex_sum, sqd, sqpc, prefac, expc, term, addesp, fjtmp, c, twoexpc, iex, jex;

	const int MO = wavy.get_nmo();
	const int nprim = wavy.get_nex();

	double* occ = new double[MO];
	double** coef = new double* [MO];
	double temp;
	int maxl, maxm, maxn;

	for (int i = 0; i < MO; i++) {
		coef[i] = new double[nprim];
		occ[i] = wavy.get_MO_occ(i);
		for (int j = 0; j < nprim; j++)
			coef[i][j] = wavy.get_MO_coef(i, j, false);
	}

	for (int iat = 0; iat < wavy.get_ncen(); iat++) 
		ESP += wavy.get_atom_charge(iat) * pow(sqrt(pow(PosGrid[0] - wavy.atoms[iat].x, 2) + pow(PosGrid[1] - wavy.atoms[iat].y, 2) + pow(PosGrid[2] - wavy.atoms[iat].z, 2)), -1);

	for (int iprim = 0; iprim < nprim; iprim++) {
		iat = wavy.get_center(iprim) - 1;
		type2vector(wavy.get_type(iprim), l_i);
		iex = wavy.get_exponent(iprim);
		for (int jprim = iprim; jprim < nprim; jprim++) {
			jat = wavy.get_center(jprim) - 1;
			type2vector(wavy.get_type(jprim), l_j);
			ex_sum = wavy.get_exponent(iprim) + wavy.get_exponent(jprim);
			jex = wavy.get_exponent(jprim);

			sqd = d2[iat][jat];
			
			prefac = 2 * 3.141592653589793 / ex_sum * exp(-iex * jex * sqd / ex_sum);
			if (prefac < 1E-10) continue;

			P[0] = (wavy.atoms[iat].x * iex + wavy.atoms[jat].x * jex) / ex_sum;
			P[1] = (wavy.atoms[iat].y * iex + wavy.atoms[jat].y * jex) / ex_sum;
			P[2] = (wavy.atoms[iat].z * iex + wavy.atoms[jat].z * jex) / ex_sum;

			Pi[0] = P[0] - wavy.atoms[iat].x;
			Pi[1] = P[1] - wavy.atoms[iat].y;
			Pi[2] = P[2] - wavy.atoms[iat].z;

			Pj[0] = P[0] - wavy.atoms[jat].x;
			Pj[1] = P[1] - wavy.atoms[jat].y;
			Pj[2] = P[2] - wavy.atoms[jat].z;
			for (int k = 0; k < 3; k++) 
				PC[k] = P[k] - PosGrid[k];

			sqpc = pow(PC[0], 2) + pow(PC[1], 2) + pow(PC[2], 2);

			expc = exp(-ex_sum * sqpc);
			MaxFn = 0;
			for (int i = 0; i < 3; i++)
				MaxFn += l_i[i] + l_j[i];
			temp = Integrate(MaxFn, ex_sum * sqpc, expc);
			Fn[MaxFn] = temp;
			twoexpc = 2 * ex_sum * sqpc;
			for (int nu = MaxFn - 1; nu >= 0; nu--)
				Fn[nu] = (expc + twoexpc * Fn[nu + 1]) / (2 * (nu + 1) - 1);

			maxl = -1;
			for (int l = 0; l <= l_i[0] + l_j[0]; l++) {
				if (l % 2 != 1)
					fjtmp = fj(l, l_i[0], l_j[0], Pi[0], Pj[0]);// *factorial[l];
				else
					fjtmp = -fj(l, l_i[0], l_j[0], Pi[0], Pj[0]);// * factorial[l];
				for (int r = 0; r <= l / 2; r++)
					for (int s = 0; s <= (l - 2 * r) / 2; s++) {
						maxl++;
						Al[maxl] = Afac(l, r, s, PC[0], ex_sum, fjtmp);
						maplrsl[maxl] = l - 2 * r - s;
					}
			}
			maxm = -1;
			for (int l = 0; l <= l_i[1] + l_j[1]; l++) {
				if (l % 2 != 1)
					fjtmp = fj(l, l_i[1], l_j[1], Pi[1], Pj[1]);// *factorial[l];
				else
					fjtmp = -fj(l, l_i[1], l_j[1], Pi[1], Pj[1]);// * factorial[l];
				for (int r = 0; r <= l / 2; r++)
					for (int s = 0; s <= (l - 2 * r) / 2; s++) {
						maxm++;
						Am[maxm] = Afac(l, r, s, PC[1], ex_sum, fjtmp);
						maplrsm[maxm] = l - 2 * r - s;
					}
			}
			maxn = -1;
			for (int l = 0; l <= l_i[2] + l_j[2]; l++) {
				if (l % 2 != 1)
					fjtmp = fj(l, l_i[2], l_j[2], Pi[2], Pj[2]);// *factorial[l];
				else
					fjtmp = -fj(l, l_i[2], l_j[2], Pi[2], Pj[2]);// * factorial[l];
				for (int r = 0; r <= l / 2; r++)
					for (int s = 0; s <= (l - 2 * r) / 2; s++) {
						maxn++;
						An[maxn] = Afac(l, r, s, PC[2], ex_sum, fjtmp);
						maplrsn[maxn] = l - 2 * r - s;
					}
			}

			term = 0.0;
			for (int l = 0; l < maxl; l++) {
				if (Al[l] == 0)
					continue;
				for (int m = 0; m < maxm; m++) {
					if (Am[m] == 0)
						continue;
					for (int n = 0; n < maxn; n++) {
						if (An[n] == 0)
							continue;
						term += Al[l] * Am[m] * An[n] * Fn[maplrsl[l] + maplrsm[m] + maplrsn[n]];
					}
				}
			}

			if (term == 0)
				continue;

			if (iprim != jprim)
				term *= 2.0;

			term *= prefac;
			addesp = 0.0;
			for (int mo = 0; mo < MO; mo++)
				addesp += occ[mo] * coef[mo][iprim] * coef[mo][jprim];
			ESP -= addesp * term;
		}
	}
	return ESP;
};

double get_lambda_1(double *a) {
	vector<double> bw, zw;
	double c, g, gapq, h, s, t, tau, term, termp, termq, theta, thresh, w;
	int run = 0;
	double eig1, eig2, eig3;
	double p1 = a[1] * a[1] + a[2] * a[2] + a[5] * a[5];
	if (p1 == 0) {
		eig1 = a[0];
		eig2 = a[4];
		eig3 = a[8];
		if ((eig1 < eig2 && eig2 < eig3) || (eig3 < eig2 && eig2 < eig1))
			return eig2;

		else if ((eig2 < eig1 && eig1 < eig3) || (eig3 < eig1 && eig1 < eig2))
			return eig1;

		else
			return eig3;
	}
	else {
		double q = (a[0] + a[4] + a[8]) / 3;
		double p2 = pow(a[0] - q, 2) + pow(a[4] - q, 2) + pow(a[8] - q, 2) + 2 * p1;
		double p = sqrt(p2 / 6);
		double B[9];
		B[0] = a[0] - q;
		B[1] = a[1];
		B[2] = a[2];
		B[3] = a[3];
		B[4] = a[4] - q;
		B[5] = a[5];
		B[6] = a[6];
		B[7] = a[7];
		B[8] = a[8] - q;
		double r = (B[0]*B[4]*B[8] 
			+ B[1]*B[5]*B[6] 
			+ B[3]*B[4]*B[7] 
			- B[0]*B[5]*B[7] 
			- B[1]*B[3]*B[8] 
			- B[2]*B[4]*B[6]) / 2;
		double phi;
		if (r <= -1)
			phi = PI / 3;
		else if (r >= 1)
			phi = 0;
		else
			phi = acos(r) / 3;
		
		eig1 = q + 2 * p * cos(phi);
		eig3 = q + 2 * p * cos(phi+2*PI/3);
		eig2 = 3 * q - eig1 - eig3;
		if ((eig1 < eig2 && eig2 < eig3) || (eig3 < eig2 && eig2 < eig1))
			return eig2;

		else if ((eig2 < eig1 && eig1 < eig3) || (eig3 < eig1 && eig1 < eig2))
			return eig1;

		else
			return eig3;
	}
};

void Calc_Spherical_Dens(
	cube& CubeSpher,
	WFN& wavy,
	int cpus,
	double radius,
	ofstream& file
)
{
#ifdef _OPENMP
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif

	time_t start;
	time(&start);

	progress_bar* progress = new progress_bar{ file, 50u, "Calculating Spherical Density" };
	const int step = max(floor(3 * CubeSpher.get_size(0) / 20), 1.0);

	vector<Thakkar> atoms;
	for (int a = 0; a < wavy.get_ncen(); a++)
		atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

#pragma omp parallel for schedule(dynamic)
	for (int i = -CubeSpher.get_size(0); i < 2 * CubeSpher.get_size(0); i++) {
		for (int j = -CubeSpher.get_size(1); j < 2 * CubeSpher.get_size(1); j++)
			for (int k = -CubeSpher.get_size(2); k < 2 * CubeSpher.get_size(2); k++) {

				double PosGrid[3];

				PosGrid[0] = i * CubeSpher.get_vector(0, 0) + j * CubeSpher.get_vector(0, 1) + k * CubeSpher.get_vector(0, 2) + CubeSpher.get_origin(0);
				PosGrid[1] = i * CubeSpher.get_vector(1, 0) + j * CubeSpher.get_vector(1, 1) + k * CubeSpher.get_vector(1, 2) + CubeSpher.get_origin(1);
				PosGrid[2] = i * CubeSpher.get_vector(2, 0) + j * CubeSpher.get_vector(2, 1) + k * CubeSpher.get_vector(2, 2) + CubeSpher.get_origin(2);

				bool skip = true;
				for (int a = 0; a < wavy.get_ncen(); a++){
					if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
						skip = false;
				}
				if (skip)
					continue;

				double dens_all = 0.0;
				double dist;
				for (int a = 0; a < wavy.get_ncen(); a++) {
					dist = sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2));
					dens_all += atoms[a].get_radial_density(dist);;
				}

				int temp_i, temp_j, temp_k;
				if (i < 0)
					temp_i = i + CubeSpher.get_size(0);
				else if (i < CubeSpher.get_size(0))
					temp_i = i;
				else
					temp_i = i - CubeSpher.get_size(0);

				if (j < 0)
					temp_j = j + CubeSpher.get_size(1);
				else if (j < CubeSpher.get_size(1))
					temp_j = j;
				else
					temp_j = j - CubeSpher.get_size(1);

				if (k < 0)
					temp_k = k + CubeSpher.get_size(2);
				else if (k < CubeSpher.get_size(2))
					temp_k = k;
				else
					temp_k = k - CubeSpher.get_size(2);

				CubeSpher.set_value(temp_i, temp_j, temp_k, CubeSpher.get_value(temp_i, temp_j, temp_k) + dens_all);
			}
		if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
			progress->write((i + CubeSpher.get_size(0)) / double(3 * CubeSpher.get_size(0)));
	}
	delete(progress);

	time_t end;
	time(&end);
	if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
	else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
	else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Static_Def(
	cube& CubeDEF,
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	double radius,
	ofstream &file
)
{
#ifdef _OPENMP
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif

	time_t start;
	time(&start);

	progress_bar* progress = new progress_bar{ file, 50u, "Calculating Deformation Density" };
	const int step = max(floor(3 * CubeDEF.get_size(0) / 20), 1.0);

	vector<Thakkar> atoms;
	for (int a = 0; a < wavy.get_ncen(); a++)
		atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

#pragma omp parallel for schedule(dynamic)
	for (int i = -CubeDEF.get_size(0); i < 2 * CubeDEF.get_size(0); i++) {
		for (int j = -CubeDEF.get_size(1); j < 2 * CubeDEF.get_size(1); j++)
			for (int k = -CubeDEF.get_size(2); k < 2 * CubeDEF.get_size(2); k++) {

				double PosGrid[3];

				PosGrid[0] = i * CubeDEF.get_vector(0, 0) + j * CubeDEF.get_vector(0, 1) + k * CubeDEF.get_vector(0, 2) + CubeDEF.get_origin(0);
				PosGrid[1] = i * CubeDEF.get_vector(1, 0) + j * CubeDEF.get_vector(1, 1) + k * CubeDEF.get_vector(1, 2) + CubeDEF.get_origin(1);
				PosGrid[2] = i * CubeDEF.get_vector(2, 0) + j * CubeDEF.get_vector(2, 1) + k * CubeDEF.get_vector(2, 2) + CubeDEF.get_origin(2);

				bool skip = true;
				for (int a = 0; a < wavy.get_ncen(); a++)
					if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
						skip = false;
				if (skip)
					continue;

				double dens_all = 0.0;
				double dist;
				for (int a = 0; a < wavy.get_ncen(); a++) {
					dist = sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2));
					dens_all += atoms[a].get_radial_density(dist);
				}

				int temp_i, temp_j, temp_k;
				if (i < 0)
					temp_i = i + CubeDEF.get_size(0);
				else if (i < CubeDEF.get_size(0))
					temp_i = i;
				else
					temp_i = i - CubeDEF.get_size(0);

				if (j < 0)
					temp_j = j + CubeDEF.get_size(1);
				else if (j < CubeDEF.get_size(1))
					temp_j = j;
				else
					temp_j = j - CubeDEF.get_size(1);

				if (k < 0)
					temp_k = k + CubeDEF.get_size(2);
				else if (k < CubeDEF.get_size(2))
					temp_k = k;
				else
					temp_k = k - CubeDEF.get_size(2);

				dens_all -= CubeRho.get_value(temp_i, temp_j, temp_k);
				CubeDEF.set_value(temp_i, temp_j, temp_k, CubeDEF.get_value(temp_i, temp_j, temp_k) - dens_all);
			}
		if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
			progress->write((i + CubeDEF.get_size(0)) / double(3 * CubeDEF.get_size(0)));
	}
	delete(progress);

	time_t end;
	time(&end);
	if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
	else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
	else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Static_Def(
	cube& CubeDEF,
	cube& CubeRho,
	cube& CubeSpher,
	WFN& wavy,
	int cpus,
	double radius,
	ofstream& file
)
{
#ifdef _OPENMP
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif

	time_t start;
	time(&start);

	progress_bar* progress = new progress_bar{ file, 50u, "Calculating Deformation Density" };
	const int step = max(floor(3 * CubeDEF.get_size(0) / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
	for (int i = -CubeDEF.get_size(0); i < 2 * CubeDEF.get_size(0); i++) {
		for (int j = -CubeDEF.get_size(1); j < 2 * CubeDEF.get_size(1); j++)
			for (int k = -CubeDEF.get_size(2); k < 2 * CubeDEF.get_size(2); k++) {

				double PosGrid[3];

				PosGrid[0] = i * CubeDEF.get_vector(0, 0) + j * CubeDEF.get_vector(0, 1) + k * CubeDEF.get_vector(0, 2) + CubeDEF.get_origin(0);
				PosGrid[1] = i * CubeDEF.get_vector(1, 0) + j * CubeDEF.get_vector(1, 1) + k * CubeDEF.get_vector(1, 2) + CubeDEF.get_origin(1);
				PosGrid[2] = i * CubeDEF.get_vector(2, 0) + j * CubeDEF.get_vector(2, 1) + k * CubeDEF.get_vector(2, 2) + CubeDEF.get_origin(2);

				bool skip = true;
				for (int a = 0; a < wavy.get_ncen(); a++)
					if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
						skip = false;
				if (skip)
					continue;

				int temp_i, temp_j, temp_k;
				if (i < 0)
					temp_i = i + CubeDEF.get_size(0);
				else if (i < CubeDEF.get_size(0))
					temp_i = i;
				else
					temp_i = i - CubeDEF.get_size(0);

				if (j < 0)
					temp_j = j + CubeDEF.get_size(1);
				else if (j < CubeDEF.get_size(1))
					temp_j = j;
				else
					temp_j = j - CubeDEF.get_size(1);

				if (k < 0)
					temp_k = k + CubeDEF.get_size(2);
				else if (k < CubeDEF.get_size(2))
					temp_k = k;
				else
					temp_k = k - CubeDEF.get_size(2);

				double rho = CubeRho.get_value(temp_i, temp_j, temp_k);
				double spher = CubeSpher.get_value(temp_i, temp_j, temp_k);
				double temp = rho - spher;
				CubeDEF.set_value(temp_i, temp_j, temp_k, CubeDEF.get_value(temp_i, temp_j, temp_k) + temp);
			}
		if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
			progress->write((i + CubeDEF.get_size(0)) / double(3 * CubeDEF.get_size(0)));
	}
	delete(progress);

	time_t end;
	time(&end);
	if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
	else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
	else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Hirshfeld(
	cube& CubeHDEF,
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	double radius,
	int ignore,
	ofstream & file
)
{
#ifdef _OPENMP
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif

	time_t start;
	time(&start);

	progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
	const int step = max(floor(3*CubeHDEF.get_size(0) / 20), 1.0);

	vector<Thakkar> atoms;
	for (int a = 0; a < wavy.get_ncen(); a++)
		atoms.push_back(Thakkar(wavy.get_atom_charge(a)));

#pragma omp parallel for schedule(dynamic)
	for (int i = -CubeHDEF.get_size(0); i < 2 * CubeHDEF.get_size(0); i++) {
		for (int j = -CubeHDEF.get_size(1); j < 2 * CubeHDEF.get_size(1); j++)
			for (int k = -CubeHDEF.get_size(2); k < 2 * CubeHDEF.get_size(2); k++) {

				double PosGrid[3];

				PosGrid[0] = i * CubeHDEF.get_vector(0, 0) + j * CubeHDEF.get_vector(0, 1) + k * CubeHDEF.get_vector(0, 2) + CubeHDEF.get_origin(0);
				PosGrid[1] = i * CubeHDEF.get_vector(1, 0) + j * CubeHDEF.get_vector(1, 1) + k * CubeHDEF.get_vector(1, 2) + CubeHDEF.get_origin(1);
				PosGrid[2] = i * CubeHDEF.get_vector(2, 0) + j * CubeHDEF.get_vector(2, 1) + k * CubeHDEF.get_vector(2, 2) + CubeHDEF.get_origin(2);

				bool skip = true;
				if (sqrt(pow(PosGrid[0] - wavy.atoms[ignore].x, 2) + pow(PosGrid[1] - wavy.atoms[ignore].y, 2) + pow(PosGrid[2] - wavy.atoms[ignore].z, 2)) < radius / 0.52)
					skip = false;
				if (skip)
					continue;

				double dens_choice = 0.0;
				double dens_all = 0.0;
				double dist, temp;
				for (int a = 0; a < wavy.get_ncen(); a++) {
					dist = sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2));
					temp = atoms[a].get_radial_density(dist);
					if (ignore == a)
						dens_choice = temp;
					dens_all += temp;
				}

				int temp_i, temp_j, temp_k;
				if (i < 0)
					temp_i = i + CubeHDEF.get_size(0);
				else if (i < CubeHDEF.get_size(0))
					temp_i = i;
				else
					temp_i = i - CubeHDEF.get_size(0);

				if (j < 0)
					temp_j = j + CubeHDEF.get_size(1);
				else if (j < CubeHDEF.get_size(1))
					temp_j = j;
				else
					temp_j = j - CubeHDEF.get_size(1);

				if (k < 0)
					temp_k = k + CubeHDEF.get_size(2);
				else if (k < CubeHDEF.get_size(2))
					temp_k = k;
				else
					temp_k = k - CubeHDEF.get_size(2);

				dens_all = dens_choice / dens_all * CubeRho.get_value(temp_i,temp_j,temp_k);
				CubeHDEF.set_value(temp_i, temp_j, temp_k, CubeHDEF.get_value(temp_i, temp_j, temp_k) + dens_all - dens_choice);
			}
		if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
			progress->write((i + CubeHDEF.get_size(0)) / double(3 * CubeHDEF.get_size(0)));
	}
	delete(progress);

	time_t end;
	time(&end);
	if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
	else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
	else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Hirshfeld(
	cube& CubeHDEF,
	cube& CubeRho,
	cube& CubeSpherical,
	WFN& wavy,
	int cpus,
	double radius,
	int ignore,
	ofstream &file
)
{
#ifdef _OPENMP
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif

	time_t start;
	time(&start);

	progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
	const int step = max(floor(3 * CubeHDEF.get_size(0) / 20), 1.0);

	Thakkar atom(wavy.get_atom_charge(ignore));

#pragma omp parallel for schedule(dynamic)
	for (int i = -CubeHDEF.get_size(0); i < 2 * CubeHDEF.get_size(0); i++) {
		for (int j = -CubeHDEF.get_size(1); j < 2 * CubeHDEF.get_size(1); j++)
			for (int k = -CubeHDEF.get_size(2); k < 2 * CubeHDEF.get_size(2); k++) {

				double PosGrid[3];

				PosGrid[0] = i * CubeHDEF.get_vector(0, 0) + j * CubeHDEF.get_vector(0, 1) + k * CubeHDEF.get_vector(0, 2) + CubeHDEF.get_origin(0);
				PosGrid[1] = i * CubeHDEF.get_vector(1, 0) + j * CubeHDEF.get_vector(1, 1) + k * CubeHDEF.get_vector(1, 2) + CubeHDEF.get_origin(1);
				PosGrid[2] = i * CubeHDEF.get_vector(2, 0) + j * CubeHDEF.get_vector(2, 1) + k * CubeHDEF.get_vector(2, 2) + CubeHDEF.get_origin(2);

				bool skip = true;
				double dist = sqrt(pow(PosGrid[0] - wavy.atoms[ignore].x, 2) + pow(PosGrid[1] - wavy.atoms[ignore].y, 2) + pow(PosGrid[2] - wavy.atoms[ignore].z, 2));
				if (dist < radius / 0.52)
					skip = false;
				if (skip)
					continue;

				double dens_choice = atom.get_radial_density(dist);
				int temp_i, temp_j, temp_k;
				if (i < 0)
					temp_i = i + CubeHDEF.get_size(0);
				else if (i < CubeHDEF.get_size(0))
					temp_i = i;
				else
					temp_i = i - CubeHDEF.get_size(0);

				if (j < 0)
					temp_j = j + CubeHDEF.get_size(1);
				else if (j < CubeHDEF.get_size(1))
					temp_j = j;
				else
					temp_j = j - CubeHDEF.get_size(1);

				if (k < 0)
					temp_k = k + CubeHDEF.get_size(2);
				else if (k < CubeHDEF.get_size(2))
					temp_k = k;
				else
					temp_k = k - CubeHDEF.get_size(2);
				CubeHDEF.set_value(temp_i, temp_j, temp_k, CubeHDEF.get_value(temp_i, temp_j, temp_k) + (dens_choice / CubeSpherical.get_value(temp_i,temp_j,temp_k) * CubeRho.get_value(temp_i,temp_j,temp_k)) - dens_choice);
			}
		if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
			progress->write((i + CubeHDEF.get_size(0)) / double(3 * CubeHDEF.get_size(0)));
	}
	delete(progress);

	time_t end;
	time(&end);
	if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
	else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
	else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Hirshfeld_atom(
	cube& CubeHirsh,
	cube& CubeRho,
	cube& CubeSpherical,
	WFN& wavy,
	int cpus,
	double radius,
	int ignore,
	ofstream& file
)
{
#ifdef _OPENMP
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif

	time_t start;
	time(&start);

	progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
	const int step = max(floor(3 * CubeHirsh.get_size(0) / 20), 1.0);

	Thakkar atom(wavy.get_atom_charge(ignore));

#pragma omp parallel for schedule(dynamic)
	for (int i = -CubeHirsh.get_size(0); i < 2 * CubeHirsh.get_size(0); i++) {
		for (int j = -CubeHirsh.get_size(1); j < 2 * CubeHirsh.get_size(1); j++)
			for (int k = -CubeHirsh.get_size(2); k < 2 * CubeHirsh.get_size(2); k++) {

				double PosGrid[3];

				PosGrid[0] = i * CubeHirsh.get_vector(0, 0) + j * CubeHirsh.get_vector(0, 1) + k * CubeHirsh.get_vector(0, 2) + CubeHirsh.get_origin(0);
				PosGrid[1] = i * CubeHirsh.get_vector(1, 0) + j * CubeHirsh.get_vector(1, 1) + k * CubeHirsh.get_vector(1, 2) + CubeHirsh.get_origin(1);
				PosGrid[2] = i * CubeHirsh.get_vector(2, 0) + j * CubeHirsh.get_vector(2, 1) + k * CubeHirsh.get_vector(2, 2) + CubeHirsh.get_origin(2);

				bool skip = true;
				double dist = sqrt(pow(PosGrid[0] - wavy.atoms[ignore].x, 2) + pow(PosGrid[1] - wavy.atoms[ignore].y, 2) + pow(PosGrid[2] - wavy.atoms[ignore].z, 2));
				if (dist < radius / 0.52)
					skip = false;
				if (skip)
					continue;

				double dens_choice = atom.get_radial_density(dist);
				int temp_i, temp_j, temp_k;
				if (i < 0)
					temp_i = i + CubeHirsh.get_size(0);
				else if (i < CubeHirsh.get_size(0))
					temp_i = i;
				else
					temp_i = i - CubeHirsh.get_size(0);

				if (j < 0)
					temp_j = j + CubeHirsh.get_size(1);
				else if (j < CubeHirsh.get_size(1))
					temp_j = j;
				else
					temp_j = j - CubeHirsh.get_size(1);

				if (k < 0)
					temp_k = k + CubeHirsh.get_size(2);
				else if (k < CubeHirsh.get_size(2))
					temp_k = k;
				else
					temp_k = k - CubeHirsh.get_size(2);
				CubeHirsh.set_value(temp_i, temp_j, temp_k, CubeHirsh.get_value(temp_i, temp_j, temp_k) + (dens_choice / CubeSpherical.get_value(temp_i, temp_j, temp_k) * CubeRho.get_value(temp_i, temp_j, temp_k)));
			}
		if (i != 0 && i % step == 0 && omp_get_thread_num() == 0)
			progress->write((i + CubeHirsh.get_size(0)) / double(3 * CubeHirsh.get_size(0)));
	}
	delete(progress);

	time_t end;
	time(&end);
	if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
	else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
	else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Rho(
	cube& CubeRho,
	WFN& wavy,
	int cpus,
	double radius,
	ofstream& file
)
{
#ifdef _OPENMP
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif

	time_t start;
	time(&start);

	progress_bar* progress = new progress_bar{ file, 50u, "Calculating Values" };
	const int step = max(floor(CubeRho.get_size(0) * 3 / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
	for (int i = -CubeRho.get_size(0); i < 2 * CubeRho.get_size(0); i++) {
		for (int j = -CubeRho.get_size(1); j < 2 * CubeRho.get_size(1); j++)
			for (int k = -CubeRho.get_size(2); k < 2 * CubeRho.get_size(2); k++) {

				double PosGrid[3],
					Rho[1];

				PosGrid[0] = i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0);
				PosGrid[1] = i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1);
				PosGrid[2] = i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2);

				bool skip = true;
				for (int a = 0; a < wavy.get_ncen(); a++)
					if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius / 0.52)
						skip = false;
				if (skip)
					continue;

				computeRho(
					PosGrid,
					wavy,
					Rho
				);

				int temp_i, temp_j, temp_k;
				if (i < 0)
					temp_i = i + CubeRho.get_size(0);
				else if (i < CubeRho.get_size(0))
					temp_i = i;
				else
					temp_i = i - CubeRho.get_size(0);

				if (j < 0)
					temp_j = j + CubeRho.get_size(1);
				else if (j < CubeRho.get_size(1))
					temp_j = j;
				else
					temp_j = j - CubeRho.get_size(1);

				if (k < 0)
					temp_k = k + CubeRho.get_size(2);
				else if (k < CubeRho.get_size(2))
					temp_k = k;
				else
					temp_k = k - CubeRho.get_size(2);

				CubeRho.set_value(temp_i, temp_j, temp_k, CubeRho.get_value(temp_i, temp_j, temp_k) + Rho[0]);
			}
		if (i != 0 && i % step == 0)
			progress->write((i + CubeRho.get_size(0)) / double(CubeRho.get_size(0) * 3));
	}
	delete(progress);

	time_t end;
	time(&end);
	if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
	else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
	else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_Prop(
	cube& CubeRho,
	cube& CubeRDG,
	cube& CubeElf,
	cube& CubeEli,
	cube& CubeLap,
	cube& CubeESP,
	WFN& wavy,
	int cpus,
	double radius,
	ofstream& file
	)
{
#ifdef _OPENMP
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif

	time_t start;
	time(&start);

	progress_bar * progress = new progress_bar{ file, 50u, "Calculating Values" };
	const int step = max(floor(CubeRho.get_size(0)*3 / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
	for (int i = -CubeRho.get_size(0); i < 2*CubeRho.get_size(0); i++) {
		for (int j = -CubeRho.get_size(1); j < 2*CubeRho.get_size(1); j++)
			for (int k = -CubeRho.get_size(2); k < 2*CubeRho.get_size(2); k++) {

				double PosGrid[3],
					Rho[1],
					Grad[1],
					Elf[1],
					Eli[1],
					Lap[1],
					Hess[9];

				PosGrid[0] = i * CubeRho.get_vector(0, 0) + j * CubeRho.get_vector(0, 1) + k * CubeRho.get_vector(0, 2) + CubeRho.get_origin(0);
				PosGrid[1] = i * CubeRho.get_vector(1, 0) + j * CubeRho.get_vector(1, 1) + k * CubeRho.get_vector(1, 2) + CubeRho.get_origin(1);
				PosGrid[2] = i * CubeRho.get_vector(2, 0) + j * CubeRho.get_vector(2, 1) + k * CubeRho.get_vector(2, 2) + CubeRho.get_origin(2);

				bool skip = true;
				for (int a = 0; a < wavy.get_ncen(); a++)
					if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius/0.52)
						skip = false;
				if (skip) 
					continue;

				if (CubeESP.get_loaded() && !CubeRDG.get_loaded())
					computeRho(
						PosGrid,
						wavy,
						Rho
					);
				
				if (CubeRDG.get_loaded() && CubeLap.get_loaded() && (CubeElf.get_loaded() || CubeEli.get_loaded()))
					computeValues(
						PosGrid,
						wavy,
						Rho,
						Grad,
						Hess,
						Elf,
						Eli,
						Lap
					);
				else if ((CubeElf.get_loaded() && CubeEli.get_loaded()) && !CubeRDG.get_loaded() && !CubeLap.get_loaded())
					computeELIELF(
						PosGrid,
						wavy,
						Elf,
						Eli
					);
				else if (CubeElf.get_loaded() && !CubeEli.get_loaded() && !CubeRDG.get_loaded() && !CubeLap.get_loaded())
					computeELF(
						PosGrid,
						wavy,
						Elf
					);
				else if (!CubeElf.get_loaded() && CubeEli.get_loaded() && !CubeRDG.get_loaded() && !CubeLap.get_loaded())
					computeELI(
						PosGrid,
						wavy,
						Eli
					);
				else if (CubeElf.get_loaded() && CubeEli.get_loaded() && CubeLap.get_loaded() && !CubeRDG.get_loaded())
					computeLapELIELF(
						PosGrid,
						wavy,
						Elf,
						Eli,
						Lap
					);
				else if (!CubeElf.get_loaded() && CubeEli.get_loaded() && CubeLap.get_loaded() && !CubeRDG.get_loaded())
					computeLapELI(
						PosGrid,
						wavy,
						Eli,
						Lap
					);
				else
					computeValues(
						PosGrid,
						wavy,
						Rho,
						Grad,
						Hess,
						Elf,
						Eli,
						Lap
					);

				if (CubeRDG.get_loaded())
					Rho[0] = get_lambda_1(Hess) < 0 ? -Rho[0] : Rho[0];

				int temp_i, temp_j, temp_k;
				if (i < 0)
					temp_i = i + CubeRho.get_size(0);
				else if (i < CubeRho.get_size(0))
					temp_i = i;
				else
					temp_i = i - CubeRho.get_size(0);

				if (j < 0)
					temp_j = j + CubeRho.get_size(1);
				else if (j < CubeRho.get_size(1))
					temp_j = j;
				else
					temp_j = j - CubeRho.get_size(1);

				if (k < 0)
					temp_k = k + CubeRho.get_size(2);
				else if (k < CubeRho.get_size(2))
					temp_k = k;
				else
					temp_k = k - CubeRho.get_size(2);

				CubeRho.set_value(temp_i, temp_j, temp_k, CubeRho.get_value(temp_i, temp_j, temp_k) + Rho[0]);
				if (CubeRDG.get_loaded()) {
					if (isnan(Grad[0])) Grad[0] = 0;
					if (isinf(Grad[0])) Grad[0] = 0;
					CubeRDG.set_value(temp_i, temp_j, temp_k, CubeRDG.get_value(temp_i, temp_j, temp_k) + Grad[0]);
				}
				if (CubeLap.get_loaded()) {
					if (isnan(Lap[0])) Lap[0] = 0;
					if (isinf(Lap[0])) Lap[0] = 0;
					CubeLap.set_value(temp_i, temp_j, temp_k, CubeLap.get_value(temp_i, temp_j, temp_k) + Lap[0]);
				}
				if (CubeElf.get_loaded()) {
					if (isnan(Elf[0])) Elf[0] = 0;
					if (isinf(Elf[0])) Elf[0] = 0;
					CubeElf.set_value(temp_i, temp_j, temp_k, CubeElf.get_value(temp_i, temp_j, temp_k) + Elf[0]);
				}
				if (CubeEli.get_loaded()) {
					if (isnan(Eli[0])) Eli[0] = 0;
					if (isinf(Eli[0])) Eli[0] = 0;
					CubeEli.set_value(temp_i, temp_j, temp_k, CubeEli.get_value(temp_i, temp_j, temp_k) + Eli[0]);
				}
			}
		if (i != 0 && i % step == 0)
			progress->write((i + CubeRho.get_size(0)) / double(CubeRho.get_size(0)*3));
	}
	delete(progress);

	time_t end;
	time(&end);
	if (difftime(end, start) < 60) file << "Time to calculate Values: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
	else if (difftime(end, start) < 3600) file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
	else file << "Time to calculate Values: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_ESP(
	cube& CubeESP,
	WFN& wavy,
	int cpus,
	double radius,
	ofstream& file
)
{
#ifdef _OPENMP
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif
	time_t start;
	time(&start);

	fill_pre();
	fill_Afac_pre();
	vector < vector <double> > d2;
	d2.resize(wavy.get_ncen());
	for (int i = 0; i < wavy.get_ncen(); i++) {
		d2[i].resize(wavy.get_ncen());
		for (int j = 0; j < wavy.get_ncen(); j++) {
			if (i == j) {
				d2[i][j] = 0;
				continue;
			}
			d2[i][j] = pow(wavy.atoms[i].x - wavy.atoms[j].x, 2) + pow(wavy.atoms[i].y - wavy.atoms[j].y, 2) + pow(wavy.atoms[i].z - wavy.atoms[j].z, 2);
		}
	}

	progress_bar * progress = new progress_bar{ file, 50u, "Calculating ESP" };
	const int step = max(floor(CubeESP.get_size(0) * 3 / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
	for (int i = -CubeESP.get_size(0); i < 2 * CubeESP.get_size(0); i++) {
		for (int j = -CubeESP.get_size(1); j < 2 * CubeESP.get_size(1); j++)
			for (int k = -CubeESP.get_size(2); k < 2 * CubeESP.get_size(2); k++) {

				double PosGrid[3];
				double ESP[1];

				PosGrid[0] = i * CubeESP.get_vector(0, 0) + j * CubeESP.get_vector(0, 1) + k * CubeESP.get_vector(0, 2) + CubeESP.get_origin(0);
				PosGrid[1] = i * CubeESP.get_vector(1, 0) + j * CubeESP.get_vector(1, 1) + k * CubeESP.get_vector(1, 2) + CubeESP.get_origin(1);
				PosGrid[2] = i * CubeESP.get_vector(2, 0) + j * CubeESP.get_vector(2, 1) + k * CubeESP.get_vector(2, 2) + CubeESP.get_origin(2);

				bool skip = true;
				for (int a = 0; a < wavy.get_ncen(); a++)
					if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius/0.52)
						skip = false;
				if (skip)
					continue;

				int temp_i, temp_j, temp_k;
				if (i < 0)
					temp_i = i + CubeESP.get_size(0);
				else if (i < CubeESP.get_size(0))
					temp_i = i;
				else
					temp_i = i - CubeESP.get_size(0);
				
				if (j < 0)
					temp_j = j + CubeESP.get_size(1);
				else if (j < CubeESP.get_size(1))
					temp_j = j;
				else
					temp_j = j - CubeESP.get_size(1);
				
				if (k < 0)
					temp_k = k + CubeESP.get_size(2);
				else if (k < CubeESP.get_size(2))
					temp_k = k;
				else
					temp_k = k - CubeESP.get_size(2);
				
				CubeESP.set_value(temp_i, temp_j, temp_k, CubeESP.get_value(temp_i, temp_j, temp_k) + computeESP(PosGrid, d2, wavy));
				//CubeESP.set_value(i, j, k, computeESP(PosGrid, d2, wavy));
			}
		if (i != 0 && i % step == 0)
			progress->write((i + CubeESP.get_size(0)) / double(CubeESP.get_size(0) * 3));
	}
	delete(progress);

	time_t end;
	time(&end);
	if (difftime(end, start) < 60) file << "Time to calculate ESP: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
	else if (difftime(end, start) < 3600) file << "Time to calculate ESP: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
	else file << "Time to calculate ESP: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

void Calc_MO(
	cube& CubeMO,
	int mo,
	WFN& wavy,
	int cpus,
	double radius,
	ofstream& file
)
{
#ifdef _OPENMP
	if (cpus != -1) {
		omp_set_num_threads(cpus);
		omp_set_dynamic(0);
		if (cpus > 1)
			omp_set_nested(1);
	}
#endif
	time_t start;
	time(&start);

	progress_bar* progress = new progress_bar{ file, 50u, "Calculating MO" };
	const int step = max(floor(CubeMO.get_size(0)*3 / 20), 1.0);

#pragma omp parallel for schedule(dynamic)
	for (int i = -CubeMO.get_size(0); i < 2 * CubeMO.get_size(0); i++) {
		for (int j = -CubeMO.get_size(1); j < 2 * CubeMO.get_size(1); j++)
			for (int k = -CubeMO.get_size(2); k < 2 * CubeMO.get_size(2); k++) {

				double PosGrid[3];
				double MO[1];

				PosGrid[0] = i * CubeMO.get_vector(0, 0) + j * CubeMO.get_vector(0, 1) + k * CubeMO.get_vector(0, 2) + CubeMO.get_origin(0);
				PosGrid[1] = i * CubeMO.get_vector(1, 0) + j * CubeMO.get_vector(1, 1) + k * CubeMO.get_vector(1, 2) + CubeMO.get_origin(1);
				PosGrid[2] = i * CubeMO.get_vector(2, 0) + j * CubeMO.get_vector(2, 1) + k * CubeMO.get_vector(2, 2) + CubeMO.get_origin(2);

				bool skip = true;
				for (int a = 0; a < wavy.get_ncen(); a++)
					if (sqrt(pow(PosGrid[0] - wavy.atoms[a].x, 2) + pow(PosGrid[1] - wavy.atoms[a].y, 2) + pow(PosGrid[2] - wavy.atoms[a].z, 2)) < radius/0.52)
						skip = false;
				if (skip)
					continue;

				computeMO(
					PosGrid,
					wavy,
					MO,
					mo
				);

				int temp_i, temp_j, temp_k;
				if (i < 0)
					temp_i = i + CubeMO.get_size(0);
				else if (i < CubeMO.get_size(0))
					temp_i = i;
				else
					temp_i = i - CubeMO.get_size(0);

				if (j < 0)
					temp_j = j + CubeMO.get_size(1);
				else if (j < CubeMO.get_size(1))
					temp_j = j;
				else
					temp_j = j - CubeMO.get_size(1);

				if (k < 0)
					temp_k = k + CubeMO.get_size(2);
				else if (k < CubeMO.get_size(2))
					temp_k = k;
				else
					temp_k = k - CubeMO.get_size(2);

				CubeMO.set_value(temp_i, temp_j, temp_k, CubeMO.get_value(temp_i, temp_j, temp_k) + MO[0]);
			}
		if (i != 0 && i % step == 0)
			progress->write((i + CubeMO.get_size(0)) / double(CubeMO.get_size(0) * 3));
	}
	delete(progress);

	time_t end;
	time(&end);
	if (difftime(end, start) < 60) file << "Time to calculate MO: " << fixed << setprecision(0) << difftime(end, start) << " s" << endl;
	else if (difftime(end, start) < 3600) file << "Time to calculate MO: " << fixed << setprecision(0) << floor(difftime(end, start) / 60) << " m " << int(floor(difftime(end, start))) % 60 << " s" << endl;
	else file << "Time to calculate MO: " << fixed << setprecision(0) << floor(difftime(end, start) / 3600) << " h " << (int(floor(difftime(end, start))) % 3600) / 60 << " m" << endl;
};

//end here
