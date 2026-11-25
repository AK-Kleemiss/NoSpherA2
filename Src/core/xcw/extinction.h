#pragma once
#include <array>
#include <cmath>
#include <string>

//Extinction corrections in cctbx's convention (cctbx/xray/extinction.h): a model returns a multiplier y on
//Fc^2, so sqrt(y)*|Fc| is compared with |F_obs|. Every model multiplies SHELX's geometry constant
//0.001*lambda^3/sin(2 theta), so the coefficient stays on SHELXL's EXTI scale with the physical constants
//(cell volume, polarisation ratio, mean path length) absorbed in it.
namespace extinction {

	enum class model { none = 0, shelx, bc_gaussian, bc_lorentzian };

	//0.001 * lambda^3 / sin(2 theta), stl = sin(theta)/lambda in Angstrom^-1; zero where the reflection is not
	//accessible at this wavelength
	inline double geometry_constant(const double wavelength, const double stl) {
		const double sin_t = wavelength * stl;
		if (sin_t <= 0.0 || sin_t >= 1.0) return 0.0;
		const double sin_2t = 2.0 * sin_t * std::sqrt(1.0 - sin_t * sin_t);
		return 0.001 * wavelength * wavelength * wavelength / sin_2t;
	}

	//cos(2 theta), which the Becker-Coppens mosaic coefficients depend on
	inline double cos_2theta(const double wavelength, const double stl) {
		const double sin_t = wavelength * stl;
		return 1.0 - 2.0 * sin_t * sin_t;
	}

	//y(t), and dy/dt when dydt is given; t = c * x * |Fc|^2, c the geometry constant, x the coefficient (for
	//the anisotropic models the quadratic form of the tensor).
	//SHELX/Zachariasen: y = (1 + t)^(-1/2), the leading term of Becker-Coppens type I.
	//Becker & Coppens (1974): y = (1 + 2t + A t^2/(1 + B t))^(-1/2), Gaussian and Lorentzian mosaic A and B
	//as in GSAS-II's GSASIIstrMath.SCExtinction.
	//D <= 0 (only for a negative coefficient) is left uncorrected rather than NaN; a backstop, as x >= 0 is kept.
	inline double correction(const model m, const double cos2t, const double t, double* dydt = nullptr) {
		if (m == model::none) {
			if (dydt) *dydt = 0.0;
			return 1.0;
		}
		double D = 0.0, dD = 0.0;
		if (m == model::shelx) {
			D = 1.0 + t;
			dD = 1.0;
		}
		else {
			const bool gaussian = (m == model::bc_gaussian);
			const double A = gaussian ? 0.58 + 0.48 * cos2t + 0.24 * cos2t * cos2t
				: 0.025 + 0.285 * cos2t;
			const double B = gaussian ? 0.02 - 0.025 * cos2t
				: (cos2t < 0.0 ? -0.45 * cos2t : 0.15 - 0.2 * (0.75 - cos2t) * (0.75 - cos2t));
			const double q = 1.0 + B * t;
			if (q <= 0.0) {
				if (dydt) *dydt = 0.0;
				return 1.0;
			}
			D = 1.0 + 2.0 * t + A * t * t / q;
			dD = 2.0 + 2.0 * A * t / q - A * t * t * B / (q * q);
		}
		if (D <= 0.0) {
			if (dydt) *dydt = 0.0;
			return 1.0;
		}
		const double y = 1.0 / std::sqrt(D);
		if (dydt) *dydt = -0.5 * y * dD / D;
		return y;
	}

	//x(h) = sum_p a_p X_p, X in Voigt order X11 X22 X33 X12 X13 X23, h_unit the normalised Cartesian
	//scattering vector. The rigorous models (Coppens & Hamilton 1970, Thornley & Nelmes 1974; XD2006 sec. 4.6.8,
	//g(D) = (D'ZD)^1/2 with D perpendicular to the diffraction plane) need the incident and diffracted beam
	//directions, which an hkl file does not carry. Averaging over the azimuth psi around h removes them:
	//<s_i s_j>_psi = (delta_ij - h_i h_j)/2 for unit s perpendicular to h, so <D'XD>_psi = (tr X - h'Xh)/2,
	//which is x for X = x*I (the isotropic limit).
	inline void aniso_coefficients(const std::array<double, 3>& h_unit, std::array<double, 6>& a) {
		for (int i = 0; i < 3; i++) a[i] = 0.5 * (1.0 - h_unit[i] * h_unit[i]);
		a[3] = -h_unit[0] * h_unit[1];
		a[4] = -h_unit[0] * h_unit[2];
		a[5] = -h_unit[1] * h_unit[2];
	}

	inline const char* name(const model m) {
		switch (m) {
		case model::shelx: return "SHELX (Zachariasen)";
		case model::bc_gaussian: return "Becker-Coppens, Gaussian mosaic";
		case model::bc_lorentzian: return "Becker-Coppens, Lorentzian mosaic";
		default: return "none";
		}
	}

	//The settings-file spelling; model::none for anything unrecognised
	inline model from_string(const std::string& s) {
		if (s == "shelx" || s == "zachariasen") return model::shelx;
		if (s == "bc_gaussian" || s == "becker_coppens_gaussian" || s == "gaussian") return model::bc_gaussian;
		if (s == "bc_lorentzian" || s == "becker_coppens_lorentzian" || s == "lorentzian") return model::bc_lorentzian;
		return model::none;
	}
}
