#pragma once
#include <array>
#include <cmath>
#include <string>

//Extinction corrections on the model F_calc^2, in cctbx's convention (see
//cctbx/xray/extinction.h): a model returns a multiplier y on Fc^2, so the amplitude that is
//compared with |F_obs| is sqrt(y)*|Fc|. Every model here multiplies the same geometry
//constant 0.001*lambda^3/sin(2 theta) SHELX uses, so the refined coefficient stays on
//SHELX's scale and the physical constants (cell volume, polarisation ratio, mean path
//length) stay absorbed in it, exactly as SHELXL's EXTI does.
namespace extinction {

	enum class model { none = 0, shelx, bc_gaussian, bc_lorentzian };

	//0.001 * lambda^3 / sin(2 theta), with sin(theta) = lambda * stl for stl = sin(theta)/lambda
	//in Angstrom^-1. Zero where the reflection is not accessible at this wavelength, which is
	//where the correction does not apply.
	inline double geometry_constant(const double wavelength, const double stl) {
		const double sin_t = wavelength * stl;
		if (sin_t <= 0.0 || sin_t >= 1.0) return 0.0;
		const double sin_2t = 2.0 * sin_t * std::sqrt(1.0 - sin_t * sin_t);
		return 0.001 * wavelength * wavelength * wavelength / sin_2t;
	}

	//cos(2 theta), which is what the Becker-Coppens mosaic coefficients depend on
	inline double cos_2theta(const double wavelength, const double stl) {
		const double sin_t = wavelength * stl;
		return 1.0 - 2.0 * sin_t * sin_t;
	}

	//y(t), and dy/dt when dydt is given. t = c * x * |Fc|^2 with c the geometry constant and
	//x the refined coefficient (for the anisotropic models, the quadratic form of the tensor).
	//  SHELX/Zachariasen:  y = (1 + t)^(-1/2), the leading term of Becker-Coppens type I.
	//  Becker-Coppens:     y = (1 + 2t + A t^2/(1 + B t))^(-1/2), with the Gaussian and
	//                      Lorentzian mosaic coefficients of Becker & Coppens (1974), taken
	//                      from GSAS-II's GSASIIstrMath.SCExtinction.
	//A non-physical shape (D <= 0, reachable only for a negative coefficient) is left
	//uncorrected rather than producing a NaN; refine_extinction keeps x >= 0 so this is a
	//backstop, not a working branch.
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

	//The anisotropic coefficient x(h) = sum_p a_p X_p, with X stored in the Voigt order
	//X11 X22 X33 X12 X13 X23 and h_unit the normalised scattering vector in Cartesian.
	//
	//The rigorous anisotropic models (Coppens & Hamilton 1970, Thornley & Nelmes 1974; XD2006
	//sec. 4.6.8 writes them as g(D) = (D'ZD)^1/2 with D perpendicular to the diffraction plane
	//and rho(N) = lambda (N'WN)^-1/2 with N in it) need the direction cosines of the incident
	//and diffracted beams, which XD demands as six extra entries per observation and which an
	//hkl file does not carry. Averaging over the azimuth psi around the scattering vector
	//removes them: <s_i s_j>_psi = (delta_ij - h_i h_j)/2 for any unit vector s perpendicular
	//to h, so <D'XD>_psi = (tr X - h'Xh)/2. For X = x*I that is exactly x, so the isotropic
	//model is this one's isotropic limit, and for merged data averaged over several azimuths
	//it is arguably the quantity the measurement actually saw.
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
