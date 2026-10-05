#include "pch.h"
#include "tuning.h"
#include "wfn_class.h"
#include "convenience.h"
#include "bondwise_analysis.h"
#include "properties.h"
#include "libCintMain.h"
#include "nos_math.h"
#include "integration_params.h"
#include "b2c.h"
#include "eli_family.h"
#include "topology.h"
#include "crystal_energies.h"
#include "spherical_density.h"
#include "citations.h"
#include "nao.h"
#include "section_log.h"
#include <algorithm>
#include <cctype>
#include <chrono>
#include <mutex>
#include <occ/core/parallel.h>
#include <occ/qm/hf.h>
#include <occ/qm/guess_kind.h>
#include <occ/qm/initial_guess.h>
#include <occ/qm/cint_interface.h>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <thread>

namespace {
	struct OhOperation {
		std::array<int, 3> permutation;
		std::array<int, 3> signs;
	};

	struct CartesianTransform {
		ivec destination;
		ivec signs;
	};

	int cartesian_shell_size(const int l) {
		return (l + 1) * (l + 2) / 2;
	}

	//Dimension of the free atom's occupied space in Aufbau order with full m degeneracy (Li 2, B-Ne 5, Fe 15).
	//An element property, so the projector rank and the index stay continuous in geometry; a partly filled shell counts in full.
	constexpr int free_atom_orbital_count(const int atomic_number) {
		//l of each shell in Aufbau order, 1s 2s 2p ... 7s 5f 6d 7p, up to Z = 118.
		//Flat array with an index loop: MSVC rejects a range-for over a constexpr int[][2] here.
		constexpr int shell_l[] = { 0, 0, 1, 0, 1, 0, 2, 1, 0, 2, 1, 0, 3, 2, 1, 0, 3, 2, 1 };
		int electrons = atomic_number;
		int dimension = 0;
		for (int i = 0; i < static_cast<int>(sizeof(shell_l) / sizeof(shell_l[0])); i++) {
			if (electrons <= 0) break;
			const int size = 2 * shell_l[i] + 1;
			dimension += size;
			electrons -= 2 * size;
		}
		return dimension;
	}

	//Orbitals an ECP core removes: complete shells in (n, then l) order until its electrons are used up
	//(60 e on Au -> 30 of 40), passing over one they cannot fill (4f under a 54-electron [Xe] core).
	//A lanthanide ECP between 46 and 60 electrons is 4f-in-core: [Kr]4d10 and 4f.
	//Without it the rank would clamp to the whole atomic block.
	constexpr int ecp_core_orbital_count(const int atomic_number, const int ecp_electrons) {
		if (atomic_number >= 57 && atomic_number <= 71 && ecp_electrons > 46 && ecp_electrons < 60)
			return ecp_core_orbital_count(atomic_number, 46) + 7;
		int electrons = ecp_electrons;
		int dimension = 0;
		for (int n = 1; electrons > 0 && n <= 7; n++)
			for (int l = 0; l < n && electrons > 0; l++) {
				const int size = 2 * l + 1;
				if (electrons < 2 * size) continue;
				dimension += size;
				electrons -= 2 * size;
			}
		return dimension;
	}

	static_assert(free_atom_orbital_count(3) == 2, "Li spans 1s and 2s");
	static_assert(free_atom_orbital_count(7) == 5, "N spans 1s, 2s and 2p - a half-filled shell in full");
	static_assert(free_atom_orbital_count(26) == 15, "Fe spans 1s..4s and 3d");
	static_assert(free_atom_orbital_count(118) == 59, "the Aufbau list reaches the last element");
	static_assert(ecp_core_orbital_count(6, 0) == 0, "no ECP removes nothing");
	static_assert(ecp_core_orbital_count(53, 28) == 14, "a 28-electron ECP covers 1s..3d");
	static_assert(free_atom_orbital_count(53) - ecp_core_orbital_count(53, 28) == 13,
		"iodine with a 28-electron ECP keeps 4s, 4p, 4d, 5s and 5p");
	static_assert(free_atom_orbital_count(79) - ecp_core_orbital_count(79, 60) == 10,
		"gold with a 60-electron ECP keeps 5s, 5p, 5d and 6s");
	static_assert(free_atom_orbital_count(55) - ecp_core_orbital_count(55, 54) == 1,
		"caesium with a 54-electron [Xe] ECP keeps 6s");
	static_assert(free_atom_orbital_count(65) - ecp_core_orbital_count(65, 54) == 5,
		"terbium with a 4f-in-core 54-electron ECP keeps 5s, 5p and 6s");

	int atomic_shell_size(const int l, const bool cartesian) {
		return cartesian ? cartesian_shell_size(l) : 2 * l + 1;
	}

	int tonto_period_number(const int atomic_number) {
		if (atomic_number < 1)
			return 0;

		int period = 1;
		int noble = 0;
		while (true) {
			const int n = (period + 2) / 2;
			noble += 2 * n * n;
			if (atomic_number <= noble)
				return period;
			period++;
		}
	}

	int tonto_column_number(const int atomic_number) {
		if (atomic_number < 1)
			return 0;

		int period = 1;
		int noble = 0;
		while (true) {
			const int n = (period + 2) / 2;
			noble += 2 * n * n;
			if (atomic_number <= noble)
				return atomic_number - (noble - 2 * n * n);
			period++;
		}
	}

	int tonto_ground_state_multiplicity(const int atomic_number) {
		const int period = tonto_period_number(atomic_number);
		const int column = tonto_column_number(atomic_number);
		switch (period) {
		case 0:
			return 1;
		case 1:
			//The first period has two columns, so He is a singlet; a 1s(1)2s(1) triplet has degenerate
			//natural orbitals and the rank-1 subspace would cut through them arbitrarily.
			return column >= 2 ? 1 : 2;
		case 2:
		case 3:
			switch (column) {
			case 8: return 1;
			case 1:
			case 3:
			case 7: return 2;
			case 2:
			case 4:
			case 6: return 3;
			case 5: return 4;
			default: return 1;
			}
		case 4:
		case 5:
			switch (column) {
			case 2:
			case 12:
			case 18: return 1;
			case 1:
			case 3:
			case 11:
			case 13:
			case 17: return 2;
			case 4:
			case 10:
			case 14:
			case 16: return 3;
			case 5:
			case 9:
			case 15: return 4;
			case 6:
			case 8: return 5;
			case 7: return 6;
			default: return 1;
			}
		case 6:
		case 7:
			switch (column) {
			case 2:
			case 16:
			case 26:
			case 32: return 1;
			case 1:
			case 3:
			case 15:
			case 17:
			case 25:
			case 27:
			case 31: return 2;
			case 4:
			case 14:
			case 18:
			case 24:
			case 28:
			case 30: return 3;
			case 5:
			case 13:
			case 19:
			case 23:
			case 29: return 4;
			case 6:
			case 12:
			case 20:
			case 22: return 5;
			case 7:
			case 11:
			case 21: return 6;
			case 8:
			case 10: return 7;
			case 9: return 8;
			default: return 1;
			}
		default:
			return 1;
		}
	}

	occ::gto::AOBasis build_occ_atomic_basis_from_wfn_atom(
		const atom &atm, const e_origin origin, const bool cartesian) {
		//The ECP core comes off the nucleus (Z - N_core, neutral pseudo-atom): occ gets no ECP potential here,
		//so declaring frozen core electrons would build a highly charged cation.
		const int effective_Z = atm.get_charge() - atm.get_ECP_electrons();
		std::vector<occ::core::Atom> occ_atoms{ { effective_Z, 0.0, 0.0, 0.0 } };
		std::vector<occ::gto::Shell> shells;
		const auto basis_set = atm.get_basis_set();
		int primitive_idx = 0;

		while (primitive_idx < static_cast<int>(basis_set.size())) {
			const int shell_id = static_cast<int>(basis_set[primitive_idx].get_shell());
			const int l = static_cast<int>(basis_set[primitive_idx].get_type()) - 1;
			err_checkf(l >= 0,
				"Encountered an invalid shell angular momentum while building an OCC atomic basis.",
				std::cout);
			vec exponents;
			vec coefficients;
			while (primitive_idx < static_cast<int>(basis_set.size()) &&
				static_cast<int>(basis_set[primitive_idx].get_shell()) == shell_id) {
				exponents.push_back(basis_set[primitive_idx].get_exponent());
				coefficients.push_back(basis_set[primitive_idx].get_coefficient());
				primitive_idx++;
			}

			occ::gto::Shell shell(l, exponents, { coefficients }, { 0.0, 0.0, 0.0 });
			shell.kind = cartesian ? occ::gto::Shell::Kind::Cartesian : occ::gto::Shell::Kind::Spherical;
			shells.push_back(shell);
		}

		occ::gto::AOBasis result(occ_atoms, shells, "wfn-atomic-basis");
		result.set_pure(!cartesian);
		//No set_ecp_electrons: the core is already off the nucleus; occ would freeze it twice.
		return result;
	}

	dMatrix2 eigen_matrix_to_dmatrix2(const occ::Mat &matrix) {
		dMatrix2 result(matrix.rows(), matrix.cols());
		for (int row = 0; row < matrix.rows(); ++row)
			for (int col = 0; col < matrix.cols(); ++col)
				result(row, col) = matrix(row, col);
		return result;
	}

	class ScopedOccLogLevel {
	public:
		explicit ScopedOccLogLevel(const spdlog::level::level_enum level) {
			auto logger = spdlog::default_logger();
			if (logger != nullptr) {
				had_logger = true;
				previous_level = logger->level();
				logger->set_level(level);
			}
		}

		~ScopedOccLogLevel() {
			if (!had_logger)
				return;
			if (auto logger = spdlog::default_logger(); logger != nullptr)
				logger->set_level(previous_level);
		}

	private:
		bool had_logger = false;
		spdlog::level::level_enum previous_level = spdlog::level::info;
	};

	//occ::parallel::nthreads is bookkeeping only; TBB is unbounded until set_num_threads() installs a
	//tbb::global_control. Open-shell free atoms have degenerate manifolds, so a nondeterministic Fock reduction
	//changes the density, hence one thread. NOS_RGBI_NO_PIN disables the pin to measure its cost; read once
	//so constructor and destructor agree.
	inline bool occ_pinning_disabled() {
		static const bool disabled = tuning("NOS_RGBI_NO_PIN") != nullptr;
		return disabled;
	}

	//Pinned means a global_control exists, not get_num_threads() == 1 (the default before any control).
	//Restore: with no control on entry, shut TBB down rather than install one at 1, or occ stays serial.
	//One global depth count: the outermost frame owns pin and restore, so the warm pass's worker frames
	//(on other threads than the enclosing one) never tear down a pool in use.
	class ScopedOccSingleThread {
	public:
		ScopedOccSingleThread() {
			if (occ_pinning_disabled())
				return;
			const std::lock_guard<std::mutex> hold(state().mutex);
			engaged = true;
			if (state().depth++ > 0)
				return;
			state().had_control = occ::parallel::get_tbb_control() != nullptr;
			state().previous = occ::parallel::get_num_threads();
			occ::parallel::set_num_threads(1);
		}
		~ScopedOccSingleThread() {
			if (!engaged)
				return;
			const std::lock_guard<std::mutex> hold(state().mutex);
			if (--state().depth > 0)
				return;
			if (state().had_control)
				occ::parallel::set_num_threads(state().previous);
			else {
				occ::parallel::shutdown_tbb();
				occ::parallel::nthreads = state().previous;
			}
		}

	private:
		bool engaged = false;
		struct State {
			std::mutex mutex;
			int depth = 0;
			bool had_control = false;
			int previous = 1;
		};
		static State &state() {
			static State s;
			return s;
		}
	};

	//Each diagnostic line is built in a local stream (own format state) and written under one lock: the
	//warm pass prints from several OpenMP threads. Callers read NOS_RGBI_DEBUG per call; tests flip it in-process.
	void rgbi_debug_line(const std::string &line) {
		static std::mutex print_mutex;
		const std::lock_guard<std::mutex> hold(print_mutex);
		std::cout << "\n" << line << std::endl;
	}

	//A free atom's density depends only on its contracted basis, electron count and conventions, not on
	//its position, so one SCF serves every atom with the same key. Compare the contraction field by field:
	//basis_set_entry::operator== also compares the owning atom index and would never match.
	struct FreeAtomKey {
		std::vector<basis_set_entry> basis;
		int charge;
		int ecp_electrons;
		e_origin origin;
		bool cartesian;
		bool operator==(const FreeAtomKey &other) const {
			if (charge != other.charge || ecp_electrons != other.ecp_electrons ||
				origin != other.origin || cartesian != other.cartesian ||
				basis.size() != other.basis.size())
				return false;
			for (size_t i = 0; i < basis.size(); i++)
				if (basis[i].get_type() != other.basis[i].get_type() ||
					basis[i].get_shell() != other.basis[i].get_shell() ||
					basis[i].get_exponent() != other.basis[i].get_exponent() ||
					basis[i].get_coefficient() != other.basis[i].get_coefficient())
					return false;
			return true;
		}
	};

	//File scope for clear_rgbi_free_atom_cache(); a linear scan suffices for a few distinct elements.
	//Locked because the warm pass fills it concurrently.
	std::vector<std::pair<FreeAtomKey, dMatrix2>> free_atom_cache;
	std::mutex free_atom_cache_mutex;

	//sf-X2C-1e scalar-relativistic correction, returned as the change to T + V in the atom's contracted basis
	//for occ's SCF as an external potential (PySCF's sfx2c1e: modified Dirac equation in the decontracted
	//basis, X from the electronic solutions, renormalise with R, contract back; p.Vp from libcint).
	//Primitives carry coefficient 1 as occ passes them, so the contraction matrix is the coefficient table;
	//S, T and V are rebuilt from it and must reproduce occ's to 1e-8.
	occ::Mat sfx2c1e_core_correction(const occ::gto::AOBasis &basis, const occ::Mat &S_ref,
		const occ::Mat &T_ref, const occ::Mat &V_ref) {
		using occ::gto::Shell;
		constexpr double c = 137.03599967994;
		const bool cart = basis.is_cartesian();

		//Unique (l, exponent) primitives. A generally contracted basis repeats its exponents across
		//shells, and duplicated primitives would make the metric singular.
		std::vector<Shell> prims;
		std::vector<std::vector<int>> prim_of(basis.size());
		for (size_t k = 0; k < basis.size(); k++) {
			const Shell &sh = basis[k];
			err_checkf(sh.num_contractions() == 1, "sf-X2C: shell with more than one contraction", std::cout);
			for (size_t p = 0; p < sh.num_primitives(); p++) {
				const double a = sh.exponents(p);
				int found = -1;
				for (int u = 0; u < static_cast<int>(prims.size()); u++)
					if (prims[u].l == sh.l && std::abs(prims[u].exponents(0) - a) <= 1e-12 * a) {
						found = u;
						break;
					}
				if (found < 0) {
					found = static_cast<int>(prims.size());
					Shell d(static_cast<int>(sh.l), { a }, { { 1.0 } }, { 0.0, 0.0, 0.0 });
					d.kind = sh.kind;
					prims.push_back(d);
				}
				prim_of[k].push_back(found);
			}
		}

		occ::qm::cint::IntegralEnvironment env(basis.atoms(), prims);
		std::vector<int> off(prims.size() + 1, 0);
		for (size_t i = 0; i < prims.size(); i++)
			off[i + 1] = off[i] + static_cast<int>(prims[i].size());
		const int n = off.back();
		auto one_e = [&](libcint::CINTIntegralFunction *fn) {
			occ::Mat M = occ::Mat::Zero(n, n);
			vec buf(64 * 64);
			for (int i = 0; i < static_cast<int>(prims.size()); i++)
				for (int j = 0; j < static_cast<int>(prims.size()); j++) {
					int shls[2] = { i, j };
					fn(buf.data(), nullptr, shls, env.atom_data_ptr(), env.num_atoms(), env.basis_data_ptr(),
						env.num_basis(), env.env_data_ptr(), nullptr, nullptr);
					const int di = off[i + 1] - off[i], dj = off[j + 1] - off[j];
					for (int b = 0; b < dj; b++)
						for (int a = 0; a < di; a++)
							M(off[i] + a, off[j] + b) = buf[a + di * b];
				}
			return M;
		};
		occ::Mat s = one_e(cart ? &libcint::int1e_ovlp_cart : &libcint::int1e_ovlp_sph);
		occ::Mat t = one_e(cart ? &libcint::int1e_kin_cart : &libcint::int1e_kin_sph);
		occ::Mat v = one_e(cart ? &libcint::int1e_nuc_cart : &libcint::int1e_nuc_sph);
		occ::Mat w = one_e(cart ? &libcint::int1e_pnucp_cart : &libcint::int1e_pnucp_sph);

		//Normalise the primitives, for the conditioning of the eigenproblem below; K absorbs the scale.
		const occ::Vec norm = s.diagonal().cwiseSqrt();
		const occ::Vec inv = norm.cwiseInverse();
		s = inv.asDiagonal() * s * inv.asDiagonal();
		t = inv.asDiagonal() * t * inv.asDiagonal();
		v = inv.asDiagonal() * v * inv.asDiagonal();
		w = inv.asDiagonal() * w * inv.asDiagonal();
		occ::Mat K = occ::Mat::Zero(n, static_cast<Eigen::Index>(basis.nbf()));
		for (size_t k = 0; k < basis.size(); k++) {
			const int first = static_cast<int>(basis.first_bf()[k]);
			const int width = static_cast<int>(basis[k].size());
			for (size_t p = 0; p < basis[k].num_primitives(); p++) {
				const int u = prim_of[k][p];
				for (int m = 0; m < width; m++)
					K(off[u] + m, first + m) += basis[k].contraction_coefficients(p, 0) * norm(off[u] + m);
			}
		}
		auto rel_dev = [](const occ::Mat &a, const occ::Mat &b) {
			return (a - b).cwiseAbs().maxCoeff() / std::max(1.0, b.cwiseAbs().maxCoeff());
		};
		const double dev = std::max({ rel_dev(K.transpose() * s * K, S_ref),
			rel_dev(K.transpose() * t * K, T_ref), rel_dev(K.transpose() * v * K, V_ref) });
		if (!(dev < 1e-8))
			throw std::runtime_error("sf-X2C: the decontracted basis does not reproduce occ's contracted "
				"integrals (max relative deviation " + std::to_string(dev) + ")");

		const int n2 = 2 * n;
		occ::Mat h4 = occ::Mat::Zero(n2, n2), m4 = occ::Mat::Zero(n2, n2);
		h4.topLeftCorner(n, n) = v;
		h4.topRightCorner(n, n) = t;
		h4.bottomLeftCorner(n, n) = t;
		h4.bottomRightCorner(n, n) = w * (0.25 / (c * c)) - t;
		m4.topLeftCorner(n, n) = s;
		m4.bottomRightCorner(n, n) = t * (0.5 / (c * c));
		Eigen::GeneralizedSelfAdjointEigenSolver<occ::Mat> dirac(h4, m4);
		if (dirac.info() != Eigen::Success)
			throw std::runtime_error("sf-X2C: the modified Dirac equation could not be solved");
		//Eigenvalues ascend, so the upper n are the electronic solutions.
		const occ::Mat cl = dirac.eigenvectors().block(0, n, n, n);
		const occ::Mat cs = dirac.eigenvectors().block(n, n, n, n);
		const occ::Mat X = cl.transpose().partialPivLu().solve(cs.transpose()).transpose();

		const occ::Mat s1 = s + X.transpose() * t * X * (0.5 / (c * c));
		const occ::Mat tx = t * X;
		const occ::Mat h1 = v + tx + tx.transpose() - X.transpose() * tx
			+ X.transpose() * w * X * (0.25 / (c * c));
		//R = S^-1/2 (S^-1/2 s1 S^-1/2)^-1/2 S^1/2, in the eigenbasis of S.
		Eigen::SelfAdjointEigenSolver<occ::Mat> es(s);
		std::vector<int> keep;
		for (int i = 0; i < n; i++)
			if (es.eigenvalues()(i) > 1e-14)
				keep.push_back(i);
		occ::Mat U(n, keep.size());
		occ::Vec ws(keep.size());
		for (size_t i = 0; i < keep.size(); i++) {
			U.col(i) = es.eigenvectors().col(keep[i]);
			ws(i) = std::sqrt(es.eigenvalues()(keep[i]));
		}
		const occ::Mat mid = ws.cwiseInverse().asDiagonal() * (U.transpose() * s1 * U) * ws.cwiseInverse().asDiagonal();
		Eigen::SelfAdjointEigenSolver<occ::Mat> em(mid);
		occ::Mat mid_isqrt = occ::Mat::Zero(mid.rows(), mid.cols());
		for (int i = 0; i < mid.rows(); i++)
			if (em.eigenvalues()(i) > 1e-14)
				mid_isqrt += em.eigenvectors().col(i) * em.eigenvectors().col(i).transpose() / std::sqrt(em.eigenvalues()(i));
		const occ::Mat R = U * (ws.cwiseInverse().asDiagonal() * mid_isqrt * ws.asDiagonal()) * U.transpose();
		const occ::Mat h_x2c = R.transpose() * h1 * R;
		occ::Mat delta = K.transpose() * (h_x2c - t - v) * K;
		return 0.5 * (delta + delta.transpose());
	}

	//Free-atom densities on disk, keyed on the in-memory key plus Hamiltonian, build_date, pin and iteration
	//cap. The full key is stored and compared on load, so a hash collision or unreadable file is a miss.
	struct FreeAtomDiskEntry {
		dMatrix2 density;
		bool converged = true;
		int iter = 0;
		double ediff_rel = 0.0, diis_error = 0.0;
	};

	std::filesystem::path free_atom_cache_dir() {
		//"off" disables the disk cache alone; NOS_RGBI_NO_FREEATOM_CACHE disables both caches.
		if (const char *d = tuning("NOS_FREEATOM_CACHE_DIR"))
			return (std::string(d).empty() || std::string(d) == "off") ? std::filesystem::path() : std::filesystem::path(d);
		if (const char *d = std::getenv("LOCALAPPDATA")) // Flawfinder: ignore - the user's cache root, used as a directory only
			return std::filesystem::path(d) / "NoSpherA2" / "free_atom_cache";
		if (const char *d = std::getenv("XDG_CACHE_HOME")) // Flawfinder: ignore - the user's cache root, used as a directory only
			return std::filesystem::path(d) / "NoSpherA2" / "free_atom_cache";
		if (const char *d = std::getenv("HOME")) // Flawfinder: ignore - the user's cache root, used as a directory only
			return std::filesystem::path(d) / ".cache" / "NoSpherA2" / "free_atom_cache";
		return {};
	}

	std::string free_atom_disk_key(const FreeAtomKey &key, const bool relativistic) {
		std::string k = relativistic ? "sfx2c1e|" : "nonrel|";
		k += build_date + "|pin=" + (occ_pinning_disabled() ? "0" : "1") + "|maxiter=";
		if (const char *cap = tuning("NOS_RGBI_FREE_ATOM_MAXITER"))
			k += cap;
		auto put = [&k](const auto value) { k.append(reinterpret_cast<const char *>(&value), sizeof(value)); };
		put(key.charge); put(key.ecp_electrons); put(static_cast<int>(key.origin)); put(key.cartesian);
		for (const auto &b : key.basis) {
			put(static_cast<int>(b.get_type())); put(static_cast<int>(b.get_shell()));
			put(b.get_exponent()); put(b.get_coefficient());
		}
		return k;
	}

	std::filesystem::path free_atom_disk_path(const std::string &disk_key, const int Z) {
		const auto dir = free_atom_cache_dir();
		if (dir.empty())
			return {};
		uint64_t h = 1469598103934665603ULL;
		for (const unsigned char ch : disk_key)
			h = (h ^ ch) * 1099511628211ULL;
		std::ostringstream name;
		name << "Z" << Z << "_" << std::hex << h << ".fad";
		return dir / name.str();
	}

	constexpr char free_atom_disk_magic[8] = { 'N', 'O', 'S', 'F', 'A', 'D', '0', '1' };

	bool load_free_atom_from_disk(const std::filesystem::path &path, const std::string &disk_key,
		FreeAtomDiskEntry &entry) {
		if (path.empty())
			return false;
		std::ifstream in(path, std::ios::binary);
		if (!in)
			return false;
		char magic[8];
		uint64_t key_size = 0;
		in.read(magic, 8);
		in.read(reinterpret_cast<char *>(&key_size), sizeof(key_size));
		if (!in || std::memcmp(magic, free_atom_disk_magic, 8) != 0 || key_size != disk_key.size())
			return false;
		std::string stored(key_size, '\0');
		in.read(stored.data(), key_size);
		if (!in || stored != disk_key)
			return false;
		int32_t conv = 0, iter = 0;
		int64_t rows = 0, cols = 0;
		in.read(reinterpret_cast<char *>(&conv), sizeof(conv));
		in.read(reinterpret_cast<char *>(&iter), sizeof(iter));
		in.read(reinterpret_cast<char *>(&entry.ediff_rel), sizeof(double));
		in.read(reinterpret_cast<char *>(&entry.diis_error), sizeof(double));
		in.read(reinterpret_cast<char *>(&rows), sizeof(rows));
		in.read(reinterpret_cast<char *>(&cols), sizeof(cols));
		if (!in || rows <= 0 || cols <= 0 || rows > 100000 || cols > 100000)
			return false;
		vec values(static_cast<size_t>(rows * cols));
		in.read(reinterpret_cast<char *>(values.data()), values.size() * sizeof(double));
		if (!in)
			return false;
		entry.density = dMatrix2(rows, cols);
		for (int64_t r = 0; r < rows; r++)
			for (int64_t col = 0; col < cols; col++)
				entry.density(r, col) = values[r * cols + col];
		entry.converged = conv != 0;
		entry.iter = iter;
		return true;
	}

	//Write to a temporary name and rename, so readers and racing writers only see complete files; failures are ignored.
	void store_free_atom_to_disk(const std::filesystem::path &path, const std::string &disk_key,
		const FreeAtomDiskEntry &entry) {
		if (path.empty())
			return;
		std::error_code ec;
		std::filesystem::create_directories(path.parent_path(), ec);
		std::ostringstream tmp_name;
		tmp_name << path.filename().string() << ".tmp."
			<< std::hash<std::thread::id>{}(std::this_thread::get_id()) << "."
			<< std::chrono::steady_clock::now().time_since_epoch().count();
		const auto tmp = path.parent_path() / tmp_name.str();
		{
			std::ofstream out(tmp, std::ios::binary);
			if (!out)
				return;
			const uint64_t key_size = disk_key.size();
			const int32_t conv = entry.converged ? 1 : 0, iter = entry.iter;
			const int64_t rows = entry.density.extent(0), cols = entry.density.extent(1);
			out.write(free_atom_disk_magic, 8);
			out.write(reinterpret_cast<const char *>(&key_size), sizeof(key_size));
			out.write(disk_key.data(), key_size);
			out.write(reinterpret_cast<const char *>(&conv), sizeof(conv));
			out.write(reinterpret_cast<const char *>(&iter), sizeof(iter));
			out.write(reinterpret_cast<const char *>(&entry.ediff_rel), sizeof(double));
			out.write(reinterpret_cast<const char *>(&entry.diis_error), sizeof(double));
			out.write(reinterpret_cast<const char *>(&rows), sizeof(rows));
			out.write(reinterpret_cast<const char *>(&cols), sizeof(cols));
			for (int64_t r = 0; r < rows; r++)
				for (int64_t col = 0; col < cols; col++) {
					const double x = entry.density(r, col);
					out.write(reinterpret_cast<const char *>(&x), sizeof(double));
				}
			if (!out) {
				out.close();
				std::filesystem::remove(tmp, ec);
				return;
			}
		}
		std::filesystem::rename(tmp, path, ec);
		if (ec)
			std::filesystem::remove(tmp, ec);
	}

	void warn_unconverged_free_atom(const std::string &label, const int Z, const size_t nbf, const int iter,
		const double ediff_rel, const double diis_error) {
		std::ostringstream line;
		line << "  WARNING: the free-atom SCF of " << label << " (Z=" << Z << ", " << nbf
			<< " functions) did not converge in " << iter << " iterations: |dE|/E=" << std::scientific
			<< std::setprecision(3) << ediff_rel << ", max|FDS-SDF|=" << diis_error
			<< ". Its density is used as the free-atom reference for every atom of this element, so"
			<< " the bond indices below inherit that residual.";
		rgbi_debug_line(line.str());
	}

	dMatrix2 compute_tonto_style_atomic_density(
		const atom &atm, const e_origin origin, const bool cartesian) {
		const FreeAtomKey key{ atm.get_basis_set(), atm.get_charge(), atm.get_ECP_electrons(),
							   origin, cartesian };
		//NOS_RGBI_NO_FREEATOM_CACHE recomputes every centre (the reference for cached numbers) and neither reads
		//nor writes the cache. Read per call so a test can flip it within one process.
		const bool cache_disabled = tuning("NOS_RGBI_NO_FREEATOM_CACHE") != nullptr;
		if (!cache_disabled) {
			const std::lock_guard<std::mutex> hold(free_atom_cache_mutex);
			for (const auto &cached : free_atom_cache)
				if (cached.first == key) {
					if (tuning("NOS_RGBI_DEBUG") != nullptr) {
						std::ostringstream line;
						line << "FREEATOM-CACHED " << atm.get_label() << " Z="
							<< atm.get_charge() - atm.get_ECP_electrons();
						rgbi_debug_line(line.str());
					}
					return cached.second;
				}
		}

		//All-electron atoms get sf-X2C-1e; an ECP carries relativity in its potential.
		//NOS_RGBI_NONREL_FREEATOM restores the non-relativistic free atom, for comparison only.
		const bool relativistic = atm.get_ECP_electrons() == 0 &&
			tuning("NOS_RGBI_NONREL_FREEATOM") == nullptr;
		const int effective_atomic_number = atm.get_charge() - atm.get_ECP_electrons();
		const std::string disk_key = cache_disabled ? std::string() : free_atom_disk_key(key, relativistic);
		const std::filesystem::path disk_path =
			cache_disabled ? std::filesystem::path() : free_atom_disk_path(disk_key, effective_atomic_number);
		if (FreeAtomDiskEntry hit; !cache_disabled && load_free_atom_from_disk(disk_path, disk_key, hit)) {
			if (tuning("NOS_RGBI_DEBUG") != nullptr)
				rgbi_debug_line("FREEATOM-DISK " + atm.get_label() + " Z=" + std::to_string(effective_atomic_number) +
					" " + disk_path.string());
			if (!hit.converged)
				warn_unconverged_free_atom(atm.get_label(), effective_atomic_number, hit.density.extent(0),
					hit.iter, hit.ediff_rel, hit.diis_error);
			const std::lock_guard<std::mutex> hold(free_atom_cache_mutex);
			free_atom_cache.push_back({ key, hit.density });
			return hit.density;
		}

		ScopedOccLogLevel quiet_occ_logs(spdlog::level::err);
		const ScopedOccSingleThread deterministic_reduction;
		const occ::gto::AOBasis basis = build_occ_atomic_basis_from_wfn_atom(atm, origin, cartesian);
		const int multiplicity = std::max(1, tonto_ground_state_multiplicity(effective_atomic_number));
		const bool restricted = multiplicity == 1 && (effective_atomic_number % 2 == 0);
		const auto spin_kind = restricted
			? occ::qm::SpinorbitalKind::Restricted
			: occ::qm::SpinorbitalKind::Unrestricted;

		//occ does not refuse more electrons than orbitals and writes past its occupied block (an ECP file read
		//without -ECP); refuse here, the caller catches this and falls back to the molecular local orbitals.
		const int n_alpha = (effective_atomic_number + multiplicity - 1) / 2;
		if (n_alpha > static_cast<int>(basis.nbf()))
			throw std::runtime_error(
				"The free-atom SCF for " + atm.get_label() + " (Z = " + std::to_string(atm.get_charge()) +
				", ECP core " + std::to_string(atm.get_ECP_electrons()) + ") would fill " +
				std::to_string(n_alpha) + " orbitals of a basis that has only " +
				std::to_string(basis.nbf()) + " functions. A wavefunction computed with an ECP has to be "
				"read as one - pass -ECP - or its cores are counted against a basis that never "
				"described them.");

		//Printed before the SCF so a crash inside occ names the atom. pinned= is 1 only with a global_control at one thread.
		if (tuning("NOS_RGBI_DEBUG") != nullptr) {
			std::ostringstream line;
			line << "FREEATOM-START " << atm.get_label() << " Z=" << effective_atomic_number
				<< " mult=" << multiplicity << (restricted ? " restricted" : " unrestricted")
				<< " nbf=" << basis.nbf() << " nsh=" << basis.size() << " n_alpha=" << n_alpha
				<< (relativistic ? " sfx2c1e" : " nonrel")
				<< " pinned=" << ((occ::parallel::get_tbb_control() != nullptr &&
					occ::parallel::get_num_threads() == 1) ? 1 : 0);
			rgbi_debug_line(line.str());
		}

		//Timed per SCF: one atom can dominate a run, so a run total divided by the SCF count is no SCF cost.
		const auto scf_t0 = std::chrono::steady_clock::now();

		occ::qm::HartreeFock hf(basis);
		occ::qm::SCF<occ::qm::HartreeFock> scf(hf, spin_kind);
		//occ's own guess is SOAD where its minimal basis covers the atom; elsewhere it nests an atomic SCF whose
		//square guess density breaks an unrestricted Fock build (heap corruption on Ce, U). Use the core guess only
		//there: forcing it everywhere converges light atoms to a different state.
		if (spin_kind == occ::qm::SpinorbitalKind::Unrestricted &&
			!occ::qm::minimal_basis_covers(basis))
			scf.set_guess_kind(occ::qm::GuessKind::Core);
		scf.set_charge_multiplicity(0, multiplicity);
		//NOS_RGBI_FREE_ATOM_MAXITER caps the iterations so the unconverged warning is testable; unset leaves occ's default.
		if (const char *cap = tuning("NOS_RGBI_FREE_ATOM_MAXITER")) {
			const int capped_iterations = std::atoi(cap);
			if (capped_iterations > 0)
				scf.maxiter = capped_iterations;
		}
		if (relativistic)
			scf.set_external_potential(sfx2c1e_core_correction(basis, hf.compute_overlap_matrix(),
				hf.compute_kinetic_matrix(), hf.compute_nuclear_attraction_matrix()), 0.0, "sfx2c1e");
		const double scf_energy = scf.compute_scf_energy();

		//occ returns the last energy without throwing when it runs out of iterations. Warn rather than refuse:
		//throwing would drop the element's ANO route, a larger change than the residual.
		if (!scf.ctx.converged)
			warn_unconverged_free_atom(atm.get_label(), effective_atomic_number, basis.nbf(), scf.iter,
				scf.ediff_rel, scf.diis_error);

		occ::qm::MolecularOrbitals mo = scf.wavefunction().mo;
		mo.update_occupied_orbitals();
		mo.update_density_matrix();
		if (tuning("NOS_RGBI_DEBUG") != nullptr) {
			std::ostringstream line;
			line << "FREEATOM " << atm.get_label() << " Z=" << effective_atomic_number
				<< " mult=" << multiplicity << (restricted ? " restricted" : " unrestricted")
				<< " nbf=" << basis.nbf() << " na=" << mo.n_alpha << " nb=" << mo.n_beta
				<< std::setprecision(14) << " E=" << scf_energy
				<< " sumD=" << mo.D.sum() << " normD=" << mo.D.norm();
			rgbi_debug_line(line.str());
		}
		//A separate line: the harness compares the FREEATOM line from " Z=" to the end, and a wall clock there
		//would make every line unique. It carries no " E=", which that test filters on.
		if (tuning("NOS_RGBI_DEBUG") != nullptr) {
			std::ostringstream line;
			line << "FREEATOM-TIME " << atm.get_label() << " Z=" << effective_atomic_number
				<< " nbf=" << basis.nbf() << " secs="
				<< std::chrono::duration<double>(std::chrono::steady_clock::now() - scf_t0).count();
			rgbi_debug_line(line.str());
		}

		dMatrix2 result = (spin_kind == occ::qm::SpinorbitalKind::Restricted)
			? eigen_matrix_to_dmatrix2(2.0 * mo.D)
			: eigen_matrix_to_dmatrix2(occ::Mat(occ::qm::block::a(mo.D) + occ::qm::block::b(mo.D)));
		//An SCF that threw is never cached; one that ran out of iterations is, so every atom of the element gets the same reference.
		if (!cache_disabled) {
			store_free_atom_to_disk(disk_path, disk_key,
				{ result, static_cast<bool>(scf.ctx.converged), static_cast<int>(scf.iter), scf.ediff_rel, scf.diis_error });
			const std::lock_guard<std::mutex> hold(free_atom_cache_mutex);
			free_atom_cache.push_back({ key, result });
		}
		return result;
	}

	//Free-atom SCFs are independent but the consuming loop is not (it mutates the overlap, a running index and
	//NAOs in order), so they can be computed up front concurrently and the serial loop finds them cached.
	//Opt-in via NOS_RGBI_PARALLEL_FREEATOM: occ's SCF is not shown re-entrant, and the result must match the
	//serial digest bit for bit. Failures are swallowed: an uncached atom is handled by the serial loop, and an
	//exception crossing the OpenMP region would terminate.
	void warm_free_atom_cache(const std::vector<atom> &ats, const e_origin origin,
		const bool cartesian) {
		//With the cache off nothing would read what this warms.
		if (tuning("NOS_RGBI_PARALLEL_FREEATOM") == nullptr ||
			tuning("NOS_RGBI_NO_FREEATOM_CACHE") != nullptr)
			return;

		//One representative per distinct cache key, i.e. exactly the set the serial loop would compute.
		std::vector<const atom *> distinct;
		std::vector<FreeAtomKey> keys;
		for (const auto &a : ats) {
			const FreeAtomKey key{ a.get_basis_set(), a.get_charge(), a.get_ECP_electrons(),
								   origin, cartesian };
			bool seen = false;
			for (const auto &k : keys)
				if (k == key) { seen = true; break; }
			if (!seen) {
				keys.push_back(key);
				distinct.push_back(&a);
			}
		}
		if (distinct.size() < 2)
			return;
		if (tuning("NOS_RGBI_DEBUG") != nullptr) {
			std::ostringstream line;
			line << "FREEATOM-WARM " << distinct.size() << " distinct of " << ats.size() << " centres";
			rgbi_debug_line(line.str());
		}

		//Pinned once here so the per-SCF guards are no-ops and no thread races on occ's thread count or pool.
		const ScopedOccSingleThread deterministic_reduction;
		//spdlog's level is global; per-SCF restores from several threads would leak occ's output.
		const ScopedOccLogLevel quiet_occ_logs(spdlog::level::err);
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < static_cast<int>(distinct.size()); i++) {
			try {
				compute_tonto_style_atomic_density(*distinct[i], origin, cartesian);
			}
			catch (...) {
			}
		}
	}

	const std::vector<OhOperation> &oh_operations() {
		static const std::vector<OhOperation> operations = [] {
			std::vector<OhOperation> result;
			result.reserve(48);
			std::array<int, 3> permutation{ 0, 1, 2 };
			do {
				for (int mask = 0; mask < 8; ++mask) {
					result.push_back({
						permutation,
						{ (mask & 1) ? -1 : 1,
						  (mask & 2) ? -1 : 1,
						  (mask & 4) ? -1 : 1 }
						});
				}
			} while (std::next_permutation(permutation.begin(), permutation.end()));
			return result;
			}();
		return operations;
	}

	std::vector<std::array<int, 3>> cartesian_exponents(const int l) {
		const int size = cartesian_shell_size(l);
		// constants::type_vector currently contains Cartesian components through i.
		const int first_type = l * (l + 1) * (l + 2) / 6 + 1;
		std::vector<std::array<int, 3>> result(size);
		for (int component = 0; component < size; ++component) {
			int exponent[3];
			constants::type2vector(first_type + component, exponent);
			err_checkf(exponent[0] >= 0 && exponent[0] + exponent[1] + exponent[2] == l,
				"Cartesian O_h symmetrization does not support angular momentum l=" +
				std::to_string(l) + ".", std::cout);
			result[component] = { exponent[0], exponent[1], exponent[2] };
		}
		return result;
	}

	std::vector<std::array<int, 3>> libcint_cartesian_exponents(const int l) {
		std::vector<std::array<int, 3>> result;
		result.reserve(cartesian_shell_size(l));
		for (int lx = l; lx >= 0; --lx) {
			for (int ly = l - lx; ly >= 0; --ly)
				result.push_back({ lx, ly, l - lx - ly });
		}
		return result;
	}

	CartesianTransform cartesian_transform_from_exponents(
		const std::vector<std::array<int, 3>> &exponents,
		const OhOperation &operation) {
		CartesianTransform transform{ ivec(exponents.size()), ivec(exponents.size(), 1) };

		for (int source = 0; source < static_cast<int>(exponents.size()); ++source) {
			std::array<int, 3> transformed{ 0, 0, 0 };
			int sign = 1;
			for (int axis = 0; axis < 3; ++axis) {
				transformed[operation.permutation[axis]] = exponents[source][axis];
				if (operation.signs[axis] < 0 && (exponents[source][axis] & 1))
					sign = -sign;
			}

			const auto destination = std::find(exponents.begin(), exponents.end(), transformed);
			err_checkf(destination != exponents.end(),
				"O_h operation produced an unknown Cartesian component.", std::cout);
			transform.destination[source] = static_cast<int>(destination - exponents.begin());
			transform.signs[source] = sign;
		}
		return transform;
	}

	dMatrix2 cartesian_transform_matrix(const int l, const OhOperation &operation,
		const bool libcint_order = false) {
		const auto exponents = libcint_order
			? libcint_cartesian_exponents(l)
			: cartesian_exponents(l);
		const auto sparse = cartesian_transform_from_exponents(exponents, operation);
		const int size = cartesian_shell_size(l);
		dMatrix2 result(size, size);
		std::fill(result.container().begin(), result.container().end(), 0.0);
		for (int source = 0; source < size; ++source)
			result(sparse.destination[source], source) = sparse.signs[source];
		return result;
	}

	dMatrix2 spherical_transform_matrix(const int l, const OhOperation &operation) {
		const dMatrix2 cart_to_spherical = cart2sph(l, true);
		const dMatrix2 cart_transform = cartesian_transform_matrix(l, operation, true);
		const int n_cart = cartesian_shell_size(l);
		const int n_spherical = 2 * l + 1;

		// Columns of C contain the real spherical functions in the Cartesian basis.
		// Solve T_cart C = C T_sph using the left pseudoinverse of C.  This keeps
		// libcint's real-spherical ordering and phase convention authoritative.
		dMatrix2 gram(n_spherical, n_spherical);
		for (int i = 0; i < n_spherical; ++i)
			for (int j = 0; j < n_spherical; ++j)
				for (int cart = 0; cart < n_cart; ++cart)
					gram(i, j) += cart_to_spherical(cart, i) * cart_to_spherical(cart, j);
		const dMatrix2 inverse_gram = LAPACKE_invert(gram, 1E-12);

		dMatrix2 transformed_cart(n_cart, n_spherical);
		for (int i = 0; i < n_cart; ++i)
			for (int j = 0; j < n_spherical; ++j)
				for (int cart = 0; cart < n_cart; ++cart)
					transformed_cart(i, j) += cart_transform(i, cart) * cart_to_spherical(cart, j);

		dMatrix2 result(n_spherical, n_spherical);
		for (int i = 0; i < n_spherical; ++i) {
			for (int j = 0; j < n_spherical; ++j) {
				for (int k = 0; k < n_spherical; ++k) {
					for (int cart = 0; cart < n_cart; ++cart) {
						result(i, j) += inverse_gram(i, k) *
							cart_to_spherical(cart, k) * transformed_cart(cart, j);
					}
				}
				if (std::abs(result(i, j)) < 1E-12)
					result(i, j) = 0.0;
			}
		}
		return result;
	}
} // namespace

//Lets the harness count SCFs per run and a long-lived process free the matrices.
void clear_rgbi_free_atom_cache() {
	const std::lock_guard<std::mutex> hold(free_atom_cache_mutex);
	free_atom_cache.clear();
}

std::string rgbi_supported_input_phrase(const std::string &refused_extension) {
	std::string ext = refused_extension;
	std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return (char)std::tolower(c); });
	if (!ext.empty() && ext.front() != '.')
		ext.insert(ext.begin(), '.');
	std::vector<std::string> works;
	for (const char *w : { ".gbw", ".molden" })
		if (ext != w)
			works.push_back(w);
	std::string phrase;
	for (size_t i = 0; i < works.size(); i++)
		phrase += (i == 0 ? "a " : (i + 1 == works.size() ? " or a " : ", a ")) + works[i];
	//both formats cannot be the refused one at once, so the list is never empty
	return phrase;
}

int highest_shell_angular_momentum(const WFN &wavy) {
	int highest = -1;
	for (const atom &a : wavy.get_atoms())
		for (const basis_set_entry &bf : a.get_basis_set())
			highest = std::max(highest, static_cast<int>(bf.get_type()) - 1);
	return highest;
}

void symmetrize_atomic_matrix_oh(dMatrix2 &matrix, const ivec &shell_angular_momenta,
	const bool spherical) {
	const bool square = matrix.extent(0) == matrix.extent(1);
	err_checkf(square, "Cannot O_h-symmetrize a non-square atomic matrix.", std::cout);
	if (!square)
		return;

	ivec shell_offsets(shell_angular_momenta.size() + 1, 0);
	for (int shell = 0; shell < static_cast<int>(shell_angular_momenta.size()); ++shell) {
		const int l = shell_angular_momenta[shell];
		const bool supported = l >= 0 && l <= 5;
		err_checkf(supported,
			"O_h symmetrization supports shells from s through h.", std::cout);
		if (!supported)
			return;
		shell_offsets[shell + 1] = shell_offsets[shell] +
			(spherical ? 2 * l + 1 : cartesian_shell_size(l));
	}

	const bool correct_size = shell_offsets.back() == static_cast<int>(matrix.extent(0));
	err_checkf(correct_size,
		"Atomic matrix size does not match its Cartesian shell description.", std::cout);
	if (!correct_size)
		return;

	dMatrix2 symmetrized(matrix.extent(0), matrix.extent(1));
	std::fill(symmetrized.container().begin(), symmetrized.container().end(), 0.0);

	const auto &operations = oh_operations();
	for (const auto &operation : operations) {
		if (spherical) {
			std::vector<dMatrix2> transforms;
			transforms.reserve(shell_angular_momenta.size());
			for (const int l : shell_angular_momenta)
				transforms.push_back(spherical_transform_matrix(l, operation));

			for (int shell_a = 0; shell_a < static_cast<int>(shell_angular_momenta.size()); ++shell_a) {
				const int first_a = shell_offsets[shell_a];
				const int size_a = 2 * shell_angular_momenta[shell_a] + 1;
				for (int shell_b = 0; shell_b < static_cast<int>(shell_angular_momenta.size()); ++shell_b) {
					const int first_b = shell_offsets[shell_b];
					const int size_b = 2 * shell_angular_momenta[shell_b] + 1;
					for (int transformed_a = 0; transformed_a < size_a; ++transformed_a) {
						for (int transformed_b = 0; transformed_b < size_b; ++transformed_b) {
							for (int a = 0; a < size_a; ++a) {
								for (int b = 0; b < size_b; ++b) {
									symmetrized(first_a + transformed_a, first_b + transformed_b) +=
										transforms[shell_a](transformed_a, a) *
										matrix(first_a + a, first_b + b) *
										transforms[shell_b](transformed_b, b);
								}
							}
						}
					}
				}
			}
			continue;
		}

		std::vector<CartesianTransform> transforms;
		transforms.reserve(shell_angular_momenta.size());
		for (const int l : shell_angular_momenta)
			transforms.push_back(cartesian_transform_from_exponents(
				cartesian_exponents(l), operation));

		for (int shell_a = 0; shell_a < static_cast<int>(shell_angular_momenta.size()); ++shell_a) {
			const int first_a = shell_offsets[shell_a];
			const int size_a = cartesian_shell_size(shell_angular_momenta[shell_a]);
			const auto &transform_a = transforms[shell_a];
			for (int shell_b = 0; shell_b < static_cast<int>(shell_angular_momenta.size()); ++shell_b) {
				const int first_b = shell_offsets[shell_b];
				const int size_b = cartesian_shell_size(shell_angular_momenta[shell_b]);
				const auto &transform_b = transforms[shell_b];
				for (int a = 0; a < size_a; ++a) {
					const int transformed_a = first_a + transform_a.destination[a];
					for (int b = 0; b < size_b; ++b) {
						const int transformed_b = first_b + transform_b.destination[b];
						symmetrized(transformed_a, transformed_b) +=
							transform_a.signs[a] * transform_b.signs[b] *
							matrix(first_a + a, first_b + b);
					}
				}
			}
		}
	}

	const double inverse_order = 1.0 / static_cast<double>(operations.size());
	for (auto &value : symmetrized.container())
		value *= inverse_order;
	matrix = std::move(symmetrized);
}

//Average over the full rotation group. By Schur's lemma the only invariant between two shells of equal l
//is a multiple of the identity on their 2l+1 components, and none exists between different l, so the
//result keeps no direction (O_h keeps e_g/t_2g and above). Reads only block diagonals, so it is
//independent of m order and phase convention and has no l limit.
void spherically_average_atomic_matrix(dMatrix2 &matrix, const ivec &shell_angular_momenta) {
	const bool square = matrix.extent(0) == matrix.extent(1);
	err_checkf(square, "Cannot spherically average a non-square atomic matrix.", std::cout);
	if (!square)
		return;

	ivec shell_offsets(shell_angular_momenta.size() + 1, 0);
	for (int shell = 0; shell < static_cast<int>(shell_angular_momenta.size()); ++shell) {
		const int l = shell_angular_momenta[shell];
		err_checkf(l >= 0, "Cannot spherically average a shell of negative angular momentum.", std::cout);
		if (l < 0)
			return;
		shell_offsets[shell + 1] = shell_offsets[shell] + 2 * l + 1;
	}

	const bool correct_size = shell_offsets.back() == static_cast<int>(matrix.extent(0));
	err_checkf(correct_size,
		"Atomic matrix size does not match its spherical shell description.", std::cout);
	if (!correct_size)
		return;

	dMatrix2 averaged(matrix.extent(0), matrix.extent(1));
	std::fill(averaged.container().begin(), averaged.container().end(), 0.0);
	for (int shell_a = 0; shell_a < static_cast<int>(shell_angular_momenta.size()); ++shell_a) {
		const int l = shell_angular_momenta[shell_a];
		const int size = 2 * l + 1;
		for (int shell_b = 0; shell_b < static_cast<int>(shell_angular_momenta.size()); ++shell_b) {
			//different l: no invariant, block dropped; equal l keeps the coupling (2p-3p)
			if (shell_angular_momenta[shell_b] != l)
				continue;
			double sum = 0.0;
			for (int m = 0; m < size; ++m)
				sum += matrix(shell_offsets[shell_a] + m, shell_offsets[shell_b] + m);
			const double mean = sum / static_cast<double>(size);
			for (int m = 0; m < size; ++m)
				averaged(shell_offsets[shell_a] + m, shell_offsets[shell_b] + m) = mean;
		}
	}
	matrix = std::move(averaged);
}

void print_dmatrix2(const dMatrix2 &EVC2, const std::string name) {
	std::cout << std::endl << name << ":\n";
	for (int i = 0; i < EVC2.extent(0); i++) {
		for (int j = 0; j < EVC2.extent(1); j++)
			std::cout << std::setw(14) << std::setprecision(8) << std::fixed << EVC2(i, j) << " ";
		std::cout << "\n";
	}
}


int compute_dens(WFN &wavy, bool debug, int *np, double *origin, double *gvector, double *incr, std::string &outname, bool rho, bool rdg, bool eli, bool lap) {
	properties_options opts;
	opts.lap = lap;
	opts.eli = eli;
	opts.rdg = rdg;
	opts.rho = rho;
	opts.resolution = *incr;
	opts.NbSteps = { np[0], np[1], np[2] };

	//indexed by cube_type: Rho, RDG, Elf, Eli, Lap, then the eleven properties this run never asks for
	std::vector<cube> cubes = { {opts.NbSteps, wavy.get_ncen(), opts.rho || opts.rdg || opts.eli || opts.lap},
		{opts.NbSteps, wavy.get_ncen(), opts.rdg},
		{},
		{opts.NbSteps, wavy.get_ncen(), opts.eli},
		{opts.NbSteps, wavy.get_ncen(), opts.lap},
		{}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {} };
	cubes[cube_type::Rho].give_parent_wfn(wavy);
	cubes[cube_type::RDG].give_parent_wfn(wavy);
	cubes[cube_type::Eli].give_parent_wfn(wavy);
	cubes[cube_type::Lap].give_parent_wfn(wavy);

	std::string Oname = outname;
	ivec ntyp;
	for (int i = 0; i < 3; i++) {
		if (debug) {
			std::cout << "gvector before: ";
			for (int j = 0; j < 3; j++) std::cout << gvector[j + 3 * i] << " ";
			std::cout << "\n";
		}
		for (int j = 0; j < 3; j++)
			gvector[j + 3 * i] = incr[i] * gvector[j + 3 * i];
		if (debug) {
			for (int j = 0; j < 3; j++) std::cout << gvector[j + 3 * i] << " ";
			std::cout << "\n";
		}
	}
	for (int i = 0; i < 3; i++) {
		cubes[cube_type::Rho].set_origin(i, origin[i]);
		cubes[cube_type::RDG].set_origin(i, origin[i]);
		cubes[cube_type::Eli].set_origin(i, origin[i]);
		cubes[cube_type::Lap].set_origin(i, origin[i]);
		for (int j = 0; j < 3; j++) {
			cubes[cube_type::Rho].set_vector(i, j, gvector[j + 3 * i]);
			cubes[cube_type::RDG].set_vector(i, j, gvector[j + 3 * i]);
			cubes[cube_type::Eli].set_vector(i, j, gvector[j + 3 * i]);
			cubes[cube_type::Lap].set_vector(i, j, gvector[j + 3 * i]);
		}
	}
	cubes[cube_type::Rho].set_path(Oname + "_rho.cube");
	cubes[cube_type::RDG].set_path(Oname + "_rdg.cube");
	cubes[cube_type::Eli].set_path(Oname + "_eli.cube");
	cubes[cube_type::Lap].set_path(Oname + "_lap.cube");

	//opt.NbAtoms[0]=wavy.get_ncen();
	std::cout << "\n   .      ___________________________________________________________      .\n";
	std::cout << "  *.                                                                      *.\n";
	std::cout << "  *.  Wavefunction              : " << std::setw(20) << wavy.get_path().filename() << " / " << std::setw(5) << wavy.get_ncen() << " atoms      *.\n";
	std::cout << "  *.  OutPut filename Prefix    : " << std::setw(40) << Oname << "*.\n";
	std::cout << "  *.                                                                      *.\n";
	std::cout << "  *.  gridBox Min               : " << std::setw(11) << std::setprecision(6) << origin[0] << std::setw(12) << origin[1] << origin[2] << "     *.\n";
	std::cout << "  *.  gridBox Max               : " << std::setw(11) << std::setprecision(6) << (origin[0] + gvector[0] * np[0] + gvector[3] * np[1] + gvector[6] * np[2]) << std::setw(12) << (origin[1] + gvector[1] * np[0] + gvector[4] * np[1] + gvector[7] * np[2]) << (origin[2] + gvector[2] * np[0] + gvector[5] * np[1] + gvector[8] * np[2]) << "     *.\n";
	std::cout << "  *.  Increments(bohr)          : " << std::setw(11) << std::setprecision(6) << array_length(d3{ gvector[0], gvector[1], gvector[2] }) << std::setw(12) << array_length(d3{ gvector[3], gvector[4], gvector[5] }) << std::setw(12) << array_length(d3{ gvector[6], gvector[7], gvector[8] }) << "     *.\n";
	std::cout << "  *.  NbSteps                   : " << std::setw(11) << opts.NbSteps[0] << std::setw(12) << opts.NbSteps[1] << std::setw(12) << opts.NbSteps[2] << "     *.\n";
	std::cout << "  *.                                                                      *.\n";
	std::cout << "  *.  Number of primitives      :     " << std::setw(5) << wavy.get_nex() << "                               *.\n";
	std::cout << "  *.  Number of MOs             :       " << std::setw(3) << wavy.get_nmo() << "                               *.\n";

	//Calc_Prop leaves the Rho cube untouched unless RDG is requested (it hands back sign(lambda2)*rho then)
	if (rho && !rdg) Calc_Rho(cubes[cube_type::Rho], wavy, 20.0, std::cout, false);
	if (rdg || eli || lap) Calc_Prop(cubes, wavy, 20.0, std::cout, false, false);

	std::cout << "  *.                                                                      *.\n";
	std::cout << "  *.  Writing .cube files ...                                             *.\n";
	if (rho && !rdg) {
		cubes[cube_type::Rho].set_path(Oname + "_rho.cube");
		cubes[cube_type::Rho].write_file(true, true);
	}
	if (rdg) {
		cubes[cube_type::Rho].set_path(Oname + "_signed_rho.cube");
		cubes[cube_type::Rho].write_file(true);
		cubes[cube_type::Rho].set_path(Oname + "_rho.cube");
		cubes[cube_type::Rho].write_file(true, true);
		cubes[cube_type::RDG].write_file(true);
	}
	if (eli) cubes[cube_type::Eli].write_file(true);
	if (lap) cubes[cube_type::Lap].write_file(true);

	std::cout << "  *                                                                       *\n";
	std::cout << "          ___________________________________________________________\n";
	return 0;
};

bond do_bonds(WFN &wavy,
	int mode_sel, bool mode_leng, bool mode_res,
	double res[], bool cub, double boxsize[],
	int atom1, int atom2, int atom3,
	const bool &debug, const bool &bohr,
	int runnumber,
	bool rho, bool rdg, bool eli, bool lap) {
	bond results{ "","","","",false,false,false };
	std::string line("");
	std::vector <std::string> label;
	label.resize(3);
	double coords1[3], coords2[3], coords3[3];
	int na = 0;
	double z[3], x[3], y[3], help[3], size[3], s2[3];
	int np[3];
	//bohr: the user's lengths (box lengths when mode_leng, spacings when !mode_res) are Angstrom; everything below is bohr
	double res_b[3], boxsize_b[3];
	for (int r = 0; r < 3; r++)
	{
		res_b[r] = (bohr && !mode_res) ? constants::ang2bohr(res[r]) : res[r];
		boxsize_b[r] = (bohr && mode_leng) ? constants::ang2bohr(boxsize[r]) : boxsize[r];
	}
	na = wavy.get_ncen();
	if (atom1 <= 0 || atom2 <= 0 || atom3 <= 0 || mode_sel <= 0 || mode_sel > 4 || atom1 > na || atom2 > na || atom3 > na || atom1 == atom2 || atom2 == atom3 || atom1 == atom3)
	{
		std::cout << "Invalid selections of atoms or mode_sel, please try again!\n";
		return(results);
	}
	if (debug) std::cout << "No. of Atoms selected: " << atom1 << atom2 << atom2 << "\n";
	for (int i = 0; i < 3; i++)
		coords1[i] = wavy.get_atom_coordinate(atom1 - 1, i);
	for (int i = 0; i < 3; i++)
		coords2[i] = wavy.get_atom_coordinate(atom2 - 1, i);
	for (int i = 0; i < 3; i++)
		coords3[i] = wavy.get_atom_coordinate(atom3 - 1, i);
	label[0] = wavy.get_atom_label(atom1 - 1);
	label[1] = wavy.get_atom_label(atom2 - 1);
	label[2] = wavy.get_atom_label(atom3 - 1);
	if (debug)
	{
		std::cout << "The Atoms found corresponding to your selection are:\n";
		std::cout << label[0] << " " << coords1[0] << " " << coords1[1] << " " << coords1[2] << "\n";
		std::cout << label[1] << " " << coords2[0] << " " << coords2[1] << " " << coords2[2] << "\n";
		std::cout << label[2] << " " << coords3[0] << " " << coords3[1] << " " << coords3[2] << "\n";
	}
	double znorm = sqrt((coords2[0] - coords1[0]) * (coords2[0] - coords1[0]) + (coords2[1] - coords1[1]) * (coords2[1] - coords1[1]) + (coords2[2] - coords1[2]) * (coords2[2] - coords1[2]));
	z[0] = (coords2[0] - coords1[0]) / znorm;
	z[1] = (coords2[1] - coords1[1]) / znorm;
	z[2] = (coords2[2] - coords1[2]) / znorm;
	help[0] = coords3[0] - coords2[0];
	help[1] = coords3[1] - coords2[1];
	help[2] = coords3[2] - coords2[2];
	y[0] = z[1] * help[2] - z[2] * help[1];
	y[1] = z[2] * help[0] - z[0] * help[2];
	y[2] = z[0] * help[1] - z[1] * help[0];
	double ynorm = sqrt(y[0] * y[0] + y[1] * y[1] + y[2] * y[2]);
	y[0] = y[0] / ynorm;
	y[1] = y[1] / ynorm;
	y[2] = y[2] / ynorm;
	x[0] = (z[1] * y[2] - z[2] * y[1]);
	x[1] = (z[2] * y[0] - z[0] * y[2]);
	x[2] = (z[0] * y[1] - z[1] * y[0]);
	z[0] = -z[0];
	z[1] = -z[1];
	z[2] = -z[2];

	double hnorm = 0.0;
	double ringhelp[3] = {};
	if (mode_sel == 2 || mode_sel == 1) hnorm = znorm;
	else if (mode_sel == 3) hnorm = sqrt((coords3[0] - coords1[0]) * (coords3[0] - coords1[0]) + (coords3[1] - coords1[1]) * (coords3[1] - coords1[1]) + (coords3[2] - coords1[2]) * (coords3[2] - coords1[2]));
	else if (mode_sel == 4)
	{
		for (int r = 0; r < 3; r++) ringhelp[r] = (coords3[r] + coords1[r] + coords2[r]) / 3;
		hnorm = (sqrt((coords3[0] - ringhelp[0]) * (coords3[0] - ringhelp[0]) + (coords3[1] - ringhelp[1]) * (coords3[1] - ringhelp[1]) + (coords3[2] - ringhelp[2]) * (coords3[2] - ringhelp[2]))
			+ sqrt((coords1[0] - ringhelp[0]) * (coords1[0] - ringhelp[0]) + (coords1[1] - ringhelp[1]) * (coords1[1] - ringhelp[1]) + (coords1[2] - ringhelp[2]) * (coords1[2] - ringhelp[2]))
			+ sqrt((coords2[0] - ringhelp[0]) * (coords2[0] - ringhelp[0]) + (coords2[1] - ringhelp[1]) * (coords2[1] - ringhelp[1]) + (coords2[2] - ringhelp[2]) * (coords2[2] - ringhelp[2]))) / 3;
	}
	if (debug)
	{
		std::cout << "Your three vectors are:\n";
		std::cout << "X= " << x[0] << " " << x[1] << " " << x[2] << "\n";
		std::cout << "Y= " << y[0] << " " << y[1] << " " << y[2] << "\n";
		std::cout << "Z= " << z[0] << " " << z[1] << " " << z[2] << "\n";
		std::cout << "From here on mode_leng and cub are used\n";
	}
	for (int r = 0; r < 3; r++)
	{
		if (res[r] < 0)
		{
			std::cout << "Wrong input in res! Try again!\n";
			return(results);
		}
		if (boxsize[r] < 0)
		{
			std::cout << "Wrong input for box scaling! Try again!\n";
			return(results);
		}
		else if (boxsize[r] >= 50)
		{
			std::cout << "Come on, be realistic!\nTry Again!\n";
			return(results);
		}
	}
	if (!mode_res)
	{
		if (debug) std::cout << "Determining resolution\n";
		if (mode_leng)
		{
			if (debug) std::cout << "mres=true; mleng=true\n";
			switch (mode_sel)
			{
			case 1:
				size[0] = ceil(9 * znorm) / 10;
				break;
			case 2:
				size[0] = ceil(15 * znorm) / 10;
				break;
			case 3:
				size[0] = ceil(15 * hnorm) / 10;
				break;
			case 4:
				size[0] = ceil(30 * hnorm) / 10;
				break;
			}
			for (int r = 0; r < 3; r++)
			{
				if (boxsize_b[r] != 0) size[r] = boxsize_b[r];
				else if (r > 0 || cub == true) size[r] = size[0];
				s2[r] = size[r] / 2;
			}
		}
		else
		{
			if (debug) std::cout << "mres=true; mleng=false\n";
			for (int r = 0; r < 3; r++)
			{
				if (cub && r > 0)
				{
					if (debug) std::cout << "I'm making things cubic!\n";
					s2[r] = s2[0];
					size[r] = size[0];
					continue;
				}
				switch (mode_sel)
				{
				case 1:
				case 2:
					size[r] = znorm;
					break;
				case 3:
				case 4:
					size[r] = hnorm;
					break;
				}
				if (boxsize[r] != 0) size[r] = size[r] * boxsize[r];
				else {
					switch (mode_sel)
					{
					case 1:
						size[r] = ceil(9 * znorm) / 10;
						break;
					case 2:
						size[r] = ceil(15 * znorm) / 10;
						break;
					case 3:
						size[r] = ceil(15 * hnorm) / 10;
						break;
					case 4:
						size[r] = ceil(30 * hnorm) / 10;
						break;
					}
				}
				s2[r] = size[r] / 2;
			}
		}
		for (int r = 0; r < 3; r++) np[r] = (int)round((2 * s2[r]) / res_b[cub ? 0 : r]) + 1;
	}
	else
	{
		if (debug) std::cout << "Boxsize is used\n";
		for (int r = 0; r < 3; r++)
		{
			if (cub) np[r] = (int)res[0];
			else np[r] = (int)res[r];
		}
		if (mode_leng)
		{
			if (debug) std::cout << "mode_leng=true; using boxsize\n";
			switch (mode_sel)
			{
			case 1:
				size[0] = ceil(9 * znorm) / 10;
				break;
			case 2:
				size[0] = ceil(15 * znorm) / 10;
				break;
			case 3:
				size[0] = ceil(15 * hnorm) / 10;
				break;
			case 4:
				size[0] = ceil(30 * hnorm) / 10;
				break;
			}
			for (int r = 0; r < 3; r++)
			{
				if (debug) std::cout << r + 1 << ". Axis:";
				if (boxsize_b[r] != 0) size[r] = boxsize_b[r];
				else if (r > 0 || cub == true) size[r] = size[0];
				s2[r] = size[r] / 2;
			}
		}
		else
		{
			if (debug) std::cout << "Enter multiplicator for length (z,y,x)\n";
			for (int r = 0; r < 3; r++)
			{
				if (cub == 1 && r > 0)
				{
					if (debug) std::cout << "I'm making things cubic!\n";
					s2[r] = s2[0];
					size[r] = size[0];
					continue;
				}
				if (debug) std::cout << r + 1 << ". Axis:";
				switch (mode_sel)
				{
				case 1:
				case 2:
					size[r] = znorm;
					break;
				case 3:
				case 4:
					size[r] = hnorm;
					break;
				}
				if (boxsize[r] != 0) size[r] = size[r] * boxsize[r];
				else {
					switch (mode_sel)
					{
					case 1:
						size[r] = ceil(9 * znorm) / 10;
						break;
					case 2:
						size[r] = ceil(15 * znorm) / 10;
						break;
					case 3:
						size[r] = ceil(15 * hnorm) / 10;
						break;
					case 4:
						size[r] = ceil(30 * hnorm) / 10;
						break;
					}
				}
				s2[r] = size[r] / 2;
			}
		}
	}
	//grid spacing: the requested one when res is a spacing, box length / points otherwise
	double incr[3];
	for (int i = 0; i < 3; i++)
		incr[i] = mode_res ? size[i] / np[i] : res_b[cub ? 0 : i];
	//origin = centre - half the grid extent (np - 1 steps per axis), so the grid points are symmetric about the centre
	double o[3] = {};
	if (debug) std::cout << "This is the origin:\n";
	for (int d = 0; d < 3; d++)
	{
		switch (mode_sel)
		{
		case 1:
			o[d] = coords1[d];
			break;
		case 2:
			o[d] = (coords1[d] + coords2[d]) / 2;
			break;
		case 3:
			o[d] = (coords3[d] + coords1[d]) / 2;
			break;
		case 4:
			o[d] = ringhelp[d];
			break;
		}
		o[d] -= 0.5 * ((np[0] - 1) * incr[0] * x[d] + (np[1] - 1) * incr[1] * y[d] + (np[2] - 1) * incr[2] * z[d]);
		if (debug) std::cout << o[d] << "\n";
	}
	std::string outname = { "" };
	outname += wavy.get_path().generic_string();
	outname += "_";
	outname += label[0];
	outname += std::to_string(atom1);
	outname += "_";
	outname += label[1];
	outname += std::to_string(atom2);
	outname += "_";
	outname += label[2];
	outname += std::to_string(atom3);
	outname += "_";
	outname += std::to_string(runnumber);
	double v[9];
	for (int i = 0; i < 3; i++) {
		v[i] = x[i];
		v[i + 3] = y[i];
		v[i + 6] = z[i];
	}
	if (compute_dens(wavy, debug, np, o, v, incr, outname, rho, rdg, eli, lap) == 0) {
		results.success = true;
		results.filename = outname;
	}
	return(results);
}

int autobonds(bool debug, WFN &wavy, const std::filesystem::path &inputfile, const bool &bohr) {
	char inputFile[2048] = "";
	if (!exists(inputfile))
	{
		std::cout << "No input file specified! I am going to make an example for you in input.example!\n";
		std::ofstream example("input.example");
		example << "!COMMENT LINES START WITH ! AND CAN ONLY BE ABOVE THE FIRST SWITCHES AND NUMBERS!\n!First row contains the following switches: rho (dens), Reduced Density Gradient (RDG), ELI-d (eli) and Laplacian of rho (lap)\n!Following rows each contain a bond you want to investigate.\n!The key to read these numbers is:\n!orientation_selection(1-4)\n!      1=atom1 centered\n!      2=bond atom1 atom2\n!      3=bond atom1 and atom3\n!      4=ring centroid of the three atoms\n!length selection(0/1)\n!      0=box will contain multiplicator of bondlength between atom1 and atom2\n!      1=box will contain length in angstrom (?)\n!resolution selection(0/1)\n!      0=res will contain number of gridpoints\n!      1=res will contain distance between gridpoints\n!res1 res2 res3\n!      based on selection above resolution in x,y,z direction (either gridpoints or distance between points)\n!cube selection (0/1)\n!      0= all selections are considered\n!      1= all selections in x-direction will be applied in the y and z direction, as well, making it a cube\n!box1 box2 box3\n!      based on selection above size of the box/cube (either multiplicator of bondlength or absolute length)\n!atom1 atom2 atom3\n!      the atoms in your wavefunction file (counting from 1) to be used as references\n! BELOW THE INPUT SECTION STARTS\n!rho, rdg, eli, lap\n!mode_sel(int) mode_leng(bool) mode_res(bool) res1(double) res2(double) res3(double) cube(bool) boxsize1(float) boxsize2(float) boxsize3(float) atom1(int) atom2(int) atom3(int)\n 1    1    1    1\n            2             1               1              20.0         20.0         20.0         0          5.0             5.1             5.2             1          2          3\n";
		example.close();
		return 0;
	}
	std::ifstream input(inputfile.c_str());
	if (!input.good())
	{
		std::cout << inputFile << " does not exist or is not readable!\n";
		return 0;
	}
	input.seekg(0);
	std::string line("");
	getline_universal(input, line);
	std::string comment("!");
	while (line.compare(0, 1, comment) == 0) getline_universal(input, line);
	int rho, rdg, eli, lap;
	if (line.length() < 10) return 0;
	else
	{
		std::istringstream iss(line);
		if (!(iss >> rho >> rdg >> eli >> lap)) {
			std::cerr << "Error parsing line for rho, rdg, eli, and lap values." << std::endl;
			return 0; // Handle the error appropriately
		}
	}
	int errorcount = 0, runnumber = 0;
	do {
		getline_universal(input, line);
		if (line.length() < 10) continue;
		runnumber++;
		int sel = 0, leng = 0, cube = 0, mres = 0, a1 = 0, a2 = 0, a3 = 0;
		double res[3], box[3];
		bool bleng = false, bres = false, bcube = false;
		if (line.length() > 1)
		{
			std::istringstream iss(line);
			if (!(iss >> sel >> leng >> mres >> res[0] >> res[1] >> res[2] >> cube >> box[0] >> box[1] >> box[2] >> a1 >> a2 >> a3)) {
				std::cerr << "Error parsing line for bond parameters.\n";
				continue; // Skip this line and move to the next
			}
		}
		if (leng == 1) { bleng = true; if (debug) std::cout << "leng=true\n"; }
		else { bleng = false; if (debug) std::cout << "leng=false\n"; }
		if (cube == 1) { bcube = true; if (debug) std::cout << "cube=true\n"; }
		else { bcube = false; if (debug) std::cout << "cube=false\n"; }
		if (mres == 1) { bres = true; if (debug) std::cout << "mres=true\n"; }
		else { bres = false; if (debug) std::cout << "mres=false\n"; }
		if (debug) std::cout << "running calculations for line " << runnumber << " of the input file:\n\n";
		bond work{ "","","","",false,false,false };
		work = do_bonds(wavy, sel, bleng, bres, res, bcube, box, a1, a2, a3, debug, bohr, runnumber, rho == 1, rdg == 1, eli == 1, lap == 1);
		if (!work.success)
		{
			std::cout << "!!!!!!!!!!!!problem somewhere during the calculations, see messages above!!!!!\n";
			errorcount++;
		}
		else {
			if (rho == 1) wavy.push_back_cube(work.filename + "_rho.cube", false, false);
			if (rdg == 1) wavy.push_back_cube(work.filename + "_rdg.cube", false, false);
			if (eli == 1) wavy.push_back_cube(work.filename + "_eli.cube", false, false);
			if (lap == 1) wavy.push_back_cube(work.filename + "_lap.cube", false, false);
			if (rho == 1 && rdg == 1) wavy.push_back_cube(work.filename + "_signed_rho.cube", false, false);
		}
	} while (!input.eof());
	std::cout << "\n  *   Finished all calculations! " << runnumber - errorcount << " out of " << runnumber << " were successful!            *\n";
	return 1;
}

std::vector<std::pair<int, int>> get_bonded_atom_pairs(const WFN &wavy) {
	std::vector<std::pair<int, int>> bonds;
	for (int i = 0; i < wavy.get_ncen(); i++)
	{
		for (int j = i + 1; j < wavy.get_ncen(); j++)
		{
			const double distance = array_length(wavy.get_atom_pos(i), wavy.get_atom_pos(j));
			const double svdW = constants::ang2bohr(constants::covalent_radii[wavy.get_atom_charge(i)] + constants::covalent_radii[wavy.get_atom_charge(j)]);
			if (distance < 1.25 * svdW)
			{
#ifdef NSA2DEBUG
				std::cout << "Bond between " << i << " (" << wavy.get_atom_charge(i) << ") and " << j << " (" << wavy.get_atom_charge(j) << ") with distance " << distance << " and svdW " << svdW << "\n";
#endif // NSA2DEBUG
				bonds.push_back(std::make_pair(i, j));
			}
		}
	}
	return bonds;
}


//Assuming square matrices
vec change_basis_sq(const vec &in, const vec &transformation, int size) {

	// new = transform^T * in * tranform
	vec result(static_cast<size_t>(size) * size);
	vec temp(static_cast<size_t>(size) * size);
	// first we do temp = t^T * i
	cblas_dgemm(CblasRowMajor,
		CblasTrans, CblasNoTrans,
		size, size, size,
		1.0,
		transformation.data(), size,
		in.data(), size,
		0.0,
		temp.data(), size);

	//Then we do res = temp * t
	cblas_dgemm(CblasRowMajor,
		CblasNoTrans, CblasNoTrans,
		size, size, size,
		1.0,
		temp.data(), size,
		transformation.data(), size,
		0.0,
		result.data(), size);
	return result;
}

//For non-square transformation matrices
//performs res = trans * in * trans^T in case of forward
//and res = trans^T * int * trans in case of !forward
dMatrix2 change_basis_general(const dMatrix2 &in, const dMatrix2 &transformation, bool forward = true) {
	//Checks should be handled by dot itself
	//int a, b, c, d;
	//c = transformation.extent(1); //cols of t
	//d = transformation.extent(0); //rows of t
	//a = in.extent(1); //cols of in
	//b = in.extent(0); //rows of in
	//err_checkf(a == b, "Input matrix must be square.", std::cout);
	//err_checkf(c == b, "Incompatible matrix dimensions for basis change.", std::cout);
	if (!forward) {
		dMatrix2 temp = dot<dMatrix2>(transformation, in, false, false); // temp = t * in
		dMatrix2 result = dot<dMatrix2>(temp, transformation, false, true); // res = temp * t^T
		return result;
	}
	else {
		dMatrix2 temp = dot<dMatrix2>(transformation, in, true, false); // temp = t^T * in
		dMatrix2 result = dot<dMatrix2>(temp, transformation, false, false); // res = temp * t
		return result;
	}
}

/**
 * Calculates Atomic Natural Orbitals for a specific atom.
 *
 * @param D_full      Pointer to the full Density Matrix (N_basis x N_basis)
 * @param S_full      Pointer to the full Overlap Matrix (N_basis x N_basis)
 * @param full_stride The leading dimension of the full matrices (usually N_basis)
 * @param atom_indices A vector containing the indices (0-based) of the basis functions for this atom
 * @return NAOResult containing sorted occupancies and coefficients
 */
Roby_information::NAOResult Roby_information::calculateAtomicNAO(const dMatrix2 &D_full,
	const dMatrix2 &S_full,
	const ivec &atom_indices,
	const ivec &shell_angular_momenta,
	const bool spherical,
	const double occupancy_cutoff,
	const int leading_orbitals_to_skip,
	const bool EVs,
	const int keep_orbitals) {

	err_checkf(D_full.extent(0) == D_full.extent(1), "Density matrix D must be square.", std::cout);
	err_checkf(S_full.extent(0) == S_full.extent(1), "Overlap matrix S must be square.", std::cout);
	err_checkf(D_full.extent(0) == S_full.extent(0), "Density and Overlap matrices must be of the same size.", std::cout);

	const int n = static_cast<int>(atom_indices.size());

	// 1. Memory Allocation for Submatrices
	// Using flat std::vectors to ensure contiguous memory for MKL
	vec D_sub(static_cast<size_t>(n) * n, 0.0); // called P in tonto
	vec S_sub(static_cast<size_t>(n) * n, 0.0); // called S in tonto
	get_submatrices(D_full, S_full, D_sub, S_sub, atom_indices);

	// Tonto averages the atomic density over all rotations; in a real-spherical basis that is exact. O_h keeps
	// axis-aligned anisotropy and is only the fallback for a Cartesian basis, whose shells mix l.
	// ponytail: a Cartesian shell of order l spans l, l-2, ...; exact averaging needs that decomposition first.
	if (!shell_angular_momenta.empty()) {
		dMatrix2 atomic_density = reshape<dMatrix2>(D_sub, Shape2D(n, n));
		if (spherical)
			spherically_average_atomic_matrix(atomic_density, shell_angular_momenta);
		else
			symmetrize_atomic_matrix_oh(atomic_density, shell_angular_momenta, spherical);
		D_sub = atomic_density.container();
	}
	vec Rho(static_cast<size_t>(n) * n);        // To store target density

	vec V = S_sub;
	vec W(n);
	// make V = Sqrt(S)
	const vec Temp = mat_sqrt(V, W);

#ifdef NSA2DEBUG
	print_dmatrix2(reshape<dMatrix2>(Temp, Shape2D(n, n)), "projection matrix V");
#endif

	const vec X = change_basis_sq(D_sub, Temp, n);
#ifdef NSA2DEBUG
	print_dmatrix2(reshape<dMatrix2>(X, Shape2D(n, n)), "projected density X");
#endif

	vec occu(n, 0);
	vec P = X;
#ifdef NSA2DEBUG
	err_checkf(isSymmetricViaEigenvalues<vec>(P, n), "Transformed matrix not symmetric!", std::cout);
#endif
	make_Eigenvalues(P, occu);
#ifdef NSA2DEBUG
	std::cout << "Eigenvalues of projected density P:\n";
	for (int i = 0; i < n; i++) {
		std::cout << std::setw(14) << std::setprecision(8) << std::fixed << W[i] << " ";
	}
	print_dmatrix2(reshape<dMatrix2>(P, Shape2D(n, n)), "Projected density P");
#endif

	for (int i = 0; i < n; i++)
		W[i] = abs(W[i]) < 1E-10 ? 0.0 : 1.0 / W[i];

	vec Temp2(static_cast<size_t>(n) * n, 0.0);
	double *T;
	int in, jn;
	for (int i = 0; i < n; i++) {
		in = i * n;
		for (int j = 0; j < n; j++) {
			jn = j * n;
			T = &Temp2[in + j];
			for (int k = 0; k < n; k++)
				*T += V[in + k] * W[k] * V[jn + k];
		}
	}

	V.clear();
	W.clear();

#ifdef NSA2DEBUG
	print_dmatrix2(reshape<dMatrix2>(Temp2, Shape2D(n, n)), "Back projection S^-0.5");
#endif

	cblas_dgemm(CblasRowMajor,
		CblasNoTrans, CblasNoTrans,
		n, n, n,
		1.0,
		Temp2.data(), n,
		P.data(), n,
		0.0,
		Rho.data(), n);

#ifdef NSA2DEBUG
	print_dmatrix2(reshape<dMatrix2>(Rho, Shape2D(n, n)), "resulting NAO");
#endif

	NAOResult result;
	result.eigenvalues = occu;
	result.eigenvectors = Rho; // Rho now contains eigenvectors
	result.matrix_elements = atom_indices;
	result.sub_DM = D_sub;
	result.sub_OM = S_sub;

	// Create an index vector to sort descending
	ivec idx(n);
	std::iota(idx.begin(), idx.end(), 0);

	//Deterministic tie-break on the index: +0.0/-0.0 ties would otherwise reorder run to run and change the kept subspace.
	std::sort(idx.begin(), idx.end(), [&](int i1, int i2) {
		const double a = result.eigenvalues[i1], b = result.eigenvalues[i2];
		if (a != b)
			return a > b;
		return i1 < i2;
		});

	// Reorder based on sorted indices
	vec sorted_evals;
	vec sorted_evecs;
	vec omitted_evals;
	vec omitted_evecs;
	sorted_evecs.reserve(static_cast<size_t>(n) * n);
	omitted_evecs.reserve(static_cast<size_t>(n) * n);


	const int skip_orbitals = std::clamp(leading_orbitals_to_skip, 0, n);
	//keep_orbitals >= 0 fixes the rank (the most occupied; eigenvalues are sorted) and ignores the threshold,
	//capped at the occupied eigenvectors: the free-atom null space is degenerate and padding into it is arbitrary.
	constexpr double null_occupation = 1E-8;
	int occupied_eigenvectors = 0;
	for (int i = skip_orbitals; i < n; i++)
		if (result.eigenvalues[idx[i]] > null_occupation)
			occupied_eigenvectors++;
	int keep = std::min(std::clamp(keep_orbitals, 0, std::max(0, n - skip_orbitals)),
		keep_orbitals >= 0 ? occupied_eigenvectors : n);
	//A rank boundary inside a degenerate set makes the projector arbitrary and breaks the molecule's symmetry,
	//so the rank extends to the end of the set. Compare with the last kept value, not the neighbour, so a slow
	//drift is not taken for one degenerate set.
	if (keep_orbitals >= 0 && keep > 0 && skip_orbitals + keep < n) {
		const double last = result.eigenvalues[idx[skip_orbitals + keep - 1]];
		const int rank_asked = keep;
		//relative, because these occupations run from 1E-8 to above 12 in the same list
		constexpr double degenerate_window = 1E-6;
		while (skip_orbitals + keep < n) {
			const double next = result.eigenvalues[idx[skip_orbitals + keep]];
			if (next <= null_occupation)
				break;
			if (std::abs(last - next) > degenerate_window * std::max(1.0, std::abs(last)))
				break;
			keep++;
		}
		if (keep != rank_asked)
			std::cout << "\n  NOTE: the atomic subspace of rank " << rank_asked << " would have cut through a "
			<< "degenerate occupation (" << std::setprecision(8) << last << "), which would have made this "
			<< "atom's projector depend on an arbitrary choice of directions inside that set and its bonds "
			<< "to symmetry-equivalent partners come out different. Rank extended to " << keep
			<< " so the set is kept whole.\n";
	}
	for (int i = 0; i < n; i++) {
		int original_idx = idx[i];
		const bool omit_orbital = i < skip_orbitals ||
			(keep_orbitals >= 0
				? i >= skip_orbitals + keep
				: (occupancy_cutoff >= 0.0 && result.eigenvalues[original_idx] < occupancy_cutoff));
		vec &target_evals = omit_orbital ? omitted_evals : sorted_evals;
		vec &target_evecs = omit_orbital ? omitted_evecs : sorted_evecs;

		target_evals.emplace_back(result.eigenvalues[original_idx]);

		for (int row = 0; row < n; row++) {
			target_evecs.emplace_back(result.eigenvectors[row * n + original_idx]);
		}
	}

	result.eigenvalues = sorted_evals;
	result.eigenvectors = sorted_evecs;
	result.omitted_eigenvalues = omitted_evals;
	result.omitted_eigenvectors = omitted_evecs;

	//Printed after the split: shows where the rank boundary falls.
	if (EVs) {
		std::cout << "Occupations of the projected density P, rank " << sorted_evals.size()
			<< " of " << n << ":\n";
		for (size_t i = 0; i < sorted_evals.size(); i++)
			std::cout << std::setw(14) << std::setprecision(8) << std::fixed << sorted_evals[i] << "  kept\n";
		for (size_t i = 0; i < omitted_evals.size(); i++)
			std::cout << std::setw(14) << std::setprecision(8) << std::fixed << omitted_evals[i] << "  omitted\n";
	}

	return result;
}

/**
 * Computes Final NAOs by symmetrically orthogonalizing the Pre-NAOs.
 *
 * @param C_PNAO    Global matrix of Pre-NAOs (N x N). Columns are orbitals.
 * @param S_AO      Original Overlap Matrix in AO basis (N x N).
 * @param n         Number of basis functions (N).
 * @return          Matrix of final NAOs (N x N), Columns are orbitals.
 */

 /* CURRENTLY NOT IN USE
 vec orthogonalizePNAOs(const vec& C_PNAO,
	 const vec& S_AO,
	 int n) {

	 vec S_PNAO(n * n);
	 vec Temp(n * n); // Intermediate buffer
	 vec C_NAO(n * n); // Result

	 // 1. Compute Overlap in PNAO basis: S_PNAO = C_PNAO^T * S_AO * C_PNAO

	 // Step A: Temp = S_AO * C_PNAO
	 // S_AO is symmetric.
	 cblas_dsymm(CblasRowMajor,
		 CblasLeft, CblasUpper,
		 n, n,
		 1.0, S_AO.data(), n,
		 C_PNAO.data(), n,
		 0.0, Temp.data(), n);

	 // Step B: S_PNAO = C_PNAO^T * Temp
	 // C_PNAO is not symmetric, so we use dgemm.
	 // Transpose the first matrix (C_PNAO^T).
	 cblas_dgemm(CblasRowMajor,
		 CblasTrans, CblasNoTrans,
		 n, n, n,
		 1.0, C_PNAO.data(), n,
		 Temp.data(), n,
		 0.0, S_PNAO.data(), n);

	 // 2. Compute S_PNAO^(-1/2) using Eigendecomposition
	 // S_PNAO is symmetric (and positive definite).

	 vec W(n); // Eigenvalues
	 // We can overwrite S_PNAO with eigenvectors to save memory,
	 // but let's keep it clear. Copy S_PNAO to 'U' (Eigenvectors).
	 vec U = S_PNAO;

	 // LAPACKE_dsyevd: Computes all eigenvalues and eigenvectors
	 err_checkf(LAPACKE_dsyevd(LAPACK_ROW_MAJOR, 'V', 'U', n, U.data(), n, W.data()) == 0, "Eigenvalue computation failed.", std::cout);

	 // Construct S^-1/2 = U * Lambda^(-1/2) * U^T
	 std::fill(S_PNAO.begin(), S_PNAO.end(), 0.0); // Reuse S_PNAO to store S^-1/2

	 for (int k = 0; k < n; k++) {
		 double scale = 1.0 / std::sqrt(W[k]);
		 // Add contribution of k-th eigenvector: scale * (v_k * v_k^T)
		 // v_k is the k-th ROW of U.
		 const double* v_k = &U[k * n];

		 // This is a rank-1 update (dger), but we can just sum manually or loop
		 // Since we need the full matrix for the next multiplication
		 for (int i = 0; i < n; i++) {
			 for (int j = 0; j < n; j++) {
				 S_PNAO[i * n + j] += scale * v_k[i] * v_k[j];
			 }
		 }
	 }

	 // 3. Transform PNAOs to NAOs
	 // C_NAO = C_PNAO * S^(-1/2)
	 // C_PNAO (Columns are PNAOs) * S_inv_sqrt (transformation matrix)

	 cblas_dgemm(CblasRowMajor,
		 CblasNoTrans, CblasNoTrans,
		 n, n, n,
		 1.0, C_PNAO.data(), n,
		 S_PNAO.data(), n, // This is now S^-1/2
		 0.0,
		 C_NAO.data(), n);

	 return C_NAO;
 }
 */

double Roby_information::projection_matrix_and_expectation(const ivec &indices, const ivec &eigvals, const ivec &eigvecs, dMatrix2 *given_NAO, dMatrix2 *proj_out) {
	const int n = indices.size();
	//vec D_Sub(n * n, 0.0);
	vec S_Sub(static_cast<size_t>(n) * n, 0.0);
	get_submatrix(overlap_matrix, S_Sub, indices);
	dMatrix2 S = reshape<dMatrix2>(S_Sub, Shape2D(n, n));
	int atom = -1;
	//dMatrix2 D = reshape<dMatrix2>(D_Sub, Shape2D(n, n));

	auto zero_projection = [&]() {
		dMatrix2 zero(n, n);
		for (int r = 0; r < n; ++r)
			for (int c = 0; c < n; ++c)
				zero(r, c) = 0.0;
		if (atom >= 0 && proj_out == nullptr) {
			projection_matrices.push_back(zero);
			overlap_matrices.push_back(S);
		}
		if (proj_out != nullptr)
			*proj_out = zero;
		return 0.0;
		};

	dMatrix2 NAOs;
	// When an explicit NAO matrix is provided together with row/col selectors,
	// use it directly — this takes priority over the full-system shortcut below.
	if (given_NAO != nullptr && !eigvals.empty() && !eigvecs.empty()) {
		const int n1 = eigvals.size();
		const int n2 = eigvecs.size();
		vec NAO_sub(static_cast<size_t>(n1) * n2);
		get_submatrix(*given_NAO, NAO_sub, eigvals, eigvecs);
		NAOs = reshape<dMatrix2>(NAO_sub, Shape2D(n1, n2));
	}
	else if (given_NAO != nullptr && eigvals.empty() && eigvecs.empty()) {
		NAOs = *given_NAO;
	}
	else if (given_NAO == nullptr && eigvals.empty() && eigvecs.empty()
		&& n == static_cast<int>(overlap_matrix.extent(0))) {
		err_checkf(static_cast<int>(total_NAOs.extent(1)) == n,
			"RGBI total NAO matrix column count (" + std::to_string(total_NAOs.extent(1)) +
			") does not match the basis-function count (" + std::to_string(n) + ").",
			std::cout);
		NAOs = total_NAOs;
	}
	else {
		//TODO: assign subspace NAOs from NAOResults for a given atom
		for (auto NAO : this->NAOs) {
			// if matrix_elemnts of this NAO are identical to indices, resahpe NAO.eigenvectors to the correct shape
			if (NAO.matrix_elements == indices) {
				NAOs = reshape<dMatrix2>(NAO.eigenvectors, Shape2D(NAO.eigenvalues.size(), n));
				atom = NAO.atom_index;
				break;
			}
		}
	}
	//For atom groups we can use submatrices of total_NAOs
	if (NAOs.size() == 0) {
		const int n1 = eigvals.size();
		const int n2 = eigvecs.size();
		if (n1 == 0 || n2 == 0)
			return zero_projection();
		vec NAO_sub(static_cast<size_t>(n1) * n2);
		if (given_NAO == nullptr)
			given_NAO = &total_NAOs;
		get_submatrix(*given_NAO, NAO_sub, eigvals, eigvecs);
		NAOs = reshape<dMatrix2>(NAO_sub, Shape2D(n1, n2));
	}
#ifdef NSA2DEBUG
	print_dmatrix2(NAOs, "W in projection matrix making");
#endif

	auto X = change_basis_general(S, transpose(NAOs));
#ifdef NSA2DEBUG
	print_dmatrix2(X, "new Basis");
#endif

	PinvRank rank{};
	auto Y = LAPACKE_invert(X, pinv_cutoff, &rank);
	last_pinv_n = rank.n;
	last_pinv_kept = rank.kept;
	last_pinv_smallest_kept = rank.smallest_kept;
	last_pinv_largest_dropped = rank.largest_dropped;
#ifdef NSA2DEBUG
	print_dmatrix2(Y, "Pseudo inverse of Y");
#endif

	//making the projection matrix
	X = change_basis_general(Y, transpose(NAOs), false);
#ifdef NSA2DEBUG
	print_dmatrix2(X, "Backtransformed X");
#endif

	if (atom >= 0 && proj_out == nullptr) {
		projection_matrices.push_back(X);
		overlap_matrices.push_back(S);
	}
	if (proj_out != nullptr)
		*proj_out = X;

	dMatrix2 S_rect = get_rectangle(overlap_matrix, indices);
#ifdef NSA2DEBUG
	print_dmatrix2(S_rect, "S_rect:");
#endif

	//overlap projection
	auto W = change_basis_general(X, S_rect);
#ifdef NSA2DEBUG
	print_dmatrix2(W, "Overlap transformed W");
#endif

	const double expect = trace_product<double>(W, density_matrix);
	return expect;
}

double Roby_information::Roby_population_analysis(const ivec atoms) {
	ivec bf_indices;
	if (atoms.size() == 0) {
		const int n_basis_functions = static_cast<int>(overlap_matrix.extent(0));
		bf_indices.reserve(n_basis_functions);
		for (int index = 0; index < n_basis_functions; ++index)
			bf_indices.push_back(index);
	}
	else {
		bf_indices = atoms;
	}
	double P = projection_matrix_and_expectation(bf_indices);
	return P;
}

void Roby_information::computeAllAtomicNAOs(WFN &wavy, const bool symmetrize, const bool use_ano_basis, const bool EVs, const bool legacy_occupancy_cutoff) {
	const int N_atoms = wavy.get_ncen();
	const std::vector<atom> ats = wavy.get_atoms();
	NAOs.reserve(N_atoms);
	ano_fallback_atoms.clear();

	density_matrix = wavy.get_dm();
	//Every index below reads DM; the .fchk readers keep the density only in UT_DensityMatrix and leave DM empty.
	err_checkf(density_matrix.extent(0) > 0 && density_matrix.extent(1) == density_matrix.extent(0),
		"RGBI needs the density matrix over the contracted basis, and " + wavy.get_path().filename().string() +
		" carries none: its reader stores the density in triangular form only. Use " +
		rgbi_supported_input_phrase(wavy.get_path().extension().string()) + " of the same calculation.",
		std::cout);

	//Refuse before the overlap and free-atom SCFs are paid for. Only the Cartesian route has the O_h transforms'
	//l <= 5 ceiling; symmetrize_atomic_matrix_oh() keeps its own check for other callers.
	if (symmetrize && wavy.get_d_f_switch()) {
		const int highest = highest_shell_angular_momentum(wavy);
		err_checkf(highest <= 5,
			"RGBI's atomic O_h symmetrization supports shells from s through h, and the Cartesian basis of " +
			wavy.get_path().filename().string() + " carries l = " + std::to_string(highest) +
			" (" + std::string(1, "spdfghiklm"[std::min(highest, 9)]) + " shells). A spherical basis has no "
			"such limit, because it is averaged over all rotations instead. -rgbi_no_sym analyses this file "
			"without the averaging, but expect it to be slow: on the 670-function i-shell file in the test "
			"set it produced no output in 1800 s.",
			std::cout);
	//A Cartesian shell mixes l, so the reference falls back to the O_h average, which keeps axis-aligned
	//anisotropy: results depend on the molecule's orientation. -rgbi_no_sym or a spherical basis avoids it.
		std::cout << "Warning: " << wavy.get_path().filename().string() << " uses a Cartesian basis, so the "
			"atomic reference is averaged over O_h instead of over all rotations and depends on the "
			"orientation of the molecule in the file.\n";
	}

	//On the ANO route over a spherical basis the free-atom reference is averaged regardless of -rgbi_no_sym;
	//the switch still applies to atoms that fall back to the molecular density.
	if (!symmetrize && use_ano_basis && !wavy.get_d_f_switch())
		std::cout << "Note: -rgbi_no_sym does not change the atomic reference on the ANO route, because "
			"that reference is a free atom and a free atom is spherically symmetric. It still applies to "
			"any atom whose ANO reference falls back to the molecular density.\n";

	if (wavy.get_d_f_switch()) {
		Int_Params basis(wavy);
		vec S_full;
		compute2C<Overlap2C_CRT>(basis, S_full);
		overlap_matrix = reshape<dMatrix2>(S_full, Shape2D(density_matrix.extent(0), density_matrix.extent(1)));
	}
	else {
		//ORCA-convention densities (gbw, or molden from one) flip the sign of |m| >= 3 components; ao_overlap
		//matches it, a plain Overlap2C_SPH would be the wrong metric for f and higher shells.
		overlap_matrix = ao_overlap(wavy);
	}

#ifdef NSA2DEBUG
	print_dmatrix2(overlap_matrix, "Overlap matrix");
#endif

	//err_checkf()

	//Legacy rule only: thresholded occupations move with geometry, so the projector rank and every index
	//step when one crosses. By default the rank is fixed per element (free_atom_orbital_count).
	const double occupancy_cutoff = use_ano_basis ? 1.0 / 14.0 : 1.0 / 6.0;

	if (use_ano_basis)
		warm_free_atom_cache(ats, wavy.get_origin(), wavy.get_d_f_switch());

	int last_index = 0;
	ivec2 indices(wavy.get_ncen());
	for (auto &a : ats) {
		indices[a.get_nr() - 1].reserve(density_matrix.extent(0) / N_atoms); // Rough estimate
		ivec shell_angular_momenta;
		int current_shell = -1;
		int nr_indices = 0;
		std::vector<basis_set_entry> basis_set = a.get_basis_set();
		for (auto &bf : basis_set) {
			if (bf.get_shell() != current_shell) {
				current_shell++;
				const int l = static_cast<int>(bf.get_type()) - 1;
				err_checkf(l >= 0,
					"Encountered an invalid shell angular momentum while building RGBI atomic NAOs.",
					std::cout);
				shell_angular_momenta.push_back(l);
				nr_indices = atomic_shell_size(l, wavy.get_d_f_switch());
				if (wavy.get_origin() == e_origin::tonto && bf.get_type() == 3) {
					// 2 - 4; 3 - 6; 5 - 6
					swap_rows_cols_symm(overlap_matrix, last_index + 1, last_index + 3);
					swap_rows_cols_symm(overlap_matrix, last_index + 2, last_index + 5);
					swap_rows_cols_symm(overlap_matrix, last_index + 4, last_index + 5);
				}

				//if (wavy.get_origin() != e_origin::tonto && wavy.get_origin() != e_origin::OCC)
				//    bf.get_type() == 1 ? nr_indices = 1 : (bf.get_type() == 2 ? nr_indices = 3 : (bf.get_type() == 3 ? nr_indices = 5 : nr_indices = 7));
				//else {
				//
				//    bf.get_type() == 1 ? nr_indices = 1 : (bf.get_type() == 2 ? nr_indices = 3 : (bf.get_type() == 3 ? nr_indices = 6 : nr_indices = 10));
				//    if (bf.get_type() == 3 && wavy.get_origin() == e_origin::tonto) {
				//        //2 - 4; 3 - 6; 5 - 6
				//        swap_rows_cols_symm(overlap_matrix, last_index + 1, last_index + 3);
				//        swap_rows_cols_symm(overlap_matrix, last_index + 2, last_index + 5);
				//        swap_rows_cols_symm(overlap_matrix, last_index + 4, last_index + 5);
				//    }
				//}
				for (int i = 0; i < nr_indices; i++) {
					indices[a.get_nr() - 1].push_back(last_index);
					last_index++;
				}
			}
		}

		//Shells adding up to more rows than the density matrix would index past its end; the shell description is at fault.
		err_checkf(last_index <= static_cast<int>(density_matrix.extent(0)),
			"The basis of atom " + std::to_string(a.get_nr()) + " (" + a.get_label() + ") describes " +
			std::to_string(last_index) + " basis functions by its shells, more than the " +
			std::to_string(density_matrix.extent(0)) + " the density matrix has. RGBI cannot index a "
			"matrix it has been given a wrong shell layout for.",
			std::cout);

		// The GBW reader converts ORCA components to PySCF/libcint order and
		// groups an atom's shells by increasing angular momentum.  Mirror that
		// layout here so each symmetry block describes the corresponding DM rows.
		if (wavy.get_origin() == e_origin::gbw)
			std::stable_sort(shell_angular_momenta.begin(), shell_angular_momenta.end());

		const bool spherical = !wavy.get_d_f_switch();
		const int keep_orbitals = legacy_occupancy_cutoff
			? -1
			: free_atom_orbital_count(a.get_charge()) - ecp_core_orbital_count(a.get_charge(), a.get_ECP_electrons());

		auto make_molecular_fallback = [&]() {
			auto fallback = calculateAtomicNAO(density_matrix, overlap_matrix,
				indices[a.get_nr() - 1],
				symmetrize ? shell_angular_momenta : ivec{},
				spherical,
				occupancy_cutoff,
				0,
				EVs,
				keep_orbitals);
			fallback.atom_index = a.get_nr() - 1;
			return fallback;
		};
		if (use_ano_basis) {
			try {
				const dMatrix2 atomic_density =
					compute_tonto_style_atomic_density(a, wavy.get_origin(), wavy.get_d_f_switch());
				const int n_local = static_cast<int>(indices[a.get_nr() - 1].size());
				if (static_cast<int>(atomic_density.extent(0)) != n_local ||
					static_cast<int>(atomic_density.extent(1)) != n_local) {
					throw std::runtime_error(
						"Atomic OCC density size does not match the loaded WFN basis layout.");
				}
				ivec local_indices(indices[a.get_nr() - 1].size());
				std::iota(local_indices.begin(), local_indices.end(), 0);
				vec S_sub(static_cast<size_t>(n_local) * n_local, 0.0);
				get_submatrix(overlap_matrix, S_sub, indices[a.get_nr() - 1]);
				const dMatrix2 atomic_overlap = reshape<dMatrix2>(S_sub, Shape2D(n_local, n_local));
				//Same fixed count as the molecular route, so both routes pick the same subspace. The matrix is a free
				//atom's density, spherical by symmetry; an open-shell single-determinant SCF breaks that and would carry the
				//file's axes into every bond of the atom, so a spherical basis is always averaged, whatever symmetrize says.
				//A Cartesian basis would fall back to O_h (l <= 5, orientation dependent), so there it follows symmetrize.
				auto ano = calculateAtomicNAO(atomic_density, atomic_overlap,
					local_indices,
					(symmetrize || spherical) ? shell_angular_momenta : ivec{},
					spherical,
					occupancy_cutoff,
					0,
					EVs,
					keep_orbitals);
				ano.sub_OM = S_sub;
				ano.sub_DM = atomic_density.container();
				ano.matrix_elements = indices[a.get_nr() - 1];
				ano.atom_index = a.get_nr() - 1;
				NAOs.emplace_back(std::move(ano));
			}
			catch (const std::exception &e) {
				std::cout << "Warning: RGBI ANO state build failed for atom "
					<< a.get_nr() - 1 << " (" << a.get_label()
					<< "). Falling back to molecular local orbitals. Reason: "
					<< e.what() << '\n';
				ano_fallback_atoms.push_back(a.get_nr() - 1);
				NAOs.emplace_back(make_molecular_fallback());
			}
			catch (...) {
				std::cout << "Warning: RGBI ANO state build failed for atom "
					<< a.get_nr() - 1 << " (" << a.get_label()
					<< "). Falling back to molecular local orbitals.\n";
				ano_fallback_atoms.push_back(a.get_nr() - 1);
				NAOs.emplace_back(make_molecular_fallback());
			}
		}
		else {
			NAOs.emplace_back(make_molecular_fallback());
		}
		NAOs.back().atom_index = a.get_nr() - 1;
	}

	//Fewer indices than rows leaves basis functions in no atom's subspace and silently lowers every Roby index.
	err_checkf(last_index == static_cast<int>(density_matrix.extent(0)),
		"The atoms' shells account for " + std::to_string(last_index) + " basis functions but the "
		"density matrix has " + std::to_string(density_matrix.extent(0)) + ". RGBI would leave the "
		"difference in no atom's subspace; the shell description of this wavefunction is incomplete.",
		std::cout);
#ifdef NSA2DEBUG
	print_dmatrix2(overlap_matrix, "Overlap matrix repaired");
#endif
}

ivec Roby_information::find_eigenvalue_pairs(const vec &eigvals, const double tolerance) {
	const int n = eigvals.size();
	ivec pairs(n, -1);
	for (int i = 0; i < n; i++) {
		if (pairs[i] >= 0)
			continue; // already paired
		for (int j = 0; j < n; j++) {
			if (pairs[j] >= 0)
				continue; // already paired
			if (abs(abs(eigvals[j]) - 1.0) < tolerance)
				continue; // skip values close ot +/- one, as they reside only on one atom
			if (abs(eigvals[i] + eigvals[j]) < tolerance) {
				pairs[i] = j;
				pairs[j] = i;
				break;
			}
		}
		if (pairs[i] == -1)
			pairs[i] = i;
	}
	return pairs;
}

void Roby_information::transform_Ionic_eigenvectors_to_Ionic_orbitals(
	dMatrix2 &EVC,
	const vec &eigvals,
	const ivec &pairs,
	const int index_a,
	const int index_b,
	const ivec &pair_matrix_indices)
{
	double fp, fm, fa, fb, s, c, s2;
	const int n_ab = EVC.extent(0);
	const int n_eigvals = eigvals.size();
	const int n_a = projection_matrices[index_a].extent(0);
	const int n_b = projection_matrices[index_b].extent(0);
	err_checkf(n_ab == n_a + n_b, "Inconsitent size in projection matrices?!", std::cout);
	err_checkf(n_eigvals == EVC.extent(1), "Inconsitency between EVC and eigvals", std::cout);

	dMatrix1 A(n_a), B(n_b);
	dMatrix1 EVC_column(n_ab);
	dMatrix2 PAS(n_a, n_ab), PBS(n_b, n_ab);

	const ivec *indices_a = nullptr;
	const ivec *indices_b = nullptr;
	for (const auto &NAO : NAOs) {
		if (NAO.atom_index == index_a)
			indices_a = &NAO.matrix_elements;
		if (NAO.atom_index == index_b)
			indices_b = &NAO.matrix_elements;
	}
	err_checkf(indices_a != nullptr, "No NAO data found for atom " + std::to_string(index_a), std::cout);
	err_checkf(indices_b != nullptr, "No NAO data found for atom " + std::to_string(index_b), std::cout);

	vec Sub_overlap(static_cast<size_t>(n_a) * n_ab);
	get_submatrix(overlap_matrix, Sub_overlap, *indices_a, pair_matrix_indices);
	dMatrix2 Sa = reshape<dMatrix2>(Sub_overlap, Shape2D(n_a, n_ab));
	PAS = dot<dMatrix2>(projection_matrices[index_a], Sa, false, false);
	Sub_overlap.clear(); Sub_overlap.resize(static_cast<size_t>(n_b) * n_ab);
	get_submatrix(overlap_matrix, Sub_overlap, *indices_b, pair_matrix_indices);
	dMatrix2 Sb = reshape<dMatrix2>(Sub_overlap, Shape2D(n_b, n_ab));
	PBS = dot<dMatrix2>(projection_matrices[index_b], Sb, false, false);
#ifdef NSA2DEBUG
	print_dmatrix2(PAS, "PAS");
	print_dmatrix2(PBS, "PBS");
#endif // NSA2DEBUG


	for (int i = 0; i < n_eigvals; i++) {
		if (pairs[i] < 0) continue;
		if (pairs[i] == i) continue;
		if (eigvals[i] < eigvals[pairs[i]]) continue;
#ifdef NSA2DEBUG
		std::cout << "Doing i=" << i << std::endl;
#endif
		for (int a = 0; a < n_ab; a++)
			EVC_column(a) = EVC(a, i);

#ifdef NSA2DEBUG
		print_dmatrix2(reshape<dMatrix2>(EVC_column, Shape2D(n_ab, 1)), "slice used");
#endif

		s = eigvals[i];
		s2 = s * s;
		if (abs(s2 - 1.0) < 1E-8) c = 0.0;
		else c = sqrt(1.0 - s2);

		fm = sqrt(1 - c) / s;
		fp = sqrt(1 + c) / s;
		fa = 0.5 * ((fm + fp) + c * (fm - fp));
		fb = 0.5 * (c * (fm + fp) + (fm - fp));

		//A and B outlive the loop, so project every pair even when the scale factor is 1.
		A = dot_BLAS<dMatrix1, dMatrix2>(PAS, EVC_column, false);
		if (abs(fa - 1.0) > 1E-8)
			for (int a = 0; a < n_a; a++)
				A(a) /= fa;
		B = dot_BLAS<dMatrix1, dMatrix2>(PBS, EVC_column, false);
		if (abs(fb - 1.0) > 1E-8)
			for (int b = 0; b < n_b; b++)
				B(b) /= fb;

#ifdef NSA2DEBUG
		std::cout << "fa: " << fa << std::endl << "fb: " << fb << std::endl;
		print_dmatrix2(reshape<dMatrix2>(A, Shape2D(n_a, 1)), "A");
		print_dmatrix2(reshape<dMatrix2>(B, Shape2D(n_b, 1)), "B");
#endif

		//build antibonding state in pairs[i]
		fa = 0.5 * (fm - fp);
		fb = 0.5 * (fm + fp);

		for (int a = 0; a < n_a; a++) {
			EVC(a, pairs[i]) = fa * A(a);
		}
		for (int b = n_a; b < n_ab; b++) {
			EVC(b, pairs[i]) = fb * B(b - n_a);
		}
	}
}

std::map<char, dMatrix2> Roby_information::make_covalent_from_ionic(
	const dMatrix2 &theta_I,
	const vec &eigvals,
	const ivec &pairs,
	bool EVs) {
	std::map<char, dMatrix2> res; // A = angle, V = eigen_value, T = Theta_vector
	const int size = eigvals.size();
	res.emplace('A', dMatrix2(size, 1));
	res.emplace('V', dMatrix2(size, 1));
	res.emplace('T', dMatrix2(theta_I.extent(0), theta_I.extent(1)));

	for (int i = 0; i < size; i++) {
		res['A'](i, 0) = 90.0;
		res['V'](i, 0) = 0.0;
	}

	if (EVs) {
		std::cout << "Eigenvalues of projected density P (in covalent representation):\n";
		for (int i = 0; i < size; i++) {
			std::cout << std::setw(14) << std::setprecision(8) << std::fixed << eigvals[i] << "\n";
		}
	}

	for (int val = 0; val < size; val++) {
		if (pairs[val] < 0) continue;
		if (pairs[val] == val) continue;
		const double s = eigvals[val];
		if (s < eigvals[pairs[val]]) continue;

		for (int i = 0; i < theta_I.extent(0); i++) {
			res['T'](i, val) = constants::INV_SQRT2 * (theta_I(i, val) + theta_I(i, pairs[val]));
			res['T'](i, pairs[val]) = constants::INV_SQRT2 * (theta_I(i, val) - theta_I(i, pairs[val]));
		}

		const double c = sqrt(1.0 - s * s);
		res['V'](val, 0) = c;
		res['V'](pairs[val], 0) = -c;

		res['A'](val, 0) = atan2(s, c) * constants::INV_PI_180;
	}

	if (EVs) {
		std::cout << "Eigenvalues of projected density P (after transformation) and matching Theta angles:\n";
		for (int i = 0; i < size; i++) {
			std::cout << std::setw(14) << std::setprecision(8) << std::fixed << res['V'](i, 0);
			for(int j = 0; j < res['T'].extent(1); j++) {
				std::cout << std::setw(14) << std::setprecision(8) << std::fixed << res['T'](i, j) << "\n";
			}
		}
	}

	return res;
}

std::string Roby_information::make_theta_info(const WFN &wavy, const std::pair<int, int> &bond,
	const vec &eigvals, const ivec &pairs, const dMatrix2 &angles,
	const vec &covalent_populations, const vec &ionic_populations) const {
	const auto &atoms = wavy.get_atoms();
	std::ostringstream out;
	out << "\nRoby-Gould theta subspaces for "
		<< atoms[bond.first].get_label() << " (A) - " << atoms[bond.second].get_label() << " (B)\n";
	out << " Pair    theta/degrees       C+       C-      Cov.       I+       I-      Ion.     Total\n";
	out << "----------------------------------------------------------------------------------------------\n";
	for (int i = 0; i < static_cast<int>(eigvals.size()); ++i) {
		if (pairs[i] < 0 || eigvals[i] < eigvals[pairs[i]])
			continue;
		const int pair = pairs[i];
		const double c_plus = covalent_populations[i];
		const double c_minus = covalent_populations[pair];
		const double covalent = pair == i ? 0.0 : 0.5 * (c_plus - c_minus);
		double i_plus = ionic_populations[i], i_minus = ionic_populations[pair], ionic = 0.0;
		if (pair == i) {
			if (eigvals[i] < 0.0) {
				i_minus = i_plus;
				i_plus = 0.0;
				ionic = -0.5 * i_minus;
			}
			else {
				i_minus = 0.0;
				ionic = 0.5 * i_plus;
			}
		}
		else
			ionic = 0.5 * (i_plus - i_minus);
		const double total = std::sqrt(covalent * covalent + ionic * ionic);
		out << std::fixed << std::setprecision(3)
			<< std::setw(4) << i + 1 << "," << std::setw(3) << pair + 1
			<< std::setw(16) << angles(i, 0)
			<< std::setw(9) << c_plus << std::setw(9) << c_minus
			<< std::setw(10) << covalent << std::setw(9) << i_plus
			<< std::setw(9) << i_minus << std::setw(10) << ionic
			<< std::setw(10) << total << '\n';
	}
	out << "----------------------------------------------------------------------------------------------\n";
	return out.str();
}

// Assembles the block-projected PAS matrix for a group: each atom's projection matrix
// is multiplied by the rectangular overlap block (atom bf rows × bond bf cols) and
// the results are stacked vertically.  Result shape: n_bf_G × n_ab.
dMatrix2 Roby_information::build_group_PAS(const ivec &group_atom_indices, const ivec &bond_bf_indices, const dMatrix2 &P_G) {
	const int n_ab = static_cast<int>(bond_bf_indices.size());
	// Collect group BF indices in the same sorted-atom order as P_G was built
	ivec group_bf;
	for (int ai : group_atom_indices)
		for (const auto &NAO : NAOs)
			if (NAO.atom_index == ai)
				for (int idx : NAO.matrix_elements)
					group_bf.push_back(idx);
	const int n_G = static_cast<int>(group_bf.size());
	err_checkf(n_G > 0, "RGBI group PAS has no basis functions for the requested atom group.", std::cout);
	err_checkf(n_ab > 0, "RGBI group PAS has no pair basis functions.", std::cout);
	vec sub_S(static_cast<size_t>(n_G) * n_ab, 0.0);
	get_submatrix(overlap_matrix, sub_S, group_bf, bond_bf_indices);
	dMatrix2 S_G = reshape<dMatrix2>(sub_S, Shape2D(n_G, n_ab));
	// PAS = P_G (n_G × n_G) × S_G (n_G × n_ab) → (n_G × n_ab)
	return dot<dMatrix2>(P_G, S_G, false, false);
}

// Group-aware version of transform_Ionic_eigenvectors_to_Ionic_orbitals.
// Builds PAS/PBS from block-assembled group projection matrices, then applies
// the same bonding/antibonding rotation as the single-atom version.
void Roby_information::transform_group_Ionic_orbitals(
	dMatrix2 &EVC,
	const vec &eigvals,
	const ivec &pairs,
	const ivec &group_a_atoms,
	const ivec &group_b_atoms,
	const ivec &bond_bf_indices,
	const dMatrix2 &P_GA,
	const dMatrix2 &P_GB)
{
	const int n_ab = EVC.extent(0);
	const int n_eigvals = static_cast<int>(eigvals.size());

	dMatrix2 PAS = build_group_PAS(group_a_atoms, bond_bf_indices, P_GA);
	dMatrix2 PBS = build_group_PAS(group_b_atoms, bond_bf_indices, P_GB);

	const int n_a = static_cast<int>(PAS.extent(0));
	const int n_b = static_cast<int>(PBS.extent(0));
	err_checkf(n_ab == n_a + n_b, "Inconsistent group sizes in transform_group_Ionic_orbitals", std::cout);

	dMatrix1 A(n_a), B(n_b);
	dMatrix1 EVC_column(n_ab);

	for (int i = 0; i < n_eigvals; i++) {
		if (pairs[i] < 0) continue;
		if (pairs[i] == i) continue;
		if (eigvals[i] < eigvals[pairs[i]]) continue;

		for (int a = 0; a < n_ab; a++)
			EVC_column(a) = EVC(a, i);

		const double s = eigvals[i];
		const double s2 = s * s;
		const double c = (abs(s2 - 1.0) < 1E-8) ? 0.0 : sqrt(1.0 - s2);

		const double fm = sqrt(1.0 - c) / s;
		const double fp = sqrt(1.0 + c) / s;
		double fa = 0.5 * ((fm + fp) + c * (fm - fp));
		double fb = 0.5 * (c * (fm + fp) + (fm - fp));

		//A and B outlive the loop body: project every pair, even when the scale factor is 1.
		A = dot_BLAS<dMatrix1, dMatrix2>(PAS, EVC_column, false);
		if (abs(fa - 1.0) > 1E-8)
			for (int a = 0; a < n_a; a++)
				A(a) /= fa;
		B = dot_BLAS<dMatrix1, dMatrix2>(PBS, EVC_column, false);
		if (abs(fb - 1.0) > 1E-8)
			for (int b = 0; b < n_b; b++)
				B(b) /= fb;

		// Build antibonding partner in column pairs[i]
		fa = 0.5 * (fm - fp);
		fb = 0.5 * (fm + fp);

		for (int a = 0; a < n_a; a++)
			EVC(a, pairs[i]) = fa * A(a);
		for (int b = n_a; b < n_ab; b++)
			EVC(b, pairs[i]) = fb * B(b - n_a);
	}
}

void Roby_information::computeGroupAnalysis(const ivec2 &group_defs, const vec &atom_pops, const ivec &atom_charges, const bool EVs) {
	RGBI_groups.clear();
	const int n_groups = static_cast<int>(group_defs.size());

	// Helper: build NAO row/col index sets for a list of sorted atom indices
	auto collect_nao_indices = [&](const ivec &atoms, ivec &eigenvals_out, ivec &eigenvecs_out) {
		eigenvals_out.clear();
		eigenvecs_out.clear();
		int start_val = 0;
		for (const auto &NAO : NAOs) {
			bool wanted = std::find(atoms.begin(), atoms.end(), NAO.atom_index) != atoms.end();
			if (wanted) {
				for (int i = 0; i < static_cast<int>(NAO.eigenvalues.size()); i++)
					eigenvals_out.push_back(start_val + i);
				for (int idx : NAO.matrix_elements)
					eigenvecs_out.push_back(idx);
			}
			start_val += static_cast<int>(NAO.eigenvalues.size());
		}
		};

	// Helper: collect BF indices for a list of atom indices in the given order.
	// Caller must pass atoms in sorted order so BF ordering matches the ionic-operator layout.
	// Scans NAOs by atom_index rather than using direct array indexing, so the result is
	// correct regardless of the order atoms were processed in computeAllAtomicNAOs.
	auto collect_bf_indices = [&](const ivec &atoms, ivec &bf_out) {
		bf_out.clear();
		for (int ai : atoms)
			for (const auto &NAO : NAOs)
				if (NAO.atom_index == ai)
					for (int idx : NAO.matrix_elements)
						bf_out.push_back(idx);
		};

	auto atom_list_string = [](const ivec &atoms) {
		std::string result;
		for (int i = 0; i < static_cast<int>(atoms.size()); ++i) {
			if (i > 0)
				result += ",";
			result += std::to_string(atoms[i]);
		}
		return result;
		};

	// Pre-compute per-group data: sorted atoms, BF indices, NAO indices, projection, population
	struct GroupData {
		ivec sorted_atoms;
		ivec bf;
		ivec nao_evals;
		ivec nao_evecs;
		dMatrix2 P;
		double pop = 0.0;
		std::string elem_list; // e.g. "N,H,H,H"
	};
	std::vector<GroupData> gdata(n_groups);
	for (int gi = 0; gi < n_groups; gi++) {
		gdata[gi].sorted_atoms = group_defs[gi];
		std::sort(gdata[gi].sorted_atoms.begin(), gdata[gi].sorted_atoms.end());
		collect_bf_indices(gdata[gi].sorted_atoms, gdata[gi].bf);
		collect_nao_indices(gdata[gi].sorted_atoms, gdata[gi].nao_evals, gdata[gi].nao_evecs);
		err_checkf(!gdata[gi].bf.empty(),
			"RGBI group G" + std::to_string(gi) + " has no basis functions for atoms " +
			atom_list_string(gdata[gi].sorted_atoms) + ".", std::cout);
		gdata[gi].pop = projection_matrix_and_expectation(
			gdata[gi].bf, gdata[gi].nao_evals, gdata[gi].nao_evecs, nullptr, &gdata[gi].P);
		{
			std::map<int, int> charge_count;
			for (int ai : gdata[gi].sorted_atoms)
				charge_count[atom_charges[ai]]++;
			std::vector<std::pair<int, int>> by_charge(charge_count.begin(), charge_count.end());
			std::sort(by_charge.begin(), by_charge.end(), [](const auto &a, const auto &b) {
				return a.first > b.first; // heaviest element first
				});
			for (const auto &[charge, count] : by_charge) {
				gdata[gi].elem_list += constants::atnr2letter(charge);
				if (count > 1) gdata[gi].elem_list += std::to_string(count);
			}
		}
	}

	// ---- Header ----
	std::cout << "\n\nRoby-Gould Bond Indices (RGBI) - Group Analysis\n";
	std::cout << "----------------------------------------------\n";
	std::cout << "Groups defined:\n";
	for (int g = 0; g < n_groups; g++) {
		std::cout << "  G" << g << ": atoms ";
		for (int k = 0; k < static_cast<int>(gdata[g].sorted_atoms.size()); k++) {
			if (k > 0) std::cout << ", ";
			const int ai = gdata[g].sorted_atoms[k];
			std::cout << ai << "(" << constants::atnr2letter(atom_charges[ai]) << ")";
		}
		std::cout << "\n";
	}
	std::cout << "----------------------------------------------\n";

	// ---- Group populations ----
	std::cout << "\nGroup populations:\n";
	for (int g = 0; g < n_groups; g++) {
		const std::string label = "G" + std::to_string(g) + " (" + gdata[g].elem_list + ")";
		std::cout << "  Population of " << std::left << std::setw(24) << label
			<< std::right << std::fixed << std::setprecision(5) << gdata[g].pop << "\n";
	}
	std::cout << "----------------------------------------------\n";

	// ---- Pair analysis ----
	for (int gi = 0; gi < n_groups; gi++) {
		for (int gj = gi + 1; gj < n_groups; gj++) {
			const ivec &ga_sorted = gdata[gi].sorted_atoms;
			const ivec &gb_sorted = gdata[gj].sorted_atoms;
			const ivec &ga_bf = gdata[gi].bf;
			const ivec &gb_bf = gdata[gj].bf;
			const dMatrix2 &P_GA = gdata[gi].P;
			const dMatrix2 &P_GB = gdata[gj].P;
			const double pop_ga = gdata[gi].pop;
			const double pop_gb = gdata[gj].pop;

			// bond_bf: ga's BFs first, gb's second — order MUST match ionic operator layout
			ivec bond_bf;
			bond_bf.insert(bond_bf.end(), ga_bf.begin(), ga_bf.end());
			bond_bf.insert(bond_bf.end(), gb_bf.begin(), gb_bf.end());

			// Pair population
			ivec bond_evals = gdata[gi].nao_evals, bond_evecs = gdata[gi].nao_evecs;
			bond_evals.insert(bond_evals.end(), gdata[gj].nao_evals.begin(), gdata[gj].nao_evals.end());
			bond_evecs.insert(bond_evecs.end(), gdata[gj].nao_evecs.begin(), gdata[gj].nao_evecs.end());
			err_checkf(!bond_bf.empty(),
				"RGBI group pair G" + std::to_string(gi) + "-G" + std::to_string(gj) +
				" has no basis functions.", std::cout);
			const bool pair_spans_full_basis =
				static_cast<int>(bond_bf.size()) == static_cast<int>(overlap_matrix.extent(0));
			// Full-basis group pairs use the same projection path as the total RGBI population.
			// This preserves the requested group BF ordering used by the established references.
			const double pair_pop = pair_spans_full_basis
				? projection_matrix_and_expectation(bond_bf)
				: projection_matrix_and_expectation(bond_bf, bond_evals, bond_evecs);

			// Build ionic operator [P_GA | 0; 0 | -P_GB] using dense group projections
			const int n_a = static_cast<int>(ga_bf.size());
			const int n_b = static_cast<int>(gb_bf.size());
			const int n_ab = n_a + n_b;
			dMatrix2 Ionic_Operator(n_ab, n_ab);
			for (int r = 0; r < n_ab; r++)
				for (int c = 0; c < n_ab; c++)
					Ionic_Operator(r, c) = 0.0;
			for (int r = 0; r < n_a; r++)
				for (int c = 0; c < n_a; c++)
					Ionic_Operator(r, c) = P_GA(r, c);
			for (int r = 0; r < n_b; r++)
				for (int c = 0; c < n_b; c++)
					Ionic_Operator(n_a + r, n_a + c) = -P_GB(r, c);

			// Sub-overlap → sqrt → solve eigenproblem
			vec S_Sub(static_cast<size_t>(n_ab) * n_ab, 0.0);
			get_submatrix(overlap_matrix, S_Sub, bond_bf);
			vec V = S_Sub;
			vec W(n_ab);
			const vec Temp = mat_sqrt(V, W);
			dMatrix2 A = reshape<dMatrix2>(Temp, Shape2D(n_ab, n_ab));
			dMatrix2 SI = LAPACKE_invert(A);

			auto X = change_basis_general(Ionic_Operator, transpose(A), true);
			vec ionic_eigenvals(X.extent(0));
			make_Eigenvalues(X.container(), ionic_eigenvals);

			auto EVC = dot<dMatrix2>(SI, X);

			// Prune near-zero eigenvalues
			ivec non_zero;
			for (int i = 0; i < static_cast<int>(ionic_eigenvals.size()); i++)
				if (abs(ionic_eigenvals[i]) > 1E-5)
					non_zero.push_back(i);

			vec pruned_eigvals;
			for (int idx : non_zero)
				pruned_eigvals.push_back(ionic_eigenvals[idx]);

			EVC = transpose(EVC);
			auto EVC2 = transpose(get_rectangle(EVC, non_zero));
			EVC.container().clear();

			const int n0 = static_cast<int>(pruned_eigvals.size());
			auto pairs = find_eigenvalue_pairs(pruned_eigvals);

			transform_group_Ionic_orbitals(EVC2, pruned_eigvals, pairs, ga_sorted, gb_sorted, bond_bf, P_GA, P_GB);

			auto covalent_info = make_covalent_from_ionic(EVC2, pruned_eigvals, pairs, EVs);

			// Per-orbital populations
			vec covalent_popul(n0), ionic_popul(n0);
			ivec vals;
			for (int i = 0; i < static_cast<int>(EVC2.extent(0)); i++)
				vals.emplace_back(i);
			for (int i = 0; i < n0; i++) {
				auto temp = transpose(covalent_info['T']);
				covalent_popul[i] = projection_matrix_and_expectation(bond_bf, { i }, vals, &temp);
				temp = transpose(EVC2);
				ionic_popul[i] = projection_matrix_and_expectation(bond_bf, { i }, vals, &temp);
			}

			const double zero_angle_cutoff = 1E-2 * constants::INV_PI_180;
			group_bond_index_result result;
			result.group_index_first = gi;
			result.group_index_second = gj;
			result.atoms_first = ga_sorted;
			result.atoms_second = gb_sorted;
			result.pair_population = pair_pop;
			result.population_first = pop_ga;
			result.population_second = pop_gb;
			result.covalent = 0.0;
			result.ionic = 0.0;

			// Lone-pair orbitals localized almost entirely within one group have
			// eigvals near ±1.  In the group basis the inter-group contamination
			// can push these to ~0.996 rather than exactly 1, so they evade the
			// angle cutoff (84.6° lies inside its [0.573°, 89.427°] window) and corrupt
			// the ionic index with large contributions of the wrong sign.  Exclude any
			// pair whose positive eigval exceeds this threshold.
			// TODO: the atom-pair loop lacks this guard, so the two paths disagree on lone pairs.
			const double lone_pair_eigval_threshold = 0.99;
			for (int i = 0; i < n0; i++) {
				if (covalent_info['A'](i, 0) < zero_angle_cutoff || covalent_info['A'](i, 0) > 90.0 - zero_angle_cutoff)
					continue;
				if (pruned_eigvals[i] < pruned_eigvals[pairs[i]])
					continue;
				if (pruned_eigvals[i] > lone_pair_eigval_threshold)
					continue;
				if (pairs[i] != i) {
					result.covalent += 0.5 * (covalent_popul[i] - covalent_popul[pairs[i]]);
					result.ionic += 0.5 * (ionic_popul[i] - ionic_popul[pairs[i]]);
				}
			}

			const double b2 = result.covalent * result.covalent + result.ionic * result.ionic;
			result.percent_covalent_Pyth = (b2 > 1E-12) ? 100.0 * (result.covalent * result.covalent / b2) : 0.0;
			const double b = sqrt(b2);
			result.total = b;
			result.percent_covalent_Arakai = (b > 1E-12) ? 200.0 * abs(asin(result.covalent / b)) / constants::PI : 0.0;

			RGBI_groups.push_back(result);
		}
	}

	// ---- Results table ----
	const std::string sep(90, '-');
	std::cout << "\n" << std::left << std::setw(14) << "Group" << std::right
		<< std::setw(8) << "n_G1"
		<< std::setw(8) << "n_G2"
		<< std::setw(8) << "n_G1G2"
		<< std::setw(8) << "s_G1G2"
		<< std::setw(8) << "Cov."
		<< std::setw(8) << "Ion."
		<< std::setw(8) << "Tot."
		<< std::setw(8) << "Pyth."
		<< std::setw(8) << "Arak."
		<< "\n" << sep << "\n";
	for (const auto &res : RGBI_groups) {
		const std::string label = "G" + std::to_string(res.group_index_first)
			+ " - G" + std::to_string(res.group_index_second);
		std::cout << std::left << std::setw(14) << label << std::right
			<< std::fixed << std::setprecision(3)
			<< std::setw(8) << res.population_first
			<< std::setw(8) << res.population_second
			<< std::setw(8) << res.pair_population
			<< std::setw(8) << (res.population_first + res.population_second - res.pair_population)
			<< std::setw(8) << res.covalent
			<< std::setw(8) << res.ionic
			<< std::setw(8) << res.total
			<< std::setw(8) << res.percent_covalent_Pyth
			<< std::setw(8) << res.percent_covalent_Arakai
			<< "\n";
	}
	std::cout << sep << "\n";
}

Roby_information::Roby_information(WFN &wavy, const ivec3 &group_sets, const bool symmetrize, const bool use_ano_basis, const bool EVs, const bool theta_info, const bool legacy_occupancy_cutoff) {
	//The tables set three or four decimals on cout; restore so later output keeps its precision.
	const ostream_format_guard restore_cout_format(std::cout);
	//The pinv cutoff decides the rank of near-singular metrics: settable to tell physics from threshold, and
	//announced so a non-default run is not mistaken for a default one.
	if (const char *env = tuning("NOS_RGBI_PINV_CUTOFF")) {
		try {
			const double v = std::stod(env);
			err_checkf(v > 0.0, "NOS_RGBI_PINV_CUTOFF must be positive, got '" + std::string(env) + "'.", std::cout);
			pinv_cutoff = v;
			std::cout << "NOS_RGBI_PINV_CUTOFF is set: RGBI pseudo-inverses cut singular values below "
				<< std::scientific << std::setprecision(3) << pinv_cutoff << " instead of the default 1.000e-05\n"
				<< std::defaultfloat;
		}
		catch (const std::invalid_argument &) {
			err_checkf(false, "NOS_RGBI_PINV_CUTOFF is not a number: '" + std::string(env) + "'.", std::cout);
		}
	}
	auto bonds = get_bonded_atom_pairs(wavy);
	//A .wfn has no shell structure, so every atom's basis is empty; OCC would segfault (uncatchable) in
	//gensqrtinv on the ANO route. One check for both routes.
	for (int a = 0; a < wavy.get_ncen(); a++)
		err_checkf(!wavy.get_atom(a).get_basis_set().empty(),
			"RGBI needs the basis set of every atom, and " + wavy.get_path().filename().string() +
			" carries none for atom " + std::to_string(a + 1) + " (" +
			constants::atnr2letter(wavy.get_atom(a).get_charge()) + "). A file that lists primitives "
			"by centre without their shell structure - .wfn, .ffn and .wfx all do - leaves every atom's "
			"basis empty; run RGBI on " + rgbi_supported_input_phrase(wavy.get_path().extension().string()) +
			" instead.", std::cout);
	citations::cite(citations::Method::RGBI, std::cout);
	const char *orbital_label = use_ano_basis ? "ANOs" : "NAOs";
	std::cout << "Calculating " << orbital_label << " for all atoms...                 " << std::flush;
	if (legacy_occupancy_cutoff)
		std::cout << "\n  (legacy occupancy cutoff: atomic subspace ranks follow the occupation numbers)\n";
	computeAllAtomicNAOs(wavy, symmetrize, use_ano_basis, EVs, legacy_occupancy_cutoff);
	std::cout << " ...done!" << std::endl;
	if (theta_info)
		std::cout << "RGBI theta-subspace reports enabled." << std::endl;
	if (use_ano_basis) {
		if (ano_fallback_atoms.empty()) {
			std::cout << "ANO fallback summary: no atom-level fallbacks were needed." << std::endl;
		}
		else {
			std::cout << "ANO fallback summary: used molecular local-orbital fallback for "
				<< ano_fallback_atoms.size() << "/" << wavy.get_ncen() << " atoms";
			std::cout << " (";
			for (int i = 0; i < static_cast<int>(ano_fallback_atoms.size()); ++i) {
				if (i > 0)
					std::cout << ", ";
				std::cout << ano_fallback_atoms[i];
			}
			std::cout << ")." << std::endl;
		}
	}
	Shape2D NAOs_size;
	NAOs_size.cols = static_cast<int>(overlap_matrix.extent(0));
	for (size_t atom_idx = 0; atom_idx < NAOs.size(); atom_idx++) {
#ifdef NSA2DEBUG
		std::cout << "Atom " << atom_idx + 1 << " NAO Occupancies:\n";
#endif
		for (size_t i = 0; i < NAOs[atom_idx].eigenvalues.size(); i++) {
#ifdef NSA2DEBUG
			std::cout << "  NAO " << i + 1 << ": " << std::setprecision(6) << NAOs[atom_idx].eigenvalues[i] << "\n";
#endif
			NAOs_size.rows++;
		}
#ifdef NSA2DEBUG
		std::cout << "\n";
#endif
	}
	//Create a matrix of size n_NAOs x n_bf and fill with subblocks of the evecs
	vec NAO_matrix(NAOs_size.cols * NAOs_size.rows, 0.0);
	total_NAOs = reshape<dMatrix2>(NAO_matrix, Shape2D(NAOs_size.rows, NAOs_size.cols));
	Shape2D temp = { 0,0 };
	for (auto NAO : NAOs) {
		const int ME = NAO.matrix_elements.size();
		const int NEV = NAO.eigenvalues.size();
		for (int i = 0; i < NEV; i++) {
			const int row = temp.rows + i;
			for (int j = 0; j < ME; j++) {
				const int index_eigenvector = i * ME + j;
				const int col = NAO.matrix_elements[j];
				total_NAOs(row, col) = NAO.eigenvectors[index_eigenvector];
			}
		}
		temp.rows += NAO.eigenvalues.size();
	}
	dMatrix2 total_NAOs_with_omitted;
	if (use_ano_basis) {
		Shape2D full_NAOs_size;
		full_NAOs_size.cols = NAOs_size.cols;
		for (const auto &NAO : NAOs) {
			full_NAOs_size.rows += static_cast<int>(
				NAO.eigenvalues.size() + NAO.omitted_eigenvalues.size());
		}
		vec full_NAO_matrix(full_NAOs_size.cols * full_NAOs_size.rows, 0.0);
		total_NAOs_with_omitted = reshape<dMatrix2>(full_NAO_matrix, full_NAOs_size);
		Shape2D full_temp = { 0,0 };
		for (const auto &NAO : NAOs) {
			const int ME = NAO.matrix_elements.size();
			const int NEV = NAO.eigenvalues.size();
			for (int i = 0; i < NEV; i++) {
				const int row = full_temp.rows + i;
				for (int j = 0; j < ME; j++) {
					const int index_eigenvector = i * ME + j;
					const int col = NAO.matrix_elements[j];
					total_NAOs_with_omitted(row, col) = NAO.eigenvectors[index_eigenvector];
				}
			}
			full_temp.rows += NEV;
			const int omitted_NEV = NAO.omitted_eigenvalues.size();
			for (int i = 0; i < omitted_NEV; i++) {
				const int row = full_temp.rows + i;
				for (int j = 0; j < ME; j++) {
					const int index_eigenvector = i * ME + j;
					const int col = NAO.matrix_elements[j];
					total_NAOs_with_omitted(row, col) = NAO.omitted_eigenvectors[index_eigenvector];
				}
			}
			full_temp.rows += omitted_NEV;
		}
	}
#ifdef NSA2DEBUG
	print_dmatrix2(total_NAOs, "Global NAO Matrix");
#endif
	// total_NAOs has n_NAOs rows (each row is an individual NAO in the basis set)
	// and n_bfs columns, where each value corresponds to a basis function from here on
	std::cout << "Calculating Population for all atoms...           " << std::flush;
	const double all_atom_population = Roby_population_analysis({});
	double all_atom_population_with_omitted = all_atom_population;
	if (use_ano_basis && total_NAOs_with_omitted.size() > 0) {
		ivec all_basis_indices(static_cast<int>(overlap_matrix.extent(0)));
		std::iota(all_basis_indices.begin(), all_basis_indices.end(), 0);
		dMatrix2 full_P;
		all_atom_population_with_omitted = projection_matrix_and_expectation(
			all_basis_indices, {}, {}, &total_NAOs_with_omitted, &full_P);
	}
	std::cout << " ...done." << std::endl;
	std::cout << std::endl << "Total Population: " << all_atom_population << "\n\n";
	vec atom_pops(NAOs.size(), 0.0);
	vec omitted_atom_pops(NAOs.size(), 0.0);
	projection_matrices.clear();
	projection_matrices.resize(NAOs.size());
	overlap_matrices.clear();
	overlap_matrices.resize(NAOs.size());
#ifndef NSA2DEBUG
	ProgressBar *pb = new ProgressBar(NAOs.size(), 40, "-", " ", "Calculating Atomic Populations");
#endif
	for (auto NAO : NAOs) {
		dMatrix2 P;
		atom_pops[NAO.atom_index] = projection_matrix_and_expectation(NAO.matrix_elements, {}, {}, nullptr, &P);
		projection_matrices[NAO.atom_index] = P;
		if (use_ano_basis && !NAO.omitted_eigenvalues.empty()) {
			vec full_eigenvectors = NAO.eigenvectors;
			full_eigenvectors.insert(full_eigenvectors.end(),
				NAO.omitted_eigenvectors.begin(), NAO.omitted_eigenvectors.end());
			dMatrix2 full_NAOs = reshape<dMatrix2>(
				full_eigenvectors,
				Shape2D(static_cast<int>(NAO.eigenvalues.size() + NAO.omitted_eigenvalues.size()),
					static_cast<int>(NAO.matrix_elements.size())));
			dMatrix2 full_P;
			const double full_atom_population =
				projection_matrix_and_expectation(NAO.matrix_elements, {}, {}, &full_NAOs, &full_P);
			omitted_atom_pops[NAO.atom_index] = full_atom_population - atom_pops[NAO.atom_index];
		}
#ifndef NSA2DEBUG
		pb->update();
#endif
	}
#ifndef NSA2DEBUG
	delete pb;
#endif
	for (int i = 0; i < atom_pops.size(); i++) {
		std::cout << "Population of atom " << i << ": " << atom_pops[i] << std::endl;
	}
	if (use_ano_basis) {
		for (int i = 0; i < omitted_atom_pops.size(); i++) {
			if (std::abs(omitted_atom_pops[i]) > 1E-8) {
				std::cout << "Atomic projector population outside selected ANO cutoff of atom " << i << ": "
					<< omitted_atom_pops[i] << "\n";
			}
		}
	}

#ifndef NSA2DEBUG
	pb = new ProgressBar(bonds.size(), 40, "-", " ", "Calculating Bond Populations");
#endif
	std::string theta_reports;
	//now perform bond analysis for all bonded atoms
	for (auto bond : bonds) {
		//std::cout << std::endl << "---------------------------- Atom Pair: " << bond.first << " " << bond.second << " ----------------------\n";
		ivec bond_indices, bond_eigenvecs, bond_eigenvals;
		//gather basis function indices for both atoms
		for (const auto &NAO : NAOs) {
			if (NAO.atom_index == bond.first || NAO.atom_index == bond.second)
				for (auto idx : NAO.matrix_elements)
					bond_indices.push_back(idx);
		}
		//just in case: sort bond indices, easy since each atom's indices are already assumed sorted and each basis function belongs to only one atom once
		std::sort(bond_indices.begin(), bond_indices.end());

		//now determine the bond atom NAO indices
		int start_val = 0;
		for (auto NAO : NAOs) {
			if (NAO.atom_index == bond.first || NAO.atom_index == bond.second) {
				const int n1 = NAO.eigenvalues.size();
				for (int i = 0; i < n1; i++) {
					bond_eigenvals.push_back(start_val + i);
				}
				for (int idx : NAO.matrix_elements)
					bond_eigenvecs.push_back(idx);
			}
			start_val += NAO.eigenvalues.size();
		}

		//calcualte population using data from both atoms
		const double bond_population = projection_matrix_and_expectation(bond_indices, bond_eigenvals, bond_eigenvecs);
		//s_AB = n_A + n_B - n_AB cancels, amplifying relative errors here into s_AB and the theta decomposition:
		//flag the bond if the cutoff, not a spectral gap, set the pair metric's rank.
		if ((last_pinv_kept > 0 && last_pinv_smallest_kept < 10.0 * pinv_cutoff)
			|| last_pinv_largest_dropped > 0.1 * pinv_cutoff) {
			std::ostringstream w;
			w << "  " << wavy.get_atoms()[bond.first].get_label() << " - " << wavy.get_atoms()[bond.second].get_label()
				<< ": kept " << last_pinv_kept << " of " << last_pinv_n
				<< " singular values, smallest kept " << std::scientific << std::setprecision(3) << last_pinv_smallest_kept
				<< ", largest dropped " << last_pinv_largest_dropped;
			pinv_warnings.push_back(w.str());
		}
		//atom_pair_populations(bond.first, bond.second) = bond_population;
		//atom_pair_populations(bond.second, bond.first) = bond_population;
		std::cout << "Bond population between atom " << bond.first + 1 << " and atom " << bond.second + 1 << ": " << bond_population << "\n";
		const int size_ion1 = projection_matrices[bond.first].extent(0) + projection_matrices[bond.second].extent(0),
			size_ion2 = projection_matrices[bond.first].extent(1) + projection_matrices[bond.second].extent(1),
			size1_1 = projection_matrices[bond.first].extent(0),
			size1_2 = projection_matrices[bond.first].extent(1),
			size2_1 = projection_matrices[bond.second].extent(0),
			size2_2 = projection_matrices[bond.second].extent(1);
		dMatrix2 Ionic_Operator(size_ion1, size_ion2);
		for (int i = 0; i < size1_1; i++) {
			for (int j = 0; j < size1_2; j++) {
				Ionic_Operator(i, j) = projection_matrices[bond.first](i, j);
			}
		}
		for (int i = 0; i < size2_1; i++) {
			for (int j = 0; j < size2_2; j++) {
				Ionic_Operator(i + size1_1, j + size1_2) = -projection_matrices[bond.second](i, j);
			}
		}
#ifdef NSA2DEBUG
		print_dmatrix2(Ionic_Operator, "Ionic Operator");
#endif
		// get matching suboverlap matrix
		const int n = bond_indices.size();
		vec S_Sub(static_cast<size_t>(n) * n, 0.0);
		get_submatrix(overlap_matrix, S_Sub, bond_indices);
		//dMatrix2 S = reshape<dMatrix2>(S_Sub, Shape2D(n, n));

		vec V = S_Sub;
		vec W(n);
		//Rank decisions on the pair overlap and its square root: routed through pinv_cutoff and reported, like the atomic ones.
		PinvRank sqrt_rank{}, inverse_rank{};
		// make V = Sqrt(S)
		const vec Temp = mat_sqrt(V, W, pinv_cutoff, &sqrt_rank);

		dMatrix2 A = reshape<dMatrix2>(Temp, Shape2D(n, n));
#ifdef NSA2DEBUG
		print_dmatrix2(A, "Overlap Sqrt SH");
#endif

		dMatrix2 SI = LAPACKE_invert(A, pinv_cutoff, &inverse_rank);
		//Shared collector: a marginal pair metric goes into the one warning block after the table, not a line per bond.
		if (sqrt_rank.marginal(pinv_cutoff) || inverse_rank.marginal(pinv_cutoff)) {
			std::ostringstream w;
			w << "  " << wavy.get_atoms()[bond.first].get_label() << " - " << wavy.get_atoms()[bond.second].get_label()
				<< ": pair overlap kept " << sqrt_rank.kept << " of " << sqrt_rank.n << " eigenvalues (smallest kept "
				<< std::scientific << std::setprecision(3) << sqrt_rank.smallest_kept << ", largest dropped "
				<< sqrt_rank.largest_dropped << "), its square root kept " << inverse_rank.kept << " of "
				<< inverse_rank.n << " (smallest kept " << inverse_rank.smallest_kept << ", largest dropped "
				<< inverse_rank.largest_dropped << ")";
			pinv_warnings.push_back(w.str());
		}

#ifdef NSA2DEBUG
		print_dmatrix2(SI, "Overlap Pseudo Inverse");
#endif

		auto X = change_basis_general(Ionic_Operator, transpose(A), true);
#ifdef NSA2DEBUG
		print_dmatrix2(X, "Overlap Eigenproblem");
#endif

		// solve symmetric eigenproblem of X
		vec ionic_eigenvals(X.extent(0));
		make_Eigenvalues(X.container(), ionic_eigenvals);

#ifdef NSA2DEBUG
		std::cout << "Ionic eigenvalues between atom " << bond.first + 1 << " and atom " << bond.second + 1 << ":\n";
		for (size_t i = 0; i < ionic_eigenvals.size(); i++) {
			std::cout << std::setw(5) << i + 1 << ": " << std::setw(10) << std::setprecision(6) << ionic_eigenvals[i] << "\n";
		}
		std::cout << "\n";
#endif

		auto EVC = dot<dMatrix2>(SI, X);
#ifdef NSA2DEBUG
		print_dmatrix2(EVC, "theta_I");
#endif

		ivec non_zero_indices;
		for (int i = 0; i < ionic_eigenvals.size(); i++) {
			if (abs(ionic_eigenvals[i]) > 1E-5)
				non_zero_indices.push_back(i);
		}
		const int n0 = non_zero_indices.size();
#ifdef NSA2DEBUG
		std::cout << "Non zero eigenvalues:\n";
		for (size_t i = 0; i < n0; i++) {
			std::cout << std::setw(5) << i + 1 << ": " << std::setw(10) << std::setprecision(6) << non_zero_indices[i] << "\n";
		}
		std::cout << "\n";
#endif

		vec pruned_eigvals;
		for (int nzv = 0; nzv < n0; nzv++)
			pruned_eigvals.push_back(ionic_eigenvals[non_zero_indices[nzv]]);
		EVC = transpose(EVC);
		auto EVC2 = transpose(get_rectangle(EVC, non_zero_indices));
#ifdef NSA2DEBUG
		print_dmatrix2(EVC2, "theta_I after pruning:");
#endif
		EVC.container().clear();

		auto pairs = find_eigenvalue_pairs(pruned_eigvals);
#ifdef NSA2DEBUG
		std::cout << "Pairs:\n";
		for (size_t i = 0; i < n0; i++) {
			std::cout << std::setw(3) << i << ": " << std::setw(5) << std::setprecision(6) << pairs[i] << "\n";
		}
		std::cout << "\n";
#endif

		transform_Ionic_eigenvectors_to_Ionic_orbitals(EVC2, pruned_eigvals, pairs, bond.first, bond.second, bond_indices);
		if (EVs) {
			print_dmatrix2(EVC2, "theta_I after Ionic");
		}
		//EVC2 is theta_Ionic, now make theta_Covalent
		auto covalent_info = make_covalent_from_ionic(EVC2, pruned_eigvals, pairs, EVs);
		if (EVs) {
			print_dmatrix2(covalent_info['A'], "theta_Angles");
			print_dmatrix2(covalent_info['V'], "eigen_C");
			print_dmatrix2(covalent_info['T'], "theta_C");
		}

		vec covalent_popul(n0), ionic_popul(n0);
		ivec vals;
		for (int i = 0; i < EVC2.extent(0); i++)
			vals.emplace_back(i);
		for (int i = 0; i < n0; i++) {
			//make covalent populations
			auto temp = transpose(covalent_info['T']);
			covalent_popul[i] = projection_matrix_and_expectation(bond_indices, { i }, vals, &(temp));
			temp = transpose(EVC2);
			ionic_popul[i] = projection_matrix_and_expectation(bond_indices, { i }, vals, &(temp));
		}
		if (theta_info)
			theta_reports += make_theta_info(wavy, bond, pruned_eigvals, pairs, covalent_info['A'], covalent_popul, ionic_popul);
#ifdef NSA2DEBUG
		std::cout << "Covalent populations:\n";
		for (int i = 0; i < n0; i++)
			std::cout << "\t" << i << ":\t" << covalent_popul[i] << "\n";
		std::cout << "Ionic populations:\n";
		for (int i = 0; i < n0; i++)
			std::cout << "\t" << i << ":\t" << ionic_popul[i] << "\n";
#endif

		const double zero_angle_cutoff = 1E-2 * constants::INV_PI_180;

		vec cov_index(n0, 0.0), ion_index(n0, 0.0);
		double b_ab = 0;
		bond_index_result results{};
		for (int i = 0; i < n0; i++) {
			if (covalent_info['A'](i, 0) < zero_angle_cutoff || covalent_info['A'](i, 0) > 90 - zero_angle_cutoff)
				continue; // skip lone pairs

			if (pruned_eigvals[i] < pruned_eigvals[pairs[i]])
				continue; // skip antibonding

			if (pairs[i] != i) {
				cov_index[i] = 0.5 * (covalent_popul[i] - covalent_popul[pairs[i]]);
				ion_index[i] = 0.5 * (ionic_popul[i] - ionic_popul[pairs[i]]);
				results.covalent += cov_index[i];
				results.ionic += ion_index[i];
			}
			else if (pruned_eigvals[i] > 0.0)
				ion_index[i] = 0.5 * ionic_popul[i];
			else if (pruned_eigvals[i] < 0.0)
				ion_index[i] = -0.5 * ionic_popul[i];
		}

		b_ab = results.covalent * results.covalent + results.ionic * results.ionic;
		results.percent_covalent_Pyth = 100 * (results.covalent * results.covalent / b_ab);
		b_ab = sqrt(b_ab);
		results.percent_covalent_Arakai = 200 * abs(asin(results.covalent / b_ab)) / constants::PI;

		int el_a = wavy.get_atom_charge(bond.first);
		int el_b = wavy.get_atom_charge(bond.second);
		if (el_a > el_b) {
			results.atom_indices = bond;
			results.atom_element_nr = std::make_pair(el_a, el_b);
			results.total = b_ab;
			results.pair_population = bond_population;
			results.population_first = atom_pops[bond.first];
			results.population_second = atom_pops[bond.second];
		}
		else {
			results.atom_indices = std::make_pair(bond.second, bond.first);
			results.atom_element_nr = std::make_pair(el_b, el_a);
			results.total = b_ab;
			results.pair_population = bond_population;
			results.population_first = atom_pops[bond.second];
			results.population_second = atom_pops[bond.first];
			results.ionic = -results.ionic;
		}
		RGBI.push_back(results);
#ifndef NSA2DEBUG
		pb->update();
#endif
	}
#ifndef NSA2DEBUG
	delete pb;
#endif
	if (theta_info)
		std::cout << theta_reports << std::endl;
	//get_nr_electrons() counts the nuclei; an ECP wavefunction does not describe the core.
	const double number_of_electrons =
		wavy.get_nr_electrons() - static_cast<double>(wavy.get_nr_ECP_electrons());
	//A difference of two sums of order N: below the cancellation floor print zero (the sign is round-off);
	//a real cutoff is far above it, so a genuinely negative value still shows.
	double omitted_population = use_ano_basis
		? all_atom_population_with_omitted - all_atom_population
		: 0.0;
	if (std::abs(omitted_population) < 1e-9 * std::max(1.0, all_atom_population))
		omitted_population = 0.0;

	const section_log::tee rgbi_tee("rgbi", "Bond populations, covalent and ionic indices");
	std::cout << "\nRoby-Gould Bond Indices (RGBI) Analysis\n----------------------------------------------\n";
	std::cout << "Number of electrons in system:         " << number_of_electrons;
	if (wavy.get_nr_ECP_electrons() > 0)
		std::cout << "  (" << wavy.get_nr_ECP_electrons() << " core electrons sit in ECPs and are "
		"not described by this wavefunction)";
	std::cout << "\n";
	std::cout << "Number of electrons in Roby Analysis:  " << all_atom_population << "\n";
	if (use_ano_basis) {
		std::cout << "Number of cutoff electrons:            " << omitted_population << "\n";
		std::cout << "Number after cutoff projection:        " << all_atom_population_with_omitted << "\n";
	}
	std::cout << "Percentage of electrons accounted for: " << std::setprecision(4) << (100.0 * all_atom_population / number_of_electrons) << " %\n";
	if (use_ano_basis) {
		std::cout << "Percentage after cutoff projection:    "
			<< std::setprecision(4)
			<< (100.0 * all_atom_population_with_omitted / number_of_electrons)
			<< " %\n";
	}
	std::cout << "----------------------------------------------\n";

	std::cout << "\n\nAtom Nr        Els  "
		<< std::setw(8) << "n_A"
		<< std::setw(8) << "n_B"
		<< std::setw(8) << "n_AB"
		<< std::setw(8) << "s_AB"
		<< std::setw(8) << "Cov."
		<< std::setw(8) << "Ion."
		<< std::setw(8) << "Tot."
		<< std::setw(8) << "Pyth."
		<< std::setw(8) << "Arak."
		<< std::endl;
	std::cout << "--------------------------------------------------------------------------------------------\n";

	//Sort results by Element of the heavier atom of the bonds and then within each group by the heavier of the second atom
	std::sort(RGBI.begin(), RGBI.end(), [](const bond_index_result &a, const bond_index_result &b) {
		if (a.atom_element_nr.first != b.atom_element_nr.first) {
			return a.atom_element_nr.first > b.atom_element_nr.first;
		}
		return a.atom_element_nr.second > b.atom_element_nr.second;
		});

	for (auto res : RGBI) {
		std::cout << std::setw(4) << res.atom_indices.first << " -" << std::setw(4) << res.atom_indices.second << "  "
			<< std::setw(3) << constants::atnr2letter(res.atom_element_nr.first) << " -" << std::setw(3) << constants::atnr2letter(res.atom_element_nr.second)
			<< std::fixed << std::setprecision(3) << std::setw(8) << res.population_first
			<< std::fixed << std::setprecision(3) << std::setw(8) << res.population_second
			<< std::fixed << std::setprecision(3) << std::setw(8) << res.pair_population
			<< std::fixed << std::setprecision(3) << std::setw(8) << res.population_first + res.population_second - res.pair_population
			<< std::fixed << std::setprecision(3) << std::setw(8) << res.covalent
			<< std::fixed << std::setprecision(3) << std::setw(8) << res.ionic
			<< std::fixed << std::setprecision(3) << std::setw(8) << res.total
			<< std::fixed << std::setprecision(3) << std::setw(8) << res.percent_covalent_Pyth
			<< std::fixed << std::setprecision(3) << std::setw(8) << res.percent_covalent_Arakai << "\n";
	}
	std::cout << "--------------------------------------------------------------------------------------------\n";

	//A rank set by the threshold rather than a spectral gap can split symmetry-equivalent bonds, so say so.
	if (!pinv_warnings.empty()) {
		std::cout << "\nWARNING: " << pinv_warnings.size() << " of " << bonds.size()
			<< " pair subspaces had their rank fixed by the " << std::scientific << std::setprecision(1) << pinv_cutoff
			<< " singular-value cutoff rather than by a gap in the spectrum, so a bond that symmetry makes "
			"identical to another can land on the other side of it and the numbers above then differ for no "
			"physical reason. NOS_RGBI_PINV_CUTOFF moves the cutoff to test that.\n";
		for (const auto &w : pinv_warnings)
			std::cout << w << "\n";
		std::cout << std::defaultfloat;
	}

	if (!group_sets.empty()) {
		const int N_atoms = wavy.get_ncen();
		ivec atom_charges(N_atoms);
		for (int i = 0; i < N_atoms; i++)
			atom_charges[i] = wavy.get_atom_charge(i);
		for (const auto &group_defs : group_sets)
			computeGroupAnalysis(group_defs, atom_pops, atom_charges, EVs);
	}
}

void bondwise_laplacian_plots(std::filesystem::path &wfn_name)
{
	WFN wavy(wfn_name);
	wavy.delete_unoccupied_MOs();

	err_checkf(wavy.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);

	int points = 1001;

	for (int i = 0; i < wavy.get_ncen(); i++)
	{
		for (int j = i + 1; j < wavy.get_ncen(); j++)
		{
			std::filesystem::path path = std::filesystem::current_path();
			double distance = sqrt(pow(wavy.get_atom_coordinate(i, 0) - wavy.get_atom_coordinate(j, 0), 2) + pow(wavy.get_atom_coordinate(i, 1) - wavy.get_atom_coordinate(j, 1), 2) + pow(wavy.get_atom_coordinate(i, 2) - wavy.get_atom_coordinate(j, 2), 2));
			double svdW = constants::ang2bohr(constants::covalent_radii[wavy.get_atom_charge(i)] + constants::covalent_radii[wavy.get_atom_charge(j)]);
			if (distance < 1.35 * svdW)
			{
				std::cout << "Bond between " << i << " (" << wavy.get_atom_charge(i) << ") and " << j << " (" << wavy.get_atom_charge(j) << ") with distance " << distance << " and svdW " << svdW << "\n";
				const vec bond_vec = { (wavy.get_atom_coordinate(j, 0) - wavy.get_atom_coordinate(i, 0)) / points, (wavy.get_atom_coordinate(j, 1) - wavy.get_atom_coordinate(i, 1)) / points, (wavy.get_atom_coordinate(j, 2) - wavy.get_atom_coordinate(i, 2)) / points };
				const double dr = distance / points;
				vec lapl(points, 0.0);
				const vec pos = { wavy.get_atom_coordinate(i, 0), wavy.get_atom_coordinate(i, 1), wavy.get_atom_coordinate(i, 2) };
#pragma omp parallel for schedule(dynamic)
				for (int k = 0; k < points; k++)
				{
					d3 t_pos = { pos[0], pos[1], pos[2] };
					t_pos[0] += k * bond_vec[0];
					t_pos[1] += k * bond_vec[1];
					t_pos[2] += k * bond_vec[2];
					lapl[k] = wavy.computeLap(t_pos);
				}
				std::filesystem::path outname(wfn_name.string() + "_bondwise_laplacian_" + std::to_string(i) + "_" + std::to_string(j) + ".dat");
				path = path / std::filesystem::path(outname);
				std::ofstream result(path, std::ios::out);
				for (int k = 0; k < points; k++)
				{
					result << std::setw(10) << std::scientific << std::setprecision(6) << dr * k << " " << std::setw(10) << std::scientific << std::setprecision(6) << lapl[k] << "\n";
				}
				result.close();
			}
			else
			{
				std::cout << "No bond between " << i << " and " << j << " with distance " << distance << " and svdW " << svdW << "\n";
			}
		}
	}
}

//The fitted density in place of the orbitals for rho and the QTAIM basins when the options
//ask for it: -ri_fit <basis> before the analysis flag on a wavefunction, or -SALTED
//<model-dir> on an xyz. rho and its gradient then come from one loop over the auxiliary
//functions, which is what makes 0.05 A grids affordable. nullptr for the orbital route.
//ELI-D keeps the orbitals: the density-only estimate (PC07 alpha, the -eli cubes) is a
//constant wherever alpha is switched off, i.e. over the bonds and lone pairs, and has no
//basins to find there.
static std::unique_ptr<Gaussian_Molecule> fitted_source(const WFN &wavy, options &opt, density_field &field, std::ostream &log)
{
	if (!Gaussian_Molecule::requested(wavy, opt)) {
		err_checkf(wavy.get_nmo() > 0, "No orbitals in " + wavy.get_path().string() + "; give -SALTED <model-dir> before the analysis flag to analyse a predicted density", log);
		return nullptr;
	}
	auto fit = std::make_unique<Gaussian_Molecule>(wavy, opt);
	log << "Density source: " << (wavy.get_nmo() == 0 ? "SALTED prediction from " + opt.salted_model_dir.string() : "RI fit of " + std::to_string(wavy.get_nmo()) + " orbitals")
		<< ", " << fit->n_aux() << " auxiliary functions.\n"
		<< "  rho and its gradient come from one pass over the auxiliary functions, N_aux per point instead of\n"
		<< "  N_basis x N_MO, so grids of 0.05 A and finer are affordable. The QTAIM basins and populations below\n"
		<< "  are those of the fitted density: it matches the orbital density to ~1e-3 e/bohr^3 at bond critical\n"
		<< "  points, and its tail (rho < 0.01) is not reliable. ELI-D "
		<< (wavy.get_nmo() == 0 ? "needs orbitals and is skipped." : "is taken from the orbitals.") << std::endl;
	const Gaussian_Molecule *f = fit.get();
	field.rho = [f](const d3 &p) { return f->rho(p); };
	field.grad = [f](const d3 &p, d3 &g) { double lap; f->values(p, g, lap); };
	return fit;
}

//Kohout's spin-resolved ELI-D (eli_family.h). For an open shell the spin-summed ELI-D is neither spin's pair
//function; the alpha-alpha, beta-beta and triplet members are, each with its own maxima and basins.
//Members come from the orbitals (eli_variants_for), never from the stated multiplicity.
static void spin_eli_analysis(const WFN &l_w, const options &opt, const std::vector<atom> &atoms, const double shell_dist, const double shell_tol)
{
	using eli_family::Member;
	std::string warn;
	const std::vector<Member> members = eli_family::eli_variants_for(l_w, &warn);
	auto has = [&](const Member m) { return std::find(members.begin(), members.end(), m) != members.end(); };
	if (!opt.spin_eli) { std::cout << "\nSpin-resolved ELI-D: off (-no_spin_eli)." << std::endl; return; }
	const eli_family::SpinSplit how = eli_family::spin_split(l_w);
	const bool unrestricted = how == eli_family::SpinSplit::unrestricted;
	double N[2]{ 0.0, 0.0 };
	int norb[2]{ 0, 0 };
	bool open_restricted = false;
	for (int mo = 0; mo < l_w.get_nmo(); mo++) {
		const double occ = l_w.get_MO_occ(mo);
		if (occ == 0.0) continue;
		double n[2];
		eli_family::mo_spin_occupations(occ, l_w.get_MO_op(mo), how, n);
		for (int s = 0; s < 2; s++) if (n[s] != 0.0) { N[s] += n[s]; norb[s]++; }
		if (how == eli_family::SpinSplit::halves && std::abs(occ - 2.0) > 1e-6) open_restricted = true;
	}
	if (!has(Member::eli_d_bb)) {
		std::cout << "\nSpin-resolved ELI-D: skipped, "
			<< (open_restricted ? "fractional occupations in one MO set without spin labels (natural orbitals) - alpha and beta densities are not resolved"
				: !unrestricted ? "restricted wavefunction - ELI-D(alpha-alpha) = ELI-D(beta-beta) is the ELI-D above and the triplet member a constant multiple of it"
				: N[1] < 1e-8 ? "the beta orbital set holds no electrons"
				: "N_alpha = N_beta and the beta orbitals are the alpha ones - not spin-polarised")
			<< "." << std::endl;
		return;
	}
	if (opt.basin_cube) { std::cout << "\nSpin-resolved ELI-D: skipped, it runs on the analytic field only and -basin_cube was given." << std::endl; return; }
	citations::cite(citations::Method::ELIFamily, std::cout);
	const double tf = eli_family::triplet_density_factor(l_w);
	std::cout << "\nSpin-resolved ELI-D (Kohout): N_alpha = " << std::fixed << std::setprecision(4) << N[0] << ", N_beta = " << N[1]
		<< (how == eli_family::SpinSplit::restricted_open ? " (restricted open shell: singly occupied MOs alpha, doubly occupied one of each)"
			: std::abs(N[0] - N[1]) < 1e-8 ? " (broken symmetry: N_alpha = N_beta, alpha and beta orbitals differ)" : "")
		<< ". The ELI-D above is the spin-summed field, not a pair function of either spin." << std::endl;
	struct member_basins { std::vector<d4> maxima; svec labels; vec pop, vol; vec2 spin; double outside = 0.0; bool done = false; };
	const char *names[3] = { "alpha-alpha", "beta-beta", "triplet" };
	member_basins res[3];
	const int nf = has(Member::eli_d_triplet) ? 3 : 2;
	for (int f = 0; f < nf; f++) {
		//A channel carried by one orbital has g_s = 0 everywhere: no field, only round-off basins
		if (f < 2 && norb[f] < 2) {
			std::cout << "\nELI-D " << names[f] << ": skipped, the " << (f ? "beta" : "alpha") << " electrons sit in a single orbital, so g_s = 0 and the field is undefined." << std::endl;
			continue;
		}
		basin_stage_timer clock;
		const eli_spin_field eval = [&l_w, f, tf](const d3 &p, double &y, d3 &g, double *aux) { l_w.computeELISpinGrad(p, f, tf, y, g, aux); };
		member_basins &r = res[f];
		const std::vector<d4> all = analytic_eli_maxima(l_w, opt.debug, &eval);
		r.maxima = all;
		cubei none;
		ivec core_map, shell_map;
		unify_core_basins(none, r.maxima, atoms, &core_map);
		unify_shell_basins(none, r.maxima, &shell_map, shell_dist, shell_tol, &atoms, &l_w, &eval);
		if (!shell_map.empty())
			for (size_t b = 1; b < core_map.size(); b++) core_map[b] = shell_map[core_map[b]];
		ivec edge_map;
		if (unify_boundary_basins(r.maxima, l_w, eval, &edge_map) > 0)
			for (size_t b = 1; b < core_map.size(); b++) core_map[b] = edge_map[core_map[b]];
		r.labels = assign_labels_to_basins(r.maxima, atoms, opt.debug, 1);
		clock.lap(std::string("ELI-D ") + names[f] + " maxima");
		r.pop = integrate_basins_on_atomic_grids(nullptr, nullptr, all, l_w, opt.accuracy, true, r.vol, r.outside, nullptr, nullptr, opt.basin_grid, nullptr, nullptr, &core_map, &eval, &r.spin);
		r.done = true;
		clock.lap(std::string("ELI-D ") + names[f] + " basins");
		double tot = 0.0, ts[2]{ 0.0, 0.0 };
		for (size_t b = 0; b < r.pop.size(); b++) {
			tot += r.pop[b];
			ts[0] += r.spin[b][0];
			ts[1] += r.spin[b][1];
		}
		//ponytail: rim seeds find empty surface attractors of the bounded field; they count in the totals but are not listed
		int empty = 0;
		for (size_t b = r.pop.size(); b-- > 0;)
			if (r.pop[b] < 5e-5) {
				r.pop.erase(r.pop.begin() + b); r.vol.erase(r.vol.begin() + b); r.spin.erase(r.spin.begin() + b);
				r.maxima.erase(r.maxima.begin() + b); r.labels.erase(r.labels.begin() + b);
				empty++;
			}
		//all.size() also counts the maxima unified into cores/shells and the empty rim attractors: say which
		const size_t merged = all.size() - r.maxima.size() - empty;
		std::cout << "\nELI-D " << names[f] << " Analysis (atomic quadrature grids), " << all.size() << " maxima";
		if (merged) std::cout << ", " << merged << " unified into core/shell basins";
		if (empty) std::cout << ", " << empty << " empty surface attractors (below 5e-5 e, not listed)";
		std::cout << ", " << r.maxima.size() << " basins";
		std::cout << ":\n"
			<< "  basin  label               electrons    N_alpha     N_beta       spin  <ELI-q>      volume         maximum        x          y          z\n";
		for (size_t b = 0; b < r.pop.size(); b++) {
			std::cout << std::setw(7) << b + 1 << "  " << std::left << std::setw(18) << r.labels[b] << std::right << std::fixed << std::setprecision(4)
				<< std::setw(11) << r.pop[b] << std::setw(11) << r.spin[b][0] << std::setw(11) << r.spin[b][1] << std::setw(11) << r.spin[b][0] - r.spin[b][1];
			//<ELI-q_s> = integral of rho_s Y_q over the basin / N_s; the triplet has no ELI-q partner
			//below 1e-3 e of the channel the average is a tail ratio of two vanishing integrals
			if (f < 2 && r.spin[b][f] > 1e-3) std::cout << std::setw(9) << r.spin[b][2] / r.spin[b][f];
			else std::cout << std::setw(9) << "-";
			std::cout << std::setw(12) << r.vol[b] << std::setw(16) << r.maxima[b][3]
				<< std::setprecision(3) << std::setw(11) << r.maxima[b][0] << std::setw(11) << r.maxima[b][1] << std::setw(11) << r.maxima[b][2] << "\n";
		}
		std::cout << std::setprecision(4) << "  total in basins: " << tot << "   N_alpha " << ts[0] << " of " << N[0] << ", N_beta " << ts[1] << " of " << N[1]
			<< "   outside every basin: " << r.outside << std::endl;
	}
	//aa vs bb: only basins with the same label pair. A label held once pairs outright (a unified core or shell
	//basin's maximum is an arbitrary point on its shell); a repeated label pairs by mutual nearest maxima within
	//a bohr. An alpha basin left without a partner is where the unpaired alpha electrons sit
	if (res[0].done && res[1].done) {
		const svec &LA = res[0].labels, &LB = res[1].labels;
		auto nearest = [](const d4 &m, const std::vector<d4> &in, const svec &lin, const std::string &l, double &d) {
			int best = -1;
			d = 1e30;
			for (size_t k = 0; k < in.size(); k++) {
				if (lin[k] != l) continue;
				const double q = std::sqrt(std::pow(m[0] - in[k][0], 2) + std::pow(m[1] - in[k][1], 2) + std::pow(m[2] - in[k][2], 2));
				if (q < d) { d = q; best = static_cast<int>(k); }
			}
			return best;
		};
		const std::vector<d4> &A = res[0].maxima, &B = res[1].maxima;
		std::vector<char> b_used(B.size(), 0);
		std::cout << "\nELI-D alpha-alpha <-> beta-beta basins (same label; a repeated label by mutual nearest maxima within 1 bohr):\n"
			<< "  aa basin                bb basin              distance  N_alpha(aa)  N_beta(bb)   N_a - N_b\n";
		for (size_t a = 0; a < A.size(); a++) {
			double d, back;
			const int b = nearest(A[a], B, LB, LA[a], d);
			const bool unique = std::count(LA.begin(), LA.end(), LA[a]) == 1 && std::count(LB.begin(), LB.end(), LA[a]) == 1;
			const bool pair = b >= 0 && (unique || (d < 1.0 && nearest(B[b], A, LA, LA[a], back) == static_cast<int>(a)));
			std::cout << std::setw(4) << a + 1 << " " << std::left << std::setw(18) << res[0].labels[a] << std::right;
			if (pair) {
				b_used[b] = 1;
				std::cout << std::setw(5) << b + 1 << " " << std::left << std::setw(18) << res[1].labels[b] << std::right << std::fixed << std::setprecision(3) << std::setw(8) << d
					<< std::setprecision(4) << std::setw(13) << res[0].spin[a][0] << std::setw(12) << res[1].spin[b][1] << std::setw(12) << res[0].spin[a][0] - res[1].spin[b][1] << "\n";
			}
			else
				std::cout << std::setw(5) << "-" << " " << std::left << std::setw(18) << "(no partner)" << std::right << std::setw(8) << "-"
					<< std::fixed << std::setprecision(4) << std::setw(13) << res[0].spin[a][0] << std::setw(12) << "-" << std::setw(12) << res[0].spin[a][0] << "\n";
		}
		for (size_t b = 0; b < B.size(); b++)
			if (!b_used[b])
				std::cout << std::setw(4) << "-" << " " << std::left << std::setw(18) << "(no partner)" << std::right << std::setw(5) << b + 1 << " " << std::left << std::setw(18) << res[1].labels[b] << std::right
					<< std::setw(8) << "-" << std::setw(13) << "-" << std::fixed << std::setprecision(4) << std::setw(12) << res[1].spin[b][1] << std::setw(12) << -res[1].spin[b][1] << "\n";
	}
	std::cout << "\nELIA (antiparallel pairs): not computed. For a single determinant the on-top pair density is exactly"
		" rho_alpha * rho_beta, so ELIA carries no pair information; it needs a correlated 2-matrix (-eli_family prints the details)." << std::endl;
	if (!warn.empty()) std::cout << "  WARNING: " << warn << std::endl;
}

void ELI_analysis(const WFN &wavy, options &opt) {
	err_checkf(wavy.get_ncen() != 0, "No Atoms in the wavefunction, this will not work!! ABORTING!!", std::cout);
	std::cout << "Analysing ELI basins in the wavefunction..." << std::endl;
	citations::cite(citations::Method::ELID, std::cout);
	density_field field;
	const std::unique_ptr<Gaussian_Molecule> fit = fitted_source(wavy, opt, field, std::cout);
	const density_field *fld = fit ? &field : nullptr;
	//A fitted source has no orbitals, so ELI comes from the density alone.
	if (fld)
		citations::cite(citations::Method::ELIOrbitalFree, std::cout);

	const double radius = opt.properties.radius;
	const double grid_spacing = opt.properties.resolution;
	properties_options prop_opt = opt.properties;
	WFN l_w = wavy;
	l_w.delete_unoccupied_MOs();
	//g = rho tau - |grad rho|^2 / 4 vanishes identically when one orbital carries the density: ELI-D is 0/0 and
	//its basins are round-off (DGrid shatters the same way). That is the definition, not the integrator, so warn.
	if (l_w.get_nmo() < 2)
		std::cout << "WARNING: this wavefunction has a single occupied orbital, so ELI-D's pair"
			" density g = rho*tau - |grad rho|^2/4 is identically zero and the field is undefined."
			" The basins below are round-off structure, not chemistry - expect thousands of them and"
			" do not quote their populations. The QTAIM basins are unaffected." << std::endl;
	//By default both basin sets are found and integrated on the analytic field (critical-point search,
	//analytic_eli_maxima, quadrature points walked up the field). The cube is built only for -basin_cube,
	//a fitted density (its attractors are not the orbitals' critical points) and a wavefunction without orbitals.
	const bool stream_qtaim = !opt.basin_cube && !fld && l_w.get_nmo() > 0;
	const bool stream_eli = !opt.basin_cube;
	const bool need_cube = !stream_qtaim || !stream_eli;
	const std::vector<atom> atoms = wavy.get_atoms();
	cube rho, eli_cube;
	basin_stage_timer T;
	if (need_cube) {
		readxyzMinMax_fromWFN(wavy, prop_opt);
		rho = cube(prop_opt.NbSteps, l_w.get_ncen(), true);
		eli_cube = cube(prop_opt.NbSteps, l_w.get_ncen(), true);
		rho.give_parent_wfn(l_w);
		eli_cube.give_parent_wfn(l_w);
		std::cout << "Calcualting grid of size " << prop_opt.NbSteps[0] << " x " << prop_opt.NbSteps[1] << " x " << prop_opt.NbSteps[2] << "..." << std::endl;
		std::cout << "Number of points: " << prop_opt.NbSteps[0] * prop_opt.NbSteps[1] * prop_opt.NbSteps[2] << std::endl;
		std::cout << "Grid parameters:\n";
		std::cout << "  Min: (" << prop_opt.MinMax[0] << ", " << prop_opt.MinMax[1] << ", " << prop_opt.MinMax[2] << ")\n";
		std::cout << "  Max: (" << prop_opt.MinMax[3] << ", " << prop_opt.MinMax[4] << ", " << prop_opt.MinMax[5] << ")\n";
		vec stepsizes{ (prop_opt.MinMax[3] - prop_opt.MinMax[0]) / prop_opt.NbSteps[0],
					  (prop_opt.MinMax[4] - prop_opt.MinMax[1]) / prop_opt.NbSteps[1],
					  (prop_opt.MinMax[5] - prop_opt.MinMax[2]) / prop_opt.NbSteps[2] };
		for (int i = 0; i < 3; i++) {
			rho.set_origin(i, prop_opt.MinMax[i]);
			rho.set_vector(i, i, stepsizes[i]);
			eli_cube.set_origin(i, prop_opt.MinMax[i]);
			eli_cube.set_vector(i, i, stepsizes[i]);
		}
		rho.calc_dv();
		eli_cube.calc_dv();
		Calc_RhoEli(rho, eli_cube, l_w, radius, fld);
		T.lap("rho and ELI-D cube");
	}

	// print table of atom positions
	std::cout << "Atom positions:\n";
	for (int a = 0; a < atoms.size(); a++) {
		std::cout << "  Atom " << a << ": (" << atoms[a].get_coordinate(0) << ", " << atoms[a].get_coordinate(1) << ", " << atoms[a].get_coordinate(2) << ")\n";
	}

	//An ECP took the core electrons out of the density. The QTAIM basins get them back from
	//Thakkar's spherical core densities, the fill the Hirshfeld grids and the scattering
	//factors apply: the nucleus is a cusp again and its basin holds the atom's full count.
	//ELI-D stays on the valence density the wavefunction has.
	ivec ecp_atoms;
	std::vector<Thakkar> ecp_cores;
	double ecp_electrons = 0.0;
	for (int a = 0; a < l_w.get_ncen(); a++)
		if (l_w.get_atom_ECP_electrons(a) > 0) {
			const int mode = l_w.get_ECP_mode() > 0 ? l_w.get_ECP_mode() : 1;
			ecp_atoms.push_back(a);
			ecp_cores.emplace_back(l_w.get_atom_charge(a), mode);
			ecp_electrons += l_w.get_atom_ECP_electrons(a);
		}
	std::function<double(const d3&)> core_density = [&](const d3 &p) {
		double s = 0.0;
		for (size_t k = 0; k < ecp_atoms.size(); k++) {
			const d3 ap = l_w.get_atom_pos(ecp_atoms[k]);
			const double d = std::sqrt(std::pow(p[0] - ap[0], 2) + std::pow(p[1] - ap[1], 2) + std::pow(p[2] - ap[2], 2));
			s += ecp_cores[k].get_core_density(d, l_w.get_atom_ECP_electrons(ecp_atoms[k]));
		}
		return s;
	};
	std::function<void(const d3&, d3&)> core_gradient = [&](const d3 &p, d3 &g) {
		const double h = 1e-4;
		for (int k = 0; k < 3; k++) {
			d3 a = p, b = p;
			a[k] += h; b[k] -= h;
			g[k] = (core_density(a) - core_density(b)) / (2.0 * h);
		}
	};
	const bool fill_cores = !ecp_atoms.empty();
	if (fill_cores) {
		std::cout << "ECP cores of " << ecp_atoms.size() << " atoms filled with Thakkar densities: " << std::fixed << std::setprecision(1) << ecp_electrons << " electrons added for the QTAIM basins" << std::endl;
		for (int x = 0; x < rho.get_size(0); x++)
			for (int y = 0; y < rho.get_size(1); y++)
				for (int z = 0; z < rho.get_size(2); z++)
					if (rho.get_value(x, y, z) > 0.0) rho.set_value(x, y, z, rho.get_value(x, y, z) + core_density(rho.get_pos(x, y, z)));
	}
	if (need_cube && opt.debug) {
		rho.set_path("rho.cube");
		eli_cube.set_path("eli.cube");
		rho.write_file(true);
		eli_cube.write_file(true);
	}

	//The critical points are refined on the orbitals (V, G and K need them) even when the
	//basins follow the fitted density; a SALTED prediction has none, so they are skipped. -topology's Newton
	//search on the analytic Hessian, checked against Poincare-Hopf: no grid decides which points exist
	std::vector<critical_point> density_critical_points;
	if (l_w.get_nmo() > 0) {
		//-topology's default 1.3 x (r_cov + r_cov) bond seeding misses hydrogen bonds and long contacts while
		//Poincare-Hopf still closes; 2.5 reaches O-H...O and C-H...O, and seeding and fragment count share it
		topology::options topo_opt;
		topo_opt.bond_scale = 2.5;
		const topology::result top = topology::analyze_topology(l_w, topology::nuclei_of(l_w), topo_opt);
		for (const topology::cp &p : top.points) {
			if (p.kind != topology::cp_kind::attractor && p.density < basin_density_cutoff) continue;
			//a saddle within a bohr of an ECP-cored nucleus sits in the pseudo-density's core hole and is no bond
			if (p.kind != topology::cp_kind::attractor && p.nearest_nucleus_distance < 1.0
				&& std::find(top.coreless_nuclei.begin(), top.coreless_nuclei.end(), p.nearest_nucleus) != top.coreless_nuclei.end()) continue;
			const critical_point_seed seed{ i3{ 0, 0, 0 }, p.position, p.density, p.gradient_norm, p.kind == topology::cp_kind::attractor && !p.is_nna, p.nearest_nucleus };
			density_critical_points.push_back(evaluate_critical_point(seed, p.position, l_w, p.iterations, true));
		}
		const section_log::tee cp_tee("qtaim", "Critical points of the electron density");
		std::cout << "Density critical points from the analytic field: NCP " << top.n_attractor << ", BCP " << top.n_bond << ", RCP "
			<< top.n_ring << ", CCP " << top.n_cage << "; Poincare-Hopf sum " << top.sum << " against " << top.target
			<< (top.complete ? " (COMPLETE)" : " (INCOMPLETE: " + top.diagnosis + ")") << std::endl;
	}
	else std::cout << "No orbitals: critical points (Hessian, V, G, K) need a wavefunction and are skipped." << std::endl;
	T.lap("density critical points");
	//Core shells make critical points of their own and an ECP atom a whole sphere of them,
	//none of which says anything about bonding and none of which any two machines find at
	//the same spots; only the nuclear attractor survives inside an atom's core radius. Sorted
	//by type, density and position so the listing reads the same everywhere.
	{
		std::vector<critical_point> kept;
		for (const critical_point &cp : density_critical_points) {
			bool core = false, nuclear = false;
			for (int a = 0; a < l_w.get_ncen(); a++) {
				const d3 apos = l_w.get_atom_pos(a);
				const double d2 = std::pow(cp.position[0] - apos[0], 2) + std::pow(cp.position[1] - apos[1], 2) + std::pow(cp.position[2] - apos[2], 2);
				if (d2 < 0.01) nuclear = true;
				else if (d2 < std::pow(core_shell_radius(l_w.get_atom_charge(a)), 2)) core = true;
			}
			if (!core || nuclear) kept.push_back(cp);
		}
		std::sort(kept.begin(), kept.end(), [](const critical_point &a, const critical_point &b) {
			if (a.type != b.type) return a.type < b.type;
			if (std::abs(a.density - b.density) > 1e-6 * std::max(1.0, std::abs(a.density))) return a.density > b.density;
			for (int k = 0; k < 3; k++)
				if (std::abs(a.position[k] - b.position[k]) > 1e-4) return a.position[k] < b.position[k];
			return false;
		});
		density_critical_points.swap(kept);
	}
	{
	const section_log::tee cp_tee("qtaim");
	std::cout << "Density Critical Points";
	if (!density_critical_points.empty())
		std::cout << " (" << density_critical_points.size() << " found)";
	std::cout << ":\n";
	if (density_critical_points.empty()) {
		std::cout << "  No critical points found.\n";
	}
	else {
		constexpr int nw = 15; // numeric field width — wide enough for large core densities
		for (size_t i = 0; i < density_critical_points.size(); i++) {
			const critical_point &cp = density_critical_points[i];

			std::string ownership_info;
			if (cp.type == "attractor" || cp.type == "bond") {
				double min_dist1 = std::numeric_limits<double>::max();
				double min_dist2 = std::numeric_limits<double>::max();
				int atom_index1 = -1;
				int atom_index2 = -1;
				for (int a = 0; a < l_w.get_ncen(); a++) {
					const d3 apos = l_w.get_atom_pos(a);
					const double dx = cp.position[0] - apos[0];
					const double dy = cp.position[1] - apos[1];
					const double dz = cp.position[2] - apos[2];
					const double dist2 = dx * dx + dy * dy + dz * dz;
					if (dist2 < min_dist1) {
						min_dist2 = min_dist1;
						atom_index2 = atom_index1;
						min_dist1 = dist2;
						atom_index1 = a;
					}
					else if (dist2 < min_dist2) {
						min_dist2 = dist2;
						atom_index2 = a;
					}
				}

				if (cp.type == "attractor" && atom_index1 >= 0) {
					ownership_info = "   " + l_w.get_atom_label(atom_index1) + std::to_string(atom_index1);
				}
				else if (cp.type == "bond" && atom_index1 >= 0 && atom_index2 >= 0) {
					ownership_info = "   "
						+ l_w.get_atom_label(atom_index1) + std::to_string(atom_index1)
						+ "-"
						+ l_w.get_atom_label(atom_index2) + std::to_string(atom_index2)
						+ " bond";
				}
			}

			std::cout << "\n  CP " << std::right << std::setw(3) << i + 1
				<< "  [" << std::left << std::setw(12) << cp.type << "]  "
				<< (cp.converged ? "converged" : "NOT converged")
				<< ownership_info << "\n";
			std::cout << std::fixed << std::setprecision(4) << std::right;
			std::cout << "    Position  :"
				<< std::setw(nw) << cp.position[0]
				<< std::setw(nw) << cp.position[1]
				<< std::setw(nw) << cp.position[2] << "\n";
			std::cout << "    Rho       :" << std::setw(nw) << cp.density << "\n";
			std::cout << "    GradRho   :"
				<< std::setw(nw) << cp.gradient[0]
				<< std::setw(nw) << cp.gradient[1]
				<< std::setw(nw) << cp.gradient[2]
				<< "   |GradRho| :" << std::setw(nw) << cp.gradient_norm << "\n";
			std::cout << "    HessRho_EigVals:"
				<< std::scientific << std::setprecision(4)
				<< std::setw(nw) << cp.hessian_eigenvalues[0]
				<< std::setw(nw) << cp.hessian_eigenvalues[1]
				<< std::setw(nw) << cp.hessian_eigenvalues[2]
				<< std::fixed << std::setprecision(4) << "\n";
			//The eigenvectors of a degenerate pair are any two in their plane and their signs
			//are free; both differ from machine to machine, so they are for -debug
			if (opt.debug) {
				std::cout << "    HessRho_EigVecs v1:"
					<< std::setw(nw) << cp.hessian_eigenvectors[0][0]
					<< std::setw(nw) << cp.hessian_eigenvectors[0][1]
					<< std::setw(nw) << cp.hessian_eigenvectors[0][2] << "\n";
				std::cout << "                    v2:"
					<< std::setw(nw) << cp.hessian_eigenvectors[1][0]
					<< std::setw(nw) << cp.hessian_eigenvectors[1][1]
					<< std::setw(nw) << cp.hessian_eigenvectors[1][2] << "\n";
				std::cout << "                    v3:"
					<< std::setw(nw) << cp.hessian_eigenvectors[2][0]
					<< std::setw(nw) << cp.hessian_eigenvectors[2][1]
					<< std::setw(nw) << cp.hessian_eigenvectors[2][2] << "\n";
			}
			std::cout << "    DelSqRho  :" << std::setw(nw) << cp.laplacian << "\n";
			if (std::isfinite(cp.ellipticity))
				std::cout << "    Bond Ellipticity:" << std::setw(nw) << cp.ellipticity << "\n";
			if (std::isfinite(cp.virial_field) || std::isfinite(cp.kinetic_lagrangian) || std::isfinite(cp.kinetic_hamiltonian) || std::isfinite(cp.lagrangian_density)) {
				std::cout << "   V         :" << std::scientific << std::setprecision(4) << std::setw(nw) << cp.virial_field
					<< "   G         :" << std::fixed << std::setprecision(4) << std::setw(nw) << cp.kinetic_lagrangian
					<< "   K         :" << std::scientific << std::setprecision(4) << std::setw(nw) << cp.kinetic_hamiltonian
					<< "   L         :" << std::scientific << std::setprecision(4) << std::setw(nw) << cp.lagrangian_density
					<< std::fixed << std::setprecision(4) << "\n";
			}
		}
		std::cout << "\n";
	}
	}

	//Every nucleus is a maximum of the density, whatever the grid says
	std::vector<d3> nuclei;
	for (const atom &a : atoms) nuclei.push_back(a.get_pos());
	//A fitted density oscillates about zero in the tail, and every positive island out there is
	//a maximum with no neighbour to merge into (2000 of them on epoxide at 0.05 A, and the
	//persistence merge is quadratic in their number), so the search stops where the fit is no
	//longer trusted; its ~1e-3 error at a bond critical point also leaves a bump there that
	//5e-3 persistence keeps, hence the looser merge
	const double floor = basin_density_cutoff, persistence = fld ? 2e-2 : 5e-3;
	//Density attractors come from the critical-point search already run, and the quadrature walks the field.
	//A fitted density keeps the cube: those critical points are the orbitals', not the fit's attractors.
	std::pair<cubei, std::vector<d4>> qtaim_results;
	if (stream_qtaim) {
		qtaim_results.second = streaming_density_attractors(l_w, density_critical_points, fill_cores ? &core_density : nullptr, fill_cores ? &core_gradient : nullptr, opt.debug);
		std::cout << "QTAIM attractors from the analytic field: " << qtaim_results.second.size() << " (" << l_w.get_ncen() << " nuclei, " << qtaim_results.second.size() - l_w.get_ncen() << " non-nuclear)" << std::endl;
	}
	else
		qtaim_results = topological_cube_analysis(&rho, atoms, opt.debug, true, floor, 1e-10, radius, persistence, &nuclei, &l_w, fill_cores ? &core_density : nullptr, fill_cores ? &core_gradient : nullptr, fld);
	T.lap("QTAIM attractors");
	svec labels = assign_labels_to_basins(qtaim_results.second, atoms, opt.debug);

	//Two integrations of the density over each basin set: the voxel sum, which is what the cube
	//resolution buys, and the atom-centred quadrature grids with the boundary decided by the
	//field itself, which is the number to compare with AIMAll and DGrid. The ELI-D basins follow
	//the orbitals' ELI-D and integrate the orbital density; only the QTAIM set uses the fit
	//A streaming ELI-D walks to every core shell's own maximum while the report keeps one core basin per atom,
	//so the integrator gets the unmerged maxima and the map to the basins
	std::vector<d4> eli_maxima_all;
	ivec eli_core_map;
	auto report = [&](const char *title, const std::pair<cubei, std::vector<d4>> &res, svec &lab, const bool eli, const bool stream) {
		//The voxel sum needs the basin cube; streaming has none
		if (!stream) {
			std::cout << "\n" << title << " (voxel sum):\n";
			integrate_values_in_basins(&rho, &(res.first), lab, opt.debug);
			T.lap(std::string(title) + " voxel sum");
		}
		vec vol;
		double outside = 0.0;
		//The overlap matrices cost one triangle per basin per thread; past the limit the delocalization indices are dropped, not the run
		basin_overlaps ovl;
		const bool orbitals = !eli && !fld && l_w.get_nmo() > 0;
		const size_t aom_bytes = orbitals ? (size_t)l_w.get_nmo() * (l_w.get_nmo() + 1) / 2 * res.second.size() * omp_get_max_threads() * sizeof(double) : 0;
		const bool want_aom = orbitals && aom_bytes < (size_t)512 * 1024 * 1024;
		if (orbitals && !want_aom)
			std::cout << "  Delocalization indices skipped: " << l_w.get_nmo() << " orbitals over " << res.second.size()
				<< " basins on " << omp_get_max_threads() << " threads would need " << aom_bytes / (1024 * 1024) << " MB of overlap matrices.\n";
		const bool mapped = eli && stream && !eli_core_map.empty();
		const vec pop = integrate_basins_on_atomic_grids(stream ? nullptr : &rho, stream ? nullptr : &(res.first), mapped ? eli_maxima_all : res.second, l_w, opt.accuracy, eli, vol, outside, fill_cores && !eli ? &core_density : nullptr, fill_cores && !eli ? &core_gradient : nullptr, opt.basin_grid, eli ? nullptr : fld, want_aom ? &ovl : nullptr, mapped ? &eli_core_map : nullptr);
		const section_log::tee basin_tee(eli ? "eli" : "qtaim", eli ? "ELI-D basins: populations, volumes and maxima"
			: "QTAIM basins: populations, charges, volumes and maxima");
		std::cout << "\n" << title << " (atomic quadrature grids):\n";
		if (!eli) citations::cite(citations::Method::QTAIM, std::cout);
		//The maximum column is 16 wide: it carries rho at the attractor, ~3e8 at a uranium nucleus
		std::cout << "  basin  label               electrons" << (eli ? "" : "     charge") << "      volume         maximum        x          y          z\n";
		double total = 0.0;
		for (size_t b = 0; b < pop.size(); b++) {
			total += pop[b];
			std::cout << std::setw(7) << b + 1 << "  " << std::left << std::setw(18) << lab[b] << std::right << std::fixed
				<< std::setprecision(4) << std::setw(11) << pop[b];
			if (!eli) {
				//The label names the atom the maximum sits on; a non-nuclear attractor has no charge
				double Z = -1.0;
				for (int a = 0; a < l_w.get_ncen(); a++)
					if (lab[b] == l_w.get_atom_label(a) + std::to_string(a)) Z = l_w.get_atom_charge(a);
				if (Z < 0) std::cout << std::setw(11) << "-";
				else std::cout << std::setw(11) << Z - pop[b];
			}
			std::cout << std::setw(12) << vol[b] << std::setw(16) << res.second[b][3]
				<< std::setprecision(3) << std::setw(11) << res.second[b][0] << std::setw(11) << res.second[b][1] << std::setw(11) << res.second[b][2] << "\n";
		}
		std::cout << "  total in basins: " << std::setprecision(4) << total << "   outside every basin: " << outside << "\n";
		//A metal's core and its outer core shell (eli_core_radius), each summed; a core far from its
		//closed-shell count (eli_core_electrons, less what an ECP took) is flagged
		if (eli)
			for (int a = 0; a < l_w.get_ncen(); a++) {
				const std::string name = l_w.get_atom_label(a) + std::to_string(a);
				double core = 0.0, shell = 0.0;
				int n_core = 0, n_shell = 0;
				for (size_t b = 0; b < pop.size(); b++) {
					if (lab[b] == name + " core") { core += pop[b]; n_core++; }
					else if (lab[b] == name + " shell") { shell += pop[b]; n_shell++; }
				}
				const int expected = std::max(0, eli_core_electrons(static_cast<int>(l_w.get_atom_charge(a))) - static_cast<int>(l_w.get_atom_ECP_electrons(a)));
				if (n_shell) std::cout << "  " << name << ": core " << std::setprecision(4) << core << " e (closed shells " << expected << "), outer core shell " << shell << " e in "
					<< n_shell << " basins, together " << core + shell << " e\n";
				if (n_core && std::abs(core - expected) > std::max(1.0, 0.15 * expected))
					std::cout << "  WARNING: " << name << " core holds " << std::setprecision(2) << core << " e, its closed shells " << expected
					<< ": the core boundary sits in the wrong shell minimum\n";
			}
		if (want_aom) {
			std::cout.flush();
			section_log::heading("qtaim", "Localization and delocalization indices");
			report_delocalization(l_w, ovl, lab, std::cout);
		}
	};
	if (l_w.get_nmo() == 0) {
		report("QTAIM Analysis", qtaim_results, labels, false, stream_qtaim);
		std::cout << "\nNo orbitals: the ELI-D basins are skipped." << std::endl;
		return;
	}

	std::pair<cubei, std::vector<d4>> eli_results;
	if (stream_eli) {
		//Gradient ascent on computeELIGrad from shells of seeds around every atom: no grid, and no Hessian (see below)
		eli_results.second = analytic_eli_maxima(l_w, opt.debug);
		eli_maxima_all = eli_results.second;
		T.lap("ELI-D maxima");
		std::cout << "ELI-D maxima from the analytic field: " << eli_results.second.size() << std::endl;
	}
	else {
	//ELI-D is undefined in the density tail, so the cube ends at the density isosurface.
	for (int x = 0; x < eli_cube.get_size(0); x++)
		for (int y = 0; y < eli_cube.get_size(1); y++)
			for (int z = 0; z < eli_cube.get_size(2); z++)
				if (rho.get_value(x, y, z) < basin_density_cutoff) eli_cube.set_value(x, y, z, 0.0);
	T.lap("ELI-D tail crop");
	//The only attractors found on the grid here (the QTAIM set is seeded from the nuclei), so for ELI-D the
	//resolution is an accuracy parameter, not a performance knob: coarser grids shift core populations and retype basins.
	if (grid_spacing > 0.05)
		std::cout << "WARNING: the ELI-D attractors are searched on the " << grid_spacing
			<< " A grid. Coarser than 0.05 A this basin set is less reliable: core populations have"
			" been seen to drift by 0.3 electrons and basins to appear or be retyped as lone pairs by"
			" 0.2 A. The QTAIM basins are unaffected." << std::endl;
	//The persistence merge absorbs a basin into its highest neighbour across their highest shared saddle. In a flat
	//valence shell every saddle is about as deep as the one to the core, so a loose threshold chains the shell into
	//the core basin. 3e-4 leaves only grid noise to it; unify_shell_basins merges shell maxima by distance (1.2 bohr)
	//and cannot chain into a core.
	//TODO: Na and Al cores come out too small; core_shell_radius does not reach the outer core shell.
	eli_results = topological_cube_analysis(&eli_cube, atoms, opt.debug, false, 0.0, 1e-10, radius, 3e-4);
	T.lap("ELI-D cube topology");
	}
	//There is no analytic ELI-D Hessian, only computeELIGrad, so maxima are not Newton-tested. Basin boundaries
	//come from walking every quadrature point up computeELIGrad, at any electron count.
	std::cout << "ELI-D basin boundaries from " << (stream_eli ? "the analytic field (the default)" : "the cube (-basin_cube)") << "." << std::endl;
	//The shells of a heavy atom's core structure ELI-D into several basins each; one core
	//basin per atom is what a bonding analysis wants, and what DGrid's ELIDcore gives
	const int core_merged = unify_core_basins(eli_results.first, eli_results.second, atoms, stream_eli ? &eli_core_map : nullptr);
	if (core_merged) std::cout << "Unified " << core_merged << " core-shell basins into their atoms' cores, " << eli_results.second.size() << " ELI-D basins remain." << std::endl;
	//Outside the cores the same sphere of maxima needs unify_shell_basins. NOS_ELI_SHELL_DIST / _TOL tune it;
	//distance 0 turns the merge off.
	double shell_dist = 1.2, shell_tol = 0.05;
	if (const char *e = tuning("NOS_ELI_SHELL_DIST")) { const double v = std::atof(e); if (v >= 0.0 && v < 10.0) shell_dist = v; }
	if (const char *e = tuning("NOS_ELI_SHELL_TOL")) { const double v = std::atof(e); if (v >= 0.0 && v < 1.0) shell_tol = v; }
	ivec eli_shell_map;
	const int shell_merged = unify_shell_basins(eli_results.first, eli_results.second, stream_eli ? &eli_shell_map : nullptr, shell_dist, shell_tol, &atoms, &l_w);
	if (shell_merged) std::cout << "Unified " << shell_merged << " shattered shell basins, " << eli_results.second.size() << " ELI-D basins remain." << std::endl;
	//The integrator walks to one of eli_maxima_all and then reads eli_core_map, so the second merge
	//has to be composed into that map rather than replacing it
	if (!eli_shell_map.empty())
		for (size_t b = 1; b < eli_core_map.size(); b++) eli_core_map[b] = eli_shell_map[eli_core_map[b]];
	svec eli_labels = assign_labels_to_basins(eli_results.second, atoms, opt.debug, 1);
	report("QTAIM Analysis", qtaim_results, labels, false, stream_qtaim);
	report("ELI-D Analysis", eli_results, eli_labels, true, stream_eli);
	const section_log::tee spin_tee("eli", "Spin-resolved ELI-D");
	spin_eli_analysis(l_w, opt, atoms, shell_dist, shell_tol);
}

// ---------------------------------------------------------------------------
// QTAIM_ELI_mask
// ---------------------------------------------------------------------------

void QTAIM_ELI_mask(
	cube& rho,
	cube& eli,
	WFN& parent_wfn,
	const std::vector<atom>& atoms,
	const ivec& selected_indices,
	double background_value,
	const std::filesystem::path& output_path,
	bool debug,
	std::ostream& log
) {
	err_checkf(!atoms.empty(), "QTAIM_ELI_mask: no atoms available.", log);
	err_checkf(rho.get_size(0) == eli.get_size(0) &&
			   rho.get_size(1) == eli.get_size(1) &&
			   rho.get_size(2) == eli.get_size(2),
			   "QTAIM_ELI_mask: rho and eli grids have different dimensions.", log);

	// Validate atom indices
	for (int idx : selected_indices) {
		err_checkf(idx >= 0 && idx < (int)atoms.size(),
			"QTAIM_ELI_mask: atom index " + std::to_string(idx) +
			" is out of range (0.." + std::to_string((int)atoms.size() - 1) + ").", log);
	}

	log << "Running QTAIM topological analysis on density grid..." << std::endl;
	auto [basin_cube, maxima] = topological_cube_analysis(&rho, atoms, debug, false, basin_density_cutoff);
	svec labels = assign_labels_to_basins(maxima, atoms, debug, 0);

	// Map selected atom indices → set of 1-based basin IDs
	// Label format: element + atom_index, e.g. "C0", "H3"
	std::set<int> selected_basins;
	for (int b = 0; b < (int)labels.size(); b++) {
		const std::string& lbl = labels[b];
		// Extract trailing digit sequence (the atom index)
		auto it = std::find_if(lbl.rbegin(), lbl.rend(),
							   [](char c) { return !std::isdigit(static_cast<unsigned char>(c)); });
		if (it == lbl.rend()) continue; // entire string is digits — skip
		const std::string num_str(it.base(), lbl.end());
		if (num_str.empty()) continue;
		int atom_idx = std::stoi(num_str);
		if (std::find(selected_indices.begin(), selected_indices.end(), atom_idx)
				!= selected_indices.end()) {
			selected_basins.insert(b + 1); // basin IDs in cubei are 1-indexed
			if (debug)
				log << "  Basin " << (b + 1) << " (\"" << lbl << "\") selected.\n";
		}
	}

	if (selected_basins.empty()) {
		log << "Warning: no QTAIM basins matched the requested atom indices. "
			   "Output will contain only the background value.\n";
	}

	const int nx = eli.get_size(0);
	const int ny = eli.get_size(1);
	const int nz = eli.get_size(2);

	// Pass 1: find tight bounding box of selected voxels
	int xmin = nx, xmax = -1, ymin = ny, ymax = -1, zmin = nz, zmax = -1;
	for (int x = 0; x < nx; x++)
		for (int y = 0; y < ny; y++)
			for (int z = 0; z < nz; z++) {
				if (selected_basins.count(basin_cube.get_value(x, y, z))) {
					xmin = std::min(xmin, x); xmax = std::max(xmax, x);
					ymin = std::min(ymin, y); ymax = std::max(ymax, y);
					zmin = std::min(zmin, z); zmax = std::max(zmax, z);
				}
			}

	if (xmax < 0) {
		// No selected voxels — write a 1×1×1 cube with background value
		xmin = xmax = 0; ymin = ymax = 0; zmin = zmax = 0;
	}

	log << "Bounding box of selected basins: ["
		<< xmin << ".." << xmax << "] x ["
		<< ymin << ".." << ymax << "] x ["
		<< zmin << ".." << zmax << "]\n";

	// Pass 2: build shrunk cube
	const std::array<int, 3> new_size = {xmax - xmin + 1, ymax - ymin + 1, zmax - zmin + 1};
	cube shrunk(new_size, (int)atoms.size(), true);
	shrunk.give_parent_wfn(parent_wfn);

	const auto new_origin = eli.get_pos(xmin, ymin, zmin);
	for (int i = 0; i < 3; i++) {
		shrunk.set_origin(i, new_origin[i]);
		for (int j = 0; j < 3; j++)
			shrunk.set_vector(i, j, eli.get_vector(i, j));
	}
	shrunk.calc_dv();

	for (int x = xmin; x <= xmax; x++)
		for (int y = ymin; y <= ymax; y++)
			for (int z = zmin; z <= zmax; z++) {
				double val = selected_basins.count(basin_cube.get_value(x, y, z))
							 ? eli.get_value(x, y, z)
							 : background_value;
				shrunk.set_value(x - xmin, y - ymin, z - zmin, val);
			}

	shrunk.set_comment1("QTAIM-masked ELI");
	shrunk.set_comment2("Selected atoms: " + [&]() {
		std::string s;
		for (int i = 0; i < (int)selected_indices.size(); i++) {
			if (i) s += ',';
			s += std::to_string(selected_indices[i]);
		}
		return s;
	}());
	shrunk.set_path(output_path);

	log << "Writing masked ELI cube to " << output_path.string() << " ..." << std::endl;
	err_checkf(shrunk.write_file(true), "QTAIM_ELI_mask: failed to write output cube.", log);
	log << "Done.\n";
}

// ---------------------------------------------------------------------------
// run_QTAIM_ELI_mask  (dispatch: cube-files mode vs WFN mode)
// ---------------------------------------------------------------------------

void run_QTAIM_ELI_mask(
	const std::filesystem::path& rho_or_wfn,
	const std::filesystem::path& eli_path,
	const ivec& selected_indices,
	double background_value,
	options& opt,
	std::ostream& log
) {
	if (!eli_path.empty()) {
		// ---- Cube-files mode ----
		log << "Reading density cube: " << rho_or_wfn.string() << std::endl;
		WFN dummy;
		cube rho(rho_or_wfn, true, dummy, log);
		log << "Reading ELI cube:     " << eli_path.string() << std::endl;
		cube eli(eli_path, true, dummy, log);
		rho.give_parent_wfn(dummy);
		eli.give_parent_wfn(dummy);

		const std::vector<atom> atoms = rho.get_parent_wfn_atoms();
		const std::filesystem::path out =
			eli_path.parent_path() / (eli_path.stem().string() + "_qtaim_masked.cube");

		QTAIM_ELI_mask(rho, eli, dummy, atoms, selected_indices, background_value, out, opt.debug, log);
	} else {
		// ---- WFN mode: compute rho and eli cubes first ----
		log << "Loading wavefunction: " << rho_or_wfn.string() << std::endl;
		WFN wavy(rho_or_wfn, opt.debug);
		wavy.delete_unoccupied_MOs();
		density_field field;
		const std::unique_ptr<Gaussian_Molecule> fit = fitted_source(wavy, opt, field, log);
		err_checkf(wavy.get_nmo() > 0, "ELI-D needs orbitals; -qtaim_eli cannot run on a SALTED density", log);

		properties_options prop_opt = opt.properties;
		readxyzMinMax_fromWFN(wavy, prop_opt);

		log << "Calculating density and ELI grid ("
			<< prop_opt.NbSteps[0] << " x " << prop_opt.NbSteps[1] << " x " << prop_opt.NbSteps[2]
			<< ") ..." << std::endl;

		cube rho(prop_opt.NbSteps, wavy.get_ncen(), true);
		cube eli(prop_opt.NbSteps, wavy.get_ncen(), true);
		rho.give_parent_wfn(wavy);
		eli.give_parent_wfn(wavy);

		const vec stepsizes{
			(prop_opt.MinMax[3] - prop_opt.MinMax[0]) / prop_opt.NbSteps[0],
			(prop_opt.MinMax[4] - prop_opt.MinMax[1]) / prop_opt.NbSteps[1],
			(prop_opt.MinMax[5] - prop_opt.MinMax[2]) / prop_opt.NbSteps[2]
		};
		for (int i = 0; i < 3; i++) {
			rho.set_origin(i, prop_opt.MinMax[i]);
			rho.set_vector(i, i, stepsizes[i]);
			eli.set_origin(i, prop_opt.MinMax[i]);
			eli.set_vector(i, i, stepsizes[i]);
		}
		rho.calc_dv();
		eli.calc_dv();

		Calc_RhoEli(rho, eli, wavy, prop_opt.radius, fit ? &field : nullptr);

		const std::vector<atom> atoms = wavy.get_atoms();
		const std::filesystem::path out =
			rho_or_wfn.parent_path() / "eli_qtaim_masked.cube";

		QTAIM_ELI_mask(rho, eli, wavy, atoms, selected_indices, background_value, out, opt.debug, log);
	}
}
