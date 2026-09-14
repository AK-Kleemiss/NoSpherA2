#pragma once
#include <occ/qm/hf.h>
#include <functional>
#include <memory>
#include <ostream>
#include <vector>

//Two-electron integrals (ab|cd) over the basis-function pairs a >= b that survive OCC's
//shell-pair screen and a Schwarz screen, packed over the 8-fold symmetry. The kept pairs are
//numbered in (a, b) order and the integrals are the lower triangle over that numbering, slot
//of k >= l at k(k+1)/2 + l. The pairs of one first index c are consecutive, first_pair()[c]
//to first_pair()[c + 1], so the part of a row over c is a contiguous segment whose second
//indices are pair_b(). With every pair kept the layout is the plain packed one and each
//segment is the run d = 0..c. Built once per run when the integrals fit the budget; a Fock
//build then contracts them instead of recomputing every quartet per iteration.
class stored_eri {
public:
	//OCC's quartet threshold: a pair whose Schwarz bound times the largest one falls below it
	//is dropped from the store, and a stored segment whose bound times the density difference
	//falls below it is skipped in a screened contraction
	static constexpr double threshold = 1e-12;
	using jk_fn = std::function<void(const occ::Mat&, occ::Mat&, occ::Mat&)>;

	//Computes and packs the integrals; false, holding nothing, when they exceed budget_bytes
	//(0 for no budget) or when hf's Fock build is not the quartet one
	bool build(const occ::qm::HartreeFock& hf, size_t budget_bytes, std::ostream& log);
	void clear();
	explicit operator bool() const { return static_cast<bool>(v_); }

	//J_ab = sum_cd (ab|cd) D_cd and K_ab = sum_cd (ac|bd) D_cd. With screen set D is a density
	//difference and every segment whose Schwarz bound times the largest difference element it
	//touches stays below the threshold is skipped, the way OCC's direct build skips shell
	//quartets on the density's shell-block norms
	void JK(const occ::Mat& D, occ::Mat& J, occ::Mat& K, bool screen = false) const;
	//OCC's two-electron Fock part for its half-scaled densities: 2J(D) - K(D) restricted, and
	//per spin 2J(Da + Db) - 2K(Ds) unrestricted. A jk of the caller's replaces JK(), for a
	//device that holds the integrals
	occ::Mat fock(const occ::qm::MolecularOrbitals& mo, bool screen = false, const jk_fn& jk = {}) const;

	int nbf() const { return nbf_; }
	int npairs() const { return static_cast<int>(pa_.size()); }
	bool dense() const { return npairs() == nbf_ * (nbf_ + 1) / 2; }
	size_t bytes() const { return sizeof(double) * (size_t)npairs() * (npairs() + 1) / 2; }
	const double* data() const { return v_.get(); }
	const std::vector<int>& pair_a() const { return pa_; }
	const std::vector<int>& pair_b() const { return pb_; }
	const std::vector<int>& first_pair() const { return first_; }
	//Kept-pair number of the packed pair a(a+1)/2 + b, -1 for a dropped pair
	const std::vector<int>& pair_index() const { return idx_; }
	//Schwarz bound sqrt((ab|ab)) of each kept pair
	const std::vector<double>& schwarz() const { return q_; }

private:
	int nbf_ = 0;
	std::unique_ptr<double[]> v_;
	std::vector<int> pa_, pb_, first_, idx_;
	//Per kept pair, per first index (the largest of its segment) and overall
	std::vector<double> q_, qseg_;
	double qmax_ = 0.0;
};
