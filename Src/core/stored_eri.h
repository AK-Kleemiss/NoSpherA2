#pragma once
#include <occ/qm/hf.h>
#include <functional>
#include <memory>
#include <ostream>
#include <vector>

//(ab|cd) over the pairs a >= b surviving OCC's shell-pair screen and a Schwarz screen, packed over the
//8-fold symmetry: kept pairs numbered in (a, b) order, integrals the lower triangle over that numbering,
//k >= l at k(k+1)/2 + l. The pairs of first index c run first_pair()[c] .. first_pair()[c + 1], so a
//row's part over c is a contiguous segment with second indices pair_b(); with every pair kept, d = 0..c.
//Built once when the integrals fit the budget; Fock builds then contract them instead of recomputing.
class stored_eri {
public:
	//OCC's quartet threshold: drops a pair whose Schwarz bound times the largest falls below it, and
	//in a screened contraction skips a segment whose bound times the density difference does
	static constexpr double threshold = 1e-12;
	using jk_fn = std::function<void(const occ::Mat&, occ::Mat&, occ::Mat&)>;

	//False, holding nothing, over budget_bytes (0: no budget) or when hf's Fock build is not the quartet one
	bool build(const occ::qm::HartreeFock& hf, size_t budget_bytes, std::ostream& log);
	void clear();
	explicit operator bool() const { return static_cast<bool>(v_); }

	//J_ab = sum_cd (ab|cd) D_cd, K_ab = sum_cd (ac|bd) D_cd. With screen, D is a density difference and a
	//segment whose Schwarz bound times the largest D element it touches is below threshold is skipped
	void JK(const occ::Mat& D, occ::Mat& J, occ::Mat& K, bool screen = false) const;
	//Two-electron Fock part for OCC's half-scaled densities: 2J(D) - K(D) restricted, per spin
	//2J(Da + Db) - 2K(Ds) unrestricted; jk replaces JK(), for a device holding the integrals
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
