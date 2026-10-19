#include "pch.h"
#include "ccsd_density.h"
#include "convenience.h"

#include <occ/core/diis.h>
#include <occ/core/util.h>
#include <occ/driver/method_parser.h>
#include <occ/main/occ_scf.h>
#include <occ/qm/cc/ccsd.h>
#include <occ/qm/cc/integrals.h>
#include <occ/qm/fitting_basis.h>
#include <Eigen/Dense>
#include <toml++/toml.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

#include <array>
#include <iomanip>
#include <ostream>

namespace
{
	//Dense row-major tensor of rank <= 6 with numpy's einsum and transpose, so that the Lambda and density
	//equations below read term for term like PySCF's ccsd_lambda.py / ccsd_rdm.py and their (T) versions
	//(restricted, chemist notation, t2 and l2 as ijab, t3-like arrays as ijkabc) they were ported from.
	//Every contraction is a transpose and one GEMM.
	struct Tn
	{
		std::vector<int> d;
		std::vector<double> v;
		Tn() = default;
		explicit Tn(const std::vector<int> &dims) : d(dims), v(count(dims), 0.0) {}
		static size_t count(const std::vector<int> &dims)
		{
			size_t n = 1;
			for (const int x : dims) n *= static_cast<size_t>(x);
			return n;
		}
		void same_shape(const Tn &o) const
		{
			if (o.d != d) throw std::runtime_error("ccsd_density: tensor shapes differ");
		}
		Tn &operator+=(const Tn &o) { same_shape(o); for (size_t i = 0; i < v.size(); i++) v[i] += o.v[i]; return *this; }
		Tn &operator-=(const Tn &o) { same_shape(o); for (size_t i = 0; i < v.size(); i++) v[i] -= o.v[i]; return *this; }
		Tn &operator*=(const double s) { for (double &x : v) x *= s; return *this; }
		double &at(const int i, const int j) { return v[static_cast<size_t>(i) * d[1] + j]; }
		double at(const int i, const int j) const { return v[static_cast<size_t>(i) * d[1] + j]; }
		double norm() const
		{
			double s = 0.0;
			for (const double x : v) s += x * x;
			return std::sqrt(s);
		}
	};
	Tn operator+(Tn a, const Tn &b) { a += b; return a; }
	Tn operator-(Tn a, const Tn &b) { a -= b; return a; }
	Tn operator*(Tn a, const double s) { a *= s; return a; }
	Tn divided(Tn a, const Tn &d) { a.same_shape(d); for (size_t i = 0; i < a.v.size(); i++) a.v[i] /= d.v[i]; return a; }
	double dot(const Tn &a, const Tn &b)
	{
		a.same_shape(b);
		double s = 0.0;
		for (size_t i = 0; i < a.v.size(); i++) s += a.v[i] * b.v[i];
		return s;
	}

	//numpy a.transpose(perm): out.d[k] = a.d[perm[k]]
	Tn T(const Tn &a, std::vector<int> perm)
	{
		constexpr int R = 6;
		const int r = static_cast<int>(a.d.size());
		std::vector<int> ad(a.d);
		ad.resize(R, 1);
		for (int k = r; k < R; k++) perm.push_back(k);
		std::array<size_t, R> stride{};
		size_t s = 1;
		for (int k = R - 1; k >= 0; k--) { stride[k] = s; s *= static_cast<size_t>(ad[k]); }
		std::vector<int> od(r);
		for (int k = 0; k < r; k++) od[k] = a.d[perm[k]];
		Tn out(od);
		std::array<int, R> n{};
		std::array<size_t, R> st{};
		for (int k = 0; k < R; k++) { n[k] = ad[perm[k]]; st[k] = stride[perm[k]]; }
		size_t w = 0;
		for (int x0 = 0; x0 < n[0]; x0++)
			for (int x1 = 0; x1 < n[1]; x1++)
				for (int x2 = 0; x2 < n[2]; x2++)
					for (int x3 = 0; x3 < n[3]; x3++)
						for (int x4 = 0; x4 < n[4]; x4++)
						{
							const size_t base = x0 * st[0] + x1 * st[1] + x2 * st[2] + x3 * st[3] + x4 * st[4];
							for (int x5 = 0; x5 < n[5]; x5++) out.v[w++] = a.v[base + x5 * st[5]];
						}
		return out;
	}

	//numpy einsum of two operands: a label in both inputs and not in the output is summed; a label may not
	//repeat within one operand, nor be shared by both inputs and the output (no batch or trace terms here)
	Tn einsum(const std::string &spec, const Tn &a, const Tn &b)
	{
		const size_t comma = spec.find(','), arrow = spec.find("->");
		const std::string la = spec.substr(0, comma), lb = spec.substr(comma + 1, arrow - comma - 1), lo = spec.substr(arrow + 2);
		if (la.size() != a.d.size() || lb.size() != b.d.size())
			throw std::runtime_error("ccsd_density: einsum '" + spec + "' does not match the operand ranks");
		std::string con, fa, fb;
		for (const char c : la) (lb.find(c) != std::string::npos && lo.find(c) == std::string::npos ? con : fa) += c;
		for (const char c : lb) if (con.find(c) == std::string::npos) fb += c;
		for (const char c : fa)
			if (fb.find(c) != std::string::npos) throw std::runtime_error("ccsd_density: einsum '" + spec + "' has a batch label");
		std::vector<int> pa, pb, cd;
		size_t M = 1, K = 1, N = 1;
		for (const char c : fa) { pa.push_back(static_cast<int>(la.find(c))); cd.push_back(a.d[la.find(c)]); M *= a.d[la.find(c)]; }
		for (const char c : con)
		{
			if (a.d[la.find(c)] != b.d[lb.find(c)]) throw std::runtime_error("ccsd_density: einsum '" + spec + "' sums over unequal dimensions");
			pa.push_back(static_cast<int>(la.find(c)));
			pb.push_back(static_cast<int>(lb.find(c)));
			K *= a.d[la.find(c)];
		}
		for (const char c : fb) { pb.push_back(static_cast<int>(lb.find(c))); cd.push_back(b.d[lb.find(c)]); N *= b.d[lb.find(c)]; }
		const Tn at = T(a, pa), bt = T(b, pb);
		Tn c(cd);
		using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
		Eigen::Map<RowMat>(c.v.data(), M, N).noalias() =
			Eigen::Map<const RowMat>(at.v.data(), M, K) * Eigen::Map<const RowMat>(bt.v.data(), K, N);
		const std::string lc = fa + fb;
		std::vector<int> po;
		for (const char x : lo) po.push_back(static_cast<int>(lc.find(x)));
		return T(c, po);
	}

	Tn from_tensor(const Eigen::Tensor<double, 2> &t)
	{
		Tn out({ static_cast<int>(t.dimension(0)), static_cast<int>(t.dimension(1)) });
		for (int i = 0; i < out.d[0]; i++)
			for (int j = 0; j < out.d[1]; j++) out.at(i, j) = t(i, j);
		return out;
	}

	Tn from_tensor(const Eigen::Tensor<double, 4> &t)
	{
		Tn out({ static_cast<int>(t.dimension(0)), static_cast<int>(t.dimension(1)), static_cast<int>(t.dimension(2)), static_cast<int>(t.dimension(3)) });
		size_t w = 0;
		for (int i = 0; i < out.d[0]; i++)
			for (int j = 0; j < out.d[1]; j++)
				for (int k = 0; k < out.d[2]; k++)
					for (int l = 0; l < out.d[3]; l++) out.v[w++] = t(i, j, k, l);
		return out;
	}

	Eigen::Tensor<double, 4> to_tensor(const Tn &a)
	{
		Eigen::Tensor<double, 4> t(a.d[0], a.d[1], a.d[2], a.d[3]);
		size_t w = 0;
		for (int i = 0; i < a.d[0]; i++)
			for (int j = 0; j < a.d[1]; j++)
				for (int k = 0; k < a.d[2]; k++)
					for (int l = 0; l < a.d[3]; l++) t(i, j, k, l) = a.v[w++];
		return t;
	}

	Tn block(const occ::Mat &m, const int r0, const int c0, const int nr, const int nc)
	{
		Tn out({ nr, nc });
		for (int i = 0; i < nr; i++)
			for (int j = 0; j < nc; j++) out.at(i, j) = m(r0 + i, c0 + j);
		return out;
	}

	//occ's integral blocks in the row-major form the equations use, plus the vvvv ladder
	struct Eris
	{
		int o = 0, v = 0;
		Tn foo, fov, fvv, oooo, ovoo, oovv, ovov, ovvo, ovvv;
		occ::Vec eo, ev;
		std::function<Eigen::Tensor<double, 4>(const Eigen::Tensor<double, 4> &)> ladder;
	};

	Eris convert(const occ::qm::cc::CCIntegrals &e)
	{
		Eris r;
		r.o = e.nocc;
		r.v = e.nvir;
		r.eo = e.mo_energy.head(e.nocc);
		r.ev = e.mo_energy.tail(e.nvir);
		occ::Mat f = e.fock;
		if (f.rows() != e.nocc + e.nvir) f = e.mo_energy.asDiagonal();
		r.foo = block(f, 0, 0, r.o, r.o);
		r.fov = block(f, 0, r.o, r.o, r.v);
		r.fvv = block(f, r.o, r.o, r.v, r.v);
		r.oooo = from_tensor(e.oooo);
		r.ovoo = from_tensor(e.ovoo);
		r.oovv = from_tensor(e.oovv);
		r.ovov = from_tensor(e.ovov);
		r.ovvo = from_tensor(e.ovvo);
		r.ovvv = from_tensor(e.ovvv);
		r.ladder = e.ladder;
		return r;
	}

	struct Imds
	{
		Tn w1, w2, w3, w4, woooo, wvooo, wVOov, wvOOv, wvvov;
	};

	//PySCF ccsd_lambda.make_intermediates, unblocked
	Imds make_intermediates(const Tn &t1, const Tn &t2, const Eris &e)
	{
		const Tn &ovvv = e.ovvv, &ovoo = e.ovoo, &ovvo = e.ovvo, &oovv = e.oovv, &oooo = e.oooo, &fov = e.fov;
		Tn w1 = e.fvv - einsum("ja,jb->ba", fov, t1);
		Tn w2 = e.foo + einsum("ib,jb->ij", fov, t1);
		Tn w3 = einsum("kc,jkbc->bj", fov, t2) * 2.0 + T(fov, { 1, 0 });
		w3 -= einsum("kc,kjbc->bj", fov, t2);
		w3 += einsum("cb,jc->bj", einsum("kc,kb->cb", fov, t1), t1);
		Tn w4 = fov;

		const Tn vvov = T(ovvv, { 2, 3, 0, 1 });
		w1 += einsum("jcba,jc->ba", ovvv, t1 * 2.0);
		w1 -= einsum("jabc,jc->ba", ovvv, t1);
		Tn theta = t2 * 2.0 - T(t2, { 1, 0, 2, 3 });
		w3 += einsum("jkcd,kdcb->bj", theta, ovvv);
		Tn wVOov = einsum("jbcd,kd->bjkc", ovvv, t1);
		Tn wvOOv = einsum("cbjd,kd->cjkb", vvov, t1 * -1.0);
		const Tn tau = t2 + einsum("ia,jb->ijab", t1, t1);
		Tn wvooo = einsum("cibd,jkbd->ckij", T(vvov, { 0, 2, 1, 3 }) * 2.0 - T(vvov, { 0, 2, 3, 1 }), tau);

		Tn wvvov = einsum("jabd,jkcd->abkc", ovvv, t2) * -1.5;
		wvvov += T(vvov, { 0, 3, 2, 1 }) * 2.0;
		wvvov -= vvov;
		{
			const Tn vackb = einsum("acjd,kjbd->ackb", vvov * 2.0 - T(ovvv, { 1, 2, 0, 3 }), t2 * 2.0 - T(t2, { 0, 1, 3, 2 }));
			wvvov += T(vackb, { 0, 3, 2, 1 });
			wvvov -= vackb * 0.5;
		}

		w2 += einsum("kbij,kb->ij", ovoo, t1) * 2.0;
		w2 -= einsum("ibkj,kb->ij", ovoo, t1);
		theta = T(t2, { 1, 0, 2, 3 }) * 2.0 - t2;
		w3 -= einsum("lckj,klcb->bj", ovoo, theta);

		Tn woooo = einsum("lc,jcik->ijkl", t1, ovoo);
		woooo += T(woooo, { 1, 0, 3, 2 });
		wvOOv += einsum("lbjk,lc->bjkc", ovoo, t1);
		wVOov -= einsum("jbkl,lc->bjkc", ovoo, t1);
		wvooo += T(ovoo, { 1, 3, 2, 0 }) * 2.0;
		wvooo -= T(ovoo, { 1, 0, 2, 3 });
		wvooo -= einsum("klbc,iblj->ckij", t2, ovoo * 1.5);
		const Tn g2ovoo = ovoo * 2.0 - T(ovoo, { 2, 1, 0, 3 });
		{
			const Tn vcjik = einsum("jlcb,lbki->cjki", t2 * 2.0 - T(t2, { 1, 0, 2, 3 }), g2ovoo);
			wvooo += T(vcjik, { 0, 3, 2, 1 });
			wvooo -= vcjik * 0.5;
		}

		const Tn voov = T(ovvo, { 1, 0, 3, 2 });
		woooo += einsum("cijd,klcd->ijkl", voov, tau);
		const Tn g2voov = voov * 2.0 - T(voov, { 0, 2, 1, 3 });
		const Tn tmpw4 = einsum("ckld,ld->kc", g2voov, t1);
		w1 -= einsum("ckja,kjcb->ba", g2voov, t2);
		w1 -= einsum("ja,jb->ba", tmpw4, t1);
		w2 += einsum("jkbc,bikc->ij", t2, g2voov);
		w2 += einsum("ib,jb->ij", tmpw4, t1);
		w3 += einsum("bc,jc->bj", einsum("kb,kc->bc", t1, tmpw4), t1);
		w4 += tmpw4;

		wvOOv += einsum("bljk,lc->bjkc", einsum("bljd,kd->bljk", voov, t1), t1);
		wVOov -= einsum("bjlk,lc->bjkc", einsum("bjld,kd->bjlk", voov, t1), t1);
		Tn VOov = einsum("bjld,klcd->bjkc", g2voov, t2);
		VOov -= einsum("bjld,kldc->bjkc", voov, t2);
		VOov += voov;
		Tn vOOv = einsum("bljd,kldc->bjkc", voov, t2);
		vOOv -= T(oovv, { 2, 1, 0, 3 });
		wVOov += VOov;
		wvOOv += vOOv;

		const Tn ov1 = vOOv * 2.0 + VOov, ov2 = VOov * 2.0 + vOOv;
		wvooo -= einsum("jb,bikc->ckij", t1, ov1);
		wvooo += einsum("kb,bijc->ckij", t1, ov2);
		w3 += einsum("ckjb,kc->bj", ov2, t1);
		wvvov += einsum("ajkc,jb->abkc", ov1, t1);
		wvvov -= einsum("ajkb,jc->abkc", ov2, t1);
		wvvov += einsum("laki,klbc->abic", g2ovoo, tau);

		woooo += T(oooo, { 0, 2, 1, 3 });
		w3 += einsum("bc,jc->bj", w1, t1);
		w3 -= einsum("kj,kb->bj", w2, t1);
		return { w1, w2, w3, w4, woooo, wvooo, wVOov, wvOOv, wvvov };
	}

	//PySCF ccsd_lambda.update_lambda, unblocked: the next (l1, l2) from the current ones
	void update_lambda(const Tn &t1, const Tn &t2, const Tn &l1, const Tn &l2, const Eris &e, const Imds &m, Tn &l1new, Tn &l2new)
	{
		const int o = e.o, v = e.v;
		const Tn theta = t2 * 2.0 - T(t2, { 0, 1, 3, 2 });
		const Tn mba = einsum("klca,klcb->ba", l2, theta);
		const Tn mij = einsum("ikcd,jkcd->ij", l2, theta);
		const Tn mba1 = einsum("jc,jb->bc", l1, t1) + mba;
		const Tn mij1 = einsum("kb,jb->kj", l1, t1) + mij;
		Tn mia1 = t1 + einsum("kc,jkbc->jb", l1, t2) * 2.0;
		mia1 -= einsum("kc,jkcb->jb", l1, t2);
		mia1 -= einsum("jk,kb->jb", einsum("ja,ka->jk", t1, l1), t1);
		mia1 -= einsum("bd,jd->jb", mba, t1);
		mia1 -= einsum("lj,lb->jb", mij, t1);

		l2new = from_tensor(e.ladder(to_tensor(l2))); //sum_cd (ac|bd) l2(i,j,c,d)
		l1new = einsum("ijab,jb->ia", l2new, t1) * 2.0;
		l1new -= einsum("jiab,jb->ia", l2new, t1);
		l2new *= 0.5; //the (ia)<->(jb) symmetrisation at the end doubles it back

		Tn w1 = m.w1, w2 = m.w2;
		for (int a = 0; a < v; a++) w1.at(a, a) -= e.ev(a);
		for (int i = 0; i < o; i++) w2.at(i, i) -= e.eo(i);

		l1new += e.fov;
		l1new += einsum("ib,ba->ia", l1, w1);
		l1new -= einsum("ja,ij->ia", l1, w2);
		l1new -= einsum("ik,ka->ia", mij, m.w4);
		l1new -= einsum("ca,ic->ia", mba, m.w4);
		l1new += einsum("ijab,bj->ia", l2, m.w3) * 2.0;
		l1new -= einsum("ijba,bj->ia", l2, m.w3);

		l2new += einsum("ia,jb->ijab", l1, m.w4);
		l2new += einsum("jibc,ca->jiba", l2, w1);
		l2new -= einsum("jk,kiba->jiba", w2, l2);

		l1new -= einsum("iajk,kj->ia", e.ovoo, mij1) * 2.0;
		l1new += einsum("jaik,kj->ia", e.ovoo, mij1);
		l2new -= einsum("jbki,ka->jiba", e.ovoo, l1);

		const Tn tau = t2 + einsum("ia,jb->ijab", t1, t1);
		const Tn l2tau = einsum("ijcd,klcd->ijkl", l2, tau);
		const Tn l2t1 = einsum("jidc,kc->ijkd", l2, t1);

		l1new -= einsum("jb,jiab->ia", l1, e.oovv);
		l1new += einsum("iabc,bc->ia", e.ovvv, mba1) * 2.0;
		l1new -= einsum("ibca,bc->ia", e.ovvv, mba1);
		l2new += einsum("jbac,ic->jiba", e.ovvv, l1);
		{
			const Tn m4 = einsum("ijkd,kadb->ijab", l2t1, e.ovvv);
			l2new -= m4;
			l1new -= einsum("ijab,jb->ia", m4, t1) * 2.0;
			l1new -= einsum("ijab,ia->jb", m4, t1) * 2.0;
			l1new += einsum("jiab,jb->ia", m4, t1);
			l1new += einsum("jiab,ia->jb", m4, t1);
		}

		const Tn voov = T(e.ovvo, { 1, 0, 3, 2 });
		l1new += einsum("jb,aijb->ia", l1, voov) * 2.0;
		l2new += T(voov, { 1, 2, 0, 3 }) * 0.5;
		l2new -= einsum("bjic,ca->jiba", voov, mba1);
		l2new -= einsum("bjka,ik->jiba", voov, mij1);
		l1new += einsum("aijb,jb->ia", voov, mia1) * 2.0;
		l1new -= einsum("bija,jb->ia", voov, mia1);
		{
			const Tn m4 = einsum("ijkl,aklb->ijab", l2tau, voov);
			l2new += m4 * 0.5;
			l1new += einsum("ijab,jb->ia", m4, t1) * 2.0;
			l1new -= einsum("ijba,jb->ia", m4, t1);
		}

		l1new -= einsum("ckij,jkca->ia", m.wvooo, l2);
		l1new += einsum("abkc,kibc->ia", m.wvvov, l2);
		l2new += einsum("iakc,bjkc->jiba", T(l2, { 0, 2, 1, 3 }) - T(l2, { 0, 3, 1, 2 }) * 0.5, m.wVOov * 2.0 + m.wvOOv);
		{
			const Tn tmp = einsum("jkca,bikc->jiba", l2, m.wvOOv);
			l2new += tmp;
			l2new += T(tmp, { 1, 0, 2, 3 }) * 0.5;
		}
		{
			const Tn m3 = einsum("ijkl,klab->ijab", m.woooo, l2);
			l2new += m3 * 0.5;
			l1new += einsum("ijab,jb->ia", m3, t1) * 2.0;
			l1new -= einsum("ijba,jb->ia", m3, t1);
		}

		for (int i = 0; i < o; i++)
			for (int a = 0; a < v; a++) l1new.at(i, a) /= e.eo(i) - e.ev(a);
		l2new += T(l2new, { 1, 0, 3, 2 });
		size_t w = 0;
		for (int i = 0; i < o; i++)
			for (int j = 0; j < o; j++)
				for (int a = 0; a < v; a++)
					for (int b = 0; b < v; b++) l2new.v[w++] /= e.eo(i) - e.ev(a) + e.eo(j) - e.ev(b);
	}

	//PySCF ccsd_t_lambda_slow / ccsd_t_rdm_slow: the triples terms that make the Lambda equations and the
	//density those of CCSD(T). w, v are the p6-symmetrised t3-like arrays over ijkabc, already divided by D.
	struct Triples
	{
		Tn l1_t, l2_t;  //constant terms added to every Lambda update
		Tn goo, gvv, dvo;
		double e_t = 0.0;
	};

	Tn p6(const Tn &t)
	{
		const Tn t1 = t + T(t, { 0, 2, 1, 3, 5, 4 });
		return t1 + T(t1, { 1, 0, 2, 4, 3, 5 }) + T(t1, { 1, 2, 0, 4, 5, 3 });
	}

	Tn r6(const Tn &w)
	{
		return w * 4.0 + T(w, { 0, 1, 2, 4, 5, 3 }) + T(w, { 0, 1, 2, 5, 3, 4 })
			- T(w, { 0, 1, 2, 5, 4, 3 }) * 2.0 - T(w, { 0, 1, 2, 3, 5, 4 }) * 2.0 - T(w, { 0, 1, 2, 4, 3, 5 }) * 2.0;
	}

	//r6 combined with the 2/3 F_ijab + 1/3 F_jiab form the RCCSD Lambda equations are solved in
	Tn as_r6(const Tn &m) { return m * 2.0 - T(m, { 0, 1, 2, 5, 4, 3 }) - T(m, { 0, 1, 2, 3, 5, 4 }); }

	//shortcut: w, v and their combinations are full o^3 v^3 arrays (PySCF's *_slow form); batch over virtual
	//triples as occ's ccsd_t does before this is wanted past ~10 correlated occupied and ~100 virtual orbitals
	Triples triples(const Tn &t1, const Tn &t2, const Eris &e)
	{
		const int o = e.o, v = e.v;
		Tn eia({ o, v });
		for (int i = 0; i < o; i++)
			for (int a = 0; a < v; a++) eia.at(i, a) = e.eo(i) - e.ev(a);
		Tn d3({ o, o, o, v, v, v });
		size_t n = 0;
		for (int i = 0; i < o; i++)
			for (int j = 0; j < o; j++)
				for (int k = 0; k < o; k++)
					for (int a = 0; a < v; a++)
						for (int b = 0; b < v; b++)
							for (int c = 0; c < v; c++) d3.v[n++] = eia.at(i, a) + eia.at(j, b) + eia.at(k, c);
		const Tn fvo = T(e.fov, { 1, 0 });
		const Tn w = p6(divided(einsum("iafb,kjcf->ijkabc", e.ovvv, t2) - einsum("iajm,mkbc->ijkabc", e.ovoo, t2), d3));
		const Tn vv = p6(divided(einsum("iajb,kc->ijkabc", e.ovov, t1) + einsum("ck,ijab->ijkabc", fvo, t2), d3));

		Triples r;
		const Tn rw = r6(w);
		r.l1_t = divided(einsum("jbkc,ijkabc->ia", e.ovov, rw), eia) * 0.5;
		{
			const Tn m = as_r6(w * 2.0 + vv * 0.5);
			Tn joovv = einsum("kfbe,ijkaef->ijab", e.ovvv, m) - einsum("ncmj,imnabc->ijab", e.ovoo, m);
			joovv += T(joovv, { 1, 0, 3, 2 });
			joovv += einsum("kc,ijkabc->ijab", e.fov, as_r6(w));
			Tn dijab({ o, o, v, v });
			n = 0;
			for (int i = 0; i < o; i++)
				for (int j = 0; j < o; j++)
					for (int a = 0; a < v; a++)
						for (int b = 0; b < v; b++) dijab.v[n++] = eia.at(i, a) + eia.at(j, b);
			r.l2_t = divided(joovv, dijab);
		}
		const Tn wv = w + vv * 0.5;
		//the response form (PySCF for_grad=True): full oo/vv blocks, not only their diagonals
		r.goo = einsum("iklabc,jklabc->ij", wv, rw);
		r.gvv = einsum("ijkacd,ijkbcd->ab", wv, rw);
		r.dvo = einsum("ijab,ijkabc->ck", t2, rw) * 0.5;
		//E(T) from the same arrays, 1/3 sum w D r6(w + v/2) (the 1/3 measured against PySCF's ccsd_t, exact to
		//10 digits): checks w and v against occ's and ORCA's (T)
		Tn wd = w;
		for (size_t q = 0; q < wd.v.size(); q++) wd.v[q] *= d3.v[q];
		r.e_t = dot(wd, r6(wv)) / 3.0;
		return r;
	}

	//PySCF ccsd_rdm._gamma1_intermediates: the CCSD blocks before symmetrisation
	struct Gamma1
	{
		Tn doo, dvv, dvo;
	};

	Gamma1 gamma1(const Tn &t1, const Tn &t2, const Tn &l1, const Tn &l2)
	{
		Gamma1 g;
		g.doo = einsum("ja,ia->ij", t1, l1) * -1.0;
		g.dvv = einsum("ia,ib->ab", t1, l1);
		g.dvo = T(t1, { 1, 0 }) - einsum("im,ma->ai", einsum("ie,me->im", t1, l1), t1);
		const Tn theta = t2 * 2.0 - T(t2, { 0, 1, 3, 2 });
		g.doo -= einsum("jkab,ikab->ij", theta, l2);
		g.dvv += einsum("jica,jicb->ab", theta, l2);
		g.dvo += einsum("imae,me->ai", theta, l1);
		g.dvo -= einsum("mi,ma->ai", einsum("mnef,inef->mi", l2, theta), t1);
		g.dvo -= einsum("ie,ae->ai", t1, einsum("mnaf,mnef->ea", l2, theta));
		return g;
	}

	//PySCF ccsd_rdm._make_rdm1: the spin-traced density in the correlated MOs, reference occupation included;
	//symmetric, so no question of <p+q> against <q+p>
	occ::Mat rdm1(const Gamma1 &g, const Tn &l1)
	{
		const int o = l1.d[0], v = l1.d[1];
		occ::Mat dm = occ::Mat::Zero(o + v, o + v);
		for (int i = 0; i < o; i++)
		{
			for (int j = 0; j < o; j++) dm(i, j) = g.doo.at(i, j) + g.doo.at(j, i);
			dm(i, i) += 2.0;
			for (int a = 0; a < v; a++) dm(i, o + a) = dm(o + a, i) = l1.at(i, a) + g.dvo.at(a, i);
		}
		for (int a = 0; a < v; a++)
			for (int b = 0; b < v; b++) dm(o + a, o + b) = g.dvv.at(a, b) + g.dvv.at(b, a);
		return dm;
	}

	occ::Mat flatten(const Tn &l1, const Tn &l2)
	{
		occ::Mat x(l1.v.size() + l2.v.size(), 1);
		std::copy(l1.v.begin(), l1.v.end(), x.data());
		std::copy(l2.v.begin(), l2.v.end(), x.data() + l1.v.size());
		return x;
	}

	void unflatten(const occ::Mat &x, Tn &l1, Tn &l2)
	{
		std::copy(x.data(), x.data() + l1.v.size(), l1.v.begin());
		std::copy(x.data() + l1.v.size(), x.data() + x.size(), l2.v.begin());
	}
}

CcsdDensityResult ccsd_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log)
{
	err_checkf(scf.mo.kind == occ::qm::SpinorbitalKind::Restricted && scf.mo.n_alpha == scf.mo.n_beta,
		"The CCSD density needs a closed-shell (restricted) reference; open-shell Lambda equations are not implemented", log);
	CcsdDensityResult result;
	const int nocc_total = static_cast<int>(scf.mo.n_alpha);
	const int nmo = static_cast<int>(scf.mo.C.cols());
	//as occ's CCSD driver: -1 is the chemical core, and at least one occupied orbital stays correlated
	result.n_frozen = std::max(0, std::min(opts.n_frozen < 0 ? occ::qm::cc::num_frozen_core(scf.basis) : opts.n_frozen, nocc_total - 1));
	const size_t budget = static_cast<size_t>(opts.max_memory_gb * 1024.0 * 1024.0 * 1024.0);

	std::string backend = occ::util::to_lower_copy(opts.backend);
	if (backend.empty() || backend == "auto" || backend == "conventional") backend = "exact";
	if (backend == "ri") backend = "df";
	err_checkf(backend == "exact" || backend == "df", "CCSD density backend '" + opts.backend + "' is not supported (exact | df)", log);
	occ::qm::cc::CCIntegrals cc_eris;
	std::string aux_name;
	if (backend == "exact")
		cc_eris = occ::qm::cc::exact_eris(scf.basis, scf.mo, result.n_frozen, budget);
	else
	{
		aux_name = opts.aux_basis.empty() ? occ::qm::resolve_fitting_basis(scf.basis.name(), occ::qm::FittingKind::Correlation) : opts.aux_basis;
		auto aux = occ::gto::AOBasis::load(scf.basis.atoms(), aux_name);
		aux.set_pure(scf.basis.is_pure());
		cc_eris = occ::qm::cc::df_eris(scf.basis, aux, scf.mo, result.n_frozen, budget);
	}
	err_checkf(result.n_frozen + cc_eris.nocc + cc_eris.nvir == nmo, "CCSD density: the correlated space does not cover the MOs", log);

	//the density is linear in T and Lambda, so T is converged further than an energy would need
	occ::qm::cc::CCSDOptions copts;
	copts.max_cycle = opts.max_cycle;
	copts.tol = opts.tol * 1e-2;
	const occ::qm::cc::CCSDResult cc = occ::qm::cc::ccsd(cc_eris, copts);
	err_checkf(cc.converged, "CCSD did not converge in " + std::to_string(opts.max_cycle) + " iterations", log);
	result.e_corr = cc.e_corr;

	const Eris e = convert(cc_eris);
	const Tn t1 = from_tensor(cc.t1), t2 = from_tensor(cc.t2);
	const Imds m = make_intermediates(t1, t2, e);
	Triples tr;
	if (opts.triples)
	{
		//about eight o^3 v^3 arrays are alive at once in triples()
		const double gb = 8.0 * 8.0 * std::pow(double(e.o), 3) * std::pow(double(e.v), 3) / (1024.0 * 1024.0 * 1024.0);
		err_checkf(gb <= opts.max_memory_gb, "The CCSD(T) density holds full o^3 v^3 arrays, about " + std::to_string(gb) +
			" GB here, over the " + std::to_string(opts.max_memory_gb) + " GB budget (ccsd_max_memory_gb)", log);
		tr = triples(t1, t2, e);
		result.e_triples = tr.e_t;
	}
	Tn l1 = t1, l2 = t2, l1new, l2new;
	occ::core::diis::DIIS diis(2, 8);
	bool converged = false;
	for (int it = 1; it <= opts.max_cycle && !converged; it++)
	{
		update_lambda(t1, t2, l1, l2, e, m, l1new, l2new);
		if (opts.triples)
		{
			l1new += tr.l1_t;
			l2new += tr.l2_t;
		}
		occ::Mat x = flatten(l1new, l2new);
		occ::Mat err = x - flatten(l1, l2);
		const double change = err.norm();
		diis.extrapolate(x, err);
		unflatten(x, l1, l2);
		result.lambda_iterations = it;
		converged = change < opts.tol;
	}
	err_checkf(converged, "The CCSD Lambda equations did not converge in " + std::to_string(opts.max_cycle) + " iterations", log);
	result.t1_norm = t1.norm();
	result.t2_norm = t2.norm();
	result.l1_norm = l1.norm();
	result.l2_norm = l2.norm();

	//frozen orbitals stay doubly occupied and uncoupled
	const int nf = result.n_frozen;
	occ::Mat dm = occ::Mat::Zero(nmo, nmo);
	for (int k = 0; k < nf; k++) dm(k, k) = 2.0;
	Gamma1 g = gamma1(t1, t2, l1, l2);
	if (opts.triples)
	{
		g.doo -= tr.goo * 0.5;
		g.dvv += tr.gvv * 0.5;
		g.dvo += tr.dvo;
	}
	dm.block(nf, nf, nmo - nf, nmo - nf) = rdm1(g, l1);
	Eigen::SelfAdjointEigenSolver<occ::Mat> es(dm);
	err_checkf(es.info() == Eigen::Success, "CCSD density: the MO density could not be diagonalised", log);
	result.occupations = es.eigenvalues().reverse();
	const occ::Mat U = es.eigenvectors().rowwise().reverse();

	occ::qm::Wavefunction &nat = result.natural;
	nat = scf;
	nat.mo.C = scf.mo.C * U;
	nat.mo.occupation = 0.5 * result.occupations;
	nat.mo.energies = occ::Vec::Zero(nmo); //natural orbitals have occupations, not energies
	nat.mo.D = 0.5 * scf.mo.C * dm * scf.mo.C.transpose();
	const std::string name = opts.triples ? "CCSD(T)" : "CCSD";
	nat.method = name;
	nat.energy.total = scf.energy.total + cc.e_corr + result.e_triples;

	log << "\n" << name << " density (unrelaxed, Lambda equations), " << (backend == "exact" ? std::string("exact integrals") : "DF integrals, " + aux_name) << "\n"
		<< "  correlated: " << cc_eris.nocc << " occupied, " << cc_eris.nvir << " virtual; frozen core: " << nf << " orbitals\n"
		<< std::fixed << std::setprecision(10)
		<< "  E(CCSD correlation) = " << cc.e_corr << "  (" << cc.iterations << " iterations)\n";
	if (opts.triples) log << "  E(T) = " << result.e_triples << "\n";
	log << "  Lambda equations converged in " << result.lambda_iterations << " iterations\n"
		<< std::setprecision(6) << "  |t1| = " << result.t1_norm << "  |t2| = " << result.t2_norm
		<< "  |l1| = " << result.l1_norm << "  |l2| = " << result.l2_norm << "\n"
		<< "  natural occupations:";
	for (int k = 0; k < std::min(nmo, 8); k++) log << " " << result.occupations(k);
	log << (nmo > 8 ? " ..." : "") << "\n  sum of occupations: " << std::setprecision(8) << result.occupations.sum()
		<< ", smallest " << std::scientific << std::setprecision(3) << result.occupations(nmo - 1) << std::defaultfloat << std::endl;
	return result;
}

CcsdDensityResult ccsd_density_from_occ_input(const std::string &toml_path, occ::io::OccInput config, std::ostream &log)
{
	const auto spec = occ::driver::parse_method_string(config.method.name);
	err_checkf(spec.kind == occ::driver::MethodKind::CCSD || spec.kind == occ::driver::MethodKind::CCSD_T,
		"Method '" + config.method.name + "' has no density here; CCSD and CCSD(T) do (method = \"ccsd\" or \"ccsd(t)\")", log);
	const auto input = toml::parse_file(toml_path);  //the view below points into it, so it must outlive the reads
	const auto scf_table = input["scf"];
	CcsdDensityOptions opts;
	opts.triples = spec.kind == occ::driver::MethodKind::CCSD_T;
	if (!spec.backend.empty()) opts.backend = spec.backend;  //"df-ccsd", as occ reads it
	opts.n_frozen = scf_table["ccsd_frozen_core"].value_or(opts.n_frozen);
	opts.backend = scf_table["ccsd_backend"].value_or(opts.backend);
	opts.max_memory_gb = scf_table["ccsd_max_memory_gb"].value_or(opts.max_memory_gb);
	opts.aux_basis = config.basis.ri_basis;
	config.method.name = "rhf";
	//occ's automatic DF/COSX would shift the reference by ~3e-5 Eh in E(HF); exact CCSD wants exact J/K, as ORCA's NoRI
	const std::string backend = occ::util::to_lower_copy(opts.backend);
	if (backend != "df" && backend != "ri") config.method.ri_policy = occ::io::RIPolicy::None;
	return ccsd_natural_orbitals(occ::main::run_scf_external(config, false), opts, log);
}

occ::qm::Wavefunction wavefunction_from_occ_input(const std::string &toml_path, const occ::io::OccInput &config, std::ostream &log)
{
	const auto kind = occ::driver::parse_method_string(config.method.name).kind;
	if (kind != occ::driver::MethodKind::CCSD && kind != occ::driver::MethodKind::CCSD_T) return occ::main::run_scf_external(config, true);
	return ccsd_density_from_occ_input(toml_path, config, log).natural;
}
