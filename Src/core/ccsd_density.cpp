#include "pch.h"
#include "ccsd_density.h"
#include "convenience.h"

#include <occ/core/diis.h>
#include <occ/core/log.h>
#include <occ/core/util.h>
#include <occ/driver/method_parser.h>
#include <occ/main/occ_scf.h>
#include <occ/qm/cc/ccsd.h>
#include <occ/qm/cc/integrals.h>
#include <occ/qm/fitting_basis.h>
#include <occ/qm/integral_engine.h>
#include <occ/qm/integral_engine_df.h>
#include <occ/qm/correlation/df_integrals.h>
#include <occ/core/units.h>
#include <Eigen/Dense>
#include <toml++/toml.hpp>
#include <unsupported/Eigen/CXX11/Tensor>

#include <array>
#include <chrono>
#include <functional>
#include <iomanip>
#include <ostream>

namespace
{
	//dense row-major tensor (rank <= 6) with numpy's einsum and transpose: the equations below are PySCF's
	//ccsd_lambda / ccsd_rdm and their (T) versions term for term (chemist notation, t2 as ijab, t3 as ijkabc)
	struct Tn
	{
		ivec d;
		vec v;
		Tn() = default;
		explicit Tn(const ivec &dims) : d(dims), v(count(dims), 0.0) {}
		static size_t count(const ivec &dims)
		{
			size_t n = 1;
			for (int k = 0; k < static_cast<int>(dims.size()); k++) n *= static_cast<size_t>(dims[k]);
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

	//numpy a.transpose(perm): out.d[k] = a.d[perm[k]]. Padded at the front, so the innermost loop runs over the
	//last real output axis, a plain copy when that axis stays last
	Tn T(const Tn &a, const ivec &perm)
	{
		constexpr int R = 6;
		const int r = static_cast<int>(a.d.size()), pad = R - r;
		std::array<int, R> ad{};
		std::array<size_t, R> stride{};
		for (int k = 0; k < R; k++) ad[k] = k < pad ? 1 : a.d[k - pad];
		size_t s = 1;
		for (int k = R - 1; k >= 0; k--) { stride[k] = s; s *= static_cast<size_t>(ad[k]); }
		ivec od(r);
		for (int k = 0; k < r; k++) od[k] = a.d[perm[k]];
		Tn out(od);
		std::array<int, R> n{};
		std::array<size_t, R> st{};
		for (int k = 0; k < R; k++)
		{
			const int src = k < pad ? k : pad + perm[k - pad];
			n[k] = ad[src];
			st[k] = stride[src];
		}
		const long long outer = static_cast<long long>(n[0]) * n[1] * n[2] * n[3] * n[4];
#pragma omp parallel for if (out.v.size() > (size_t(1) << 18))
		for (long long q = 0; q < outer; q++)
		{
			size_t rem = static_cast<size_t>(q), base = 0;
			for (int k = 4; k >= 0; k--)
			{
				base += (rem % n[k]) * st[k];
				rem /= n[k];
			}
			double *dst = out.v.data() + static_cast<size_t>(q) * n[5];
			const double *src = a.v.data() + base;
			if (st[5] == 1)
				std::copy_n(src, n[5], dst);
			else
				for (int x5 = 0; x5 < n[5]; x5++) dst[x5] = src[x5 * st[5]];
		}
		return out;
	}

	bool identity(const ivec &perm)
	{
		for (size_t k = 0; k < perm.size(); k++)
			if (perm[k] != static_cast<int>(k)) return false;
		return true;
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
		ivec pa, pb, cd;
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
		Tn at, bt;
		if (!identity(pa)) at = T(a, pa);
		if (!identity(pb)) bt = T(b, pb);
		const double *ap = at.v.empty() ? a.v.data() : at.v.data(), *bp = bt.v.empty() ? b.v.data() : bt.v.data();
		Tn c(cd);
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, static_cast<int>(M), static_cast<int>(N), static_cast<int>(K),
			1.0, ap, static_cast<int>(K), bp, static_cast<int>(N), 0.0, c.v.data(), static_cast<int>(N));
		const std::string lc = fa + fb;
		ivec po;
		for (const char x : lo) po.push_back(static_cast<int>(lc.find(x)));
		return identity(po) ? c : T(c, po);
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
	//density those of CCSD(T), evaluated one occupied triple at a time. W[ijk](abc) is PySCF's p6-symmetrised
	//t3-like array w (likewise V for v), already divided by D, as a v^3 block; nothing o^3 v^3 is ever formed.
	struct Triples
	{
		Tn l1_t, l2_t;  //constant terms added to every Lambda update
		Tn goo, gvv, dvo;
		Tn gvv_c;  //the vv block's share of the occupied directions in proj (local clusters)
		double e_t = 0.0;
	};

	using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	using MapR = Eigen::Map<RowMat>;
	using CMapR = Eigen::Map<const RowMat>;

	Triples triples(const Tn &t1, const Tn &t2, const Eris &e, const occ::Mat *proj = nullptr)
	{
		const int o = e.o, v = e.v;
		const size_t vv = size_t(v) * v, v3 = vv * v;
		const Tn L1 = T(e.ovvv, { 0, 1, 3, 2 });  //[i][a][b][f] = (ia|fb)
		const Tn L2 = T(e.ovvv, { 0, 3, 1, 2 });  //[k][e][f][b] = (kf|be)
		const Tn t2T = T(t2, { 0, 1, 3, 2 });     //[k][j][f][c] = t2(kj,cf)
		const Tn t2_t = T(t2, { 1, 0, 2, 3 });    //[k][m][b][c] = t2(mk,bc)
		const Tn ovoo_t = T(e.ovoo, { 0, 2, 1, 3 }), ovov_t = T(e.ovov, { 0, 2, 1, 3 });  //[i][j][a][m], [i][j][a][b]
		auto at = [](const Tn &a, const size_t block, const size_t size) { return a.v.data() + block * size; };
		auto ix = [v](const int a, const int b, const int c) { return (size_t(a) * v + b) * v + c; };

		//PySCF's un-symmetrised w for one ordered triple: (iafb,kjcf) - (iajm,mkbc), into out as abc
		const int iv = v, ivv = static_cast<int>(vv);
		auto gemm = [](const bool ta, const bool tb, const int m, const int n, const int k, const double alpha, const double *a, const int lda,
			const double *b, const int ldb, const double beta, double *c, const int ldc) {
			cblas_dgemm(CblasRowMajor, ta ? CblasTrans : CblasNoTrans, tb ? CblasTrans : CblasNoTrans, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
		};
		auto base_w = [&](const int i, const int j, const int k, double *out) {
			gemm(false, false, ivv, iv, iv, 1.0, at(L1, i, v3), iv, at(t2T, size_t(k) * o + j, vv), iv, 0.0, out, iv);
			gemm(false, false, iv, ivv, o, -1.0, at(ovoo_t, size_t(i) * o + j, size_t(v) * o), o, at(t2_t, k, size_t(o) * vv), ivv, 1.0, out, ivv);
		};

		Triples r;
		r.l1_t = Tn({ o, v });
		r.goo = Tn({ o, o });
		r.gvv = Tn({ v, v });
		r.dvo = Tn({ v, o });
		if (proj) r.gvv_c = Tn({ v, v });
		Tn joovv({ o, o, v, v });
		std::vector<Tn> gvv_i(o, Tn({ v, v })), dvo_i(o, Tn({ v, o }));
		vec e_i(o, 0.0);
		//W + V/2 and r6(W) of every i of the current (j, k): the oo block pairs W[ijk] with W[i'jk]
		RowMat WV(o, v3), RW(o, v3);
		const int nt = omp_get_max_threads();
		vec2 buf(size_t(nt) * 8, vec(v3));

		//everything one block (i, j, k) contributes, given its W and V; wv, rw receive W + V/2 and r6(W)
		auto block = [&](const int i, const int j, const int k, const double *W, const double *V, double *wv, double *rw, double *m, double *arw) {
			const double eijk = e.eo(i) + e.eo(j) + e.eo(k);
			double et = 0.0;
			for (int a = 0; a < v; a++)
				for (int b = 0; b < v; b++)
					for (int c = 0; c < v; c++)
					{
						const size_t x = ix(a, b, c), cab = ix(c, a, b), bca = ix(b, c, a), cba = ix(c, b, a), acb = ix(a, c, b), bac = ix(b, a, c);
						//r6 and as_r6 act on abc only: r6(w) = 4w + w(cab) + w(bca) - 2 [w(cba) + w(acb) + w(bac)]
						rw[x] = 4.0 * W[x] + W[cab] + W[bca] - 2.0 * (W[cba] + W[acb] + W[bac]);
						wv[x] = W[x] + 0.5 * V[x];
						m[x] = 2.0 * (2.0 * W[x] + 0.5 * V[x]) - (2.0 * W[cba] + 0.5 * V[cba]) - (2.0 * W[acb] + 0.5 * V[acb]);
						arw[x] = 2.0 * W[x] - W[cba] - W[acb];
						//E(T) = 1/3 sum W D r6(W + V/2) = 1/3 sum D r6(W) (W + V/2), r6 self-adjoint and D symmetric in abc
						et += (eijk - e.ev(a) - e.ev(b) - e.ev(c)) * rw[x] * wv[x];
					}
			e_i[i] += et / 3.0;
			//l1_t[i] and joovv[i] are this i's alone; the vv, vo parts are per i and summed after the loop
			MapR(r.l1_t.v.data() + size_t(i) * v, 1, v).noalias() += (CMapR(rw, v, vv) * CMapR(at(ovov_t, size_t(j) * o + k, vv), vv, 1)).transpose();
			gemm(false, false, iv, iv, ivv, 1.0, m, ivv, at(L2, k, v3), iv, 1.0, joovv.v.data() + (size_t(i) * o + j) * vv, iv);
			MapR(joovv.v.data() + (size_t(i) * o + j) * vv, vv, 1).noalias() += CMapR(arw, vv, v) * CMapR(e.fov.v.data() + size_t(k) * v, v, 1);
			//joovv[i, J, ab] -= sum_c (kc|jJ) m[ab, c]
			gemm(true, true, o, ivv, iv, -1.0, at(ovoo_t, size_t(k) * o + j, size_t(v) * o), o, m, iv, 1.0, joovv.v.data() + size_t(i) * o * vv, ivv);
			gemm(false, true, iv, iv, ivv, 1.0, wv, ivv, rw, ivv, 1.0, gvv_i[i].v.data(), iv);
			const Eigen::VectorXd dv = CMapR(rw, vv, v).transpose() * CMapR(at(t2, size_t(i) * o + j, vv), vv, 1);
			for (int c = 0; c < v; c++) dvo_i[i].at(c, k) += 0.5 * dv(c);
		};

		//block (i, k, j) is block (i, j, k) with b and c swapped (joint permutation symmetry): same oo share, hence
		//the weight 2; as_r6 does not commute with the swap, so its other terms come from the swapped W and V
		for (int j = 0; j < o; j++)
			for (int k = 0; k <= j; k++)
			{
#pragma omp parallel for schedule(dynamic)
				for (int i = 0; i < o; i++)
				{
					vec *B = &buf[size_t(omp_get_thread_num()) * 8];
					double *b0 = B[0].data(), *b1 = B[1].data(), *b2 = B[2].data(), *b3 = B[3].data(), *b4 = B[4].data(), *b5 = B[5].data();
					double *W = B[6].data(), *V = B[7].data();
					base_w(i, j, k, b0);
					base_w(i, k, j, b1);
					base_w(j, i, k, b2);
					base_w(j, k, i, b3);
					base_w(k, i, j, b4);
					base_w(k, j, i, b5);
					//p6 and the division by D in one pass; v's six permutations coincide in pairs, so three terms twice
					const double eijk = e.eo(i) + e.eo(j) + e.eo(k);
					const double *Xij = at(ovov_t, size_t(i) * o + j, vv), *Xik = at(ovov_t, size_t(i) * o + k, vv), *Xjk = at(ovov_t, size_t(j) * o + k, vv);
					const double *Tij = at(t2, size_t(i) * o + j, vv), *Tik = at(t2, size_t(i) * o + k, vv), *Tjk = at(t2, size_t(j) * o + k, vv);
					const double *t1i = at(t1, i, v), *t1j = at(t1, j, v), *t1k = at(t1, k, v);
					const double *fi = at(e.fov, i, v), *fj = at(e.fov, j, v), *fk = at(e.fov, k, v);
					for (int a = 0; a < v; a++)
						for (int b = 0; b < v; b++)
							for (int c = 0; c < v; c++)
							{
								const size_t x = ix(a, b, c);
								const double d = eijk - e.ev(a) - e.ev(b) - e.ev(c);
								W[x] = (b0[x] + b1[ix(a, c, b)] + b2[ix(b, a, c)] + b3[ix(b, c, a)] + b4[ix(c, a, b)] + b5[ix(c, b, a)]) / d;
								const size_t ab = size_t(a) * v + b, ac = size_t(a) * v + c, bc = size_t(b) * v + c;
								V[x] = 2.0 * (Xij[ab] * t1k[c] + fk[c] * Tij[ab] + Xik[ac] * t1j[b] + fj[b] * Tik[ac] + Xjk[bc] * t1i[a] + fi[a] * Tjk[bc]) / d;
							}
					block(i, j, k, W, V, WV.row(i).data(), RW.row(i).data(), b0, b1);
					if (j == k) continue;
					double *Wm = b2, *Vm = b3;
					for (int a = 0; a < v; a++)
						for (int b = 0; b < v; b++)
							for (int c = 0; c < v; c++)
							{
								Wm[ix(a, b, c)] = W[ix(a, c, b)];
								Vm[ix(a, b, c)] = V[ix(a, c, b)];
							}
					block(i, k, j, Wm, Vm, b4, b5, b0, b1);
				}
				//the vv block attributed by the block's first occupied index, projected onto proj's directions
				if (proj)
					for (int c = 0; c < proj->cols(); c++)
					{
						const Eigen::VectorXd wp = WV.transpose() * proj->col(c), rp = RW.transpose() * proj->col(c);
						gemm(false, true, iv, iv, ivv, j == k ? 1.0 : 2.0, wp.data(), ivv, rp.data(), ivv, 1.0, r.gvv_c.v.data(), iv);
					}
				//the response form (PySCF for_grad=True): the full oo block, not only its diagonal
				cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, o, o, static_cast<int>(v3), j == k ? 1.0 : 2.0,
					WV.data(), static_cast<int>(v3), RW.data(), static_cast<int>(v3), 1.0, r.goo.v.data(), o);
			}
		for (int i = 0; i < o; i++)
		{
			r.gvv += gvv_i[i];
			r.dvo += dvo_i[i];
			r.e_t += e_i[i];
			for (int a = 0; a < v; a++) r.l1_t.at(i, a) *= 0.5 / (e.eo(i) - e.ev(a));
		}
		joovv += T(joovv, { 1, 0, 3, 2 });
		r.l2_t = joovv;
		size_t n = 0;
		for (int i = 0; i < o; i++)
			for (int j = 0; j < o; j++)
				for (int a = 0; a < v; a++)
					for (int b = 0; b < v; b++) r.l2_t.v[n++] /= e.eo(i) - e.ev(a) + e.eo(j) - e.ev(b);
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

	//CCSD, the triples terms and the Lambda equations on prepared integrals. proj (o x n, optional) holds the
	//occupied directions whose share of the triples vv block is formed too (a local cluster's central orbitals)
	struct CcSolution
	{
		Tn t1, t2, l1, l2;
		Triples tr;
		double e_corr = 0.0;
		int cc_iterations = 0, lambda_iterations = 0;
	};

	CcSolution solve_cc(const occ::qm::cc::CCIntegrals &cc_eris, const CcsdDensityOptions &opts, std::ostream &log,
		const std::function<void(const char *)> &mark, const occ::Mat *proj = nullptr)
	{
		CcSolution s;
		occ::qm::cc::CCSDOptions copts;
		copts.max_cycle = opts.max_cycle;
		copts.tol = opts.tol * 1e-2;
		const occ::qm::cc::CCSDResult cc = occ::qm::cc::ccsd(cc_eris, copts);
		err_checkf(cc.converged, "CCSD did not converge in " + std::to_string(opts.max_cycle) + " iterations", log);
		s.e_corr = cc.e_corr;
		s.cc_iterations = cc.iterations;
		mark("CCSD amplitudes");
		const Eris e = convert(cc_eris);
		s.t1 = from_tensor(cc.t1);
		s.t2 = from_tensor(cc.t2);
		const Tn &t1 = s.t1, &t2 = s.t2;
		const Imds m = make_intermediates(t1, t2, e);
		mark("Lambda intermediates");
		if (opts.triples)
		{
			//W + V/2 and r6(W) for every i of one (j, k), two re-laid-out ovvv, eight v^3 buffers per thread
			const double v3 = std::pow(double(e.v), 3);
			const double gb = 8.0 * (4.0 * e.o * v3 + 8.0 * v3 * omp_get_max_threads()) / (1024.0 * 1024.0 * 1024.0);
			err_checkf(gb <= opts.max_memory_gb, "The CCSD(T) density needs about " + std::to_string(gb) +
				" GB here, over the " + std::to_string(opts.max_memory_gb) + " GB budget (ccsd_max_memory_gb)", log);
			s.tr = triples(t1, t2, e, proj);
			mark("triples terms");
		}
		s.l1 = t1;
		s.l2 = t2;
		Tn &l1 = s.l1, &l2 = s.l2;
		Tn l1new, l2new;
		occ::core::diis::DIIS diis(2, 8);
		bool converged = false;
		for (int it = 1; it <= opts.max_cycle && !converged; it++)
		{
			update_lambda(t1, t2, l1, l2, e, m, l1new, l2new);
			if (opts.triples)
			{
				l1new += s.tr.l1_t;
				l2new += s.tr.l2_t;
			}
			occ::Mat x = flatten(l1new, l2new);
			occ::Mat err = x - flatten(l1, l2);
			const double change = err.norm();
			diis.extrapolate(x, err);
			unflatten(x, l1, l2);
			s.lambda_iterations = it;
			converged = change < opts.tol;
		}
		mark("Lambda iterations");
		err_checkf(converged, "The CCSD Lambda equations did not converge in " + std::to_string(opts.max_cycle) + " iterations", log);
		return s;
	}

	//occ's integral blocks for mo (exact or DF, as opts.backend), the n_frozen lowest occupied left out
	occ::qm::cc::CCIntegrals build_eris(const occ::gto::AOBasis &basis, const occ::qm::MolecularOrbitals &mo, const int n_frozen,
		const CcsdDensityOptions &opts, std::string &backend, std::string &aux_name, std::ostream &log)
	{
		backend = occ::util::to_lower_copy(opts.backend);
		if (backend.empty() || backend == "auto" || backend == "conventional") backend = "exact";
		if (backend == "ri") backend = "df";
		err_checkf(backend == "exact" || backend == "df", "CCSD density backend '" + opts.backend + "' is not supported (exact | df)", log);
		const size_t budget = static_cast<size_t>(opts.max_memory_gb * 1024.0 * 1024.0 * 1024.0);
		if (backend == "exact") return occ::qm::cc::exact_eris(basis, mo, n_frozen, budget);
		aux_name = opts.aux_basis.empty() ? occ::qm::resolve_fitting_basis(basis.name(), occ::qm::FittingKind::Correlation) : opts.aux_basis;
		auto aux = occ::gto::AOBasis::load(basis.atoms(), aux_name);
		aux.set_pure(basis.is_pure());
		return occ::qm::cc::df_eris(basis, aux, mo, n_frozen, budget);
	}

	//natural orbitals of dm (in scf's MO basis) as an occ Wavefunction: C_NO in mo.C, n/2 in mo.occupation
	void natural_orbitals(const occ::qm::Wavefunction &scf, const occ::Mat &dm, const std::string &name, const double e_total,
		CcsdDensityResult &result, std::ostream &log)
	{
		const int nmo = static_cast<int>(dm.rows());
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
		nat.method = name;
		nat.energy.total = e_total;
	}

	//Pipek-Mezey localisation with Mulliken charges by Jacobi sweeps (Pipek & Mezey, JCP 90, 4916 (1989)): rotates
	//C's columns among themselves to maximise sum_k sum_A (Q^A_kk)^2
	occ::Mat pipek_mezey(occ::Mat C, const occ::Mat &S, const ivec &bf_atom, const int natoms, std::ostream &log)
	{
		const int nbf = static_cast<int>(C.rows()), n = static_cast<int>(C.cols());
		occ::Mat SC = S * C;
		vec qss(natoms), qtt(natoms), qst(natoms);
		double max_rot = 1.0;
		int sweep = 0;
		for (; sweep < 500 && max_rot > 1e-8; sweep++)
		{
			max_rot = 0.0;
			for (int a = 0; a < n; a++)
				for (int b = a + 1; b < n; b++)
				{
					std::fill(qss.begin(), qss.end(), 0.0);
					std::fill(qtt.begin(), qtt.end(), 0.0);
					std::fill(qst.begin(), qst.end(), 0.0);
					for (int mu = 0; mu < nbf; mu++)
					{
						const int A = bf_atom[mu];
						qss[A] += C(mu, a) * SC(mu, a);
						qtt[A] += C(mu, b) * SC(mu, b);
						qst[A] += 0.5 * (C(mu, a) * SC(mu, b) + C(mu, b) * SC(mu, a));
					}
					double x = 0.0, y = 0.0;
					for (int A = 0; A < natoms; A++)
					{
						const double d = qss[A] - qtt[A];
						x += qst[A] * qst[A] - 0.25 * d * d;
						y += qst[A] * d;
					}
					if (x * x + y * y < 1e-30) continue;
					//cos 4g = -x / r, sin 4g = y / r maximises the functional for this pair
					const double g = 0.25 * std::atan2(y, -x);
					if (std::abs(g) < 1e-14) continue;
					max_rot = std::max(max_rot, std::abs(g));
					const double c = std::cos(g), sn = std::sin(g);
					const occ::Vec ca = C.col(a), cb = C.col(b), sa = SC.col(a), sb = SC.col(b);
					C.col(a) = c * ca + sn * cb;
					C.col(b) = c * cb - sn * ca;
					SC.col(a) = c * sa + sn * sb;
					SC.col(b) = c * sb - sn * sa;
				}
		}
		if (max_rot > 1e-8) log << "  Pipek-Mezey: largest rotation still " << max_rot << " after " << sweep << " sweeps" << std::endl;
		return C;
	}

	Tn as_tn(const occ::Mat &R)
	{
		Tn r({ static_cast<int>(R.rows()), static_cast<int>(R.cols()) });
		for (int k = 0; k < R.rows(); k++)
			for (int i = 0; i < R.cols(); i++) r.at(k, i) = R(k, i);
		return r;
	}

	//one cluster's half-density H (H + H^T: the correlation part of the 1-RDM) in [cluster LMOs | virtuals],
	//R(k, i) = <LMO k | semicanonical i>; a term belongs to the cluster of its distinguished occupied index (Li &
	//Piecuch): the free one for oo and ov, Lambda's first for vv, the triples block's first through proj
	occ::Mat cluster_half_density(const CcSolution &s, const occ::Mat &R, const ivec &central, const bool with_triples)
	{
		const int o = s.t1.d[0], v = s.t1.d[1];
		const Tn Rt = as_tn(R);
		auto occ1 = [&](const Tn &x) { return einsum("ki,ia->ka", Rt, x); };
		auto occ2 = [&](const Tn &x) { return einsum("lj,kjab->klab", Rt, einsum("ki,ijab->kjab", Rt, x)); };
		const Tn t1 = occ1(s.t1), l1 = occ1(s.l1), t2 = occ2(s.t2), l2 = occ2(s.l2);
		bvec is_c(o, false);
		for (int c = 0; c < static_cast<int>(central.size()); c++) is_c[central[c]] = true;
		const Gamma1 g = gamma1(t1, t2, l1, l2);
		Tn l1c = l1, l2c = l2;
		const size_t ovv = size_t(o) * v * v;
		for (int k = 0; k < o; k++)
			if (!is_c[k])
			{
				std::fill_n(l1c.v.begin() + size_t(k) * v, v, 0.0);
				std::fill_n(l2c.v.begin() + size_t(k) * ovv, ovv, 0.0);
			}
		const Tn theta = t2 * 2.0 - T(t2, { 0, 1, 3, 2 });
		Tn dvv = einsum("ia,ib->ab", t1, l1c) + einsum("jica,jicb->ab", theta, l2c);
		occ::Mat doo(o, o), dvo(v, o);
		for (int k = 0; k < o; k++)
		{
			for (int j = 0; j < o; j++) doo(k, j) = g.doo.at(k, j);
			for (int a = 0; a < v; a++) dvo(a, k) = g.dvo.at(a, k);
		}
		if (with_triples)
		{
			occ::Mat goo(o, o), dvo_t(v, o);
			for (int i = 0; i < o; i++)
			{
				for (int j = 0; j < o; j++) goo(i, j) = s.tr.goo.at(i, j);
				for (int a = 0; a < v; a++) dvo_t(a, i) = s.tr.dvo.at(a, i);
			}
			doo -= 0.5 * R * goo * R.transpose();
			dvo += dvo_t * R.transpose();
			dvv += s.tr.gvv_c * 0.5;
		}
		occ::Mat H = occ::Mat::Zero(o + v, o + v);
		for (int c = 0; c < static_cast<int>(central.size()); c++)
		{
			const int k = central[c];
			H.block(k, 0, 1, o) = doo.row(k);
			for (int a = 0; a < v; a++) H(k, o + a) = l1.at(k, a) + dvo(a, k);
		}
		for (int a = 0; a < v; a++)
			for (int b = 0; b < v; b++) H(o + a, o + b) = dvv.at(a, b);
		return H;
	}

	//the cluster's share of E(CCSD corr), sum_jab [2 (Ia|jb) - (Ib|ja)] tau_Ij^ab with LMO I central
	double cluster_energy(const CcSolution &s, const occ::Mat &R, const Tn &ovov, const int I)
	{
		const int o = s.t1.d[0], v = s.t1.d[1];
		const Tn Rt = as_tn(R);
		const Tn t1 = einsum("ki,ia->ka", Rt, s.t1);
		const Tn t2 = einsum("lj,kjab->klab", Rt, einsum("ki,ijab->kjab", Rt, s.t2));
		const Tn g = einsum("lj,kajb->kalb", Rt, einsum("ki,iajb->kajb", Rt, ovov));
		double e = 0.0;
		for (int j = 0; j < o; j++)
			for (int a = 0; a < v; a++)
				for (int b = 0; b < v; b++)
				{
					const double tau = t2.v[((size_t(I) * o + j) * v + a) * v + b] + t1.at(I, a) * t1.at(j, b);
					e += (2.0 * g.v[((size_t(I) * v + a) * o + j) * v + b] - g.v[((size_t(I) * v + b) * o + j) * v + a]) * tau;
				}
		return e;
	}
}

CcsdDensityResult ccsd_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log)
{
	err_checkf(scf.mo.kind == occ::qm::SpinorbitalKind::Restricted && scf.mo.n_alpha == scf.mo.n_beta,
		"The CCSD density needs a closed-shell (restricted) reference; open-shell Lambda equations are not implemented", log);
	CcsdDensityResult result;
	auto clock = std::chrono::steady_clock::now();
	auto mark = [&clock](const char *stage) {
		const auto now = std::chrono::steady_clock::now();
		occ::log::info("CCSD density: {} took {:.3f} s", stage, std::chrono::duration<double>(now - clock).count());
		clock = now;
	};
	const int nocc_total = static_cast<int>(scf.mo.n_alpha);
	const int nmo = static_cast<int>(scf.mo.C.cols());
	//as occ's CCSD driver: -1 is the chemical core, and at least one occupied orbital stays correlated
	result.n_frozen = std::max(0, std::min(opts.n_frozen < 0 ? occ::qm::cc::num_frozen_core(scf.basis) : opts.n_frozen, nocc_total - 1));
	std::string backend, aux_name;
	const occ::qm::cc::CCIntegrals cc_eris = build_eris(scf.basis, scf.mo, result.n_frozen, opts, backend, aux_name, log);
	mark("MO integrals");
	err_checkf(result.n_frozen + cc_eris.nocc + cc_eris.nvir == nmo, "CCSD density: the correlated space does not cover the MOs", log);
	const CcSolution s = solve_cc(cc_eris, opts, log, mark);
	result.e_corr = s.e_corr;
	result.e_triples = s.tr.e_t;
	result.lambda_iterations = s.lambda_iterations;
	const Tn &t1 = s.t1, &t2 = s.t2, &l1 = s.l1, &l2 = s.l2;
	const Triples &tr = s.tr;
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
	const std::string name = opts.triples ? "CCSD(T)" : "CCSD";
	natural_orbitals(scf, dm, name, scf.energy.total + s.e_corr + result.e_triples, result, log);
	mark("density and natural orbitals");
	log << "\n" << name << " density (unrelaxed, Lambda equations), " << (backend == "exact" ? std::string("exact integrals") : "DF integrals, " + aux_name) << "\n"
		<< "  correlated: " << cc_eris.nocc << " occupied, " << cc_eris.nvir << " virtual; frozen core: " << nf << " orbitals\n"
		<< std::fixed << std::setprecision(10)
		<< "  E(CCSD correlation) = " << s.e_corr << "  (" << s.cc_iterations << " iterations)\n";
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

CcsdDensityResult local_ccsd_natural_orbitals(const occ::qm::Wavefunction &scf, const CcsdDensityOptions &opts, std::ostream &log)
{
	err_checkf(scf.mo.kind == occ::qm::SpinorbitalKind::Restricted && scf.mo.n_alpha == scf.mo.n_beta,
		"The local CCSD density needs a closed-shell (restricted) reference", log);
	CcsdDensityResult result;
	auto clock = std::chrono::steady_clock::now();
	auto mark = [&clock](const char *stage) {
		const auto now = std::chrono::steady_clock::now();
		occ::log::info("local CCSD density: {} took {:.3f} s", stage, std::chrono::duration<double>(now - clock).count());
		clock = now;
	};
	const auto quiet = [](const char *) {};
	const auto &basis = scf.basis;
	const occ::Mat &C = scf.mo.C;
	const int nbf = static_cast<int>(C.rows()), nmo = static_cast<int>(C.cols()), nocc = static_cast<int>(scf.mo.n_alpha);
	result.n_frozen = std::max(0, std::min(opts.n_frozen < 0 ? occ::qm::cc::num_frozen_core(basis) : opts.n_frozen, nocc - 1));
	const int nf = result.n_frozen;
	const occ::qm::IntegralEngine engine(basis);
	const occ::Mat S = engine.one_electron_operator(occ::qm::cint::Operator::overlap);
	//the converged Fock matrix in the AO basis, exact on the MO span: F = S C e C^T S
	const occ::Mat F = S * C * scf.mo.energies.head(nmo).asDiagonal() * C.transpose() * S;
	const occ::Mat Cocc = C.leftCols(nocc);
	//projected atomic orbitals: column mu is AO mu with every occupied orbital (core included) projected out
	const occ::Mat Q = occ::Mat::Identity(nbf, nbf) - Cocc * (Cocc.transpose() * S);
	const auto &bf_atom = basis.bf_to_atom();
	const int natoms = static_cast<int>(basis.atoms().size());
	const occ::Mat L = pipek_mezey(C.middleCols(nf, nocc - nf), S, bf_atom, natoms, log);
	const int n = static_cast<int>(L.cols());
	const occ::Mat SL = S * L;
	occ::Mat pop = occ::Mat::Zero(natoms, n);  //Mulliken population of LMO k on atom A
	for (int mu = 0; mu < nbf; mu++)
		for (int k = 0; k < n; k++) pop(bf_atom[mu], k) += L(mu, k) * SL(mu, k);
	occ::Mat centre = occ::Mat::Zero(3, n);
	for (int k = 0; k < n; k++)
	{
		for (int A = 0; A < natoms; A++)
			centre.col(k) += pop(A, k) * occ::Vec3(basis.atoms()[A].x, basis.atoms()[A].y, basis.atoms()[A].z);
		centre.col(k) /= pop.col(k).sum();
	}
	mark("localisation");
	const double r_cluster = opts.cluster_radius * occ::units::ANGSTROM_TO_BOHR;
	occ::Mat Dhalf = occ::Mat::Zero(nbf, nbf);
	double e_corr = 0.0, e_t = 0.0;
	long long sum_o = 0, sum_v = 0, sum_cand = 0, sum_dom = 0;
	int max_o = 0, max_v = 0;
	std::string backend = occ::util::to_lower_copy(opts.backend), aux_name;
	const bool use_df = backend == "df" || backend == "ri";
	//one DF store for the molecule: every cluster's MP2 screening and, with the df backend, its CC integrals
	const std::string fit = opts.aux_basis.empty() ? occ::qm::resolve_fitting_basis(basis.name(), occ::qm::FittingKind::Correlation) : opts.aux_basis;
	occ::gto::AOBasis aux = occ::gto::AOBasis::load(basis.atoms(), fit);
	aux.set_kind(basis.kind());
	occ::qm::IntegralEngineDF df_engine(basis.atoms(), basis.shells(), aux.shells());
	occ::qm::DFIntegrals df(df_engine, static_cast<size_t>(opts.max_memory_gb * 1024.0 * 1024.0 * 1024.0));
	const occ::Mat FL = F * L;
	occ::Vec fdiag(n);
	for (int k = 0; k < n; k++) fdiag(k) = L.col(k).dot(FL.col(k));
	mark("DF store");
	//per LMO I: MP2 pairs (I, J) over the candidates within r_cluster, in their PAO domain; pairs above
	//pair_threshold select I's LMOs, and their MP2 virtual density is kept with the domain's virtuals
	ivec2 sel(n), doms(n);
	std::vector<occ::Mat> Cdom(n), Dvir(n);
	for (int I = 0; I < n; I++)
	{
		ivec cand{ I };
		for (int J = 0; J < n; J++)
			if (J != I && (centre.col(I) - centre.col(J)).norm() < r_cluster) cand.push_back(J);
		const int nc = static_cast<int>(cand.size());
		ivec dom;
		for (int mu = 0; mu < nbf; mu++)
			for (int c = 0; c < nc; c++)
				if (pop(bf_atom[mu], cand[c]) > opts.domain_threshold)
				{
					dom.push_back(mu);
					break;
				}
		occ::Mat P(nbf, dom.size());
		for (size_t k = 0; k < dom.size(); k++) P.col(k) = Q.col(dom[k]);
		//orthonormal virtuals from the domain's PAOs, the occupied directions (overlap eigenvalue 0) dropped
		Eigen::SelfAdjointEigenSolver<occ::Mat> sp(P.transpose() * S * P);
		ivec keep;
		for (int k = 0; k < sp.eigenvalues().size(); k++)
			if (sp.eigenvalues()(k) > 1e-6) keep.push_back(k);
		const int nd = static_cast<int>(keep.size());
		occ::Mat Cd(nbf, nd);
		for (int k = 0; k < nd; k++) Cd.col(k) = P * sp.eigenvectors().col(keep[k]) / std::sqrt(sp.eigenvalues()(keep[k]));
		Eigen::SelfAdjointEigenSolver<occ::Mat> fd(Cd.transpose() * F * Cd);
		const occ::Mat Cds = Cd * fd.eigenvectors();
		const occ::Vec &ed = fd.eigenvalues();
		//semicanonical MP2 pairs (I, J) with F_II + F_JJ as the occupied energies, the prescreening of DLPNO
		occ::Mat Lc(nbf, nc);
		for (int k = 0; k < nc; k++) Lc.col(k) = L.col(cand[k]);
		const occ::Mat B = df.build_b_tilde(Lc, Cds);  //row k*nd + a: (J_k a|P)
		ivec cl{ I };
		occ::Mat Dv = occ::Mat::Zero(nd, nd);
		for (int k = 0; k < nc; k++)
		{
			const occ::Mat K = B.topRows(nd) * B.middleRows(static_cast<Eigen::Index>(k) * nd, nd).transpose();  //(Ia|Jb)
			occ::Mat t(nd, nd);
			double e_pair = 0.0;
			for (int b = 0; b < nd; b++)
				for (int a = 0; a < nd; a++)
				{
					t(a, b) = K(a, b) / (fdiag(I) + fdiag(cand[k]) - ed(a) - ed(b));
					e_pair += t(a, b) * (2.0 * K(a, b) - K(b, a));
				}
			if (k > 0 && !(std::abs(e_pair) > opts.pair_threshold)) continue;
			if (k > 0) cl.push_back(cand[k]);
			//the MP2 virtual density of pairs (I, J) and (J, I): sum_c t_ac (2 t_bc - t_cb), t_JI = t_IJ^T
			Dv += t * (2.0 * t - t.transpose()).transpose();
			if (k > 0) Dv += t.transpose() * (2.0 * t.transpose() - t).transpose();
		}
		sel[I] = cl;
		doms[I] = dom;
		Cdom[I] = Cds;
		Dvir[I] = 0.5 * (Dv + Dv.transpose());
		sum_cand += nc;
		sum_dom += nd;
	}
	mark("MP2 pair screening");
	//merging: an LMO whose selection lies inside a larger one's is one more central LMO of that cluster
	ivec order(n), owner(n, -1);
	for (int I = 0; I < n; I++) order[I] = I;
	std::stable_sort(order.begin(), order.end(), [&](const int a, const int b) { return sel[a].size() > sel[b].size(); });
	ivec2 groups;
	for (int x = 0; x < n; x++)
	{
		const int I = order[x];
		if (owner[I] >= 0) continue;
		bvec in(n, false);
		for (int k = 0; k < static_cast<int>(sel[I].size()); k++) in[sel[I][k]] = true;
		ivec g{ I };
		owner[I] = static_cast<int>(groups.size());
		for (int y = x + 1; y < n; y++)
		{
			const int J = order[y];
			if (owner[J] >= 0) continue;
			bool inside = true;
			for (int k = 0; k < static_cast<int>(sel[J].size()) && inside; k++) inside = in[sel[J][k]];
			if (!inside) continue;
			g.push_back(J);
			owner[J] = static_cast<int>(groups.size());
		}
		groups.push_back(g);
	}
	//per cluster: the leader's LMOs, and the natural orbitals above lno_threshold of the members' summed MP2 virtual
	//densities in their joint PAO domain (local natural orbitals)
	for (int gi = 0; gi < static_cast<int>(groups.size()); gi++)
	{
		const ivec &g = groups[gi], &cl = sel[g[0]];
		bvec use(nbf, false);
		for (int m = 0; m < static_cast<int>(g.size()); m++)
			for (int k = 0; k < static_cast<int>(doms[g[m]].size()); k++) use[doms[g[m]][k]] = true;
		ivec dom;
		for (int mu = 0; mu < nbf; mu++)
			if (use[mu]) dom.push_back(mu);
		occ::Mat P(nbf, dom.size());
		for (size_t k = 0; k < dom.size(); k++) P.col(k) = Q.col(dom[k]);
		Eigen::SelfAdjointEigenSolver<occ::Mat> sp(P.transpose() * S * P);
		ivec keep;
		for (int k = 0; k < sp.eigenvalues().size(); k++)
			if (sp.eigenvalues()(k) > 1e-6) keep.push_back(k);
		const int nu = static_cast<int>(keep.size());
		occ::Mat Cu(nbf, nu);
		for (int k = 0; k < nu; k++) Cu.col(k) = P * sp.eigenvectors().col(keep[k]) / std::sqrt(sp.eigenvalues()(keep[k]));
		const occ::Mat CuS = Cu.transpose() * S;
		occ::Mat Dv = occ::Mat::Zero(nu, nu);
		for (int m = 0; m < static_cast<int>(g.size()); m++)
		{
			const occ::Mat X = CuS * Cdom[g[m]];
			Dv += X * Dvir[g[m]] * X.transpose();
		}
		Eigen::SelfAdjointEigenSolver<occ::Mat> no(Dv);
		ivec lno;
		for (int k = 0; k < nu; k++)
			if (no.eigenvalues()(k) > opts.lno_threshold) lno.push_back(k);
		const int v = static_cast<int>(lno.size());
		occ::Mat Cv(nbf, v);
		for (int k = 0; k < v; k++) Cv.col(k) = Cu * no.eigenvectors().col(lno[k]);
		const int o = static_cast<int>(cl.size());
		occ::Mat Co(nbf, o);
		for (int k = 0; k < o; k++) Co.col(k) = L.col(cl[k]);
		//semicanonical within the cluster's occupied and within its virtual space, so occ's canonical CCSD applies;
		//F_ov vanishes because the virtuals are orthogonal to every occupied orbital of the converged SCF
		Eigen::SelfAdjointEigenSolver<occ::Mat> fo(Co.transpose() * F * Co), fv(Cv.transpose() * F * Cv);
		const occ::Mat &R = fo.eigenvectors();
		const occ::Mat Cvs = Cv * fv.eigenvectors();
		occ::qm::MolecularOrbitals mo;
		mo.kind = occ::qm::SpinorbitalKind::Restricted;
		mo.n_alpha = mo.n_beta = o;
		mo.n_ao = nbf;
		mo.C.resize(nbf, o + v);
		mo.C << Co * R, Cvs;
		mo.energies.resize(o + v);
		mo.energies << fo.eigenvalues(), fv.eigenvalues();
		const occ::qm::cc::CCIntegrals eris = use_df ? occ::qm::cc::df_eris(df, mo, 0) : build_eris(basis, mo, 0, opts, backend, aux_name, log);
		ivec central(g.size());
		occ::Mat proj(o, g.size());  //the central LMOs in the semicanonical orbitals
		for (int m = 0; m < static_cast<int>(g.size()); m++)
		{
			central[m] = static_cast<int>(std::find(cl.begin(), cl.end(), g[m]) - cl.begin());
			proj.col(m) = R.row(central[m]).transpose();
		}
		const CcSolution s = solve_cc(eris, opts, log, quiet, &proj);
		const occ::Mat H = cluster_half_density(s, R, central, opts.triples);
		occ::Mat Ccl(nbf, o + v);
		Ccl << Co, Cvs;
		Dhalf += Ccl * H * Ccl.transpose();
		const Tn ovov = from_tensor(eris.ovov);
		for (int m = 0; m < static_cast<int>(g.size()); m++) e_corr += cluster_energy(s, R, ovov, central[m]);
		e_t += s.tr.e_t * g.size() / o;  //a log figure: (T) shared evenly by the cluster's LMOs
		sum_o += o;
		sum_v += v;
		max_o = std::max(max_o, o);
		max_v = std::max(max_v, v);
	}
	mark("cluster CCSD densities");
	const occ::Mat D = 2.0 * Cocc * Cocc.transpose() + Dhalf + Dhalf.transpose();
	const occ::Mat dm = C.transpose() * S * D * S * C;
	const std::string name = opts.triples ? "local CCSD(T)" : "local CCSD";
	result.e_corr = e_corr;
	result.e_triples = opts.triples ? e_t : 0.0;
	natural_orbitals(scf, dm, name, scf.energy.total + e_corr + result.e_triples, result, log);
	mark("density and natural orbitals");
	log << "\n" << name << " density (cluster-in-molecule, unrelaxed Lambda densities), " << (use_df ? "DF integrals, " + fit : std::string("exact integrals")) << "\n"
		<< "  " << n << " Pipek-Mezey LMOs, frozen core: " << nf << " orbitals; candidate radius " << opts.cluster_radius << " A, domain threshold "
		<< opts.domain_threshold << ", pair threshold " << opts.pair_threshold << " Eh, LNO threshold " << opts.lno_threshold << "\n"
		<< "  candidates " << std::fixed << std::setprecision(1) << double(sum_cand) / n << " LMOs, PAO domains " << double(sum_dom) / n << " virtuals on average\n"
		<< "  " << groups.size() << " clusters: occupied " << double(sum_o) / groups.size() << " on average (max " << max_o << "), virtual "
		<< double(sum_v) / groups.size() << " (max " << max_v << ")\n"
		<< std::setprecision(10) << "  E(CCSD correlation, sum over clusters) = " << e_corr << "\n"
		<< "  natural occupations:";
	for (int k = 0; k < std::min(nmo, 8); k++) log << " " << std::setprecision(6) << result.occupations(k);
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
	opts.local = scf_table["ccsd_local"].value_or(opts.local);
	opts.cluster_radius = scf_table["ccsd_cluster_radius"].value_or(opts.cluster_radius);
	opts.domain_threshold = scf_table["ccsd_domain_threshold"].value_or(opts.domain_threshold);
	opts.pair_threshold = scf_table["ccsd_pair_threshold"].value_or(opts.pair_threshold);
	opts.lno_threshold = scf_table["ccsd_lno_threshold"].value_or(opts.lno_threshold);
	config.method.name = "rhf";
	//exact CCSD on an exact-J/K reference, as ORCA's NoRI
	const std::string backend = occ::util::to_lower_copy(opts.backend);
	if (backend != "df" && backend != "ri") config.method.ri_policy = occ::io::RIPolicy::None;
	const occ::qm::Wavefunction scf = occ::main::run_scf_external(config, false);
	return opts.local ? local_ccsd_natural_orbitals(scf, opts, log) : ccsd_natural_orbitals(scf, opts, log);
}

occ::qm::Wavefunction wavefunction_from_occ_input(const std::string &toml_path, const occ::io::OccInput &config, std::ostream &log)
{
	const auto kind = occ::driver::parse_method_string(config.method.name).kind;
	if (kind != occ::driver::MethodKind::CCSD && kind != occ::driver::MethodKind::CCSD_T) return occ::main::run_scf_external(config, true);
	return ccsd_density_from_occ_input(toml_path, config, log).natural;
}
