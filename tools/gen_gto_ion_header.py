"""Scalar-relativistic charge states for the heavy elements the Slater ion tables do not cover.

    average  <molden> <out.npz>          spherically averaged density of one atomic state
    thakkar  <Thakkar_coefs.h> <Z> ...   electron count of the ported Thakkar neutral (self-check)
    header   <dir> <Thakkar_coefs.h> <out.h> [z_min=56]
                                         GTO_ion_delta_coefs.h from <dir>/<El>/q<q>.npz

The states come from ORCA (UHF, DKH2, SARC-DKH-TZVP plus one diffuse s/p/d
shell, lowest multiplicity per charge state; ion_states.sh) as molden files.
Koga/Thakkar is non-relativistic HF, which is wrong for a U core, so for
these elements every state - the neutral included - is the DKH2 density. It
is stored as a correction to the Thakkar neutral that NoSpherA2 already
evaluates analytically:

    T_q = rho_q(DKH2) - rho_0(Thakkar),   q = -1 (bound anions only), 0, 1 .. 4

rho_q = rho_0(Thakkar) + T_q is then the DKH2 density itself; storing the
difference keeps the steep core out of the interpolated table.
"""
import sys, os, re, math
import numpy as np
_trapz = getattr(np, "trapezoid", None) or np.trapz

N_R, R_MIN, R_MAX = 300, 1e-5, 25.0          # down to 1e-5 bohr: T_q carries the relativistic core
N_S, S_MIN, S_MAX = 481, 0.0, 6.0                # 0.0125: the 6s/7s rearrangement puts structure below s = 0.2
A0 = 0.529177210903                          # Bohr in Angstrom; k = 4 pi s a0 (bohr^-1)
R_TAB = np.exp(np.linspace(math.log(R_MIN), math.log(R_MAX), N_R))
S_TAB = np.linspace(S_MIN, S_MAX, N_S)
# Fine log grid for the quadratures: at r = 20 bohr it is 0.02 bohr apart, eight points per
# period of sin(k r) at the largest k (40 bohr^-1); 1e-7 resolves a Z = 103 1s cusp.
R_FINE = np.exp(np.linspace(math.log(1e-7), math.log(60.0), 30000))
Q_MIN, Q_MAX = -1, 4
ELEMENTS = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se "
            "Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy "
            "Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf "
            "Es Fm Md No Lr").split()
LETTERS = "spdfghi"


# ---------------------------------------------------------------- molden -> spherical average
def read_molden(path):
    """Shells [(l, exps, coefs)] of the single atom and MOs [(occ, coef vector)]."""
    shells, mos, sec, cur = [], [], None, None
    lines = open(path).read().splitlines()
    i = 0
    while i < len(lines):
        t = lines[i].strip()
        if t.startswith("["):
            sec = t.split("]")[0][1:].upper()
        elif sec == "GTO" and t:
            w = t.split()
            if w[0].lower() in LETTERS and len(w) >= 2:
                l, n = LETTERS.index(w[0].lower()), int(w[1])
                e, c = [], []
                for j in range(n):
                    x = lines[i + 1 + j].replace("D", "E").split()
                    e.append(float(x[0])); c.append(float(x[1]))
                shells.append((l, np.array(e), np.array(c)))
                i += n
        elif sec == "MO" and t:
            if t.startswith("Occup="):
                cur = [float(t.split("=")[1]), []]
                mos.append(cur)
            elif "=" not in t:
                cur[1].append(float(t.split()[1]))
        i += 1
    return shells, [(o, np.array(c)) for o, c in mos]


def gauss_moment(l, p):
    """Int_0^inf r^(2l+2) exp(-p r^2) dr."""
    return math.gamma(l + 1.5) / (2.0 * p ** (l + 1.5))


def radial_parts(shells, prim_normalised):
    """Normalised contracted radial functions (Int R^2 r^2 dr = 1) as (l, exps, weights)."""
    out = []
    for l, e, c in shells:
        w = c * np.array([gauss_moment(l, 2 * a) ** -0.5 for a in e]) if prim_normalised else c.copy()
        G = np.array([[gauss_moment(l, a + b) for b in e] for a in e])
        out.append((l, e, w / math.sqrt(w @ G @ w)))
    return out


def spherical_average(path):
    """rho(r) averaged over the sphere on R_FINE: on one centre only equal (l, m) AO pairs
    survive the angular integral, each contributing P_ij R_i R_j / (4 pi)."""
    shells, mos = read_molden(path)
    aos = [(si, m) for si, (l, _, _) in enumerate(shells) for m in range(2 * l + 1)]
    nao = len(aos)
    occ = [(o, c) for o, c in mos if o > 1e-10]
    assert all(len(c) == nao for _, c in occ), "MO length != number of pure AOs (cartesian molden?)"
    best = None
    # Molden contraction coefficients are meant for normalised primitives; the raw reading is
    # kept only as a fall-back and has to win the normalisation test to be used.
    for prim_norm in (True, False):
        rad = radial_parts(shells, prim_norm)
        S = np.zeros((nao, nao))
        for a, (sa, ma) in enumerate(aos):
            la, ea, wa = rad[sa]
            for b, (sb, mb) in enumerate(aos):
                lb, eb, wb = rad[sb]
                if la == lb and ma == mb:
                    S[a, b] = wa @ np.array([[gauss_moment(la, x + y) for y in eb] for x in ea]) @ wb
        err = max(abs(c @ S @ c - 1.0) for _, c in occ)
        if best is None or err < best[0]:
            best = (err, rad, S)
    err, rad, S = best
    assert err < 1e-5, f"{path}: MOs not normalised in either reading ({err:.2e})"
    P = sum(o * np.outer(c, c) for o, c in occ)
    rho = np.zeros_like(R_FINE)
    for l in sorted({s[0] for s in rad}):
        idx = [si for si, s in enumerate(rad) if s[0] == l]
        Rl = np.array([(R_FINE ** l) * (rad[si][2][:, None] * np.exp(-np.outer(rad[si][1], R_FINE ** 2))).sum(0)
                       for si in idx])
        for m in range(2 * l + 1):
            cols = [aos.index((si, m)) for si in idx]
            rho += np.einsum("ig,ij,jg->g", Rl, P[np.ix_(cols, cols)], Rl) / (4 * math.pi)
    return rho, float(np.sum(P * S)), float(integrate(rho)), err


# ---------------------------------------------------------------- Thakkar neutral, ported
def read_c_array(text, name):
    body = re.search(r"\b" + name + r"\[\d*\]\s*=?\s*\{(.*?)\};", text, re.S).group(1)
    body = re.sub(r"//[^\n]*", "", body)
    return [float(x) for x in body.replace("\n", " ").split(",") if x.strip()]


class ThakkarNeutral:
    """Thakkar::get_radial_density (spherical_density.cpp) in numpy: Slater orbitals
    sum c r^(n-1) exp(-z r) per occupied subshell, rho = sum occ Orb^2 / 4 pi."""

    def __init__(self, header):
        t = open(header).read()
        self.a = {k: read_c_array(t, "Thakkar_" + k) for k in ("nex", "ns", "np", "nd", "nf", "occ", "n", "z", "c")}

    def rho(self, Z, r):
        a = self.a
        if Z == 1:
            return 6.0835 * np.exp(-2.3 * r) / (4 * math.pi)
        occ, off = a["occ"], (Z - 1) * 19
        nr_ex = int(sum(a["nex"][:Z - 1]))
        nr_coef = 0                            # Spherical_Atom::previous_element_coef
        for t in range(Z - 2, -1, -1):
            for lo, hi, key in ((0, 7, "ns"), (7, 13, "np"), (13, 17, "nd"), (17, 19, "nf")):
                nr_coef += sum(int(a[key][t]) for m in range(lo, hi) if occ[t * 19 + m] != 0)
        orb = np.zeros((19, len(r)))
        for lo, hi, key in ((0, 7, "ns"), (7, 13, "np"), (13, 17, "nd"), (17, 19, "nf")):
            for _ in range(int(a[key][Z - 1])):
                prim = r ** (a["n"][nr_ex] - 1) * np.exp(-a["z"][nr_ex] * r)
                for m in range(lo, hi):
                    if occ[off + m] == 0:
                        continue
                    orb[m] += a["c"][nr_coef] * prim
                    nr_coef += 1
                nr_ex += 1
        return sum(occ[off + m] * orb[m] ** 2 for m in range(19)) / (4 * math.pi)


# ---------------------------------------------------------------- quadratures
def integrate(rho, r=R_FINE):
    """4 pi Int r^2 rho dr on a log grid (dr = r dlnr)."""
    return float(_trapz(4 * math.pi * r ** 3 * rho, np.log(r)))


def form_factor(rho, s, r=R_FINE):
    """4 pi Int r^2 rho j0(k r) dr, k = 4 pi s a0 in bohr^-1."""
    lr = np.log(r)
    return np.array([_trapz(4 * math.pi * r ** 3 * rho * np.sinc(4 * math.pi * si * A0 * r / math.pi), lr)
                     for si in s])


def cubic_at(t, x):
    """spherical_density.cpp cubic_at, vectorised: 4-point Lagrange, one-sided in the end intervals."""
    n = len(t)
    j = np.clip(x.astype(int), 1, n - 3)
    f = x - j
    cub = (-f * (f - 1) * (f - 2) / 6 * t[j - 1] + (f + 1) * (f - 1) * (f - 2) / 2 * t[j]
           - (f + 1) * f * (f - 2) / 2 * t[j + 1] + (f + 1) * f * (f - 1) / 6 * t[j + 2])
    return np.where(x >= n - 1, t[-1], cub)


# ---------------------------------------------------------------- commands
def cmd_average(molden, out):
    rho, nt, nq, err = spherical_average(molden)
    np.savez(out, rho=rho, n_trace=nt, n_quad=nq, norm_err=err)
    print(f"{molden}: N(trace) {nt:.6f}  N(quad) {nq:.6f}  MO norm err {err:.1e}")


def cmd_thakkar(header, *zs):
    th = ThakkarNeutral(header)
    for z in zs:
        print(f"Z={z}: N = {integrate(th.rho(int(z), R_FINE)):.6f}")


def cmd_header(d, thakkar_h, out, z_min="56"):
    """z_min: lightest element taken from <dir>. 56 (Ba) leaves Rb-Cs on Slater/Thakkar;
    37 would make the 5th row relativistic too."""
    th = ThakkarNeutral(thakkar_h)
    span = Q_MAX - Q_MIN + 1
    qmax, anion, offset, rho_t, f_t, log = [0] * 103, [0] * 103, [-1] * (103 * span), [], [], []
    for Z, el in enumerate(ELEMENTS, 1):
        if Z < int(z_min) or not os.path.exists(os.path.join(d, el, "q+0.npz")):
            continue
        r0 = th.rho(Z, R_FINE)
        n0 = integrate(r0)
        f0 = form_factor(r0, S_TAB) * (Z / n0)
        qs = [0]
        for q in range(1, Q_MAX + 1):          # contiguous: rho_q interpolates between neighbours
            if not os.path.exists(os.path.join(d, el, f"q{q:+d}.npz")):
                break
            qs.append(q)
        if os.path.exists(os.path.join(d, el, "q-1.npz")):   # ion_states.sh keeps bound anions only
            qs.insert(0, -1)
        for q in qs:
            rq = np.load(os.path.join(d, el, f"q{q:+d}.npz"))["rho"]
            nq = integrate(rq)
            # Exact electron counts on both sides, so f_T(0) = -q to the last digit: the
            # quadrature leaves each ~1e-6 off and the evaluator sums them.
            rq = rq * ((Z - q) / nq)
            fq = form_factor(rq, S_TAB)
            fq *= (Z - q) / fq[0]
            T = rq - r0 * (Z / n0)
            offset[(Z - 1) * span + (q - Q_MIN)] = len(rho_t)
            rho_t.append(np.interp(np.log(R_TAB), np.log(R_FINE), T))
            f_t.append(fq - f0)
            # What NoSpherA2 will see: the tables read back through the C++ interpolation.
            x = (np.log(np.clip(R_FINE, R_MIN, None)) - math.log(R_MIN)) / (math.log(R_MAX / R_MIN) / (N_R - 1))
            Ti = np.where(R_FINE <= R_MAX, cubic_at(rho_t[-1], x), 0.0)
            s_mid = S_TAB[:-1] + 0.5 * (S_TAB[1] - S_TAB[0])
            df = form_factor(rq, s_mid) - form_factor(r0 * (Z / n0), s_mid) \
                 - cubic_at(f_t[-1], (s_mid - S_MIN) / (S_TAB[1] - S_TAB[0]))
            log.append(f"{el} q={q:+d} N(GTO)={nq:.6f} N(Thakkar)={n0:.6f} T(s=0)={fq[0] - f0[0]:+.6f}"
                       f"  read back: N err {integrate(r0 * (Z / n0) + Ti) - (Z - q):+.1e}"
                       f"  max|df| at s midpoints {np.abs(df).max():.1e}")
        qmax[Z - 1] = qs[-1]
        anion[Z - 1] = int(qs[0] == -1)
    fmt = lambda a: ",\n".join("\t" + ", ".join(f"{v:.10e}" for v in a[i:i + 6]) for i in range(0, len(a), 6))
    ifmt = lambda a: ",\n".join("\t" + ", ".join(str(v) for v in a[i:i + 12]) for i in range(0, len(a), 12))
    with open(out, "w", newline="\n") as h:
        h.write("#pragma once\n// Generated by tools/gen_gto_ion_header.py, do not edit by hand.\n"
                "// Scalar-relativistic charge states for the heavy elements without a Slater ion: UHF/DKH2 atoms from ORCA\n"
                "// (SARC-DKH-TZVP plus one diffuse s/p/d shell, lowest multiplicity per state), spherically averaged.\n"
                "// T_q = rho_q(DKH2) - rho_0(Thakkar), so rho_0(Thakkar) + T_q is the DKH2 density of charge state q;\n"
                "// q = -1 only for bound anions. Normalised so that T_q integrates to exactly -q.\n"
                "// r: 300 points log-uniform over [1e-5, 25] Bohr; s = sin(theta)/lambda: 481 points over [0, 6] 1/Angstrom,\n"
                "// k = 4 pi s a0. GTOIon_qmax[Z-1] = highest cation state (0 = only the neutral, -1 = element absent),\n"
                "// GTOIon_anion[Z-1] = 1 if q = -1 exists, GTOIon_offset[(Z-1)*6 + (q+1)] = table index or -1.\n\n")
        h.write(f"const int GTOIon_n_r = {N_R};\nconst int GTOIon_n_s = {N_S};\n"
                f"const int GTOIon_q_min = {Q_MIN};\nconst int GTOIon_q_max = {Q_MAX};\n"
                f"const double GTOIon_r_min = {R_MIN:g};\nconst double GTOIon_r_max = {R_MAX:g};\n"
                f"const double GTOIon_s_min = {S_MIN:g};\nconst double GTOIon_s_max = {S_MAX:g};\n\n")
        present = [offset[(Z - 1) * span - Q_MIN] >= 0 for Z in range(1, 104)]
        h.write(f"const int GTOIon_qmax[103] = {{\n{ifmt([qmax[i] if present[i] else -1 for i in range(103)])}\n}};\n\n")
        h.write(f"const int GTOIon_anion[103] = {{\n{ifmt(anion)}\n}};\n\n")
        h.write(f"const int GTOIon_offset[{len(offset)}] = {{\n{ifmt(offset)}\n}};\n\n")
        h.write(f"const double GTOIon_rho[{len(rho_t) * N_R}] = {{\n{fmt(np.concatenate(rho_t))}\n}};\n\n")
        h.write(f"const double GTOIon_f[{len(f_t) * N_S}] = {{\n{fmt(np.concatenate(f_t))}\n}};\n")
    print("\n".join(log))
    print(f"{len(rho_t)} tables, {sum(present)} elements, {sum(anion)} anions -> {out}")


if __name__ == "__main__":
    {"average": cmd_average, "thakkar": cmd_thakkar, "header": cmd_header}[sys.argv[1]](*sys.argv[2:])
