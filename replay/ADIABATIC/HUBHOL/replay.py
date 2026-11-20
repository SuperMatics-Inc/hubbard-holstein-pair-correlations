#!/usr/bin/env python3
"""HUBHOL replay: Hubbard-Holstein extension of STAGE7R Theorem 3.3.

Exact Fraction arithmetic and outward-rounded mpmath.iv intervals (70 digits)
check every finite claim of RIGOR.md. Lines tagged NUMERICAL are floating
cross-checks (quadrature, matrix exponentials) and are not proofs.

Default exits 0; --all-mutants requires every guarded child to exit 1.
"""

import argparse
import subprocess
import sys
from fractions import Fraction as Q
from pathlib import Path

import mpmath as mp

mp.iv.dps = 70
iv = mp.iv

MUTANTS = (
    "lf-uncentered",          # Hubbard term written as U_c n_up n_dn without shifting mu
    "ueff-sign",              # atomic reference with (U_c-U)/2 instead of (U-U_c)/2
    "mott-gap-no-mu",         # Mott gap G=|U_eff|/2 ignoring mu
    "jw-gauge",               # frustrate the K22 four-cycle: no sign gauge
    "sos-c1",                 # edge SOS with 1*g^2/alpha instead of 2*g^2/alpha
    "sos-ndouble-dilute",     # dilute cost 1[n=2] instead of 1[n>=1]
    "shiba-site-only",        # Shiba flip of one site only
    "ph-for-neven",           # claim full p-h flip maps N_odd to N_even
    "cap-6-5",                # attractive/Mott cap 6/5 instead of 5/4
    "dilute-cap-5-4",         # dilute cap 5/4 instead of 3/2
    "mott-uniform-mu",        # Mott absorption |U_eff|/4 for every mu
    "dilute-L1024-source",    # dilute canonical transfer at L=1024
    "source-L512",            # attractive source disk at L=512
    "refined-q-600",          # refined rate 600*sqrt(eps) instead of 640*sqrt(eps)
    "stage4-budget",          # STAGE4 extension at beta*t<=1/256
    "phi-bound-1-8",          # Lemma 5.1 with 1/8 instead of 1/2
    "barrier-rate-quarter",   # single-site phonon cost (U/4)Delta instead of (U/2)Delta
    "edge-witness-bt-half",   # Mott-edge witness at beta*t=1/2
    "free-bond-cosh",         # free single-bond activity cosh^2-1 (one spin dropped)
    "free-witness-d1",        # free-fermion witness at distance 1
)


class Rejected(Exception):
    pass


def require(ok, name):
    if not ok:
        raise Rejected(name)


def ivq(x):
    """Exact rational -> tight interval."""
    x = Q(x)
    return iv.mpf(x.numerator) / iv.mpf(x.denominator)


def lt(a, b):
    """Certified a < b for intervals/rationals."""
    a = ivq(a) if isinstance(a, (int, Q)) else iv.mpf(a)
    b = ivq(b) if isinstance(b, (int, Q)) else iv.mpf(b)
    return a.b < b.a


def icosh(x):
    return (iv.exp(x) + iv.exp(-x)) / 2


def isinh(x):
    return (iv.exp(x) - iv.exp(-x)) / 2


# ----------------------------------------------------------------------------
# Section 1: Lang-Firsov reference (Lemma 1.2)
# ----------------------------------------------------------------------------

SITE_STATES = [(0, 0), (1, 0), (0, 1), (1, 1)]


def site_energy(ueff, mu, nu, nd):
    n = nu + nd
    return Q(ueff, 2) * (1 if n == 1 else 0) - mu * n


def check_lf_identity(mutant):
    points = [(Q(1), Q(0), Q(0)), (Q(0), Q(1), Q(0)), (Q(0), Q(0), Q(1)), (Q(0), Q(0), Q(0)),
              (Q(7, 3), Q(-5, 2), Q(1, 9)), (Q(1, 2), Q(4), Q(-3, 7))]
    # affine independence of the first four points: det of [1,U,Uc,mu] rows
    m = [[Q(1), *p] for p in points[:4]]

    def det(a):
        if len(a) == 1:
            return a[0][0]
        return sum((-1) ** j * a[0][j] * det([r[:j] + r[j + 1:] for r in a[1:]]) for j in range(len(a)))

    require(det(m) != 0, "lf-affine-points")
    for (U, Uc, mu) in points:
        for (nu, nd) in SITE_STATES:
            n = nu + nd
            if mutant == "lf-uncentered":
                hub = Uc * nu * nd
            else:
                hub = Uc * (nu - Q(1, 2)) * (nd - Q(1, 2))
            lhs = -Q(U, 2) * (n - 1) ** 2 + hub - mu * n
            ueff = (Uc - U) if mutant == "ueff-sign" else (U - Uc)
            rhs = site_energy(ueff, mu, nu, nd) + (Q(Uc, 4) - Q(U, 2))
            require(lhs == rhs, "lf-identity")
    # Remark 1.1: uncentered U_c n_up n_dn - mu' N == centered with mu = mu' - U_c/2 (+ const)
    for (U, Uc, mu1) in points:
        consts = set()
        for (nu, nd) in SITE_STATES:
            n = nu + nd
            a = Uc * nu * nd - mu1 * n
            b = Uc * (nu - Q(1, 2)) * (nd - Q(1, 2)) - (mu1 - Q(Uc, 2)) * n
            consts.add(a - b)
        require(len(consts) == 1, "uncentered-shift")
    return len(points)


# ----------------------------------------------------------------------------
# Section 2: regimes, bond structure, edge SOS
# ----------------------------------------------------------------------------

def regime(ueff, mu, mutant=None):
    """Return (name, G, cost) per Definition 2.1; cost maps (n) -> 0/1."""
    if ueff > 0:
        return "A", Q(ueff, 2), lambda n: 1 if n == 1 else 0
    a = -ueff  # |U_eff|
    if abs(mu) < Q(a, 2):
        g = Q(a, 2) if mutant == "mott-gap-no-mu" else Q(a, 2) - abs(mu)
        return "M", g, lambda n: 0 if n == 1 else 1
    if mu < -Q(a, 2):
        return "D-", -mu - Q(a, 2), lambda n: 1 if n >= 1 else 0
    if mu > Q(a, 2):
        return "D+", mu - Q(a, 2), lambda n: 1 if n <= 1 else 0
    return "edge", Q(0), None


def check_gap_property(mutant):
    samples = []
    for ueff in (Q(3), Q(1, 5), Q(0), Q(-1), Q(-6), Q(-41, 7)):
        for mu in (Q(-9), Q(-3), Q(-1, 2), Q(-1, 10), Q(0), Q(1, 3), Q(2), Q(13, 4)):
            samples.append((ueff, mu))
    counted = {"A": 0, "M": 0, "D-": 0, "D+": 0, "edge": 0}
    for ueff, mu in samples:
        name, G, cost = regime(ueff, mu, mutant)
        counted[name] += 1
        if cost is None:
            continue
        require(G > 0, "gap-positive")
        es = [site_energy(ueff, mu, nu, nd) for nu, nd in SITE_STATES]
        emin = min(es)
        for (nu, nd), e in zip(SITE_STATES, es):
            n = nu + nd
            if cost(n):
                require(e - emin >= G, f"gap-{name}")
        if name == "D-":
            require(site_energy(ueff, mu, 1, 1) - emin >= 2 * G, "gap-D-double")
        if name == "D+":
            require(site_energy(ueff, mu, 0, 0) - emin >= 2 * G, "gap-D+empty")
    require(all(counted[k] > 0 for k in ("A", "M", "D-", "D+", "edge")), "regime-coverage")
    return counted


def ann(j):
    """Jordan-Wigner annihilation matrix on 16 bond states; dict (row,col)->sign."""
    out = {}
    for s in range(16):
        if (s >> j) & 1:
            sign = (-1) ** bin(s & ((1 << j) - 1)).count("1")
            out[(s ^ (1 << j), s)] = sign
    return out


def matmul_sparse(a, b):
    out = {}
    for (i, k), x in a.items():
        for (k2, j), y in b.items():
            if k == k2:
                out[(i, j)] = out.get((i, j), 0) + x * y
    return {k: v for k, v in out.items() if v != 0}


def transpose(a):
    return {(j, i): v for (i, j), v in a.items()}


def bond_hopping(mutant):
    c = [ann(j) for j in range(4)]
    T = [[0] * 16 for _ in range(16)]
    for sg in (0, 1):
        for term in (matmul_sparse(transpose(c[sg]), c[2 + sg]), matmul_sparse(transpose(c[2 + sg]), c[sg])):
            for (i, j), v in term.items():
                T[i][j] += v
    if mutant == "jw-gauge":
        T[3][6] = -T[3][6]
        T[6][3] = -T[6][3]
    return T


def signless_B():
    """Independent construction: one spin moves across the bond."""
    B = [[0] * 16 for _ in range(16)]
    for a in range(16):
        for j in (0, 1):
            mask = (1 << j) | (1 << (j + 2))
            if bin(a & mask).count("1") == 1:
                B[a][a ^ mask] = 1
    return B


def components(B):
    seen, comps = set(), []
    for r in range(16):
        if r in seen:
            continue
        stack, comp = [r], set()
        while stack:
            a = stack.pop()
            if a in comp:
                continue
            comp.add(a)
            stack.extend(b for b in range(16) if B[a][b])
        seen |= comp
        comps.append(tuple(sorted(comp)))
    return sorted(comps)


def is_psd(M):
    """Exact symmetric LDL^T test with rational entries."""
    M = [row[:] for row in M]
    n = len(M)
    for k in range(n):
        p = M[k][k]
        if p < 0:
            return False
        if p == 0:
            if any(M[k][j] != 0 for j in range(k + 1, n)):
                return False
            continue
        for i in range(k + 1, n):
            f = M[i][k] / p
            if f:
                for j in range(k + 1, n):
                    M[i][j] -= f * M[k][j]
    return True


def site_n(s, site):
    return ((s >> (2 * site)) & 1) + ((s >> (2 * site + 1)) & 1)


COSTS = {
    "odd": lambda n: 1 if n == 1 else 0,
    "even": lambda n: 0 if n == 1 else 1,
    "occ": lambda n: 1 if n >= 1 else 0,
    "hole": lambda n: 1 if n <= 1 else 0,
}


def gamma_bond(cost):
    return [cost(site_n(s, 0)) + cost(site_n(s, 1)) for s in range(16)]


def check_bond(mutant):
    T = bond_hopping(mutant)
    B = signless_B()
    require(all(abs(T[i][j]) == B[i][j] for i in range(16) for j in range(16)), "abs-T-equals-B")
    require(all(T[i][j] == T[j][i] for i in range(16) for j in range(16)), "T-symmetric")
    comps = components(B)
    expected = [(0,), (1, 4), (2, 8), (3, 6, 9, 12), (5,), (7, 13), (10,), (11, 14), (15,)]
    require(comps == expected, "components")
    for a in (3, 12):
        for b in (6, 9):
            require(B[a][b] == 1 and B[b][a] == 1, "k22-complete")
    require(B[3][12] == 0 and B[6][9] == 0, "k22-bipartite")
    # Lemma 6.1: sign gauge
    sgn = [0] * 16
    consistent = True
    for r in range(16):
        if sgn[r]:
            continue
        sgn[r] = 1
        stack = [r]
        while stack:
            a = stack.pop()
            for b in range(16):
                if T[a][b]:
                    want = sgn[a] * (1 if T[a][b] > 0 else -1)
                    if sgn[b] == 0:
                        sgn[b] = want
                        stack.append(b)
                    elif sgn[b] != want:
                        consistent = False
    require(consistent, "jw-gauge-consistent")
    require(all(sgn[i] * T[i][j] * sgn[j] >= 0 for i in range(16) for j in range(16)), "jw-gauge-nonneg")
    # STAGE4 hypothesis ||T_e|| <= 2: 4I - T^2 PSD
    T2 = [[sum(T[i][k] * T[k][j] for k in range(16)) for j in range(16)] for i in range(16)]
    M = [[Q(4 * (i == j) - T2[i][j]) for j in range(16)] for i in range(16)]
    require(is_psd(M), "norm-T-le-2")
    return T, B, comps, sgn


def check_sos(B, comps, mutant):
    c_const = Q(1) if mutant == "sos-c1" else Q(2)
    costs = dict(COSTS)
    if mutant == "sos-ndouble-dilute":
        costs["occ"] = lambda n: 1 if n == 2 else 0
    xs = [Q(0), Q(1, 64), Q(1, 4), Q(1, 2), Q(1), Q(3), Q(250)]
    kpq = {}
    for name, cost in costs.items():
        gam = gamma_bond(cost)
        # structural certificate (proof of Theorem 2.4)
        for comp in comps:
            if len(comp) == 2:
                require(max(gam[comp[0]], gam[comp[1]]) >= 1, f"sos-pair-{name}")
            elif len(comp) == 4:
                require(gam[3] == gam[12] and gam[6] == gam[9], f"sos-k22-const-{name}")
                require(gam[3] + gam[6] >= 2, f"sos-k22-sum-{name}")
                kpq[name] = (gam[3], gam[6])
        # exact PSD at sample ratios x = g/alpha
        for x in xs:
            R = [[(Q(gam[i]) + c_const * x * x if i == j else Q(0)) - x * B[i][j] for j in range(16)]
                 for i in range(16)]
            require(is_psd(R), f"sos-psd-{name}-x={x}")
    if mutant != "sos-ndouble-dilute":
        require(kpq == {"odd": (0, 2), "even": (2, 0), "occ": (1, 2), "hole": (1, 2)}, "sos-k22-table")
    # optimality of constant 2 for N_odd: c=1, x=1 gives det -1 on the K22 + block
    gp, gq, x, c1 = 0, 2, 1, 1
    require(gp * gq + c1 * x * x * (gp + gq) + c1 * c1 * x ** 4 - 4 * x * x == -1, "sos-constant-2-needed")
    # Remark 2.5 conjugations
    shiba = 0b0010 if mutant == "shiba-site-only" else 0b1010
    perm = [s ^ shiba for s in range(16)]
    require(all(B[perm[i]][perm[j]] == B[i][j] for i in range(16) for j in range(16)), "shiba-B")
    g_odd, g_even = gamma_bond(COSTS["odd"]), gamma_bond(COSTS["even"])
    require(all(g_odd[perm[s]] == g_even[s] for s in range(16)), "shiba-odd-even")
    ph = [s ^ 0b1111 for s in range(16)]
    require(all(B[ph[i]][ph[j]] == B[i][j] for i in range(16) for j in range(16)), "ph-B")
    g_occ, g_hole = gamma_bond(COSTS["occ"]), gamma_bond(COSTS["hole"])
    require(all(g_occ[ph[s]] == g_hole[s] for s in range(16)), "ph-occ-hole")
    target = g_even if mutant == "ph-for-neven" else g_odd
    require(all(g_odd[ph[s]] == target[s] for s in range(16)), "ph-odd-invariant")
    return len(xs) * len(costs), kpq


# ----------------------------------------------------------------------------
# Section 3: caps
# ----------------------------------------------------------------------------

def check_cap_identities():
    vs = [Q(0), Q(1, 7), Q(2, 5), Q(1, 2), Q(3, 4), Q(1), Q(2), Q(-3)]
    for v in vs:
        require(5 * (1 + v * v) - 4 * (1 + v) == 5 * (v - Q(2, 5)) ** 2 + Q(1, 5), "cap-am-identity")
        require(3 * (1 + v * v) ** 2 - 2 * (1 + 2 * v + v ** 3)
                == (1 - 2 * v) ** 2 + v * v * (3 * (v - Q(1, 3)) ** 2 + Q(5, 3)), "cap-d-identity")
    for v in vs[1:5]:
        R = lambda y: (1 + y / v) / (1 + y)
        require(R(v * v) == (1 + v) / (1 + v * v), "R-at-v2")
        require(R(2 * v * v + v ** 4) == (1 + 2 * v + v ** 3) / (1 + v * v) ** 2, "R-at-D")
        for y1, y2 in ((Q(0), Q(1, 3)), (Q(1, 5), Q(4)), (Q(2), Q(9))):
            require(R(y2) - R(y1) == (y2 - y1) * (1 / v - 1) / ((1 + y1) * (1 + y2)), "R-monotone-identity")
            require(R(y2) > R(y1), "R-monotone")


def site_weights(b, m, absorb, cost):
    """z_a and z_b for dimensionless b=beta*U_eff, m=beta*mu, absorb=beta*G/2 (intervals)."""
    za = iv.mpf(0)
    zb = iv.mpf(0)
    for nu, nd in SITE_STATES:
        n = nu + nd
        be = (b / 2 if n == 1 else iv.mpf(0)) - m * n
        w = iv.exp(-be)
        za += w
        zb += w * (iv.exp(absorb) if cost(n) else 1)
    return za, zb


def check_caps(mutant):
    cap_am = Q(6, 5) if mutant == "cap-6-5" else Q(5, 4)
    cap_d = Q(5, 4) if mutant == "dilute-cap-5-4" else Q(3, 2)
    worst_a = 4 * iv.log(1 / (iv.sqrt(2) - 1))  # v = sqrt2 - 1
    worst_d = 2 * iv.log(ivq(Q(5, 2)))          # v = 2/5 at U_eff=0
    samples = []
    for b in (ivq(Q(1, 10)), ivq(1), worst_a, ivq(10), ivq(40)):
        for m in (ivq(-5), ivq(0), ivq(Q(3, 10)), ivq(7)):
            samples.append(("A", b, m))
    for b in (ivq(-1), ivq(-6), ivq(-30)):
        for frac in (Q(0), Q(3, 10), Q(-9, 10), Q(99, 100)):
            samples.append(("M", b, -b / 2 * ivq(frac)))
    for b in (ivq(0), ivq(-1), ivq(-6)):
        for g in (ivq(Q(1, 100)), worst_d, ivq(1), ivq(9)):
            samples.append(("D-", b, b / 2 - g))
            samples.append(("D+", b, -b / 2 + g))
    worst = {"A": iv.mpf(0), "M": iv.mpf(0), "D": iv.mpf(0)}
    for name, b, m in samples:
        if name == "A":
            G2, cost, cap = b / 4, COSTS["odd"], cap_am
        elif name == "M":
            G2, cost, cap = (-b / 2 - abs(m)) / 2, COSTS["even"], cap_am
        elif name == "D-":
            G2, cost, cap = (-m + b / 2) / 2, COSTS["occ"], cap_d
        else:
            G2, cost, cap = (m + b / 2) / 2, COSTS["hole"], cap_d
        za, zb = site_weights(b, m, G2, cost)
        ratio = zb / za
        v = iv.exp(-G2)
        # y bound (Lemma 3.3): costed/uncosted weight ratio, computed directly
        w1 = iv.mpf(0)
        w0 = iv.mpf(0)
        for nu, nd in SITE_STATES:
            n = nu + nd
            wt = iv.exp(-((b / 2 if n == 1 else iv.mpf(0)) - m * n))
            if cost(n):
                w1 += wt
            else:
                w0 += wt
        y = w1 / w0
        ybound = v * v if name in ("A", "M") else 2 * v * v + v ** 4
        require(y.b <= (ybound * (1 + ivq(Q(1, 10 ** 60)))).b, f"y-bound-{name}")
        require(ratio.b < ivq(cap), f"cap-{name}")
        key = "D" if name.startswith("D") else name
        if ratio.a > worst[key].a:
            worst[key] = ratio
    # attaining witnesses: (1+sqrt2)/2 in A, (1+2v+v^3)/(1+v^2)^2 at v=2/5 in D
    require(worst["A"].a > ivq(Q(6, 5)), "cap-am-attained")
    require(worst["D"].a > ivq(Q(138, 100)), "cap-d-attained")
    # Remark 3.5: uniform Mott absorption |U_eff|/4 fails outside the window
    za, zb = site_weights(ivq(-8), ivq(-10), ivq(2), COSTS["even"])
    bad = zb / za
    if mutant == "mott-uniform-mu":
        require(bad.b < ivq(cap_am), "mott-uniform-cap")
    require(bad.a > 7, "mott-uniform-witness")
    return worst, bad


# ----------------------------------------------------------------------------
# Section 4: cluster constants
# ----------------------------------------------------------------------------

def cluster(c, L):
    q = 64 * c / L
    qs = 144 * c / L
    return {
        "eps_max": Q(1, 64 * L * L),
        "q": q,
        "A": 3 * c / (1 - q),
        "site": 3 * c * q / (1 - q),
        "qs": qs,
        "src": 3 * c * qs / (1 - qs),
    }


def check_constants(mutant):
    require(lt(iv.exp(1), 3), "e<3")
    require(lt(iv.exp(ivq(Q(1, 4))), Q(4, 3)), "e^(1/4)<4/3")
    require(Q(1, 8) <= Q(1, 4), "source-penalty")
    L_am = 512 if mutant == "source-L512" else 1024
    am = cluster(Q(5, 4), L_am)
    require(am["site"] <= 1 and am["src"] <= 1, "am-budgets")
    require(am["eps_max"] == Q(1, 2 ** 26), "am-eps")
    require((am["q"], am["A"], am["site"], am["qs"], am["src"])
            == (Q(5, 64), Q(240, 59), Q(75, 236), Q(45, 256), Q(675, 844)), "am-table")
    dgc = cluster(Q(3, 2), 1024)
    require(dgc["site"] <= 1 and (dgc["q"], dgc["A"], dgc["site"]) == (Q(3, 32), Q(144, 29), Q(27, 58)), "d-gc")
    require(dgc["src"] == Q(243, 202) and dgc["src"] > 1, "d-1024-source-fails")
    L_d = 1024 if mutant == "dilute-L1024-source" else 2048
    dcan = cluster(Q(3, 2), L_d)
    require(dcan["src"] <= 1, "d-canonical-source")
    require((dcan["q"], dcan["A"], dcan["site"], dcan["qs"], dcan["src"], dcan["eps_max"])
            == (Q(3, 64), Q(288, 61), Q(27, 122), Q(27, 256), Q(243, 458), Q(1, 2 ** 28)), "d-can-table")
    # Corollary 4.3: q = 512 c sqrt(eps); at eps=2^-26, sqrt(eps)=2^-13
    coef = 600 if mutant == "refined-q-600" else 640
    require(512 * Q(5, 4) == coef, "refined-coef")
    require(Q(coef, 2 ** 13) == am["q"] == Q(5, 64), "refined-q-at-edge")
    require(512 * Q(3, 2) == 768 and Q(768, 2 ** 13) == Q(3, 32), "refined-q-dilute")
    # Theorem 4.1 attractive form: G = U_eff/2 so eps<=2^-26 <=> beta t^2/U_eff <= 2^-27
    require(Q(1, 2 ** 26) / 2 == Q(1, 2 ** 27), "attractive-threshold")
    # chain |w(E)| <= c (4c/(3L))^k with v <= k+1: check e^{k/4} c^{k+1} L^-k <= c(4c/(3L))^k at k=1..6
    for k in range(1, 7):
        lhs = iv.exp(ivq(Q(k, 4))) * ivq(Q(5, 4)) ** (k + 1) / ivq(1024) ** k
        rhs = ivq(Q(5, 4)) * (ivq(Q(4, 1) * Q(5, 4) / (3 * 1024))) ** k
        require(lhs.b <= rhs.a, "activity-chain")
        require(iv.exp(ivq(k + 1)).b <= 3 ** (k + 1), "e^v<=3^(k+1)")
    return am, dgc, dcan


def check_stage4(mutant):
    s = Q(1, 128) if mutant == "stage4-budget" else Q(1, 512)  # beta*v_e = 2*beta*t
    x = s / (1 - s)
    require((iv.exp(ivq(s)) - 1).b <= ivq(x).a, "x-bound")
    q = 16 * 3 * x
    budget = 3 * q / (1 - q)
    require(budget < 1, "stage4-budget")
    q6 = 36 * 3 * x
    b6 = 3 * q6 / (1 - q6)
    require(b6 < 1, "stage4-source-budget")
    C = 3 / (1 - q)
    require((x, q, budget, C, q6, b6) == (Q(1, 511), Q(48, 511), Q(144, 463), Q(1533, 463), Q(108, 511), Q(324, 403)),
            "stage4-table")
    return C, q


# ----------------------------------------------------------------------------
# Section 5: phonon factor
# ----------------------------------------------------------------------------

def phi_closed(U, Om, beta, times, jumps):
    """STAGE6 closed form, -log of the normalized oscillator trace (sum c = 0)."""
    alpha2 = U / (2 * Om)
    x = beta * Om / 2
    tot = iv.mpf(0)
    pairs = sorted(zip(times, jumps), key=lambda p: p[0].a)
    for i in range(len(pairs)):
        for j in range(i + 1, len(pairs)):
            ui, ci = pairs[i]
            uj, cj = pairs[j]
            tot += ci * cj * (icosh(x - Om * (uj - ui)) - icosh(x))
    return alpha2 * tot / isinh(x)


def check_phi(mutant):
    fac = Q(1, 8) if mutant == "phi-bound-1-8" else Q(1, 2)
    trajectories = [
        ([Q(3, 10), Q(9, 10), Q(7, 5), Q(11, 5)], [1, 1, -1, -1]),
        ([Q(1, 2), Q(3, 2)], [-1, 1]),
        ([Q(1, 10), Q(1, 5), Q(2), Q(29, 10)], [1, -1, -1, 1]),
        ([Q(1, 1000), Q(299, 100)], [2, -2]),
    ]
    count = 0
    for U in (ivq(Q(13, 10)), ivq(5)):
        for Om in (ivq(Q(1, 10)), ivq(1), ivq(Q(5, 2)), ivq(50)):
            beta = ivq(3)
            for times, jumps in trajectories:
                ph = phi_closed(U, Om, beta, [ivq(t) for t in times], jumps)
                require(ph.b >= 0, "phi-nonneg")
                bound = ivq(fac) * U * sum(abs(c) for c in jumps) / (Om * (1 - iv.exp(-beta * Om)))
                require(ph.b <= bound.a, "phi-lemma-5.1")
                count += 1
    # two-jump: Phi <= (U/2) Delta  (closed form); mutant: (U/4) Delta
    rate = Q(1, 4) if mutant == "barrier-rate-quarter" else Q(1, 2)
    for Om in (ivq(1), ivq(Q(1, 10)), ivq(40)):
        for delta in (Q(1, 100), Q(1, 2), Q(9, 2)):
            beta = ivq(10)
            ph = phi_closed(ivq(1), Om, beta, [ivq(0), ivq(delta)], [-1, 1])
            require(ph.b <= ivq(rate) * ivq(delta), "phi-two-jump")
            count += 1
    return count


def phi_quadrature():
    """NUMERICAL: closed form vs direct double integral of STAGE6 (2.2)."""
    mp.mp.dps = 20
    U, Om, beta = mp.mpf("1.3"), mp.mpf("2.5"), mp.mpf(3)
    times = [mp.mpf("0.3"), mp.mpf("0.9"), mp.mpf("1.4"), mp.mpf("2.2")]
    jumps = [1, 1, -1, -1]
    knots = [mp.mpf(0)] + times + [beta]
    levels = [0, 1, 2, 1, 0]

    def K(w):
        w = w % beta
        return U * Om / 2 * (mp.exp(-Om * w) + mp.exp(-Om * (beta - w))) / (1 - mp.exp(-beta * Om))

    tot = mp.mpf(0)
    for i in range(5):
        for j in range(5):
            if levels[i] == levels[j]:
                continue
            val = mp.quad(lambda u, v: K(u - v), [knots[i], knots[i + 1]], [knots[j], knots[j + 1]])
            tot += (levels[i] - levels[j]) ** 2 * val
    quad = tot / 4
    closed = phi_closed(ivq(Q(13, 10)), ivq(Q(5, 2)), ivq(3), [ivq(Q(3, 10)), ivq(Q(9, 10)), ivq(Q(7, 5)), ivq(Q(11, 5))], jumps)
    return float(quad), float(closed.mid)


# ----------------------------------------------------------------------------
# Section 6: barrier and witness numbers
# ----------------------------------------------------------------------------

def F_beta(b, r):
    return (b * r - 1 + iv.exp(-b * r)) / (r * r)


def check_barriers(mutant):
    # F_beta(r) >= beta/(2r) for beta*r >= 2
    for br in (2, 3, 10, 96):
        lhs = ivq(br) - 1 + iv.exp(-ivq(br))
        require(lhs.a >= ivq(br).b / 2, "F-lower")
    # Lemma 6.2 witness: U=t=1, U_c=1/2, beta=64 -> r=3/2
    w62 = F_beta(ivq(64), ivq(Q(3, 2))) / 4
    require(w62.a > ivq(Q(95, 9)).a and w62.a > 1, "lemma-6.2-witness")
    # Lemma 6.3: 8 e^-2 < 2 so 8e^-2 * 2^27 < 2^28; and > 1
    b63 = 8 * iv.exp(-2)
    require(b63.b < 2 and b63.a > 1, "lemma-6.3-factor")
    # Lemma 6.4: (cosh(beta t)-1)/4 at beta t = 3 exceeds 2
    bt = ivq(Q(1, 2)) if mutant == "edge-witness-bt-half" else ivq(3)
    w64 = (icosh(bt) - 1) / 4
    require(w64.a > 2, "lemma-6.4-witness")
    # Lemma 6.5(2): free single bond
    power = 2 if mutant == "free-bond-cosh" else 4
    wa = icosh(ivq(Q(51, 200))) ** power - 1
    wb = icosh(ivq(Q(61, 100))) ** power - 1
    require(wa.a > iv.exp(-2).b, "lemma-6.5-e^-2")
    require(wb.a > 1, "lemma-6.5-one")
    # per-spin identity (1+e^x)(1+e^-x)/4 = cosh^2(x/2) at x = 1.22
    x = ivq(Q(61, 50))
    lhs = (1 + iv.exp(x)) * (1 + iv.exp(-x)) / 4
    rhs = icosh(x / 2) ** 2
    require(abs(lhs - rhs).b < ivq(Q(1, 10 ** 60)), "per-spin-identity")
    return w62, w64, wa, wb


COS6 = [Q(1), Q(1, 2), Q(-1, 2), Q(-1), Q(-1, 2), Q(1, 2)]


def free_G0(r, bt):
    tot = iv.mpf(0)
    for n1 in range(6):
        for n2 in range(6):
            ck = COS6[n1] + COS6[n2]
            cr = COS6[(n1 * r[0] + n2 * r[1]) % 6]
            if cr == 0:
                continue
            tot += ivq(cr) / (1 + iv.exp(-2 * bt * ivq(ck)))
    return tot / 36


def check_free_witness(mutant):
    r, d = ((1, 0), 1) if mutant == "free-witness-d1" else ((2, 3), 5)
    require(min(r[0], 6 - r[0]) + min(r[1], 6 - r[1]) == d, "torus-distance")
    G0 = free_G0(r, ivq(1))
    require(G0.a > 0 or G0.b < 0, "G0-sign")
    sq = G0 * G0
    b4 = ivq(Q(1533, 463)) * ivq(Q(48, 511)) ** d
    b7 = ivq(Q(240, 59)) * ivq(Q(5, 64)) ** d
    require(sq.a > b4.b, "free-witness-stage4")
    require(sq.a > b7.b, "free-witness-stage7r")
    # sanity: G0 at r=0 is the density per spin, 1/2 at mu=0
    dens = free_G0((0, 0), ivq(1))
    require(abs(dens - ivq(Q(1, 2))).b < ivq(Q(1, 10 ** 50)), "half-filling-density")
    return sq, b4, b7


def numeric_bonds(T, sgn):
    """NUMERICAL: exact 16x16 bond traces at Omega=infinity vs the lemma lower bounds."""
    mp.mp.dps = 30
    Tg = mp.matrix(16, 16)
    for i in range(16):
        for j in range(16):
            Tg[i, j] = sgn[i] * T[i][j] * sgn[j]

    def w_inf(beta, t, ueff, mu):
        D = mp.matrix(16, 16)
        za = mp.mpf(0)
        for nu, nd in SITE_STATES:
            za += mp.exp(-beta * float(site_energy(Q(ueff), Q(mu), nu, nd)))
        for s in range(16):
            e = 0
            for site in (0, 1):
                nu, nd = (s >> (2 * site)) & 1, (s >> (2 * site + 1)) & 1
                e += site_energy(Q(ueff), Q(mu), nu, nd)
            D[s, s] = -beta * mp.mpf(e.numerator) / e.denominator
        Z = sum(mp.expm(D + beta * t * Tg)[i, i] for i in range(16))
        return Z / za ** 2 - 1

    out = {}
    # Lemma 6.3 at mu = 0
    for ueff in (1, -1):
        beta, t, r = mp.mpf(8), mp.mpf(1) / 2, mp.mpf(abs(ueff))
        w = w_inf(beta, t, ueff, 0)
        b = t ** 2 / 4 * (beta * r - 1 + mp.exp(-beta * r)) / r ** 2
        out[f"6.3 ueff={ueff}"] = (float(w), float(b))
    # Lemma 6.4 at the Mott edge, U_eff = -10, mu = -5, beta t = 3
    w = w_inf(3, 1, -10, -5)
    out["6.4 edge"] = (float(w), float((mp.cosh(3) - 1) / 4))
    # Lemma 6.5(2) at U_eff = 0, mu = 0, beta t = 1.22
    w = w_inf(mp.mpf("1.22"), 1, 0, 0)
    out["6.5 free"] = (float(w), float(mp.cosh(mp.mpf("0.61")) ** 4 - 1))
    return out


# ----------------------------------------------------------------------------

def certificate(mutant=None, verbose=True):
    require(mutant is None or mutant in MUTANTS, "known-mutant")
    n_lf = check_lf_identity(mutant)
    counted = check_gap_property(mutant)
    T, B, comps, sgn = check_bond(mutant)
    n_sos, kpq = check_sos(B, comps, mutant)
    check_cap_identities()
    worst, bad = check_caps(mutant)
    am, dgc, dcan = check_constants(mutant)
    C4, q4 = check_stage4(mutant)
    n_phi = check_phi(mutant)
    w62, w64, wa, wb = check_barriers(mutant)
    sq, b4, b7 = check_free_witness(mutant)
    if not verbose:
        return
    print(f"LF_IDENTITY: h_el = (U_eff/2)1[n=1] - mu n + (U_c/4-U/2) at {n_lf} parameter points (4 affinely independent)")
    print(f"REGIMES: gap property on samples {counted}")
    print(f"BOND: components {comps}; |T|=B; JW gauge consistent; ||T_e||<=2 exact")
    print(f"EDGE_SOS: {n_sos} exact 16x16 PSD checks; K22 (Gamma_P,Gamma_Q) = {kpq}; constant 2 needed")
    print(f"CAPS: worst sampled A/M ratio in [{float(worst['A'].a):.12f}, {float(worst['A'].b):.12f}] "
          f"(= (1+sqrt2)/2) < 5/4; worst sampled D ratio >= {float(worst['D'].a):.10f} < 3/2; "
          f"uniform-Mott counterexample ratio >= {float(bad.a):.6f}")
    print(f"CONSTANTS A/M: q={am['q']} A={am['A']} site={am['site']} q_s={am['qs']} src={am['src']} eps<={am['eps_max']}")
    print(f"CONSTANTS D(GC, L=1024): q={dgc['q']} A={dgc['A']} site={dgc['site']} src={dgc['src']} (>1, no transfer)")
    print(f"CONSTANTS D(L=2048): q={dcan['q']} A={dcan['A']} site={dcan['site']} q_s={dcan['qs']} src={dcan['src']} eps<={dcan['eps_max']}")
    print(f"STAGE4_EXTENSION: beta t<=1/1024 -> C={C4}, q={q4}")
    print(f"PHI: {n_phi} interval checks of Lemma 5.1 and Phi<=(U/2)Delta")
    print(f"BARRIERS: Lemma 6.2 witness w >= {float(w62.a):.10f}; Mott edge (cosh3-1)/4 >= {float(w64.a):.10f}; "
          f"free bond cosh^4-1 at 0.51: {float(wa.a):.10f} > e^-2, at 1.22: {float(wb.a):.10f} > 1")
    print(f"FREE_WITNESS: 6x6, beta t=1, r=(2,3), d=5: G0^2 >= {float(sq.a):.10e} > STAGE4 {float(b4.b):.10e} "
          f"and STAGE7R {float(b7.b):.10e} (ratios {float(sq.a) / float(b4.b):.4f}, {float(sq.a) / float(b7.b):.4f})")
    quad, closed = phi_quadrature()
    print(f"NUMERICAL Phi quadrature {quad:.15f} vs closed form {closed:.15f} (diff {abs(quad - closed):.2e})")
    for k, (w, bnd) in numeric_bonds(T, sgn).items():
        print(f"NUMERICAL bond {k}: w_inf={w:.10f} lemma lower bound={bnd:.10f} holds={w >= bnd * (1 - 1e-12)}")
    print("SCOPE: finite exact/interval witnesses; analytic theorems in RIGOR.md require independent review")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mutant", choices=MUTANTS)
    parser.add_argument("--all-mutants", action="store_true")
    args = parser.parse_args()

    if args.mutant:
        try:
            certificate(mutant=args.mutant, verbose=False)
        except Rejected as exc:
            print(f"REJECTED {args.mutant}: {exc}")
            return 1
        print(f"UNCAUGHT {args.mutant}")
        return 0

    try:
        certificate(verbose=True)
    except Rejected as exc:
        print(f"BASELINE FAILURE: {exc}")
        return 1

    if not args.all_mutants:
        return 0

    script = str(Path(__file__).resolve())
    for mutant in MUTANTS:
        child = subprocess.run(
            [sys.executable, "-B", script, "--mutant", mutant],
            capture_output=True,
            text=True,
            check=False,
        )
        print(child.stdout.strip())
        if child.stderr:
            print(child.stderr.strip())
        print(f"mutant {mutant}: exit {child.returncode}")
        if child.returncode != 1 or not child.stdout.startswith(f"REJECTED {mutant}:"):
            return 1
    print(f"CERTIFIED: all {len(MUTANTS)} mutants rejected with exit 1")
    return 0


if __name__ == "__main__":
    sys.exit(main())
