"""Independent re-implementation of the Borodzik-Friedl rank-2 obstruction.

Written from scratch for the review; imports nothing from the project.

For a 2-change sequence (Theorem 1.1) there is a hermitian 2x2 A(t) with
A(1)=diag(e1,e2), e_i in {+1,-1}, presenting the Blanchfield form.  C=A(-1) then
satisfies (BF Lemma 3.2, Lemma 3.3, proof of Theorem 3.4):
  * |det C| = D
  * sign(C) = sigma + e1 + e2          (sigma = sign(V+V^T))
  * C_ii = e_i mod 4, C_12 even        (in the basis where A(1) is diagonal)
  * l(C) ~= 2 l(V+V^T)
Modes:
  'signed' : all four conditions (the project's refined test)
  'mod2'   : only |det C|=D, C congruent to I mod 2 in some basis, l(C)~=2l
             (BF Proposition 5.3 as published, for |sigma|<4; no signature use)
"""
import itertools
import json
import sys
from fractions import Fraction

import sympy


def charpoly_signature(M):
    """Exact signature of a symmetric integer matrix: all roots real, so
    Descartes' rule of signs is exact."""
    from sympy.polys.matrices import DomainMatrix
    dM = DomainMatrix([[sympy.ZZ(v) for v in row] for row in M], (len(M), len(M)), sympy.ZZ)
    coeffs = [int(c) for c in dM.charpoly()]
    assert coeffs[-1] != 0, 'singular'

    def changes(cs):
        cs = [c for c in cs if c != 0]
        return sum(1 for u, v in zip(cs, cs[1:]) if (u > 0) != (v > 0))

    n = len(coeffs) - 1
    pos = changes(coeffs)
    neg = changes([c * (-1) ** (n - i) for i, c in enumerate(coeffs)])
    assert pos + neg == n
    return pos - neg


class LinkingForm:
    """H = Z^n / M Z^n with lambda(x,y) = x^T M^{-1} y mod 1, enumerated by BFS."""

    def __init__(self, M):
        n = len(M)
        from sympy.polys.matrices import DomainMatrix
        dM = DomainMatrix([[sympy.ZZ(v) for v in row] for row in M], (n, n), sympy.ZZ)
        det = int(dM.det())
        assert det != 0
        self.D = D = abs(det)
        inv = dM.convert_to(sympy.QQ).inv().to_Matrix()
        self.B = [[int(inv[i, j] * D) for j in range(n)] for i in range(n)]  # M^{-1} = B / D
        assert all((inv[i, j] * D).q == 1 for i in range(n) for j in range(n))
        # check M B = D I exactly
        assert all(sum(M[i][k] * self.B[k][j] for k in range(n)) == (D if i == j else 0) for i in range(n) for j in range(n))
        self.n = n
        # BFS enumeration of the group; key(x) = B x mod D is canonical on H
        Bcols = [[self.B[i][j] % D for i in range(n)] for j in range(n)]
        zero = tuple([0] * n)
        reps = {zero: zero}          # key -> representative x (mod D)
        frontier = [(zero, zero)]
        while frontier:
            nxt = []
            for kx, x in frontier:
                for i in range(n):
                    col = Bcols[i]
                    ky = tuple((kx[t] + col[t]) % D for t in range(n))
                    if ky not in reps:
                        y = list(x); y[i] = (y[i] + 1) % D; y = tuple(y)
                        reps[ky] = y
                        nxt.append((ky, y))
            frontier = nxt
        assert len(reps) == D, (len(reps), D)
        self.elems = list(reps.values())
        self.Bx = [list(kk) for kk in reps.keys()]
        self.self2 = [(2 * sum(x[i] * bx[i] for i in range(n))) % D for x, bx in zip(self.elems, self.Bx)]
        self.by_self2 = {}
        for idx, s in enumerate(self.self2):
            self.by_self2.setdefault(s, []).append(idx)
        import numpy as np
        E = np.array(self.elems, dtype=np.int64)
        BXm = np.array(self.Bx, dtype=np.int64)
        assert D < 3_000_000 // max(1, n)
        self.P = (2 * (E @ BXm.T)) % D
        assert all(self.P[i, i] == self.self2[i] for i in range(D))
        self.by_self2_np = {s: np.array(v) for s, v in self.by_self2.items()}

    def pair2(self, i, j):
        x = self.elems[i]
        by = self.Bx[j]
        return (2 * sum(x[k] * by[k] for k in range(self.n))) % self.D

    def invariant_factors(self):
        S = sympy.matrices.normalforms.smith_normal_form(sympy.Matrix(self._M), domain=sympy.ZZ)
        return S


def gl2_z4():
    out = []
    for p in itertools.product(range(4), repeat=4):
        a, b, c, d = p
        if (a * d - b * c) % 2 == 1:
            out.append(((a, b), (c, d)))
    assert len(out) == 96
    return out


GL2Z4 = gl2_z4()


def congruent_residues(a, b, c):
    """Set of (C11 mod 4, C22 mod 4, C12 mod 2) over all bases mod 4 for C=[[a,c],[c,b]]."""
    res = set()
    for (p11, p12), (p21, p22) in GL2Z4:
        # columns of P are new basis vectors u=(p11,p21), w=(p12,p22)
        def q(u1, u2, w1, w2):
            return a * u1 * w1 + c * (u1 * w2 + u2 * w1) + b * u2 * w2
        n11 = q(p11, p21, p11, p21) % 4
        n22 = q(p12, p22, p12, p22) % 4
        n12 = q(p11, p21, p12, p22) % 2
        res.add((n11, n22, n12))
    return res


def isometric(LF, a, b, c, detC):
    """l([[a,c],[c,b]]) ~= 2 l(M)?  Lemma 5.4: exists v1,v2 with
    2l(v_i,v_j) = (C^{-1})_ij mod 1.  C^{-1} = [[b,-c],[-c,a]]/detC.
    Returns (verdict, reason)."""
    D = LF.D
    s = 1 if detC > 0 else -1  # 1/detC = s/D
    t11, t22, t12 = (s * b) % D, (s * a) % D, (-s * c) % D
    L1 = LF.by_self2.get(t11, [])
    L2 = LF.by_self2.get(t22, [])
    if not L1 or not L2:
        return False, 'self-value-missing'
    import numpy as np
    sub = LF.P[np.ix_(LF.by_self2_np[t11], LF.by_self2_np[t22])]
    hits = np.argwhere(sub == t12)
    if len(hits):
        i = int(LF.by_self2_np[t11][hits[0][0]])
        j = int(LF.by_self2_np[t22][hits[0][1]])
        assert LF.pair2(i, j) == t12
        return True, (i, j)
    return False, 'cross-value'


def candidates(D, detC, bound):
    """All symmetric [[a,c],[c,b]] with ab-c^2=detC, |a|,|b|<=bound, 0<=c<=bound."""
    out = []
    for a in range(-bound, bound + 1):
        for c in range(0, bound + 1):
            if a == 0:
                if -c * c == detC:
                    for b in range(-bound, bound + 1):
                        out.append((a, b, c))
                continue
            num = detC + c * c
            if num % a == 0:
                b = num // a
                if abs(b) <= bound:
                    out.append((a, b, c))
    return out


def candidates_reduced(D, detC):
    """Lemma 5.2 representatives: 0<|a|<=|b|<=D, 0<=c<=|a|/2 (so |a|<=sqrt(4D/3)),
    plus the isotropic family [[a,r],[r,0]], r^2=D, |a|<=r, when detC=-D is -square."""
    from math import isqrt
    out = []
    amax = isqrt(4 * D // 3) + 2
    for a in range(-amax, amax + 1):
        if a == 0:
            continue
        for c in range(0, abs(a) // 2 + 1):
            num = detC + c * c
            if num % a == 0:
                b = num // a
                if abs(a) <= abs(b) <= D:
                    out.append((a, b, c))
    r = isqrt(D)
    if r * r == D and detC == -D:
        for a in range(-r, r + 1):
            out.append((a, 0, r))
    return out


def sig2(a, b, c):
    det = a * b - c * c
    if det < 0:
        return 0
    return 2 if a > 0 else -2


def analyse(V, mode='signed', bound=None, want_witness=False, reduced=False):
    n = len(V)
    M = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    sigma = charpoly_signature(M)
    LF = LinkingForm(M)
    D = LF.D
    bound = bound or D
    cache = {}
    result = {'D': D, 'sigma': sigma, 'patterns': {}}
    cand_by_det = {d: (candidates_reduced(D, d) if reduced else candidates(D, d, bound)) for d in (D, -D)}
    res_cache = {}

    def residues(a, b, c):
        k = (a % 4, b % 4, c % 4)
        if k not in res_cache:
            res_cache[k] = congruent_residues(a, b, c)
        return res_cache[k]

    def iso(a, b, c, detC):
        k = (a % D, b % D, c % D, detC > 0)
        if k not in cache:
            cache[k] = isometric(LF, a, b, c, detC)
        return cache[k]

    if mode == 'signed':
        patterns = [(1, 1), (1, -1), (-1, -1)]
    else:
        patterns = ['mod2']
    survived_any = False
    for pat in patterns:
        info = {'n_candidates': 0, 'n_residue_ok': 0, 'n_sig_ok': 0,
                'reject_reasons': {}, 'witness': None}
        if pat == 'mod2':
            dets = [D, -D]
        else:
            e1, e2 = pat
            req = sigma + e1 + e2
            if req not in (-2, 0, 2):
                info['verdict'] = 'EXCLUDED_SIGNATURE'
                result['patterns'][str(pat)] = info
                continue
            dets = [D] if req != 0 else [-D]
        for detC in dets:
            for (a, b, c) in cand_by_det[detC]:
                info['n_candidates'] += 1
                rs = residues(a, b, c)
                if pat == 'mod2':
                    ok = any(r[0] % 2 == 1 and r[1] % 2 == 1 and r[2] == 0 for r in rs)
                else:
                    ok = (e1 % 4, e2 % 4, 0) in rs
                if not ok:
                    continue
                info['n_residue_ok'] += 1
                if pat != 'mod2' and sig2(a, b, c) != req:
                    continue
                info['n_sig_ok'] += 1
                v, why = iso(a, b, c, detC)
                if v:
                    info['witness'] = {'C': [[a, c], [c, b]], 'v_idx': why}
                    break
                info['reject_reasons'][why] = info['reject_reasons'].get(why, 0) + 1
            if info['witness']:
                break
        info['verdict'] = 'SURVIVES' if info['witness'] else 'EXCLUDED_LINKING'
        if info['witness']:
            survived_any = True
        result['patterns'][str(pat)] = info
    result['obstructed'] = not survived_any
    return result


if __name__ == '__main__':
    pass
