#!/usr/bin/env python
"""Independent Owens check for non-cyclic H_1(Sigma_2) (owens_independent.mvector needs a cyclic group).

Minima: owens_independent.class_minima (exact, integer arithmetic) for Qtilde (classes = keys A x mod D, A = adjugate,
adding componentwise, spin at 0); plumbing side d(Y, g) = -m_{-Q}(g) by v2's m_Q (the validated sharp formula;
ClassMap coordinates, spin at 0).  The Owens step (forms, minima of Qtilde, isomorphisms, matching) is independent.  Group isomorphisms G(Qtilde) -> G(P) are enumerated from a basis of
G(Qtilde) found by brute force (direct-sum decomposition), with images of matching orders, checked bijective.
Forms: every (m1, m2, a) with the required parity, 0 <= a <= D//2 + 2 (not only reduced ones).  Also recomputes
Traczyk's eps and dim H1(Sigma_2; F3) from V_K(e^{i pi/3}) with exact cyclotomic arithmetic.
Usage: owens_general_group.py KNOT [KNOT ...]      Output: owens_general_group.json (merged)
"""
import itertools, json, math, os, sys, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'bernhard_jablan', 'code'))
sys.path.insert(0, os.path.join(ROOT, 'lib'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
import owens_independent as OI
from owens_traczyk_8_6_4_1 import qt_ind, gram_reduce
import numpy as np
import owens_obstruction as OO
import plumbing as P
import kinfo
from owens_u3_verify import forms_all


def group_of(keys, mods):
    elems = list(keys)
    add = lambda a, b: tuple((x + y) % m for x, y, m in zip(a, b, mods))
    zero = tuple(0 for _ in elems[0])
    def order(g):
        k, x = 1, g
        while x != zero:
            x = add(x, g); k += 1
        return k
    return elems, add, zero, order


def span(basis, orders, add, zero):
    out = {}
    for cs in itertools.product(*[range(o) for o in orders]):
        x = zero
        for c, b in zip(cs, basis):
            for _ in range(c):
                x = add(x, b)
        out[x] = cs
    return out


def find_basis(elems, add, zero, order):
    N = len(elems)
    byord = sorted(elems, key=order, reverse=True)
    for k in (1, 2, 3):
        for basis in itertools.combinations(byord[: min(len(byord), 400)], k):
            ords = [order(b) for b in basis]
            if math.prod(ords) != N:
                continue
            sp = span(basis, ords, add, zero)
            if len(sp) == N:
                return list(basis), ords
    raise RuntimeError('no basis with at most 3 generators found')


def isomorphisms(basis, ords, addA, zeroA, elemsB, addB, zeroB, orderB):
    """All isomorphisms G_A -> G_B as dicts, from images of the basis with order dividing ord(b)."""
    spanA = span(basis, ords, addA, zeroA)            # element -> coefficient tuple
    cand = [[h for h in elemsB if o % orderB(h) == 0] for o in ords]
    for imgs in itertools.product(*cand):
        phi = {}
        for x, cs in spanA.items():
            y = zeroB
            for c, h in zip(cs, imgs):
                for _ in range(c):
                    y = addB(y, h)
            phi[x] = y
        if len(set(phi.values())) == len(spanA):
            yield phi


def jones_at_omega(J):
    """V(e^{i pi/3}) in Z[w]/(w^2 - w + 1) as (a, b) meaning a + b w, w = e^{i pi/3}; then eps and d with
    V = eps (i sqrt 3)^d (Lickorish-Millett), i sqrt3 = 2w - 1."""
    pw = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]    # w^k for k mod 6 (w^2 = w - 1)
    a = b = 0
    for e, c in J.items():
        x, y = pw[e % 6]; a += c * x; b += c * y
    v = (a, b)
    mul = lambda p, q: (p[0] * q[0] - p[1] * q[1], p[0] * q[1] + p[1] * q[0] + p[1] * q[1])
    r = (1, 0)
    for d in range(0, 12):
        for eps in (1, -1):
            if v == (eps * r[0], eps * r[1]):
                return eps, d
        r = mul(r, (-1, 2))                                    # times i sqrt 3 = 2w - 1
    return None, None


if __name__ == '__main__':
    T = json.load(open('owens_u3_targets.json'))
    outf = 'owens_general_group.json'
    res = json.load(open(outf)) if os.path.exists(outf) else {}
    for n in sys.argv[1:]:
        r = T[n]; f = r['plumbing']; D = r['D']; sigma = r['sigma']
        n_even = 2 if abs(sigma) == 4 else 1
        Q = P.plumbing_matrix(f['weights'], f['edges'])
        # plumbing side: v2's exact m_Q (ClassMap coordinates add componentwise mod cm.d; spin at 0) -- the 8x8 box
        # search of class_minima is too large here; the Owens step below is computed independently
        mP, cm = OO.m_Q(-np.array(Q, dtype=object))
        assert cm.order == D
        dY = {g: -v for g, v in mP.items()}
        elemsB, addB, zeroB, orderB = group_of(dY.keys(), list(cm.d))
        fl = forms_all(D, n_even, D // 2 + 2)
        fit = {1: [], -1: []}; checked = 0
        print(n, 'D', D, 'forms', len(fl), 'plumbing classes', len(dY), flush=True)
        for g in fl:
            mQ, DQ, _ = OI.class_minima(gram_reduce(qt_ind(*g)))      # same lattice, smaller box
            elemsA, addA, zeroA, orderA = group_of(mQ.keys(), [D] * 4)
            if sorted(orderA(x) for x in elemsA) != sorted(orderB(x) for x in elemsB):
                continue                                       # non-isomorphic groups: no phi
            checked += 1
            basis, ords = find_basis(elemsA, addA, zeroA, orderA)
            for s in (1, -1):
                for phi in isomorphisms(basis, ords, addA, zeroA, elemsB, addB, zeroB, orderB):
                    if all((mQ[x] - s * dY[phi[x]]) >= 0 and (mQ[x] - s * dY[phi[x]]).denominator == 1
                           and (mQ[x] - s * dY[phi[x]]).numerator % 2 == 0 for x in mQ):
                        fit[s].append(g); break
        J = kinfo.parse_jones(n)
        eps, d3 = jones_at_omega(J if sigma > 0 else {-e: c for e, c in J.items()})
        res[n] = dict(D=D, sigma=sigma, n_even=n_even, forms=len(fl), forms_with_isomorphic_group=checked,
                      fit_plus=fit[1][:4], fit_minus=fit[-1][:4], blocked=[s for s in (1, -1) if not fit[s]],
                      traczyk_eps=eps, dim3=d3, traczyk_excludes_02=(d3 == 2 and eps == -1),
                      v2code=dict(verdict=r['verdict'], eps=r.get('traczyk_eps'), dim3=r.get('dim3')))
        print(n, res[n], flush=True)
        json.dump(res, open(outf, 'w'), indent=1, default=str)
