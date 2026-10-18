#!/usr/bin/env python
"""Casson-Walker invariant of a plumbing tree in the normalisation of Section 3.4 of the paper, for any signature:
    L(Y) = (|det Q|/12) * ( tr Q - 3 sigma(Q) + sum_v (2 - d_v) (Q^{-1})_{vv} )
(Lescop's plumbing formula; for negative-definite trees sigma = -s and this is graph_verify.L_closed).
Checks on every algebraic KnotInfo diagram (validate_plumbing.json): L(raw tree) = L(normal form), and exactly one
orientation SIGMA of the construction satisfies L = R(K) (graph_u1.R_of), the same for all knots.
Output: cw_general.json
"""
import json, os, sys, warnings
warnings.filterwarnings('ignore')
from collections import Counter
from fractions import Fraction
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from validate_plumbing import setup
from algebraic_plumbing import plumbing_from_algebraic, normal_form, plumbing_matrix


def L_general(w, edges):
    import sympy
    from identify_trees import exact_inertia
    Q = sympy.Matrix(plumbing_matrix(w, edges))
    n = Q.rows; H = abs(Q.det()); Qi = Q.inv()
    deg = Counter(v for e in edges for v in e)
    I = exact_inertia(w, edges)
    assert I['nullity'] == 0
    sig = I['pplus'] - I['pminus']
    val = H * (sum(w) - 3 * sig + sum((2 - deg[i]) * Qi[i, i] for i in range(n))) / 12
    num, den = sympy.fraction(sympy.nsimplify(val))
    return Fraction(int(num), int(den))


def one(n):
    setup()
    import graph_u1 as GU
    pd = json.loads(GU.BY[n]['pd_notation'])
    R = GU.R_of(n)
    rec = {}
    for s in (1, -1):
        w, e = plumbing_from_algebraic(pd, s)
        Lr = L_general(w, e)
        wn, en = normal_form(w, e)
        Ln = L_general(wn, en)
        rec[s] = (Lr == Ln, Lr == R)
    return n, rec


if __name__ == '__main__':
    V = json.load(open(os.path.join(HERE, 'validate_plumbing.json')))['records']
    names = [n for n, r in V.items() if r['D'] and r[str(1)]['weights'] is not None]
    out = {}
    with Pool(8) as p:
        for n, r in p.imap_unordered(one, names, chunksize=20):
            out[n] = r
    inv = Counter(r[1][0] and r[-1][0] for r in out.values())
    which = Counter(tuple(s for s in (1, -1) if r[s][1]) for r in out.values())
    print('L(raw) = L(normal form):', dict(inv))
    print('orientations with L = R:', {str(k): v for k, v in which.items()})
    json.dump({n: {str(s): list(v) for s, v in r.items()} for n, r in sorted(out.items())}, open(os.path.join(HERE, 'cw_general.json'), 'w'))
