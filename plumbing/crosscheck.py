#!/usr/bin/env python
"""Cross-checks of the star-type u = 1 obstructions (star_u1.json) that do not depend on SnapPy's peripheral
framing of covers.

(A) Sigma_2 by the orbifold route: fill the exterior with (2,0) (cone angle pi along K) and take the 2-fold cyclic
    cover; its cusp is filled with the lift of mu^2, so the cover is Sigma_2(K).  star_u1 used the (1,0) filling of
    the unfilled cyclic cover, which is not Sigma_2 for some torus knots (5_1 gives Z/5+Z/5, 10_124 gives Z/8).
    The Regina name by the orbifold route must equal the name used in star_u1 (= name of the plumbing's boundary).
(B) Casson-Walker invariant, independent of SnapPy and Regina.
    Rustamov: |H| lambda(Y) = sum_s [chi(HF_red(Y,s)) - d(Y,s)/2].  For Y = boundary of a negative-definite plumbing
    tree with <= 1 bad vertex, HF^+(-Y) is lattice homology (even degrees), so HF_red(Y) sits in odd degrees and
        L(Y) := sum_s d(Y,s) + 2 sum_s rank HF_red(Y,s) = -2 |H| lambda(Y).
    Mullins (v2 eq. casson-sum):  R(K) := -D sigma/4 + D V'(-1) / (6 V(-1)) = -2 D lambda(Sigma_2 K)
    (KnotInfo sigma and Jones).  Sigma_2(K) = s Y forces s L(Y) = R(K); s must be the orientation that the
    linking form fixed in star_u1.
(C) Calibration of (B) and of the linking-form orientation rule on Montesinos knots (KnotInfo notation, 3+ tangles,
    10-12 crossings) with v2's notation plumbings: s from Casson-Walker, s from the linking form, and for alternating
    knots s from d(spin) = -sigma/4 must agree.
Usage: crosscheck.py star | controls   (5 worker processes; results written after every knot)
Output: crosscheck.json (star), crosscheck_controls.json (controls: rank <= 9, det <= 100, up to 30 knots)
"""
import json, math, os, sys, time, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'montesinos'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
import _compat  # noqa
import numpy as np
import snappy, regina
import montesinos_u1 as MU
from hf_red import hf_red_all
import kinfo
import general as G
import plumbing as P
import d_owens as DO
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}


def sigma2_orbifold(pd, tries=6):
    E = snappy.Link(pd).exterior()
    h = None
    for t in range(tries):
        E2 = E.copy()
        if t:
            E2.randomize()
        E2.dehn_fill((2, 0))
        C = E2.covers(2, cover_type='cyclic')[0]
        h = str(C.homology())
        T = regina.Triangulation3(C.filled_triangulation()._to_string()); T.intelligentSimplify()
        rec = regina.StandardTriangulation.recognise(T)
        if rec and rec.manifold():
            return rec.manifold().name(), h
    return None, h


def R_of(name):
    J = kinfo.parse_jones(name); s = int(BY[name]['signature']); D = int(BY[name]['determinant'])
    V1 = sum(c * (-1) ** (e % 2) for e, c in J.items())
    dV1 = sum(e * c * (-1) ** ((e - 1) % 2) for e, c in J.items())
    assert abs(V1) == D
    return Fraction(-D * s, 4) + Fraction(D * dV1, 6 * V1)


def L_of(Q):
    Q = np.array(Q, dtype=np.int64)
    d, cm = MU.correction_terms(Q)
    ranks, _ = hf_red_all(Q)
    if ranks is None:
        return None, None, None
    assert set(ranks) == set(d)
    d_spin = d[tuple(0 for _ in cm.d)]          # ClassMap is centred at the spin structure (also for non-cyclic H_1)
    return sum(d.values()) + 2 * sum(ranks.values()), d_spin, sorted(ranks.values())


def link_sign(Q, seifert):
    VV = (np.array(seifert, dtype=object) + np.array(seifert, dtype=object).T).tolist()
    tl, _ = G.linking_multiset(VV)
    lv, _ = G.linking_multiset(np.array(Q, dtype=np.int64).tolist())
    tn = sorted((-x) % 1 for x in tl)
    if tl == tn:
        return 0 if lv == tl else None          # orientation not detected by the linking form
    return 1 if lv == tl else (-1 if lv == tn else None)


def cw_signs(L, R):
    return [s for s in (1, -1) if s * L == R]


def star_one(n):
    r = RES[n]
    pd = json.loads(BY[n]['pd_notation']); D = int(BY[n]['determinant'])
    nm_orb, h = sigma2_orbifold(pd)
    rec = dict(det=D, verdict=r['verdict'], name_star_u1=r['sigma2'], name_orbifold=nm_orb, H1_orbifold=h,
               names_equal=(nm_orb == r['sigma2']))
    if r.get('plumbing'):
        Q = P.plumbing_matrix(r['plumbing']['weights'], r['plumbing']['edges'])
        L, d_spin, ranks = L_of(Q)
        R = R_of(n)
        s_link = 1 if r['plumbing']['orientation'] == '+' else -1
        sc = cw_signs(L, R) if L is not None else None
        rec.update(L=str(L), R=str(R), hf_red_ranks_nonzero=[x for x in ranks if x] if ranks else None,
                   cw_signs=sc, link_sign=s_link, link_ambiguous=r['plumbing']['orientation_ambiguous'],
                   cw_consistent=(sc is not None and (s_link in sc or r['plumbing']['orientation_ambiguous']) and bool(sc)))
    return n, rec


def control_one(nm):
    r = BY[nm]; D = int(r['determinant'])
    try:
        Q, flipped = MU.plumbing_for(r['montesinos_notation'].strip(), D)
        if Q is None:
            return nm, None
        L, d_spin, ranks = L_of(Q)
        if L is None:
            return nm, None
        R = R_of(nm)
    except Exception as e:
        return nm, dict(error=repr(e)[:80])
    sc = cw_signs(L, R)
    sl = link_sign(Q, json.loads(r['seifert_matrix']))
    sig = int(r['signature'])
    ss = [s for s in (1, -1) if s * d_spin == Fraction(-sig, 4)] if r['alternating'] == 'Y' else None
    ok = bool(sc) and (sl in (0, None) or sl in sc) and (ss is None or any(x in sc for x in ss))
    return nm, dict(notation=r['montesinos_notation'], alternating=r['alternating'], D=D, rank=len(Q), L=str(L), R=str(R),
                    cw_signs=sc, link_sign=sl, spin_signs=ss, consistent=ok)


RES = json.load(open('star_u1.json'))
if __name__ == '__main__':
    from multiprocessing import Pool
    part = sys.argv[1]
    T0 = time.time()
    if part == 'star':
        outf, items, fn = 'crosscheck.json', list(RES), star_one
    else:
        ctl = []
        for nm, r in BY.items():
            m = (r.get('montesinos_notation') or '').strip()
            c = r.get('crossing_number', '')
            if m.startswith('K(') and m.count(';') >= 2 and c.isdigit() and 10 <= int(c) <= 12 and int(r['determinant']) <= 100:
                try:
                    Q, _ = MU.plumbing_for(m, int(r['determinant']))
                except Exception:
                    continue
                if Q is not None and len(Q) <= 9:
                    ctl.append(nm)
        step = max(1, len(ctl) // 30)
        outf, items, fn = 'crosscheck_controls.json', ctl[::step][:30], control_one
    out = {}
    with Pool(5) as pool:
        for n, rec in pool.imap_unordered(fn, items):
            out[n] = rec
            print(f"[{time.time()-T0:5.0f}s] {n:9s} {rec}", flush=True)
            json.dump({'star': out} if part == 'star' else out, open(outf, 'w'), indent=1, default=str)
    if part == 'star':
        print('STAR: orbifold names equal', sum(v['names_equal'] for v in out.values()), '/', len(out),
              '| Casson-Walker consistent', sum(bool(v.get('cw_consistent')) for v in out.values()), '/', sum('L' in v for v in out.values()))
    else:
        print('CONTROLS consistent', sum(bool(v and v.get('consistent')) for v in out.values()), '/', sum(1 for v in out.values() if v and 'L' in v))
