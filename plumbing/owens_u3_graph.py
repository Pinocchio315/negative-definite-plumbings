#!/usr/bin/env python
"""Owens's u = 2 obstruction with graph-manifold correction terms: new u >= 3 for non-alternating knots.

Owens, Theorem 1 (as implemented in v2's owens_obstruction): if K is unknotted by p positive and n negative crossing
changes with n = sigma/2 (KnotInfo convention, sigma(right trefoil) = -2), then Sigma_2(K) bounds Qtilde(m1,m2,a)
(candidates(D, n): reduced forms, exactly n of m1, m2 even) and m_Qtilde(g) - d(Sigma_2, phi(g)) is a non-negative
even integer for some group isomorphism phi (matching_exists).
  |sigma| = 4: u = 2 forces the (0,2) pattern of the sigma > 0 representative, n = 2.
  |sigma| = 2: (1,1) by Owens with n = 1, and (0,2) only when Traczyk excludes it (dim H1(Sigma_2;F3) = 2, eps = -1),
               exactly as in v2's obstruct_u2.
Correction terms: negative-definite plumbing found as in graph_u1 (orbifold-route Regina names, exact name match;
stars included), d = max (K^2+|V|)/4 = -m_{-Q}.
Orientation: Owens's Sigma_2 is the canonical one of the sigma > 0 representative.  Canonical Sigma_2(K) = s_cw Y with
s_cw from the closed-form Casson-Walker identity (graph_verify.py), so the expected orientation is
sign(sigma) * s_cw.  Verdicts are also computed in both orientations: 'OBSTRUCTED' means no form fits in either.
Controls (mode 'controls'): KnotInfo Montesinos knots (<= 12 crossings, v2 notation plumbing) with |sigma| in {2,4}
and exact u in {2,3}: |sigma| = 4, u = 2 knots must pass in the expected orientation; alternating ones must give
v2's Goeritz verdict (obstruct_u2).
Usage: owens_u3_graph.py controls|targets [workers]     Output: owens_u3_controls.json / owens_u3_targets.json
"""
import json, os, sys, time, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'montesinos'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
import _compat  # noqa
import numpy as np
import sympy
import montesinos_u1 as MU
import owens_obstruction as OO
import kinfo
from traczyk import mirror_jones, traczyk_excludes_00_2
import general as G
import plumbing as P
import graph_u1 as GU                      # names_orb, shape_key, R_of; also patches plumbing.regina_name
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
GU.SHAPES['star'] = ('S2',)
GU.SEARCH['star'] = dict(erange=range(-8, 2), chain_w=(-2,), max_chain=0)
v2T = json.load(open(os.path.join(ROOT, 'data', 'v2_final_u_table.json')))


def L_closed(Q):
    S = sympy.Matrix(np.array(Q, dtype=np.int64).tolist())
    n = S.rows; H = abs(S.det()); Si = S.inv()
    deg = [sum(1 for j in range(n) if j != i and S[i, j] != 0) for i in range(n)]
    val = H * (sum(S[i, i] for i in range(n)) + 3 * n + sum((2 - deg[i]) * Si[i, i] for i in range(n))) / 12
    num, den = sympy.fraction(val)
    return Fraction(int(num), int(den))


def find_plumbing_for_knot(n):
    import spherogram
    pd = json.loads(BY[n]['pd_notation']); D = int(BY[n]['determinant'])
    mpd = [list(x) for x in spherogram.Link(pd).mirror().PD_code()]
    nK, nM = GU.names_orb(pd, 16), GU.names_orb(mpd, 16)
    if not any(G.parse_name(x)[0] in GU.SHAPES for x in set(nK) | set(nM)):
        nK, nM = GU.names_orb(pd, 40), GU.names_orb(mpd, 40)
    allnames = set(nK) | set(nM)
    shapes = {}
    for nm in sorted(allnames):
        kind, pieces = G.parse_name(nm)
        if kind in GU.SHAPES and tuple(b for b, _ in pieces) == GU.SHAPES[kind]:
            shapes.setdefault((kind, GU.shape_key(kind, pieces)), nm)
    info = dict(names_sigma2=nK, names_mirror=nM, shapes=len(shapes), searches=[])
    for (kind, _), nm in shapes.items():
        f, st = G.find_plumbing(nm, D, json.loads(BY[n]['seifert_matrix']), alt_names=tuple(allnames - {nm}),
                                log=lambda *a: None, **GU.SEARCH[kind])
        info['searches'].append(dict(shape_name=nm, stats=st))
        if f is not None:
            return dict(f, kind=kind), info
    return None, info


def owens_test(Q, D, sigma, s_cw, n_even):
    """Owens's test for the representative with sigma > 0, in both orientations of Y = boundary of the plumbing."""
    Q = np.array(Q, dtype=np.int64)
    mq, cm = OO.m_Q(-np.array(Q, dtype=object))
    dY = {g: -v for g, v in mq.items()}
    forms = []
    for f in OO.candidates(D, n_even):
        Qt = OO.qtilde(*f)
        assert OO.det_int(Qt) == D
        forms.append((f, OO.m_Q(Qt)))
    res = {}
    for s in (1, -1):
        mG = {g: s * v for g, v in dY.items()}
        rs = [(f, OO.matching_exists(mQt, cmQt, mG, cm)) for f, (mQt, cmQt) in forms]
        res[s] = dict(passing=[f for f, r in rs if r is True], undecided=[f for f, r in rs if r is None])
    s_exp = (1 if sigma > 0 else -1) * s_cw if s_cw else None
    blocked = {s: not res[s]['passing'] and not res[s]['undecided'] for s in (1, -1)}
    if blocked[1] and blocked[-1]:
        verdict = 'OBSTRUCTED'
    elif s_exp is not None and blocked[s_exp]:
        verdict = 'OBSTRUCTED (expected orientation only)'
    elif s_exp is not None and res[s_exp]['passing']:
        verdict = 'PASS'
    else:
        verdict = 'UNDECIDED'
    return dict(verdict=verdict, forms=len(forms), group=list(cm.d), s_expected=s_exp, d_spin=str(dY[tuple(0 for _ in cm.d)]),
                passing_plus=res[1]['passing'][:4], passing_minus=res[-1]['passing'][:4],
                undecided=len(res[1]['undecided']) + len(res[-1]['undecided']))


def run_one(n, Q, sigma, D):
    """Apply the |sigma|-dependent test given a plumbing matrix Q for (+-)Sigma_2(K)."""
    L = L_closed(Q); R = GU.R_of(n)
    cw = [s for s in (1, -1) if s * L == R]
    s_cw = cw[0] if len(cw) == 1 else None
    rec = dict(L=str(L), R=str(R), cw_signs=cw)
    if abs(sigma) == 4:
        rec.update(owens_test(Q, D, sigma, s_cw, 2), pattern='(0,2) [n=2]')
    elif abs(sigma) == 2:
        J = kinfo.parse_jones(n)
        excl, eps, d3 = traczyk_excludes_00_2(J if sigma > 0 else mirror_jones(J))
        o = owens_test(Q, D, sigma, s_cw, 1)
        rec.update(o, pattern='(1,1) [n=1]', traczyk_excludes_02=excl, traczyk_eps=eps, dim3=d3)
        if d3 > 2:                                   # as in v2's obstruct_u2: Wendt, u >= dim H1(Sigma_2; F3)
            rec['verdict'] = 'OBSTRUCTED (Wendt)'
        elif not excl and o['verdict'].startswith('OBSTRUCTED'):
            rec['verdict'] = '(1,1) OBSTRUCTED, (0,2) open'
    else:
        rec['verdict'] = 'NOT APPLICABLE'
    return rec


def control_one(nm):
    r = BY[nm]; D = int(r['determinant']); sigma = int(r['signature'])
    t0 = time.time()
    try:
        Q, _ = MU.plumbing_for(r['montesinos_notation'].strip(), D)
        if Q is None:
            return nm, None
        rec = run_one(nm, Q, sigma, D)
        rec.update(notation=r['montesinos_notation'], alternating=r['alternating'], D=D, sigma=sigma, u=v2T.get(nm))
        if r['alternating'] == 'Y':
            jones = kinfo.parse_jones(nm)
            v = OO.obstruct_u2(json.loads(r['pd_notation']), sigma, D, jones=jones)
            rec['v2_goeritz_verdict'] = v['verdict']
    except Exception as e:
        return nm, dict(error=repr(e)[:120])
    rec['seconds'] = round(time.time() - t0, 1)
    return nm, rec


def target_one(n):
    t0 = time.time()
    D = int(BY[n]['determinant']); sigma = int(BY[n]['signature'])
    try:
        f, info = find_plumbing_for_knot(n)
        if f is None:
            return n, dict(verdict='NO PLUMBING', D=D, sigma=sigma, info=info, seconds=round(time.time() - t0))
        Q = P.plumbing_matrix(f['weights'], f['edges'])
        rec = run_one(n, Q, sigma, D)
        deg = [0] * len(f['weights'])
        for u, v in f['edges']:
            deg[u] += 1; deg[v] += 1
        rec.update(D=D, sigma=sigma, plumbing=f, bad_vertices=sum(1 for i, w in enumerate(f['weights']) if deg[i] > -w),
                   names=info, u_before=v2T.get(n))
    except Exception as e:
        return n, dict(verdict='ERROR', error=repr(e)[:160], D=D, sigma=sigma)
    rec['seconds'] = round(time.time() - t0)
    return n, rec


if __name__ == '__main__':
    from multiprocessing import Pool
    mode = sys.argv[1]; workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    if mode == 'controls':
        items = []
        for nm, r in BY.items():
            m = (r.get('montesinos_notation') or '').strip(); c = r.get('crossing_number', '')
            u = v2T.get(nm)
            if (m.startswith('K(') and c.isdigit() and int(c) <= 12 and abs(int(r['signature'] or 0)) in (2, 4)
                    and u and u[0] == u[1] and u[0] in (2, 3) and int(r['determinant']) <= 300):
                items.append(nm)
        fn, outf = control_one, 'owens_u3_controls.json'
    else:
        rows = json.load(open('reach_u23_orb.json'))
        items = [x['name'] for x in rows if x['kind'] in ('star', 'pair', 'triple') and abs(x['sigma']) in (2, 4)]
        fn, outf = target_one, 'owens_u3_targets.json'
    out = json.load(open(outf)) if os.path.exists(outf) else {}
    todo = [x for x in items if x not in out]
    print(mode, 'items', len(items), '| to do', len(todo), flush=True)
    T0 = time.time()
    with Pool(workers) as pool:
        for n, rec in pool.imap_unordered(fn, todo):
            out[n] = rec
            json.dump(out, open(outf, 'w'), indent=1, default=str)
            if rec is None:
                continue
            print(f"[{time.time()-T0:5.0f}s] {n:9s} sigma {rec.get('sigma')} D {rec.get('D')} u {rec.get('u') or rec.get('u_before')} | "
                  f"{rec.get('verdict')} | v2 {rec.get('v2_goeritz_verdict', '-')} | CW {rec.get('cw_signs')} | forms {rec.get('forms')} "
                  f"| pass+ {rec.get('passing_plus')} pass- {rec.get('passing_minus')} | {rec.get('error', '')}", flush=True)
    from collections import Counter
    print('verdicts', dict(Counter(r['verdict'] for r in out.values() if r and 'verdict' in r)))
