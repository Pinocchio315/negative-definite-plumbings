#!/usr/bin/env python
"""Verification of the u = 1 obstructions in graph_u1.json (graph-manifold Sigma_2), and a closed-form Casson-Walker
check that needs no Floer homology.

Closed form.  For Y = boundary of a negative-definite plumbing tree of spheres (s vertices, weights e_v, degrees d_v,
intersection matrix Q), with L(Y) := -2 |H| lambda(Y) (Casson normalisation; L = sum d + 2 sum rank HF_red),
    L(Y) = (|H|/12) * ( sum_v e_v + 3 s + sum_v (2 - d_v) (Q^{-1})_{vv} )
(Lescop's plumbing formula in the normalisation fixed by E8: L = 2, L(3,1): L = -1/6, Sigma(2,3,7): L = 2).
Calibration here: it must reproduce L = sum d + 2 sum rank HF_red on the 19 star plumbings (crosscheck.json,
star_u1_orb.json) and on the 30 Montesinos controls (crosscheck_controls.json, v2 notation plumbings).
Then for every obstruction: s L(Y) = R(K) = -D sigma/4 + D V'(-1)/(6 V(-1)) with s the linking-form orientation;
independent Ni-Wu re-check of the same d (d = max (K^2+|V|)/4 via m_Q); novelty.
Output: graph_verify.json
"""
import json, os, sys, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'montesinos'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
sys.path.insert(0, os.path.join(ROOT, 'bernhard_jablan', 'code'))
import _compat  # noqa
import numpy as np
import sympy
import montesinos_u1 as MU
import u1_independent as U1
import plumbing as P
import d_owens as DO
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
k916 = json.load(open(os.path.join(ROOT, 'data', 'knotinfo_2026-09-16_unknotting.json')))
v2T = json.load(open(os.path.join(ROOT, 'data', 'v2_final_u_table.json')))
JAB = {'13n_' + x for x in '128 152 365 391 492 839 842 967 976 993 1082 1084 1146 1345 1478 1524 1846 2453 2456 2612 2681 2689 2729 3251 3263 3639 3729 3811 3853 3881 4029 4381 4623'.split()}


def L_closed(Q):
    """(|H|/12)(sum e + 3s + sum (2 - deg)(Q^-1)_vv) for a plumbing matrix Q (diagonal = weights, 1 on edges)."""
    S = sympy.Matrix(np.array(Q, dtype=np.int64).tolist())
    n = S.rows; H = abs(S.det()); Si = S.inv()
    deg = [sum(1 for j in range(n) if j != i and S[i, j] != 0) for i in range(n)]
    val = sum(S[i, i] for i in range(n)) + 3 * n + sum((2 - deg[i]) * Si[i, i] for i in range(n))
    return Fraction(int(sympy.fraction(H * val / 12)[0]), int(sympy.fraction(H * val / 12)[1]))


def R_of(name):
    import kinfo
    J = kinfo.parse_jones(name); s = int(BY[name]['signature']); D = int(BY[name]['determinant'])
    V1 = sum(c * (-1) ** (e % 2) for e, c in J.items())
    dV1 = sum(e * c * (-1) ** ((e - 1) % 2) for e, c in J.items())
    assert abs(V1) == D
    return Fraction(-D * s, 4) + Fraction(D * dV1, 6 * V1)


out = dict(calibration={}, verify={})
# calibration: stars (L from d + HF_red) and Montesinos controls
cc = json.load(open('crosscheck.json'))['star']; su = json.load(open('star_u1.json'))
so = json.load(open('star_u1_orb.json'))
cal = []
for n, r in su.items():
    if r.get('plumbing') and cc[n].get('L') not in (None, 'None'):
        cal.append((n, P.plumbing_matrix(r['plumbing']['weights'], r['plumbing']['edges']), Fraction(cc[n]['L'])))
for n, r in so.items():
    if r.get('plumbing') and r.get('L') not in (None, 'None'):
        cal.append((n, P.plumbing_matrix(r['plumbing']['weights'], r['plumbing']['edges']), Fraction(r['L'])))
for n, r in json.load(open('crosscheck_controls.json')).items():
    if r and r.get('L'):
        Q, _ = MU.plumbing_for(r['notation'].strip(), r['D'])
        cal.append((n, Q, Fraction(r['L'])))
for n, Q, L in cal:
    Lc = L_closed(Q)
    out['calibration'][n] = dict(L_floer=str(L), L_closed=str(Lc), equal=(Lc == L))
ok = sum(v['equal'] for v in out['calibration'].values())
print(f'calibration: closed form = sum d + 2 sum rank HF_red on {ok} / {len(cal)} plumbings', flush=True)
for n, v in out['calibration'].items():
    if not v['equal']:
        print('   MISMATCH', n, v)

# verification of the graph-manifold obstructions
G1 = json.load(open('graph_u1.json'))
for n, r in sorted(G1.items(), key=lambda x: (x[0].split('_')[0], int(x[0].split('_')[1]))):
    if not r['verdict'].startswith('OBSTRUCTED'):
        continue
    f = r['plumbing']; D = r['det']
    Q = P.plumbing_matrix(f['weights'], f['edges'])
    Lc = L_closed(Q); R = R_of(n)
    cw = [s for s in (1, -1) if s * Lc == R]
    s_link = 1 if f['orientation'] == '+' else -1
    Dp, dZ = DO.plumbing_d(f['weights'], f['edges'])
    surv = U1.ni_wu_survivors(D, dZ) if r['verdict'] == 'OBSTRUCTED' else None
    matched = f['regina_name_of_plumbing']
    rec = dict(kind=f['kind'], det=D, bad_vertices=r['bad_vertices'], vertices=len(f['weights']), name=matched,
               name_from=('Sigma_2' if matched in r['names_sigma2'] else '') + ('/mirror' if matched in r['names_mirror'] else ''),
               L_closed=str(Lc), R=str(R), cw_signs=cw, link_sign=s_link, link_ambiguous=f['orientation_ambiguous'],
               cw_consistent=bool(cw) and (s_link in cw or f['orientation_ambiguous']),
               independent_ni_wu_survivors=None if surv is None else len(surv),
               knotinfo_2026_09_16=k916.get(n), v2_table=v2T.get(n), in_jablonowski_36=n in JAB)
    out['verify'][n] = rec
    print(f"{n:9s} {rec['kind']:6s} D {D:3d} | bad {rec['bad_vertices']} | name from {rec['name_from']:15s} | L {str(Lc):>8s} R {str(R):>8s} "
          f"| CW {cw} link {s_link}{'(amb)' if rec['link_ambiguous'] else ''} consistent {rec['cw_consistent']} "
          f"| Ni-Wu survivors {rec['independent_ni_wu_survivors']} | KI {k916.get(n, {}).get('u')} v2 {v2T.get(n)} JAB {n in JAB}", flush=True)
json.dump(out, open('graph_verify.json', 'w'), indent=1, default=str)
