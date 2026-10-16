#!/usr/bin/env python
"""star_u1.py on the knots that the orbifold-route classification (reach_u12_orb.json) recognises as stars but
reach.py did not (Regina's recognition depends on the triangulation).  Same pipeline: star plumbing search accepted
on the Regina name, v2's correction terms, graded-root HF_red, Ni-Wu u = 1 test (both orientations, all affine
identifications), then the Casson-Walker identity of crosscheck.py.
Output: star_u1_orb.json
"""
import gzip, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'montesinos'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
import _compat  # noqa
import numpy as np
import montesinos_u1 as MU
from hf_red import hf_red_all
import general as G
import plumbing as P
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
done = json.load(open('star_u1.json'))
rows = [r for r in json.load(open('reach_u12_orb.json')) if r['kind'] == 'star' and r['name'] not in done]
greene = {}
for line in gzip.open(os.path.join(ROOT, 'data', 'greene_sweep_2026-09-08.jsonl.gz'), 'rt'):
    r = json.loads(line)
    greene[r['name']] = (r.get('verdict'), r.get('unresolved_classes'), r.get('test'))


def R_of(name):
    import kinfo
    from fractions import Fraction
    J = kinfo.parse_jones(name); s = int(BY[name]['signature']); D = int(BY[name]['determinant'])
    V1 = sum(c * (-1) ** (e % 2) for e, c in J.items())
    dV1 = sum(e * c * (-1) ** ((e - 1) % 2) for e, c in J.items())
    assert abs(V1) == D
    return Fraction(-D * s, 4) + Fraction(D * dV1, 6 * V1)


out = {}
T0 = time.time()
for row in rows:
    n, nm = row['name'], row['sigma2']
    D = int(BY[n]['determinant'])
    found, st = G.find_plumbing(nm, D, json.loads(BY[n]['seifert_matrix']), erange=range(-8, 2), max_chain=0)
    rec = dict(sigma2=nm, det=D, in_v2_greene=greene.get(n), search=st)
    if found is None:
        rec['verdict'] = 'NO PLUMBING'
    else:
        Q = np.array(P.plumbing_matrix(found['weights'], found['edges']), dtype=np.int64)
        d, cm = MU.correction_terms(Q)
        adm_d = MU.u1_admissible(d, cm, D, Q)
        ranks, _ = hf_red_all(Q)
        if adm_d is None:
            verdict = 'OBSTRUCTED (H1 not cyclic)'
        elif not adm_d:
            verdict = 'OBSTRUCTED'
        elif ranks is None:
            verdict = 'PASS (ranks not computed)'
        else:
            verdict = 'PASS' if MU.u1_admissible(d, cm, D, Q, ranks_by_class=ranks) else 'OBSTRUCTED_HF'
        L = sum(d.values()) + 2 * sum(ranks.values()) if ranks is not None else None
        R = R_of(n)
        s_link = 1 if found['orientation'] == '+' else -1
        cw = [s for s in (1, -1) if L is not None and s * L == R]
        rec.update(plumbing=found, verdict=verdict, hf_red_ranks=sorted(ranks.values()) if ranks else None,
                   L=str(L), R=str(R), cw_signs=cw, link_sign=s_link,
                   cw_consistent=bool(cw) and (s_link in cw or found['orientation_ambiguous']))
    out[n] = rec
    print(f"[{time.time()-T0:5.0f}s] {n:9s} det {D:4d} | {nm[:48]:48s} | v2 Greene: {rec['in_v2_greene']} | u=1 test: {rec['verdict']}"
          f" | CW L {rec.get('L')} R {rec.get('R')} consistent {rec.get('cw_consistent')}", flush=True)
json.dump(out, open('star_u1_orb.json', 'w'), indent=1, default=str)
print('OBSTRUCTED (new u = 2 candidates):', [n for n, r in out.items() if r['verdict'].startswith('OBSTRUCTED')])
