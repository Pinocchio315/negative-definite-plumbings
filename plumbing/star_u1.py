#!/usr/bin/env python
"""u = 1 obstruction for open [1,2] knots whose Sigma_2 is a Seifert fibred space over S^2 (Regina-recognised).

Plumbing: star found by search (centre weight, legs from both normalisations of each exceptional fibre), accepted
when Regina recognises its boundary with the same (canonical) SFS name; negative definite; |det| = det K.
Then v2's validated Montesinos machinery on that star (at most one bad vertex):
  correction terms (sharp formula), reduced Floer ranks (graded root), Ni-Wu half-integral surgery test with the
  reduced-homology mapping cone (u1_admissible), both orientations and all affine identifications.
Records whether the knot was already in v2's Greene sweep.
Output: star_u1.json
"""
import gzip, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'montesinos'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
import _compat  # noqa
import montesinos_u1 as MU
from hf_red import hf_red_all
import general as G
import plumbing as P
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
rows = json.load(open('reach_u12.json'))
stars = [(n, nm) for n, k, nm in rows if k == 'star']
greene = {}
for line in gzip.open(os.path.join(ROOT, 'data', 'greene_sweep_2026-09-08.jsonl.gz'), 'rt'):
    r = json.loads(line)
    greene[r['name']] = (r.get('verdict'), r.get('unresolved_classes'), r.get('test'))
out = {}
T0 = time.time()
for n, nm in stars:
    D = int(BY[n]['determinant'])
    found, st = G.find_plumbing(nm, D, json.loads(BY[n]['seifert_matrix']), erange=range(-8, 2), max_chain=0)
    rec = dict(sigma2=nm, det=D, in_v2_greene=greene.get(n), search=st)
    if found is None:
        rec['verdict'] = 'NO PLUMBING'
    else:
        import numpy as np
        Q = np.array(P.plumbing_matrix(found['weights'], found['edges']), dtype=np.int64)
        d, cm = MU.correction_terms(Q)
        adm_d = MU.u1_admissible(d, cm, D, Q)
        ranks, _ = hf_red_all(Q) if adm_d else (None, None)
        if adm_d and ranks is not None:
            adm = MU.u1_admissible(d, cm, D, Q, ranks_by_class=ranks)
            verdict = 'PASS' if adm else 'OBSTRUCTED_HF'
        elif adm_d is None:
            verdict = 'OBSTRUCTED (H1 not cyclic)'
        else:
            verdict = 'OBSTRUCTED' if not adm_d else 'PASS (ranks not computed)'
        rec.update(plumbing=found, verdict=verdict, hf_red_ranks=sorted(ranks.values()) if ranks else None)
    out[n] = rec
    print(f"[{time.time()-T0:5.0f}s] {n:9s} det {D:4d} | {nm[:60]:60s} | v2 Greene: {rec['in_v2_greene']} | u=1 test: {rec['verdict']}", flush=True)
json.dump(out, open('star_u1.json', 'w'), indent=1, default=str)
print('OBSTRUCTED (new u = 2 candidates):', [n for n, r in out.items() if r['verdict'].startswith('OBSTRUCTED')])
