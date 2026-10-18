#!/usr/bin/env python
"""Independent checks for the 32 knots of targets.txt.
  (1) a second algebraic diagram (random search with another seed) gives the same normal forms (canonical strings) in
      both orientations; Neumann's normal form is unique, so any difference would expose an error;
  (2) SnapPy: the exterior of each algebraic diagram is isometric to the exterior of the KnotInfo diagram;
  (3) numpy eigenvalues of the normal forms agree with the exact inertia.
Output: crosscheck_u.json
"""
import json, os, sys, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from validate_plumbing import setup
from alg_diagrams import search, BY
from algebraic_plumbing import plumbing_from_algebraic, normal_form, plumbing_matrix


def one(n):
    setup()
    import numpy as np, snappy, spherogram
    from both_orientations import canon
    U = json.load(open(os.path.join(HERE, 'u_plumbing.json')))
    AD = json.load(open(os.path.join(HERE, 'alg_diagrams.json')))
    _, second = search(n, seed=7)
    rec = dict(second_crossings=second['crossings'] if second else None)
    if second:
        same = {}
        for s in ('1', '-1'):
            wn, en = normal_form(*plumbing_from_algebraic(second['pd'], int(s)))
            same[s] = canon(wn, en) == U[n][s]['canon']
        rec['second_same_normal_forms'] = same
    E0 = spherogram.Link(json.loads(BY[n]['pd_notation'])).exterior()
    iso = []
    for d in [AD[n]['pd']] + ([second['pd']] if second else []):
        E1 = spherogram.Link(d).exterior()
        ok = None
        for t in range(8):
            try:
                ok = E0.is_isometric_to(E1); break
            except Exception:
                E1.randomize()
        iso.append(ok)
    rec['isometric_to_knotinfo'] = iso
    rec['knotinfo_exterior'] = dict(solution=E0.solution_type(), volume=float(E0.volume()))
    ev = {}
    for s in ('1', '-1'):
        Q = np.array(plumbing_matrix(U[n][s]['weights'], U[n][s]['edges']), dtype=float)
        x = np.linalg.eigvalsh(Q)
        ev[s] = [int((x > 1e-9).sum()), int((x < -1e-9).sum())] == [U[n][s]['inertia']['pplus'], U[n][s]['inertia']['pminus']]
    rec['numpy_inertia_agrees'] = ev
    return n, rec


if __name__ == '__main__':
    AD = json.load(open(os.path.join(HERE, 'alg_diagrams.json')))
    out = {}
    with Pool(8) as p:
        for n, r in p.imap_unordered(one, sorted(AD)):
            out[n] = r
            print(n, r, flush=True)
    json.dump(dict(sorted(out.items())), open(os.path.join(HERE, 'crosscheck_u.json'), 'w'), indent=1)
    print('second diagram found', sum(1 for r in out.values() if r['second_crossings']),
          '| same normal forms', sum(1 for r in out.values() if r.get('second_same_normal_forms') and all(r['second_same_normal_forms'].values())),
          '| isometric', sum(1 for r in out.values() if all(r['isometric_to_knotinfo'])),
          '| numpy inertia agrees', sum(1 for r in out.values() if all(r['numpy_inertia_agrees'].values())))
