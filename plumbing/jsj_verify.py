#!/usr/bin/env sage -python
"""Certify the hyperbolic pieces found by jsj.py (Sage + SnapPy verify_hyperbolicity on the ideal triangulation of a piece
cut off by an incompressible torus).  If such a piece is hyperbolic, Sigma_2 has positive simplicial volume (Gromov:
additive under gluing along incompressible tori), so Sigma_2 is not a graph manifold and plumbing methods cannot apply.
    sage -python jsj_verify.py [workers]          Output: jsj_verify.json
"""
import sys
import os
if os.environ.get('VENV_SITE_PACKAGES'):      # the venv with snappy, regina, database_knotinfo
    sys.path.append(os.environ['VENV_SITE_PACKAGES'])
import json, time, warnings
warnings.filterwarnings('ignore')
from collections import Counter


def one(args):
    n, isosigs = args
    import sage.all  # noqa: F401
    import snappy
    for sig in isosigs:
        try:
            M = snappy.Manifold(sig)
        except Exception:
            continue
        for t in range(20):
            if 'positively' in M.solution_type():
                break
            M.randomize()
        if 'positively' not in M.solution_type():
            continue
        try:
            ok = M.verify_hyperbolicity()
            ok = bool(ok[0] if isinstance(ok, tuple) else ok)
        except Exception:
            ok = False
        if ok:
            return n, dict(certified_hyperbolic_piece=True, isosig=sig, volume=float(M.volume()), cusps=M.num_cusps())
    return n, dict(certified_hyperbolic_piece=False)


if __name__ == '__main__':
    from multiprocessing import Pool
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    J = json.load(open('jsj.json'))
    items = []
    for n, r in J.items():
        sigs = [p['ideal_isosig'] for f in r['incompressible_tori'] for p in f['pieces'] if 'positively' in (p.get('snappy_solution') or '')]
        if sigs:
            items.append((n, sigs))
    with Pool(workers) as p:
        out = dict(p.map(one, items))
    json.dump(out, open('jsj_verify.json', 'w'), indent=1)
    print('knots with a candidate hyperbolic piece:', len(items), '| certified:', Counter(r['certified_hyperbolic_piece'] for r in out.values()))
    v = Counter(round(r['volume'], 3) for r in out.values() if r['certified_hyperbolic_piece'])
    print('piece volumes (most common):', v.most_common(8))
