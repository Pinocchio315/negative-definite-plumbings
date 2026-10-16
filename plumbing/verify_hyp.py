#!/usr/bin/env sage -python
"""Certify hyperbolicity of Sigma_2 (orbifold route) with SnapPy's verified computations inside Sage (interval
arithmetic: a verified positively oriented solution of the gluing and Dehn filling equations).  A certified hyperbolic
Sigma_2 is not a graph manifold, so plumbing (lattice) methods cannot apply to it.
Run with Sage's Python (the venv's SnapPy is appended to sys.path; Sage provides the interval arithmetic):
    sage -python verify_hyp.py [workers]          Output: verify_hyp.json
Targets: the 'hyperbolic' (451) and 'unrecognised' (193) knots of recognise_more.json.
"""
import sys
import os
if os.environ.get('VENV_SITE_PACKAGES'):      # the venv with snappy, regina, database_knotinfo
    sys.path.append(os.environ['VENV_SITE_PACKAGES'])
import json, time, warnings
warnings.filterwarnings('ignore')
from collections import Counter


def one(n):
    import sage.all  # noqa: F401  (interval arithmetic for SnapPy's verify module)
    import snappy
    from database_knotinfo import link_list
    global BY
    try:
        BY
    except NameError:
        BY = {r['name']: r for r in link_list()[1:]}
    # The filled cover itself rarely carries a geometric solution; drilling one of its short closed geodesics and
    # refilling the new cusp along that geodesic's meridian gives another triangulation of the same Sigma_2 (checked by
    # H_1), on which SnapPy usually finds a positively oriented solution that verify_hyperbolicity can certify.
    E = snappy.Link(json.loads(BY[n]['pd_notation'])).exterior()
    rec = dict(certified=False)
    for t in range(4):
        E2 = E.copy()
        if t:
            E2.randomize()
        E2.dehn_fill((2, 0))
        C = E2.covers(2, cover_type='cyclic')[0]
        H = str(C.homology())
        for k in range(6):
            try:
                F = C.drill(k).filled_triangulation([0])
            except Exception:
                continue
            F.dehn_fill((1, 0))
            for s in range(8):
                if 'positively' in F.solution_type():
                    break
                F.randomize()
            if 'positively' not in F.solution_type() or str(F.homology()) != H:
                continue
            try:
                ok = F.verify_hyperbolicity()
                ok = bool(ok[0] if isinstance(ok, tuple) else ok)
            except Exception as e:
                rec['error'] = repr(e)[:120]; ok = False
            if ok:
                rec.update(certified=True, tries=t + 1, drilled=k, volume=float(F.volume()), H1=H)
                return n, rec
    return n, rec


if __name__ == '__main__':
    from multiprocessing import Pool
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    M = json.load(open('recognise_more.json'))
    targets = sorted(n for n, r in M.items() if r['new'] in ('still unrecognised',))
    olds = {n: M[n]['old'] for n in targets}
    t0 = time.time(); out = {}
    with Pool(workers) as p:
        for i, (n, rec) in enumerate(p.imap_unordered(one, targets, chunksize=4)):
            rec['old'] = olds[n]; out[n] = rec
            if (i + 1) % 50 == 0:
                print(f'[{time.time()-t0:5.0f}s] {i+1}/{len(targets)}', dict(Counter((r['old'], r['certified']) for r in out.values())), flush=True)
                json.dump(out, open('verify_hyp.json', 'w'), indent=1)
    json.dump(out, open('verify_hyp.json', 'w'), indent=1)
    print('done', len(out), dict(Counter((r['old'], r['certified']) for r in out.values())), f'{time.time()-t0:.0f}s')
