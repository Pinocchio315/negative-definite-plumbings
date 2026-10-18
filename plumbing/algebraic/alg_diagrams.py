#!/usr/bin/env python
"""Algebraic diagrams of the 32 knots whose double branched covers Regina does not recognise. Random Reidemeister moves (spherogram backtrack), then simplification
('level', 'basic', 'global'); a diagram is algebraic when bigon contraction reduces it to one vertex (1*).
The moves preserve the knot type, so every diagram found is a diagram of the same knot. Spherogram's moves use
Python's global random generator, which is seeded per knot, but the search is still not exactly reproducible, since
spherogram also iterates over sets of objects. Every algebraic diagram gives the same normal forms (uniqueness of the
Neumann normal form); a rerun found other diagrams for 30 of the 32 knots, with the same normal forms.
Input: targets.txt, the 32 knots whose double branched covers Regina does not recognise (toroidal, no piece certified
hyperbolic).  Output: alg_diagrams.json  {name: {pd, crossings, iteration, mode, seed}}
"""
import json, random, sys, time, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
from basic_polyhedron import reduce_diagram
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
ITERS = 6000


def search(n, seed=1, avoid=None):
    """An algebraic diagram of n with the fewest crossings found; diagrams for which avoid(pd) holds are skipped."""
    import spherogram
    random.seed(seed)
    K0 = spherogram.Link(json.loads(BY[n]['pd_notation']))
    found = []
    for it in range(ITERS):
        K = K0.copy()
        K.backtrack(random.randint(5, 40))
        for mode in ('level', 'basic', 'global'):
            L = K.copy()
            try:
                L.simplify(mode)
            except Exception:
                continue
            pd = [list(x) for x in L.PD_code()]
            if len(L.link_components) == 1 and reduce_diagram(pd) == 1 and not (avoid and avoid(pd)):
                found.append(dict(pd=pd, crossings=len(pd), iteration=it, mode=mode, seed=seed))
        if found and it >= 200:
            break
    if not found:
        return n, None
    return n, min(found, key=lambda d: (d['crossings'], d['iteration']))


if __name__ == '__main__':
    names = sys.argv[1:] or [l.strip() for l in open('targets.txt') if l.strip()]
    out, t0 = {}, time.time()
    with Pool(8) as p:
        for n, rec in p.imap_unordered(search, names):
            out[n] = rec
            print(f"[{time.time()-t0:5.0f}s] {n}", 'none' if rec is None else f"{rec['crossings']} crossings (it {rec['iteration']}, {rec['mode']})", flush=True)
    json.dump(dict(sorted(out.items())), open('alg_diagrams.json', 'w'), indent=1)
    print('found', sum(1 for r in out.values() if r), 'of', len(out))
