#!/usr/bin/env python
"""Recognition by exploring the Pachner graph (Regina retriangulate) for the 35 graph-manifold candidates and the one
atoroidal unknown of jsj.json: every triangulation reached within the given height (from the smallest orbifold-route
triangulation of Sigma_2, and of the mirror) is offered to StandardTriangulation.recognise; stop at the first name.
Recognised knots then go through neumann/run_fresh.one (Neumann normal form, p+ test, u = 1 or Owens test).
Usage: pachner_recognise.py [workers]        Output: pachner_recognise.json
"""
import json, os, sys, time, warnings
warnings.filterwarnings('ignore')
from collections import Counter
from multiprocessing import Pool
import jsj as JS


def explore(T, height, limit):
    import regina
    found = {}
    count = [0]
    def action(sig, tri):
        count[0] += 1
        r = regina.StandardTriangulation.recognise(tri)
        if r and r.manifold():
            found[r.manifold().name()] = 1
            return True
        return count[0] >= limit
    T.retriangulate(height, 1, action)
    return found, count[0]


def one(n):
    t0 = time.time()
    T = JS.smallest_sigma2(n, 16)
    names, visited = explore(T, 2, 60000)
    if not names:
        names, v3 = explore(T, 3, 150000); visited += v3
    return n, dict(names=sorted(names), visited=visited, seconds=round(time.time() - t0, 1))


if __name__ == '__main__':
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    J = json.load(open('jsj.json'))
    targets = sorted(n for n, r in J.items() if not r['verdict'].startswith('toroidal, a piece'))
    out = {}
    t0 = time.time()
    with Pool(workers) as p:
        for n, rec in p.imap_unordered(one, targets):
            out[n] = rec
            print(f"[{time.time()-t0:5.0f}s] {n:9s} {rec['names'][:1]} visited {rec['visited']} {rec['seconds']}s", flush=True)
            json.dump(out, open('pachner_recognise.json', 'w'), indent=1)
    import general as G
    print('recognised', sum(1 for r in out.values() if r['names']), '/', len(out),
          dict(Counter(G.parse_name(r['names'][0])[0] or 'other' for r in out.values() if r['names'])))
