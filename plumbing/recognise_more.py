#!/usr/bin/env python
"""How much would better recognition of Sigma_2 gain?  For the knots whose Sigma_2 was 'unrecognised' or 'hyperbolic' in
reach_u12_orb.json / reach_u23_orb.json:
  (1) SnapPy verified hyperbolicity (interval arithmetic) of the orbifold-route Sigma_2: certified hyperbolic => not a
      graph manifold, out of scope for plumbing methods;
  (2) for the rest, up to 60 more randomised triangulations through Regina's recognition.
Also: Conway notation (KnotInfo, <= 12 crossings) - algebraic (no '*', '.', ':') knots have graph-manifold Sigma_2.
Usage: recognise_more.py [workers]      Output: recognise_more.json
"""
import json, sys, time, warnings
warnings.filterwarnings('ignore')
from collections import Counter
from multiprocessing import Pool
import general as G
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}


def one(args):
    n, old = args
    import snappy, regina
    pd = json.loads(BY[n]['pd_notation'])
    E = snappy.Link(pd).exterior()
    rec = dict(old=old)
    verified = False
    names = Counter()
    for t in range(60):
        E2 = E.copy()
        if t:
            E2.randomize()
        E2.dehn_fill((2, 0))
        C = E2.covers(2, cover_type='cyclic')[0]
        if not verified and t < 6:
            try:
                ok = C.verify_hyperbolicity()
                verified = bool(ok[0] if isinstance(ok, tuple) else ok)
            except Exception:
                pass
        if verified:
            break
        T = regina.Triangulation3(C.filled_triangulation()._to_string()); T.intelligentSimplify()
        r = regina.StandardTriangulation.recognise(T)
        if r and r.manifold():
            names[r.manifold().name()] += 1
            if t >= 8:
                break
    kinds = Counter(G.parse_name(x)[0] or 'other' for x in names)
    rec.update(verified_hyperbolic=verified, names=dict(names), kinds=dict(kinds))
    if verified:
        rec['new'] = 'hyperbolic (verified)'
    elif names:
        rec['new'] = sorted(kinds)[0] if kinds else 'other'
    else:
        rec['new'] = 'still unrecognised'
    v = (BY[n].get('conway_notation') or '').strip()
    rec['conway'] = v
    rec['algebraic_notation'] = bool(v) and not any(ch in v for ch in '*.:')
    return n, rec


if __name__ == '__main__':
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    items = []
    for f in ('reach_u12_orb.json', 'reach_u23_orb.json'):
        for x in json.load(open(f)):
            if x['kind'] in ('unrecognised', 'hyperbolic'):
                items.append((x['name'], x['kind']))
    t0 = time.time(); out = {}
    with Pool(workers) as p:
        for n, rec in p.imap_unordered(one, items, chunksize=2):
            out[n] = rec
    json.dump(out, open('recognise_more.json', 'w'), indent=1)
    print(f'{len(out)} knots, {time.time()-t0:.0f}s')
    print('old -> new:', sorted(Counter((r['old'], r['new']) for r in out.values()).items()))
    print('algebraic Conway notation among them:', sorted(Counter((r['old'], r['new']) for r in out.values() if r['algebraic_notation']).items()))
    print('newly recognised graph-type:', sorted(n for n, r in out.items() if r['new'] in ('star', 'pair', 'triple')))
