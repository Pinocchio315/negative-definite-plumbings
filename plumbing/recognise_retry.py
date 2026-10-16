#!/usr/bin/env python
"""Second recognition pass for the Sigma_2 still unrecognised after recognise_more.py (60 tries): up to 300 randomised
triangulations of the orbifold-route cover, each simplified by Regina's intelligentSimplify and, every tenth try,
additionally by simplifyExhaustive(1) when the triangulation is small (<= 30 tetrahedra); also the mirror diagram.
Usage: recognise_retry.py [workers]     Output: recognise_retry.json
"""
import json, sys, time, warnings
warnings.filterwarnings('ignore')
from collections import Counter
from multiprocessing import Pool
import general as G
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}


def tries(pd, k):
    import snappy, regina
    E = snappy.Link(pd).exterior()
    names = Counter()
    for t in range(k):
        E2 = E.copy()
        for _ in range(1 + t % 3):
            E2.randomize()
        E2.dehn_fill((2, 0))
        C = E2.covers(2, cover_type='cyclic')[0]
        T = regina.Triangulation3(C.filled_triangulation()._to_string()); T.intelligentSimplify()
        if t % 10 == 9 and T.size() <= 30:
            T.simplifyExhaustive(1)
        r = regina.StandardTriangulation.recognise(T)
        if r and r.manifold():
            names[r.manifold().name()] += 1
            if sum(names.values()) >= 3:
                break
    return names


def one(n):
    import spherogram
    pd = json.loads(BY[n]['pd_notation'])
    names = tries(pd, 150)
    if not names:
        mpd = [list(x) for x in spherogram.Link(pd).mirror().PD_code()]
        names = tries(mpd, 150)
    kinds = Counter(G.parse_name(x)[0] or 'other' for x in names)
    return n, dict(names=dict(names), kinds=dict(kinds))


if __name__ == '__main__':
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    M = json.load(open('recognise_more.json'))
    targets = sorted(n for n, r in M.items() if r['new'] == 'still unrecognised' and r['old'] == 'unrecognised')
    t0 = time.time(); out = {}
    with Pool(workers) as p:
        for n, rec in p.imap_unordered(one, targets):
            out[n] = rec
            if rec['names']:
                print(f'[{time.time()-t0:5.0f}s] {n}: {rec["kinds"]} {list(rec["names"])[:1]}', flush=True)
    json.dump(out, open('recognise_retry.json', 'w'), indent=1)
    print('targets', len(targets), '| recognised now:', sum(1 for r in out.values() if r['names']),
          dict(Counter(k for r in out.values() for k in r['kinds'])), f'| {time.time()-t0:.0f}s')
