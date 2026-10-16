"""Sigma_2 classes (orbifold route, reach_orb.classify) for the open non-alternating knots with u in [2,3] or [2,4]:
the knots for which graph-manifold correction terms could feed Owens's u = 2 obstruction.  Output: reach_u23_orb.json"""
import json, sys, time
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0, '../lib')
from ranges import load_ranges
import reach_orb as RO

if __name__ == '__main__':
    R = load_ranges('baseline')          # open ranges before data/new_values_20261015.json
    tg = sorted(n for n, r in R.items() if r['lo'] == 2 and r['hi'] in (3, 4) and RO.BY[n]['alternating'] == 'N')
    print('targets', len(tg), Counter((R[n]['lo'], R[n]['hi']) for n in tg), flush=True)
    t0 = time.time()
    with Pool(8) as p:
        rows = p.map(RO.classify, tg, chunksize=4)
    print('Sigma_2 classes:', dict(Counter(k for _, k, _, _ in rows)), f'{time.time()-t0:.0f}s', flush=True)
    json.dump([dict(name=n, kind=k, sigma2=nm, H1=h, range=[R[n]['lo'], R[n]['hi']], sigma=int(RO.BY[n]['signature']))
               for n, k, nm, h in rows], open('reach_u23_orb.json', 'w'), indent=1)
    for kind in ('star', 'pair', 'triple'):
        ks = [(n, int(RO.BY[n]['signature'])) for n, k, _, _ in rows if k == kind]
        print(kind, len(ks), '| by |sigma|:', dict(Counter(abs(s) for _, s in ks)))
