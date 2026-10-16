"""Re-classify Sigma_2 for the 437 knots with u in [1,2] by the orbifold route ((2,0) filling, then the 2-fold cyclic
cover), since the (1,0) filling of the unfilled cover used in reach.py is not always Sigma_2.
Regina recognition is tried first (several randomisations); otherwise hyperbolic if SnapPy finds a geometric
solution of volume > 0.5.  Records H1 (must be Z/det) and the reach.py class for comparison.
Output: reach_u12_orb.json
"""
import json, sys, time, warnings
warnings.filterwarnings('ignore')
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0, '../lib')
import general as G
from ranges import load_ranges
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}


def classify(n):
    import snappy, regina
    pd = json.loads(BY[n]['pd_notation'])
    E = snappy.Link(pd).exterior()
    h = None; hyp = False
    for t in range(4):
        E2 = E.copy()
        if t:
            E2.randomize()
        E2.dehn_fill((2, 0))
        C = E2.covers(2, cover_type='cyclic')[0]
        h = str(C.homology())
        T = regina.Triangulation3(C.filled_triangulation()._to_string()); T.intelligentSimplify()
        rec = regina.StandardTriangulation.recognise(T)
        if rec and rec.manifold():
            nm = rec.manifold().name()
            return n, G.parse_name(nm)[0] or 'other-recognised', nm, h
        st = C.solution_type()
        if ('positive' in st or 'negative' in st) and C.volume() > 0.5:
            hyp = True
    return n, ('hyperbolic' if hyp else 'unrecognised'), None, h


if __name__ == '__main__':
    R = load_ranges('baseline')          # open ranges before data/new_values_20261015.json
    targets = sorted(n for n, r in R.items() if r['lo'] == 1 and r['hi'] == 2)
    old = {n: (k, nm) for n, k, nm in json.load(open('reach_u12.json'))}
    t0 = time.time(); rows = []; c = Counter(); changed = []
    with Pool(6) as pool:
        for i, (n, kind, nm, h) in enumerate(pool.imap_unordered(classify, targets, chunksize=2)):
            D = int(BY[n]['determinant'])
            h_ok = h == (f'Z/{D}' if D > 1 else '0')
            rows.append(dict(name=n, kind=kind, sigma2=nm, H1=h, H1_ok=h_ok, old_kind=old[n][0], old_name=old[n][1]))
            c[kind] += 1
            if kind != old[n][0] or (nm and nm != old[n][1]):
                changed.append(n)
                print(f'[{time.time()-t0:4.0f}s] changed {n:9s} {old[n][0]:12s} -> {kind:12s} {nm} (H1 {h}, ok {h_ok})', flush=True)
            if (i + 1) % 50 == 0:
                print(f'[{time.time()-t0:4.0f}s] {i+1}/{len(targets)}', dict(c), flush=True)
    rows.sort(key=lambda r: r['name'])
    print('orbifold route:', dict(c), '| changed vs reach.py:', len(changed), '| H1 = Z/det in all:', all(r['H1_ok'] for r in rows))
    json.dump(rows, open('reach_u12_orb.json', 'w'), indent=1)
