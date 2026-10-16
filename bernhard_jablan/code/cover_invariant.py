"""Exact separation of candidates that share Alexander polynomial and HFK with the target.
Invariant: multiset of first homology groups of all connected covers of degree <= DMAX of the
knot complement (permutation representations of pi_1 + Smith normal form; purely combinatorial).
Different multisets => non-homeomorphic complements => different knots (mirror included)."""
import json, sys, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
import snappy
from database_knotinfo import link_list
sys.path.insert(0, __import__('os').path.dirname(__file__))
by = {r['name']: r for r in link_list()[1:]}
DMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 5

def change(pd, idxs):
    n2 = 2 * len(pd); new = [list(r) for r in pd]
    for i in idxs:
        a, b, c, d = pd[i]
        assert (c - a) % n2 == 1
        new[i] = [b, c, d, a] if (d - b) % n2 == 1 else [d, a, b, c]
    return new

def inv(pd):
    M = snappy.Link(pd).exterior()
    out = {}
    for d in range(2, DMAX + 1):
        out[d] = sorted(str(C.homology()) for C in M.covers(d))
    return out

def work(job):
    ref, S = job
    pd = change(json.loads(by[ref]['pd_notation']), S)
    return ref, S, inv(pd)

if __name__ == '__main__':
    target = sys.argv[1]; comp = json.load(open(sys.argv[2]))
    T = inv(json.loads(by[target]['pd_notation']))
    jobs = [(r['ref'], c['subset']) for r in comp['references'] for c in r['candidates']
            if c['decision'].startswith('UNRESOLVED') or c['decision'].startswith('ACCEPTED')]
    with Pool(11) as p:
        res = p.map(work, jobs)
    dec = {}
    rows = []
    for (ref, S, I), (_, _) in zip(res, jobs):
        same = (I == T)
        c = next(c for r in comp['references'] if r['ref'] == ref for c in r['candidates'] if c['subset'] == S)
        status = c['decision'].split(':')[0]
        first_diff = None if same else next(d for d in T if I[d] != T[d])
        rows.append(dict(ref=ref, subset=S, prior=status, snappy=c.get('snappy_identify', [])[:1], cover_invariant_equal=same, first_differing_degree=first_diff))
    from collections import Counter
    print('target cover homology counts by degree:', {d: len(v) for d, v in T.items()})
    print(Counter((r['prior'], r['snappy'][0] if r['snappy'] else '-', r['cover_invariant_equal'], r['first_differing_degree']) for r in rows))
    json.dump(dict(target=target, dmax=DMAX, target_invariant=T, rows=rows), open(sys.argv[2].replace('completeness_', 'cover_separation_'), 'w'), indent=0)
