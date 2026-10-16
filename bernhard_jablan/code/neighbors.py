"""Step 2: BJSet = knots obtained by one crossing change in any minimal diagram (BH Prop 2.3)."""
import json, sys, re, gzip, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
from collections import defaultdict
from database_knotinfo import link_list
from enum_diagrams import change, by

def norm(s):
    s = s.split('(')[0]
    m = re.match(r'K(\d+)([an])(\d+)$', s)
    if m: return f'{m.group(1)}{m.group(2)}_{m.group(3)}'
    return s  # Rolfsen name like 7_4

def ident(pd):
    import snappy, spherogram
    L = spherogram.Link(pd)
    L.simplify('global')
    if len(L.crossings) == 0:
        return ['unknot'], 'simplify'
    try:
        ids = [norm(str(x)) for x in snappy.Link(pd).exterior().identify()]
    except Exception:
        ids = []
    if ids:
        # prefer the KnotInfo-style name
        ki = [i for i in ids if i in by]
        return (ki or ids), 'identify'
    try:
        parts = L.deconnect_sum()
    except Exception:
        parts = [L]
    if len(parts) > 1:
        names = []
        for P in parts:
            P.simplify('global')
            try:
                pid = [norm(str(x)) for x in P.exterior().identify()]
            except Exception:
                pid = []
            names.append(([i for i in pid if i in by] or pid or [f'?{len(P.crossings)}c'])[0])
        return ['#'.join(sorted(names))], 'deconnect'
    return [f'NONHYP?{len(L.crossings)}c'], 'unresolved'

def work(d):
    pd0 = json.loads(by[d['ref']]['pd_notation'])
    S = set(d['subset'])
    out = []
    for c in range(len(pd0)):
        T = sorted(S ^ {c})
        pd = change(pd0, T)
        ids, how = ident(pd)
        out.append(dict(ref=d['ref'], subset=d['subset'], crossing=c, ids=ids, how=how, pd=pd))
    return out

if __name__ == '__main__':
    D = json.load(open(sys.argv[1]))
    rows = []
    with Pool(11) as pool:
        for out in pool.imap_unordered(work, D['diagrams']):
            rows.extend(out)
    bj = defaultdict(int)
    for r in rows: bj[r['ids'][0]] += 1
    json.dump(dict(target=D['target'], n_diagrams=len(D['diagrams']), changes=rows), open(sys.argv[1].replace('diagrams_', 'bjset_'), 'w'), indent=0)
    print(D['target'], 'diagrams', len(D['diagrams']), 'single changes', len(rows), 'distinct knots', len(bj))
    for k, v in sorted(bj.items()): print('  ', k, v)
