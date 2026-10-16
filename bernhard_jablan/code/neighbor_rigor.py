"""Neighbour-side rigour for the 13n1587 BJSet.
(1) Cover-homology invariant (degree <= 4) of every prime neighbour diagram vs. its named knot.
(2) Name-free McCoy check: simplify each neighbour diagram; if it becomes a reduced alternating
    diagram, test every crossing for unknotting (HFK rank 1 <=> unknot). No unknotting crossing
    => u >= 2 by McCoy's theorem, without identifying the knot."""
import json, sys, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
import snappy, spherogram
from knot_floer_homology import pd_to_hfk
from database_knotinfo import link_list
by = {r['name']: r for r in link_list()[1:]}
DMAX = 4

def change(pd, idxs):
    n2 = 2 * len(pd); new = [list(r) for r in pd]
    for i in idxs:
        a, b, c, d = pd[i]
        assert (c - a) % n2 == 1
        new[i] = [b, c, d, a] if (d - b) % n2 == 1 else [d, a, b, c]
    return new

def cover_inv(pd):
    M = snappy.Link(pd).exterior()
    return {d: sorted(str(C.homology()) for C in M.covers(d)) for d in range(2, DMAX + 1)}

def is_alternating(pd):
    pos = {}
    for row in pd:
        for i, lab in enumerate(row):
            pos.setdefault(lab, []).append(i % 2)
    return all(sorted(v) == [0, 1] for v in pos.values())

def work(r):
    name = r['ids'][0]
    out = dict(ref=r['ref'], subset=r['subset'], crossing=r['crossing'], name=name)
    if name in by:
        out['cover_invariant_matches_named'] = cover_inv(r['pd']) == cover_inv(json.loads(by[name]['pd_notation']))
    L = spherogram.Link(r['pd']); L.simplify('global')
    spd = L.PD_code()
    if len(L.link_components) == 1 and len(spd) > 0 and is_alternating(spd):
        # PD from spherogram may not satisfy the consecutive-label convention; relabel via snappy round trip
        try:
            unk = []
            for c in range(len(spd)):
                h = pd_to_hfk(change(spd, [c]))
                unk.append(h['total_rank'] == 1)
            out['simplified_alternating_crossings'] = len(spd)
            out['mccoy_no_unknotting_crossing'] = not any(unk)
        except Exception as e:
            out['mccoy_error'] = repr(e)
    return out

if __name__ == '__main__':
    B = json.load(open(sys.argv[1]))
    rows = [r for r in B['changes'] if '#' not in r['ids'][0]]
    with Pool(11) as p:
        res = p.map(work, rows)
    from collections import Counter
    print('cover invariant vs named knot:', Counter((x['name'], x.get('cover_invariant_matches_named')) for x in res))
    print('name-free McCoy:', Counter((x['name'], x.get('simplified_alternating_crossings'), x.get('mccoy_no_unknotting_crossing'), 'err' if 'mccoy_error' in x else '') for x in res if 'simplified_alternating_crossings' in x or 'mccoy_error' in x))
    json.dump(res, open(sys.argv[2], 'w'), indent=0)
