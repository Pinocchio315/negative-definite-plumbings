"""Rigorous treatment of the composite neighbours of a target (review point 3).
For each composite single-change diagram: total signature (exact), and for each connected
summand: crossings after simplification, HFK, comparison (mirror included) with named knots.
HFK detects the trefoil (Ghiggini 2008); hyperbolic summands are also identified by SnapPy."""
import json, sys, warnings
warnings.filterwarnings('ignore')
import spherogram, snappy
from knot_floer_homology import pd_to_hfk
from database_knotinfo import link_list
sys.path.insert(0, __import__('os').path.dirname(__file__))
from bf2 import charpoly_signature
by = {r['name']: r for r in link_list()[1:]}

def key(h, mirror=False):
    r = h['ranks'].items()
    if mirror: r = [((-a, -m), v) for (a, m), v in r]
    return sorted((list(k), v) for k, v in r)

REF = {n: pd_to_hfk(json.loads(by[n]['pd_notation'])) for n in ['3_1', '4_1', '5_1', '5_2', '7_3', '7_5', '7_7']}
B = json.load(open(sys.argv[1]))
out = []
for r in B['changes']:
    if '#' not in r['ids'][0]:
        continue
    L = spherogram.Link(r['pd'])
    V = L.seifert_matrix(); n = len(V)
    sig = charpoly_signature([[V[i][j] + V[j][i] for j in range(n)] for i in range(n)])
    L2 = spherogram.Link(r['pd']); L2.simplify('global')
    parts = []
    for P in L2.deconnect_sum():
        P.simplify('global')
        pdp = P.PD_code()
        h = pd_to_hfk(pdp)
        try:
            ids = [str(x) for x in P.exterior().identify()]
        except Exception:
            ids = []
        match = [m for m, hr in REF.items() if key(h) in (key(hr), key(hr, True))]
        parts.append(dict(crossings=len(P.crossings), hfk_total_rank=h['total_rank'], seifert_genus=h['seifert_genus'],
                          hfk_matches=match, snappy_identify=ids, nontrivial=h['total_rank'] > 1))
    out.append(dict(ref=r['ref'], subset=r['subset'], crossing=r['crossing'], total_signature=sig, summands=parts))
json.dump(out, open(sys.argv[2], 'w'), indent=1)
from collections import Counter
print('composite neighbour diagrams:', len(out))
print('signatures:', Counter(o['total_signature'] for o in out))
print('summand signatures (crossings, HFK match, SnapPy):', Counter(tuple((p['crossings'], tuple(p['hfk_matches']), tuple(p['snappy_identify'][:1])) for p in o['summands']) for o in out))
print('all summands nontrivial:', all(p['nontrivial'] for o in out for p in o['summands']))
