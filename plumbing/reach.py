"""How many open knots does the graph-manifold tool reach?  Classify Sigma_2 for the 437 knots with u in [1,2]."""
import json, sys, warnings
warnings.filterwarnings('ignore')
from collections import Counter
sys.path.insert(0, '../lib')
import general as G
from ranges import load_ranges
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
R = load_ranges('baseline')          # open ranges before data/new_values_20261015.json
targets = sorted(n for n, r in R.items() if r['lo'] == 1 and r['hi'] == 2)
c = Counter(); rows = []
for n in targets:
    nm, C = G.sigma2_name(json.loads(BY[n]['pd_notation']), tries=2)
    kind = 'hyperbolic' if nm == 'HYPERBOLIC' else ('unrecognised' if nm is None else (G.parse_name(nm)[0] or 'other-recognised'))
    c[kind] += 1; rows.append((n, kind, nm))
print(dict(c)); json.dump(rows, open('reach_u12.json', 'w'), indent=1)
