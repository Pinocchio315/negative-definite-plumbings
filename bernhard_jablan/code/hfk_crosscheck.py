"""Cross-check every neighbour identification by comparing knot Floer homology with the named knot's PD (up to mirror)."""
import json, sys, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
from collections import Counter
import enum_diagrams as en
from knot_floer_homology import pd_to_hfk
def key(h, mirror=False):
    r = h['ranks'].items()
    if mirror: r = [((-a, -m), v) for (a, m), v in r]
    return tuple(sorted(r))
def work(row):
    name = row['ids'][0]
    if name not in en.by:
        return name, 'skip (composite)'
    ref = pd_to_hfk(json.loads(en.by[name]['pd_notation']))
    h = pd_to_hfk(row['pd'])
    return name, ('HFK match' if key(h) in (key(ref), key(ref, True)) else 'HFK MISMATCH')
if __name__ == '__main__':
    B = json.load(open(sys.argv[1]))
    with Pool(11) as p:
        res = p.map(work, B['changes'])
    c = Counter(res)
    for k in sorted(c): print(k, c[k])
    json.dump(sorted(c.items()), open(sys.argv[1].replace('bjset_', 'hfk_crosscheck_'), 'w'))
