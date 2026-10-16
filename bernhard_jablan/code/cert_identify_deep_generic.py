"""Deeper Regina Pachner search (height <= 3, both sides) for neighbours of TARGET whose names are in NAMES
and that cert_identify.py did not certify with Regina; updates results/cert_identify_TARGET.json."""
import json, sys, warnings
warnings.filterwarnings('ignore')
import regina
from database_knotinfo import link_list
by = {r['name']: r for r in link_list()[1:]}
target, names = sys.argv[1], sys.argv[2].split(',')
RP = f'../results/cert_identify_{target}.json'; BP = f'../results/bjset_{target}.json'
R = json.load(open(RP)); B = json.load(open(BP))
cache = {}
def sigs_of(T0, hmax):
    out = set()
    for h in range(hmax + 1):
        T = regina.Triangulation3(T0)
        if h: T.simplifyExhaustive(h, 8)
        out.add(T.isoSig())
    return out
for r in R:
    if r['name'] not in names or r['regina_cert']:
        continue
    if r['name'] not in cache:
        Tk = regina.Link.fromPD(json.loads(by[r['name']]['pd_notation'])).complement(); Tk.simplify()
        cache[r['name']] = sigs_of(Tk, 3)
    pd = next(c['pd'] for c in B['changes'] if c['ref'] == r['ref'] and c['subset'] == r['subset'] and c['crossing'] == r['crossing'])
    Tn = regina.Link.fromPD(pd).complement(); Tn.simplify()
    common = sigs_of(Tn, 3) & cache[r['name']]
    r['regina_cert'] = bool(common)
    if common:
        r['regina_height'] = 3; r['regina_sig'] = sorted(common)[0]
    print(r['name'], r['ref'], r['subset'], r['crossing'], 'regina certificate:', bool(common), flush=True)
json.dump(R, open(RP, 'w'), indent=0)
from collections import Counter
print(Counter((r['name'], r['snappy_cert'], r['regina_cert']) for r in R))
