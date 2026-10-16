"""Deeper Regina Pachner search (height <= 3, both sides) for 12n563 neighbours not certified by cert_identify.py; updates cert_identify_13n_1587.json."""
import json, sys, warnings
warnings.filterwarnings('ignore')
import regina
from database_knotinfo import link_list
by = {r['name']: r for r in link_list()[1:]}
R = json.load(open('../results/cert_identify_13n_1587.json'))
B = json.load(open('../results/bjset_13n_1587.json'))
todo = [r for r in R if r['name'] == '12n_563' and not r['regina_cert']]
pdk = json.loads(by['12n_563']['pd_notation'])
Tk = regina.Link.fromPD(pdk).complement(); Tk.simplify()
sigk = set()
for h in (0, 1, 2, 3):
    T = regina.Triangulation3(Tk)
    if h: T.simplifyExhaustive(h, 8)
    sigk.add(T.isoSig())
for r in todo:
    pd = next(c['pd'] for c in B['changes'] if c['ref'] == r['ref'] and c['subset'] == r['subset'] and c['crossing'] == r['crossing'])
    Tn = regina.Link.fromPD(pd).complement(); Tn.simplify()
    found = None
    for h in (0, 1, 2, 3):
        T = regina.Triangulation3(Tn)
        if h: T.simplifyExhaustive(h, 8)
        if T.isoSig() in sigk:
            found = (h, T.size(), T.isoSig()); break
    # fallback: direct isomorphism test of minimal found triangulations
    print(r['ref'], r['subset'], r['crossing'], 'regina certificate:', found, '| sizes', Tn.size(), Tk.size())
    r['regina_cert'] = found is not None
    if found:
        r['regina_height'], r['regina_tets'], r['regina_sig'] = found[0], (found[1], found[1]), found[2]
json.dump(R, open('../results/cert_identify_13n_1587.json', 'w'), indent=0)
