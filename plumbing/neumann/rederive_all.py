#!/usr/bin/env python
"""Re-derive every adopted graph-manifold value with the single method of the paper (Neumann normal form built from the
Regina name, precheck.one -> run_fresh pipeline), and compare with the verdicts recorded by the earlier routes.
Covers the 43 strict entries of data/new_values_20261015.json; each normal form is also triangulated (surgery on the
plumbing link) and recognised by Regina, independently of the formula conventions.
Usage: rederive_all.py [workers]          Output: rederive_all.json
"""
import json, os, sys
from multiprocessing import Pool
from collections import Counter
import run_fresh as RF
CEN = os.path.join(RF.TOOL, '..')


def names_for(n):
    T = RF.TOOL
    for f in ('star_u1.json', 'star_u1_orb.json'):
        d = json.load(open(os.path.join(T, f)))
        if n in d:
            return [d[n]['sigma2']]
    g = json.load(open(os.path.join(T, 'graph_u1.json')))
    if n in g:
        return sorted(set(g[n]['names_sigma2']) | set(g[n]['names_mirror']))
    o = json.load(open(os.path.join(T, 'owens_u3_targets.json')))
    if n in o and o[n]:
        nm = o[n].get('names') or o[n].get('info') or {}
        return sorted(set(nm.get('names_sigma2', {})) | set(nm.get('names_mirror', {})))
    fr = json.load(open(os.path.join(T, 'recognise_more_fresh.json')))
    return fr[n]['names']


def one(e):
    import spherogram
    n = e['knot']
    pd = json.loads(RF.BY[n]['pd_notation']); mpd = [list(x) for x in spherogram.Link(pd).mirror().PD_code()]
    names = sorted(set(names_for(n)) | set(RF.GU.names_orb(pd, 16)))
    _, rec = RF.one((n, names))          # adds 16 recognitions of the mirror
    pl = rec.get('plumbing') or {}
    out = dict(expected=e['lo'], source=e['source'], applicability=rec['applicability'], verdict=rec.get('verdict'),
               cw_consistent=rec.get('cw_consistent'), ni_wu=rec.get('independent_ni_wu_survivors'),
               vertices=len(pl.get('weights', [])), seconds=rec['seconds'], plumbing=pl)
    if pl:   # independent identification: triangulate the plumbed manifold (surgery on the plumbing link) and recognise it
        nm, _, _ = RF.GU.fast_regina_name(RF.P.plumbed_manifold(pl['weights'], pl['edges']), tries=12)
        out.update(plumbing_name=nm, plumbing_name_matches=nm in set(names) | set(RF.GU.names_orb(mpd, 16)))
    return n, out


if __name__ == '__main__':
    vals = [e for e in json.load(open(os.path.join(CEN, 'data', 'new_values_20261015.json')))['values'] if e['strict']]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as p:
        out = dict(p.map(one, vals))
    json.dump(out, open('rederive_all.json', 'w'), indent=1)
    ok = {n: r for n, r in out.items() if str(r['verdict']).startswith('OBSTRUCTED') and r['cw_consistent']}
    print('entries', len(out), '| reproduced (obstructed + CW consistent):', len(ok))
    print('verdicts', dict(Counter(r['verdict'] for r in out.values())), '| applicability', dict(Counter(r['applicability'] for r in out.values())))
    print('not reproduced:', {n: (r['verdict'], r['applicability']) for n, r in out.items() if n not in ok})
    print('seconds per knot: max', max(r['seconds'] for r in out.values()), '| mean', round(sum(r['seconds'] for r in out.values()) / len(out), 2))
    print('triangulated normal form recognised as Sigma_2:', sum(bool(r.get('plumbing_name_matches')) for r in out.values()), '/', len(out),
          '| unmatched:', {n: r.get('plumbing_name') for n, r in out.items() if not r.get('plumbing_name_matches')})
