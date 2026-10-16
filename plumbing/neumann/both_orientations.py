#!/usr/bin/env python
"""Both normal forms for the knots whose verdict is NOT APPLICABLE (class (2) of Section 3.1 of the paper).

Every tree accepted by precheck.realisations is in normal form (nodes of valency >= 3, all other weights <= -2) and
its plumbed manifold has a Regina name of Sigma_2.  Regina's names do not distinguish Sigma_2 from -Sigma_2, so by the
uniqueness of the normal form (Neumann 1989, Thm 1.2) each accepted tree is the normal form of Sigma_2 or of -Sigma_2,
and two different accepted trees are the normal forms of both.  This script lists the different accepted trees of
each knot and their positive indices p+; class (2) requires two trees, both with p+ >= 1.

The best tree of each orientation choice is read from the records of precheck.py (precheck_u12/u23/sigma0.json,
../classify47_normal.json) or recomputed from the recorded Regina names (run_fresh/run_retry13/run_pachner4.json keep
only p+).  precheck.one keeps one tree per orientation choice; when these coincide, all accepted trees are recomputed
from the recorded names.
Usage: both_orientations.py [workers]          Output: both_orientations.json
"""
import json, os, sys
from collections import Counter
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
import precheck as PC


def J(*p):
    return json.load(open(os.path.join(TOOL, *p)))


def canon(w, edges):
    """Canonical string of a weighted tree (minimum over all roots)."""
    adj = {i: [] for i in range(len(w))}
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    enc = lambda v, p: '(' + str(w[v]) + ''.join(sorted(enc(u, v) for u in adj[v] if u != p)) + ')'
    return min(enc(r, -1) for r in range(len(w)))


def in_normal_form(w, edges):
    deg = Counter(v for e in edges for v in e)
    return all(w[i] <= -2 for i in range(len(w)) if deg[i] <= 2)


def recorded_names(n):
    out = set()
    g = J('graph_u1.json').get(n)
    if g:
        out |= set(g.get('names_sigma2', [])) | set(g.get('names_mirror', []))
    o = J('owens_u3_targets.json').get(n)
    if o:
        out |= set(o.get('info', {}).get('names_sigma2', {})) | set(o.get('info', {}).get('names_mirror', {}))
    for f in ('reach_u12_orb.json', 'reach_u23_orb.json'):
        out |= {x['sigma2'] for x in J(f) if x['name'] == n and x.get('sigma2')}
    for f in ('recognise_more_fresh.json', 'recognise_retry.json', 'pachner_recognise.json'):
        r = J(f).get(n)
        if r and r.get('names'):
            out |= set(r['names'])
    c = J('classify47_normal.json').get(n)
    if c:
        out |= set(c.get('names', []))
    return sorted(out)


def best_trees(args):
    n, names = args
    return n, PC.one((n, names))[1]


def all_trees(n):
    names = recorded_names(n)
    D = int(PC.BY[n]['determinant'])
    trees = {}
    for nm in names:
        for orient, w, edges, nodes, pplus in PC.realisations(nm, set(names), D)[1]:
            trees.setdefault(canon(w, edges), dict(weights=w, edges=edges, pplus=pplus))
    return n, len(names), trees


if __name__ == '__main__':
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    recs, source = {}, {}
    for f in ('precheck_u12.json', 'precheck_u23.json', 'precheck_sigma0.json'):
        for n, r in J('neumann', f).items():
            if r['verdict'] == 'NOT APPLICABLE':
                recs[n], source[n] = r, f
    for n, r in J('classify47_normal.json').items():
        if r['verdict'] == 'NOT APPLICABLE':
            recs[n], source[n] = r, 'classify47_normal.json'
    redo = []
    for f, names_file in (('run_fresh.json', 'recognise_more_fresh.json'), ('run_retry13.json', 'recognise_retry.json'),
                          ('run_pachner4.json', 'pachner_recognise.json')):
        nm = J(names_file)
        for n, r in J('neumann', f).items():
            if r['applicability'] == 'NOT APPLICABLE' and n not in recs:
                redo.append((n, sorted(nm[n]['names']))); source[n] = f
    with Pool(workers) as p:
        for n, r in p.imap_unordered(best_trees, redo):
            assert r['verdict'] == 'NOT APPLICABLE', (n, r['verdict'])
            recs[n] = r
    out = {}
    for n, r in recs.items():
        trees = {}
        for v in r['normal_forms'].values():
            assert in_normal_form(v['weights'], v['edges']), n
            trees.setdefault(canon(v['weights'], v['edges']), dict(weights=v['weights'], edges=v['edges']))
        out[n] = dict(source=source[n], recomputed_all=False, trees=list(trees.values()))
    single = sorted(n for n, r in out.items() if len(r['trees']) < 2)
    print('records', len(out), '| one tree per orientation choice coincides for', len(single), flush=True)
    with Pool(workers) as p:
        for n, nn, trees in p.imap_unordered(all_trees, single):
            assert nn and all(in_normal_form(t['weights'], t['edges']) for t in trees.values()), n
            out[n].update(recomputed_all=True, names=nn, trees=list(trees.values()))
    for n, r in out.items():
        for t in r['trees']:
            t['pplus'] = PC.index(t['weights'], t['edges'])[0]
        r['distinct_trees'] = len(r['trees'])
        r['pplus'] = sorted(t['pplus'] for t in r['trees'])
    json.dump(dict(sorted(out.items())), open(os.path.join(HERE, 'both_orientations.json'), 'w'), indent=1)
    print('distinct accepted trees:', dict(Counter(r['distinct_trees'] for r in out.values())))
    print('min p+ over the trees:', dict(Counter(min(r['pplus']) for r in out.values())))
