#!/usr/bin/env python
"""Identification of every accepted tree by triangulation, independently of the conventions of formula.py, and the
exact inertia of its intersection matrix.

Trees: both trees of each knot of both_orientations.json (class (2): the normal forms of Sigma_2 and -Sigma_2), and the
negative-definite tree of each APPLICABLE record (class (3); the records run_pachner4, run_retry13, run_fresh,
rederive_all, normal_case3, precheck_sigma0 and ../classify47_normal.json, read in the order of
../paper/make_tables.py).  The plumbed manifold of each tree is triangulated as surgery on the plumbing link
(plumbing.plumbed_manifold), recognised with Regina, and its name is compared with the recorded Regina names of
Sigma_2; when it is not among them, 16 more randomised recognitions of the cover of the knot and of its mirror are
added.  Regina's names do not depend on the orientation.  The inertia (p+, p-, nullity) of Q_Gamma is computed
exactly from the characteristic polynomial (all roots are real, so Descartes' rule of signs counts them), together
with the exact determinant.
Usage: identify_trees.py [workers]          Output: identify_trees.json
"""
import json, os, sys
from collections import Counter
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, '..')
sys.path.insert(0, HERE)
sys.path.insert(0, TOOL)
from both_orientations import J, canon, recorded_names


def all_names(n):
    out = set(recorded_names(n))
    for f in ('star_u1.json', 'star_u1_orb.json'):
        r = J(f).get(n)
        if r and r.get('sigma2'):
            out.add(r['sigma2'])
    return out


def exact_inertia(w, edges):
    import sympy
    n = len(w)
    Q = sympy.zeros(n, n)
    for i, x in enumerate(w):
        Q[i, i] = x
    for a, b in edges:
        Q[a, b] = Q[b, a] = 1
    x = sympy.Symbol('x')
    c = [int(t) for t in sympy.Poly(Q.charpoly(x).as_expr(), x).all_coeffs()]      # leading coefficient first
    changes = lambda s: sum(1 for u, v in zip(s, s[1:]) if u * v < 0)
    nonzero = [t for t in c if t]
    zero = len(c) - 1 - max(i for i, t in enumerate(c) if t)                         # multiplicity of the root 0
    neg = [t * (-1) ** (len(c) - 1 - i) for i, t in enumerate(c)]                    # coefficients of p(-x)
    return dict(pplus=changes(nonzero), pminus=changes([t for t in neg if t]), nullity=zero, det=int(Q.det()))


def one(job):
    import snappy, spherogram
    import graph_u1 as GU, plumbing as P
    n, cls, trees = job
    names = all_names(n)
    out, extra = [], None
    for w, edges in trees:
        nm = None
        for tries in (12, 48):
            nm, _, _ = GU.fast_regina_name(P.plumbed_manifold(w, edges), tries=tries)
            if nm:
                break
        matched = nm in names
        if nm and not matched and extra is None:
            import precheck as PC
            pd = json.loads(PC.BY[n]['pd_notation'])
            mpd = [list(x) for x in spherogram.Link(pd).mirror().PD_code()]
            extra = set(GU.names_orb(pd, 16)) | set(GU.names_orb(mpd, 16))
        if nm and not matched and extra:
            matched = nm in extra
        out.append(dict(weights=w, edges=[list(e) for e in edges], canon=canon(w, edges), name=nm, matched=matched,
                        **exact_inertia(w, edges)))
    return n, dict(cls=cls, extra_names_used=extra is not None, trees=out)


def jobs():
    BO = J('neumann', 'both_orientations.json')
    out = {n: (n, 2, [(t['weights'], t['edges']) for t in r['trees']]) for n, r in BO.items()}
    seen = set(out)
    for f in ('run_pachner4.json', 'run_retry13.json', 'run_fresh.json', 'rederive_all.json', 'normal_case3.json'):
        for n, r in J('neumann', f).items():
            if n in seen:
                continue
            seen.add(n)
            if r['applicability'] == 'APPLICABLE':
                out[n] = (n, 3, [(r['plumbing']['weights'], r['plumbing']['edges'])])
    for src in (J('neumann', 'precheck_sigma0.json'), J('classify47_normal.json')):
        for n, r in src.items():
            if n in seen:
                continue
            seen.add(n)
            if r['verdict'] == 'APPLICABLE':
                o = [k for k, v in r['pplus'].items() if v == 0][0]
                out[n] = (n, 3, [(r['normal_forms'][o]['weights'], r['normal_forms'][o]['edges'])])
    return sorted(out.values())


if __name__ == '__main__':
    todo = jobs()
    print('knots', len(todo), dict(Counter(c for _, c, _ in todo)), '| trees', sum(len(t) for _, _, t in todo), flush=True)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as p:
        res = dict(p.imap_unordered(one, todo))
    json.dump(dict(sorted(res.items())), open(os.path.join(HERE, 'identify_trees.json'), 'w'), indent=1)
    T = [(r['cls'], t) for r in res.values() for t in r['trees']]
    print('recognised:', dict(Counter((c, t['name'] is not None) for c, t in T)))
    print('identified with Sigma_2:', dict(Counter((c, t['matched']) for c, t in T)),
          '| extra randomised names needed for', sum(r['extra_names_used'] for r in res.values()), 'knots')
    print('exact p+ by class:', dict(Counter((c, t['pplus']) for c, t in T)), '| nullity > 0:', sum(t['nullity'] > 0 for _, t in T))
