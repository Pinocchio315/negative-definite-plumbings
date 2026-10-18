#!/usr/bin/env python
"""Validation of algebraic_plumbing.py on every KnotInfo knot with at most 13 crossings whose KnotInfo diagram is
algebraic (bigon merging reaches one vertex).
  (1) |det Q| = det K for the trees of both signs;
  (2) agreement with the trees identified by triangulation (neumann/identify_trees.json): for a knot of class (2) the two
      normal forms are exactly its two identified trees, the normal forms of Sigma_2(K) and -Sigma_2(K); for a knot of
      class (3) its identified negative-definite tree is one of the two normal forms.
  (3) every returned tree is in normal form: each vertex of degree <= 2 has weight <= -2.
Output: validate_plumbing.json
"""
import json, os, sys, time, warnings
from collections import Counter
warnings.filterwarnings('ignore')
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.abspath(os.path.join(HERE, '..'))                 # plumbing/
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))           # repository root
sys.path.insert(0, HERE)
from basic_polyhedron import reduce_diagram
from algebraic_plumbing import plumbing_from_algebraic, normal_form, plumbing_matrix


def setup():
    sys.path.insert(0, TOOL)
    sys.path.insert(0, os.path.join(TOOL, 'neumann'))
    sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'montesinos'))
    sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))


_BY = {}


def knotinfo():
    if not _BY:
        from database_knotinfo import link_list
        _BY.update({r['name']: r for r in link_list()[1:]})
    return _BY


def one(n):
    setup()
    from both_orientations import canon
    import sympy
    BY = knotinfo()
    pd = json.loads(BY[n]['pd_notation'])
    if reduce_diagram([list(x) for x in pd]) != 1:
        return n, None
    rec = dict(D=int(BY[n]['determinant']))
    for sigma in (1, -1):
        w, edges = plumbing_from_algebraic(pd, sigma)
        d = abs(int(sympy.Matrix(plumbing_matrix(w, edges)).det()))
        try:
            wn, en = normal_form(w, edges)
            c = canon(wn, en)
        except ValueError as e:
            wn, en, c = None, None, str(e)
        rec[str(sigma)] = dict(det=d, canon=c, weights=wn, edges=en)
    return n, rec


if __name__ == '__main__':
    from database_knotinfo import link_list
    names = [r['name'] for r in link_list()[1:] if r['name'][0].isdigit() and r.get('pd_notation')
             and int(r['name'].split('_')[0].rstrip('an')) <= 13]
    IT = json.load(open(os.path.join(TOOL, 'neumann', 'identify_trees.json')))
    t0 = time.time(); out = {}
    with Pool(8) as p:
        for k, (n, r) in enumerate(p.imap_unordered(one, names, chunksize=20)):
            if r is not None:
                out[n] = r
            if k % 2000 == 0:
                print(f'[{time.time()-t0:5.0f}s] {k}/{len(names)} algebraic so far {len(out)}', flush=True)
    bad_det = sorted(n for n, r in out.items() if r['1']['det'] != r['D'] or r['-1']['det'] != r['D'])
    compared = {2: [], 3: []}; miss = []
    for n, r in out.items():
        if n not in IT:
            continue
        cls = IT[n]['cls']
        rec_c = {t['canon'] for t in IT[n]['trees']}
        mine = {r['1']['canon'], r['-1']['canon']}
        ok = (rec_c == mine) if cls == 2 else (len(rec_c) == 1 and rec_c <= mine)
        compared[cls].append(n)
        if not ok:
            miss.append(n)
        r['identified'] = dict(cls=cls, agrees=ok)
    def in_normal_form(t):
        deg = Counter(v for e in t['edges'] for v in e)
        return all(x <= -2 for v, x in enumerate(t['weights']) if deg[v] <= 2)
    not_nf = sorted(n for n, r in out.items() for s in ('1', '-1') if r[s]['weights'] is None or not in_normal_form(r[s]))
    summary = dict(algebraic=len(out), det_mismatch=bad_det, identified_class2=len(compared[2]),
                   identified_class3=len(compared[3]), identified_disagree=sorted(miss), not_in_normal_form=not_nf)
    json.dump(dict(summary=summary, records=dict(sorted(out.items()))), open(os.path.join(HERE, 'validate_plumbing.json'), 'w'), indent=1)
    print(json.dumps(summary, indent=1))
