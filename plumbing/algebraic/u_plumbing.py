#!/usr/bin/env python
"""Plumbing trees of Sigma_2(K) for the 32 knots of targets.txt, from their algebraic diagrams (alg_diagrams.json).
For both orientations (SIGMA = +-1): the tree, its normal form, |det| (must equal det K), the exact inertia of the
normal form, the nodes and the bad vertices (-e(v) < deg(v)), and the number of vertices at which bigon merging of the
KnotInfo diagram stops.  The orientation of Sigma_2(K) is fixed by the Casson-Walker identity L(Q) = R(K)
(cw_general.L_general, valid for indefinite trees; on all algebraic KnotInfo diagrams the orientation SIGMA = -1 of the
construction is Sigma_2(K)) and by the linking form; the conclusions do not need it, since both orientations are tested.
Output: u_plumbing.json
"""
import json, os, sys, warnings
warnings.filterwarnings('ignore')
from collections import Counter
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from validate_plumbing import setup, TOOL
from algebraic_plumbing import plumbing_from_algebraic, normal_form, plumbing_matrix
from cw_general import L_general
from basic_polyhedron import reduce_diagram


def bad_vertices(w, edges):
    deg = Counter(v for e in edges for v in e)
    return [i for i in range(len(w)) if -w[i] < deg[i]]


def one(n):
    setup()
    import sympy
    from sympy.matrices.normalforms import smith_normal_form
    import owens_u3_graph as OW, graph_u1 as GU
    from identify_trees import exact_inertia
    from both_orientations import canon
    sys.path.insert(0, os.path.join(TOOL, 'neumann'))
    from run_invert import orientation
    BY = GU.BY
    AD = json.load(open(os.path.join(HERE, 'alg_diagrams.json')))
    pd = AD[n]['pd']
    D = int(BY[n]['determinant'])
    rec = dict(D=D, signature=int(BY[n]['signature']), diagram_crossings=len(pd),
               knotinfo_diagram_vertices=reduce_diagram(json.loads(BY[n]['pd_notation'])))
    for sigma in (1, -1):
        w, edges = plumbing_from_algebraic(pd, sigma)
        Q = plumbing_matrix(w, edges)
        wn, en = normal_form(w, edges)
        Qn = plumbing_matrix(wn, en)
        deg = Counter(v for e in en for v in e)
        rec[str(sigma)] = dict(raw_vertices=len(w), det_raw=abs(int(sympy.Matrix(Q).det())),
                               weights=wn, edges=en, canon=canon(wn, en), inertia=exact_inertia(wn, en),
                               nodes=[i for i in range(len(wn)) if deg[i] >= 3], bad=bad_vertices(wn, en),
                               L=str(L_general(wn, en)),
                               smith=[int(x) for x in smith_normal_form(sympy.Matrix(Qn), domain=sympy.ZZ).diagonal() if abs(x) != 1],
                               linking_orientation=orientation(Qn, json.loads(BY[n]['seifert_matrix'])))
    R = GU.R_of(n)
    rec['R'] = str(R)
    rec['cw_sigma'] = [s for s in (1, -1) if rec[str(s)]['L'] == str(R)]
    return n, rec


if __name__ == '__main__':
    AD = json.load(open(os.path.join(HERE, 'alg_diagrams.json')))
    out = {}
    with Pool(8) as p:
        for n, r in p.imap_unordered(one, sorted(AD)):
            out[n] = r
            s = r['cw_sigma']
            line = [n, 'D', r['D'], 'sig', r['signature'], 'cw sigma', s, 'link', r['-1']['linking_orientation'], r['1']['linking_orientation']]
            for sg in ('1', '-1'):
                t = r[sg]; I = t['inertia']
                line += [f"| s={sg}: det {t['det_raw']} H1 {t['smith']} p+ {I['pplus']} p- {I['pminus']} nodes {len(t['nodes'])}"]
            print(*line, flush=True)
    json.dump(dict(sorted(out.items())), open(os.path.join(HERE, 'u_plumbing.json'), 'w'), indent=1)
