#!/usr/bin/env python
"""u = 1 obstruction for open [1,2] knots whose Sigma_2 is a graph manifold with two or three Seifert pieces
(Regina GraphPair / GraphTriple), via negative-definite plumbing trees.

Identification.  Sigma_2 by the orbifold route ((2,0) filling, 2-fold cyclic cover).  All Regina names over many
randomised triangulations of Sigma_2 and of -Sigma_2 (= Sigma_2 of the mirror diagram) are collected.  For every
recognised shape (one node per Seifert piece, legs from both normalisations of each exceptional fibre, chains between
consecutive nodes) plumbing trees are filtered by |det| = D, negative definiteness and linking form = +-(V+V^T)^{-1};
a tree is accepted when Regina recognises its boundary with one of the collected names.  Equal names are equal
manifolds; a different name for the same manifold can only lose a knot, never give a wrong identification.
Correction terms d = max (K^2+|V|)/4: Ozsvath-Szabo for <= 1 bad vertex, Ozsvath-Stipsicz-Szabo for 2, Zemke
(lattice = Heegaard Floer homology for negative-definite plumbing trees) in general.  Floer ranks (graded roots) and
the Casson-Walker identity of crosscheck.py only for <= 1 bad vertex.
u = 1 test: v2's u1_admissible (both orientations, all affine identifications), with HF_red ranks when available.
Usage: graph_u1.py [workers]          Output: graph_u1.json (written after every knot), graph_u1.log
"""
import json, os, sys, time, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'montesinos'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
import _compat  # noqa
import numpy as np
import montesinos_u1 as MU
import owens_obstruction as OO
from hf_red import hf_red_all
import general as G
import plumbing as P
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
SHAPES = {'pair': ('D', 'D'), 'triple': ('D', 'A', 'D')}
SEARCH = {'pair': dict(erange=range(-8, 2), chain_w=(-1, -2, -3, -4, -5, -6), max_chain=3),
          'triple': dict(erange=range(-6, 1), chain_w=(-1, -2, -3, -4), max_chain=2)}


def fast_regina_name(M, tries=6, seed=0):
    """plumbing.regina_name without simplifyExhaustive (which does not terminate in practice on the larger plumbing
    triangulations): SnapPy randomisations + intelligentSimplify, first recognised name."""
    import regina
    T = None
    for t in range(tries):
        M2 = M.copy()
        if t:
            M2.randomize()
        T = regina.Triangulation3(M2.filled_triangulation()._to_string()); T.intelligentSimplify()
        rec = regina.StandardTriangulation.recognise(T)
        if rec is not None and rec.manifold() is not None:
            return rec.manifold().name(), T.size(), T
    return None, T.size(), T


P.regina_name = fast_regina_name          # used by general.find_plumbing


def names_orb(pd, tries):
    import snappy, regina
    E = snappy.Link(pd).exterior()
    names = {}
    for t in range(tries):
        E2 = E.copy()
        if t:
            E2.randomize()
        E2.dehn_fill((2, 0))
        C = E2.covers(2, cover_type='cyclic')[0]
        T = regina.Triangulation3(C.filled_triangulation()._to_string()); T.intelligentSimplify()
        rec = regina.StandardTriangulation.recognise(T)
        if rec and rec.manifold():
            nm = rec.manifold().name(); names[nm] = names.get(nm, 0) + 1
    return names


def shape_key(kind, pieces):
    """Tree family of a recognised name: leg_options depends on each fibre only through (alpha, {beta, -beta} mod alpha),
    and the search is symmetric under reversing the row of nodes."""
    fib = [tuple(sorted((a, min(b % a, (-b) % a)) for a, b in f)) for _, f in pieces]
    return min(tuple(fib), tuple(reversed(fib)))


def R_of(name):
    import kinfo
    J = kinfo.parse_jones(name); s = int(BY[name]['signature']); D = int(BY[name]['determinant'])
    V1 = sum(c * (-1) ** (e % 2) for e, c in J.items())
    dV1 = sum(e * c * (-1) ** ((e - 1) % 2) for e, c in J.items())
    assert abs(V1) == D
    return Fraction(-D * s, 4) + Fraction(D * dV1, 6 * V1)


def one(n):
    import spherogram
    t0 = time.time()
    pd = json.loads(BY[n]['pd_notation']); D = int(BY[n]['determinant'])
    mpd = [list(x) for x in spherogram.Link(pd).mirror().PD_code()]
    nK, nM = names_orb(pd, 16), names_orb(mpd, 16)
    if not any(G.parse_name(x)[0] in SHAPES for x in set(nK) | set(nM)):
        nK, nM = names_orb(pd, 40), names_orb(mpd, 40)
    allnames = set(nK) | set(nM)
    shapes = {}
    for nm in sorted(allnames):
        kind, pieces = G.parse_name(nm)
        if kind in SHAPES and tuple(b for b, _ in pieces) == SHAPES[kind]:
            shapes.setdefault((kind, shape_key(kind, pieces)), nm)
    rec = dict(det=D, names_sigma2=nK, names_mirror=nM, shapes=len(shapes), searches=[])
    found = None
    for (kind, _), nm in shapes.items():
        t1 = time.time()
        f, st = G.find_plumbing(nm, D, json.loads(BY[n]['seifert_matrix']), alt_names=tuple(allnames - {nm}),
                                log=lambda *a: None, **SEARCH[kind])
        rec['searches'].append(dict(shape_name=nm, stats=st))
        print(f'    {n}: shape {nm[:70]} -> {"found" if f else "none"} {st} ({time.time()-t1:.0f}s)', flush=True)
        if f is not None:
            found = dict(f, kind=kind); break
    if found is None:
        rec['verdict'] = 'NO PLUMBING' if shapes else 'NO PAIR/TRIPLE NAME'
        rec['seconds'] = round(time.time() - t0)
        return n, rec
    w, e = found['weights'], found['edges']
    Q = np.array(P.plumbing_matrix(w, e), dtype=np.int64)
    deg = [0] * len(w)
    for u, v in e:
        deg[u] += 1; deg[v] += 1
    nbad = sum(1 for i in range(len(w)) if deg[i] > -w[i])
    if nbad <= 1:
        d, cm = MU.correction_terms(Q)
    else:
        mq, cm = OO.m_Q(-np.array(Q, dtype=object))
        d = {g: -v for g, v in mq.items()}
    adm_d = MU.u1_admissible(d, cm, D, Q)
    # graded roots (v2's hf_red) enumerate a box of prod(|w|+1) points per Spin^c class with exact sympy arithmetic:
    # computed when the d-test passes (needed for the Floer test; box*D <= 2e6) or, for the Casson-Walker check, when
    # box*D <= 3e5.  Otherwise a passing knot is reported as 'PASS (d only)'.
    box = int(np.prod([abs(x) + 1 for x in w]))
    ranks = hf_red_all(Q)[0] if nbad <= 1 and ((adm_d and box * D <= 2_000_000) or box * D <= 300_000) else None
    if adm_d is None:
        verdict = 'OBSTRUCTED (H1 not cyclic)'
    elif not adm_d:
        verdict = 'OBSTRUCTED'
    elif ranks is None:
        verdict = 'PASS (d only)'
    else:
        verdict = 'PASS' if MU.u1_admissible(d, cm, D, Q, ranks_by_class=ranks) else 'OBSTRUCTED_HF'
    rec.update(plumbing=found, bad_vertices=nbad, verdict=verdict, d_spin=str(d[tuple(0 for _ in cm.d)]), hf_box=box)
    if ranks is not None:
        L = sum(d.values()) + 2 * sum(ranks.values()); R = R_of(n)
        cw = [s for s in (1, -1) if s * L == R]
        s_link = 1 if found['orientation'] == '+' else -1
        rec.update(hf_red_ranks_nonzero=[x for x in ranks.values() if x], L=str(L), R=str(R), cw_signs=cw,
                   cw_consistent=bool(cw) and (s_link in cw or found['orientation_ambiguous']))
    rec['seconds'] = round(time.time() - t0)
    return n, rec


if __name__ == '__main__':
    from multiprocessing import Pool
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    new = {r['name']: r for r in json.load(open('reach_u12_orb.json'))}
    old = {n: k for n, k, nm in json.load(open('reach_u12.json'))}
    targets = sorted(n for n in new if new[n]['kind'] in SHAPES or (old.get(n) in SHAPES and new[n]['kind'] == 'unrecognised'))
    outf = 'graph_u1.json'
    out = json.load(open(outf)) if os.path.exists(outf) else {}
    todo = [n for n in targets if n not in out]
    print('targets', len(targets), '| to do', len(todo), flush=True)
    T0 = time.time()
    with Pool(workers) as pool:
        for n, rec in pool.imap_unordered(one, todo):
            out[n] = rec
            json.dump(out, open(outf, 'w'), indent=1, default=str)
            print(f"[{time.time()-T0:5.0f}s] {n:9s} D {rec['det']:4d} | {rec['verdict']:22s} | shapes {rec['shapes']} "
                  f"| bad {rec.get('bad_vertices')} | CW {rec.get('cw_consistent')} | {rec['seconds']}s", flush=True)
    obs = sorted(n for n, r in out.items() if r['verdict'].startswith('OBSTRUCTED'))
    print('verdicts', {v: sum(1 for r in out.values() if r['verdict'] == v) for v in sorted({r['verdict'] for r in out.values()})})
    print('OBSTRUCTED:', obs)
