#!/usr/bin/env python
"""Calibrate / validate formula.py against Regina's recognition of triangulated plumbings (plumbing.plumbed_manifold,
the route validated on stars, 11n22 and the 34 u = 2 knots).

Random negative-definite trees: pairs (two nodes with 2 legs each, chain of length 0..4) and triples (end nodes with 2
legs, centre node with 0..1 leg, chains of length 0..2), weights -1..-6 on chains, -2..-5 on legs, -1..-6 on nodes.
A sample counts when Regina recognises the triangulation (several randomisations); the formula name must be equal.
For triples every convention for the two matrices is tried and the ones matching all samples are reported.
Output: calibrate.json
"""
import itertools, json, random, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
sys.path.insert(0, '..')
import plumbing as P
import graph_u1 as GU            # fast_regina_name patch
import formula as F
import regina

PAIR_CONV = ('J', 'b', 1, 1, False)


def random_tree(rng, kind):
    w, edges = [], []
    def add(x):
        w.append(x); return len(w) - 1
    nnodes = 2 if kind == 'pair' else 3
    nodes = [add(rng.randint(-6, -1)) for _ in range(nnodes)]
    nlegs = [2] * nnodes if kind == 'pair' else [2, rng.randint(0, 1), 2]
    for v, k in zip(nodes, nlegs):
        for _ in range(k):
            prev = v
            for _ in range(rng.randint(1, 3)):
                x = add(rng.randint(-5, -2)); edges.append((prev, x)); prev = x
    for i in range(nnodes - 1):
        prev = nodes[i]
        for _ in range(rng.randint(0, 4 if kind == 'pair' else 2)):
            x = add(rng.randint(-6, -1)); edges.append((prev, x)); prev = x
        edges.append((prev, nodes[i + 1]))
    return w, edges, nodes


def recognise(w, edges, tries=10):
    M = P.plumbed_manifold(w, edges)
    nm, size, T = GU.fast_regina_name(M, tries=tries)
    return nm


def triple_names(w, edges, nodes):
    legs, chains = F.decompose(w, edges, nodes)
    out = {}
    for E, V, sb, se, t0, t1, inv0, inv1 in itertools.product(F.EDGES, F.VERTS, (1,), (1,), (False, True), (False, True), (False, True), (False, True)):
        M0 = F.chain_matrix(chains[0], E, V); M1 = F.chain_matrix(chains[1], E, V)
        def fix(M, t, inv):
            if inv:
                d = M[0][0] * M[1][1] - M[0][1] * M[1][0]
                M = [[d * M[1][1], -d * M[0][1]], [-d * M[1][0], d * M[0][0]]]
            if t:
                M = [[M[0][0], M[1][0]], [M[0][1], M[1][1]]]
            return M
        A, B = fix(M0, t0, inv0), fix(M1, t1, inv1)
        try:
            g = regina.GraphTriple(F.piece(w[nodes[0]], legs[nodes[0]], 1), F.piece(w[nodes[1]], legs[nodes[1]], 2),
                                   F.piece(w[nodes[2]], legs[nodes[2]], 1),
                                   regina.Matrix2(*A[0], *A[1]), regina.Matrix2(*B[0], *B[1]))
            out[(E, V, t0, t1, inv0, inv1)] = g.name()
        except Exception as e:
            out[(E, V, t0, t1, inv0, inv1)] = 'ERR ' + repr(e)[:40]
    return out


if __name__ == '__main__':
    rng = random.Random(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
    res = dict(pair=[], triple=[])
    t0 = time.time()
    # pairs
    while len(res['pair']) < 60 and time.time() - t0 < 1500:
        w, edges, nodes = random_tree(rng, 'pair')
        Q = P.plumbing_matrix(w, edges)
        if np.linalg.eigvalsh(Q.astype(float)).max() > -1e-9 or abs(round(np.linalg.det(Q.astype(float)))) < 2:
            continue
        nm = recognise(w, edges)
        if not nm or not nm.startswith('SFS [D'):
            continue
        legs, chains = F.decompose(w, edges, nodes)
        fn = F.pair_name(w[0], legs[0], w[1], legs[1], chains[0], PAIR_CONV)
        res['pair'].append(dict(weights=w, edges=edges, recognised=nm, formula=fn, equal=(nm == fn), chain=chains[0]))
    print('pairs: formula = recognised on', sum(r['equal'] for r in res['pair']), '/', len(res['pair']),
          '| chain lengths', sorted({len(r['chain']) for r in res['pair']}), flush=True)
    for r in res['pair']:
        if not r['equal']:
            print('  MISMATCH chain', r['chain'], '\n    rec', r['recognised'], '\n    for', r['formula'])
    # triples
    t1 = time.time(); conv_ok = None; ntr = 0
    while ntr < 30 and time.time() - t1 < 1500:
        w, edges, nodes = random_tree(rng, 'triple')
        Q = P.plumbing_matrix(w, edges)
        if np.linalg.eigvalsh(Q.astype(float)).max() > -1e-9 or abs(round(np.linalg.det(Q.astype(float)))) < 2:
            continue
        nm = recognise(w, edges)
        if not nm or nm.count('SFS') != 3:
            continue
        names = triple_names(w, edges, nodes)
        good = {c for c, x in names.items() if x == nm}
        conv_ok = good if conv_ok is None else conv_ok & good
        ntr += 1
        res['triple'].append(dict(weights=w, edges=edges, recognised=nm, matching_conventions=len(good)))
    print('triples:', ntr, 'samples | conventions matching all:', sorted(conv_ok)[:6] if conv_ok else conv_ok, flush=True)
    json.dump(dict(pair=res['pair'], triple=res['triple'], triple_conventions=sorted(conv_ok) if conv_ok else []),
              open('calibrate.json', 'w'), indent=1, default=str)
