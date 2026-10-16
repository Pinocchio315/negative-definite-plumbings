#!/usr/bin/env python
"""Inverse problem: Regina name of Sigma_2 (star / GraphPair / GraphTriple) -> negative-definite plumbing tree, using
the calibrated formula (formula.py) instead of triangulating candidates.

Normal-form search space (a minimal negative-definite plumbing has no -1 vertex of degree <= 2, since blowing one down
keeps the form negative definite): legs and chains have weights <= -2; nodes have weights <= -1.
  legs: for every exceptional fibre (alpha, beta) of a piece, the Hirzebruch-Jung leg of alpha/(beta mod alpha) or,
        for the other orientation of the manifold, of alpha/(-beta mod alpha) (one global choice);
  chains (pairs): length 0..LMAX, weights in [-WMAX, -2];  triples: each chain length 0..3;
  nodes: the determinant of a tree is affine in each weight, so for given legs and chains the node weights with
        |det| = D are solved for (the last node from the others).
A candidate is accepted when it is negative definite, |det| = D, and its formula name is one of the names Regina gives
for Sigma_2 or -Sigma_2 (equal names = homeomorphic manifolds).  Accepted plumbings are then re-checked by the
pipeline's own tests (linking form, Casson-Walker identity) before any use.
"""
import itertools, sys, time
import numpy as np
sys.path.insert(0, '..')
import formula as F
import general as G
import plumbing as P

TRIPLE_CONV = ('J', 'b', False, False, True, False)


def leg_of(alpha, beta):
    b = beta % alpha
    return [-x for x in G.cf(alpha, b)]


def build(nodes_w, legs, chains):
    """Tree from node weights (row), legs per node, chain weights between consecutive nodes; returns (w, edges, nodes)."""
    w = list(nodes_w); edges = []; nodes = list(range(len(nodes_w)))
    def add(x):
        w.append(x); return len(w) - 1
    for v, ls in zip(nodes, legs):
        for leg in ls:
            prev = v
            for x in leg:
                y = add(x); edges.append((prev, y)); prev = y
    for i, ch in enumerate(chains):
        prev = nodes[i]
        for x in ch:
            y = add(x); edges.append((prev, y)); prev = y
        edges.append((prev, nodes[i + 1]))
    return w, edges, nodes


def det_of(w, edges):
    return round(np.linalg.det(P.plumbing_matrix(w, edges).astype(float)))


def negdef(w, edges):
    return np.linalg.eigvalsh(P.plumbing_matrix(w, edges).astype(float)).max() < -1e-9


def name_of(kind, w, edges, nodes):
    legs, chains = F.decompose(w, edges, nodes)
    if kind == 'star':
        return F.star_name(w[0], legs[0])
    if kind == 'pair':
        return F.pair_name(w[nodes[0]], legs[nodes[0]], w[nodes[1]], legs[nodes[1]], chains[0], ('J', 'b', 1, 1, False))
    import regina
    E, V, t0, t1, inv0, inv1 = TRIPLE_CONV
    M0 = F.chain_matrix(chains[0], E, V); M1 = F.chain_matrix(chains[1], E, V)
    d = M0[0][0] * M0[1][1] - M0[0][1] * M0[1][0]
    A = [[d * M0[1][1], -d * M0[0][1]], [-d * M0[1][0], d * M0[0][0]]]
    g = regina.GraphTriple(F.piece(w[nodes[0]], legs[nodes[0]], 1), F.piece(w[nodes[1]], legs[nodes[1]], 2),
                           F.piece(w[nodes[2]], legs[nodes[2]], 1), regina.Matrix2(*A[0], *A[1]), regina.Matrix2(*M1[0], *M1[1]))
    return g.name()


def chains_upto(L, wmax):
    for l in range(L + 1):
        yield from itertools.product(range(-2, -wmax - 1, -1), repeat=l)


def solve_last(nodes_w_partial, legs, chains, D, emin):
    """With all node weights but the last fixed, det is affine in the last one: det = a e + b."""
    w0, e0, n0 = build(nodes_w_partial + [-1], legs, chains)
    w1, e1, n1 = build(nodes_w_partial + [-2], legs, chains)
    d0, d1 = det_of(w0, e0), det_of(w1, e1)
    a = d0 - d1                       # det(e) = a e + b with det(-1) = d0, det(-2) = d1, so det(e) = d0 + a (e + 1)
    out = []
    for target in (D, -D):
        if a == 0:
            continue
        num = target - d0
        if num % a == 0:
            e = num // a - 1
            if emin <= e <= -1:
                out.append(e)
    return out


def search(names, D, LMAX=5, WMAX=7, EMIN=-30, time_limit=600):
    """names: Regina names of Sigma_2 and -Sigma_2.  Returns list of found plumbings (w, edges, nodes, name)."""
    t0 = time.time(); found = []; stats = dict(candidates=0, det=0, negdef=0, named=0)
    shapes = {}
    for nm in names:
        kind, pieces = G.parse_name(nm)
        if kind in ('star', 'pair', 'triple'):
            shapes.setdefault((kind, tuple(tuple(sorted(f)) for _, f in pieces)), (kind, pieces))
    for kind, pieces in shapes.values():
        for orient in (1, -1):
            legs = [[leg_of(a, orient * b) for a, b in f if a > 1] for _, f in pieces]
            if kind == 'star':
                for e in range(-1, EMIN - 1, -1):
                    w, edges, nodes = build([e], legs, [])
                    stats['candidates'] += 1
                    if abs(det_of(w, edges)) != D:
                        continue
                    stats['det'] += 1
                    if not negdef(w, edges):
                        continue
                    stats['negdef'] += 1
                    nm = name_of(kind, w, edges, nodes); stats['named'] += 1
                    if nm in names:
                        found.append((w, edges, nodes, nm))
                continue
            if kind == 'pair':
                chain_sets = ((c,) for c in chains_upto(LMAX, WMAX))
                node_ranges = [range(-1, EMIN - 1, -1)]
            else:
                chain_sets = itertools.product(chains_upto(3, 6), chains_upto(3, 6))
                node_ranges = [range(-1, -16, -1), range(-1, -16, -1)]
            for chs in chain_sets:
                if time.time() - t0 > time_limit:
                    stats['timeout'] = True
                    return found, stats
                for partial in itertools.product(*node_ranges):
                    for e_last in solve_last(list(partial), legs, list(chs), D, EMIN):
                        w, edges, nodes = build(list(partial) + [e_last], legs, list(chs))
                        stats['candidates'] += 1
                        if abs(det_of(w, edges)) != D:
                            continue
                        stats['det'] += 1
                        if not negdef(w, edges):
                            continue
                        stats['negdef'] += 1
                        try:
                            nm = name_of(kind, w, edges, nodes)
                        except Exception:
                            continue
                        stats['named'] += 1
                        if nm in names:
                            found.append((w, edges, nodes, nm))
                            if len(found) >= 3:
                                stats['seconds'] = round(time.time() - t0)
                                return found, stats
    stats['seconds'] = round(time.time() - t0)
    return found, stats
