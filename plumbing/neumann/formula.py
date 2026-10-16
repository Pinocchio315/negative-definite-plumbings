"""Plumbing tree -> Regina name of its boundary by formula (no triangulation), calibrated against Regina's recognition
of triangulated plumbings (calibrate.py).

Star (one node): SFS over S^2 with a fibre (alpha, beta) per leg, alpha/beta = [b1, ..., br] (Hirzebruch-Jung continued
fraction of the leg weights -b1, ..., -br read from the node) and the obstruction fibre (1, e) for the node weight e;
SFSpace.reduce().  Pairs / triples: one bounded SFS per node (the node's legs and (1, e_node)), joined along tori; the
matrix of a chain is a product of an edge matrix E and vertex matrices V(c), conventions fixed by calibration.
"""
import itertools
from fractions import Fraction
import regina


def hj(ws):
    x = Fraction(-ws[-1])
    for w in reversed(ws[:-1]):
        x = Fraction(-w) - 1 / x
    return x.numerator, x.denominator


def adjacency(n, edges):
    adj = {i: [] for i in range(n)}
    for u, v in edges:
        adj[u].append(v); adj[v].append(u)
    return adj


def decompose(weights, edges, nodes):
    """For a tree whose nodes (vertices carrying Seifert pieces) are given in a row: legs of each node and the chain of
    weights between consecutive nodes."""
    adj = adjacency(len(weights), edges)
    def path(a, b):
        prev = {a: None}; st = [a]
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in prev:
                    prev[y] = x; st.append(y)
        p = [b]
        while p[-1] != a:
            p.append(prev[p[-1]])
        return p[::-1]
    chains = [path(nodes[i], nodes[i + 1])[1:-1] for i in range(len(nodes) - 1)]
    on_chain = {x for c in chains for x in c} | set(nodes)
    legs = {v: [] for v in nodes}
    for v in nodes:
        for s in adj[v]:
            if s in on_chain:
                continue
            leg, prev, cur = [], v, s
            while True:
                leg.append(weights[cur])
                nxt = [x for x in adj[cur] if x != prev]
                if not nxt:
                    break
                assert len(nxt) == 1, 'legs must be chains'
                prev, cur = cur, nxt[0]
            legs[v].append(leg)
    return legs, [[weights[x] for x in c] for c in chains]


def piece(e, legs, punctures, sb=1, se=1):
    S = regina.SFSpace(regina.SFSpace.bo1 if punctures else regina.SFSpace.o1, 0, punctures)
    for leg in legs:
        a, b = hj(leg)
        S.insertFibre(a, sb * b)
    S.insertFibre(1, se * e)
    return S


def mat_mul(A, B):
    return [[A[0][0] * B[0][0] + A[0][1] * B[1][0], A[0][0] * B[0][1] + A[0][1] * B[1][1]],
            [A[1][0] * B[0][0] + A[1][1] * B[1][0], A[1][0] * B[0][1] + A[1][1] * B[1][1]]]


EDGES = {'J': [[0, 1], [1, 0]], '-J': [[0, -1], [-1, 0]], 'R': [[0, 1], [-1, 0]], '-R': [[0, -1], [1, 0]]}
VERTS = {'a': lambda c: [[1, 0], [-c, -1]], 'b': lambda c: [[1, 0], [c, -1]], 'c': lambda c: [[-1, 0], [-c, 1]],
         'd': lambda c: [[-1, 0], [c, 1]]}


def chain_matrix(chain, E, V):
    M = EDGES[E]
    for c in chain:
        M = mat_mul(EDGES[E], mat_mul(VERTS[V](c), M))
    return M


def pair_name(eA, legsA, eB, legsB, chain, conv):
    E, V, sb, se, transpose = conv
    M = chain_matrix(chain, E, V)
    if transpose:
        M = [[M[0][0], M[1][0]], [M[0][1], M[1][1]]]
    g = regina.GraphPair(piece(eA, legsA, 1, sb, se), piece(eB, legsB, 1, sb, se), M[0][0], M[0][1], M[1][0], M[1][1])
    return g.name()


def star_name(e, legs):
    S = regina.SFSpace()
    for leg in legs:
        a, b = hj(leg); S.insertFibre(a, b)
    S.insertFibre(1, e)
    S.reduce()
    return S.name()


CONVENTIONS = list(itertools.product(EDGES, VERTS, (1, -1), (1, -1), (False, True)))
