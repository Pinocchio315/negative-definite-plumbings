#!/usr/bin/env python
"""Plumbing tree of the double branched cover of a knot given by an algebraic diagram (Lemma 2.7 of the paper).

Bigon contraction (basic_polyhedron.reduce_diagram) merges two vertex tangles A, B along a bigon; with the boundary
edges of each tangle listed counterclockwise, the merge is the tangle sum X + Z with X = A read from position i+2 and
Z = B read from position j (A[i], A[i+1] = B[j+1], B[j] are the bigon edges).  A tangle is stored as a rooted plumbing
tree together with a frame parity p: read from a position of parity p, its double cover is the plumbed manifold of the
tree with a regular fibre over the root removed, the lifts of the slopes 0 and infinity being a section curve and a
fibre of the boundary torus.  The spheres around the merged tangles cut the knot into crossings (rational tangles,
covered by solid tori) and hollow pieces between three spheres, covered by S^1 x (S^2 - 3 discs) with fibres of slope
infinity and sections of slope 0 (Bonahon-Siebenmann, Prop. A.1 and Section A.3).  Then
  crossing          : one vertex of weight SIGMA (a half-twisted band), parity 0;
  rotation by 90deg : a new root of weight 0 joined to the old root (parity changes);
  sum X + Z         : the roots are fused and their weights added;
  closure           : the numerator closure of the frame of the tree is the closed plumbing.
SIGMA = +1 or -1 gives the two orientations; which one is Sigma_2(K) is decided afterwards (Casson-Walker identity).
"""
import json, sys
from collections import Counter

SIGMA = 1


class Graph:
    def __init__(self):
        self.w, self.adj, self.n = {}, {}, 0

    def add(self, weight):
        v = self.n; self.n += 1
        self.w[v] = weight; self.adj[v] = set()
        return v

    def edge(self, a, b):
        assert a != b and b not in self.adj[a]
        self.adj[a].add(b); self.adj[b].add(a)

    def fuse(self, a, b):
        """Fuse vertex b into a (weights add); in a tree the two are never adjacent here."""
        assert a != b and b not in self.adj[a]
        self.w[a] += self.w.pop(b)
        for c in self.adj.pop(b):
            self.adj[c].discard(b); self.adj[c].add(a); self.adj[a].add(c)
        return a

    def export(self):
        idx = {v: i for i, v in enumerate(sorted(self.w))}
        w = [self.w[v] for v in sorted(self.w)]
        edges = sorted({tuple(sorted((idx[a], idx[b]))) for a in self.adj for b in self.adj[a]})
        return w, [list(e) for e in edges]


def plumbing_from_algebraic(pd, sigma=SIGMA):
    """Returns (weights, edges) of a plumbing tree whose boundary is +-Sigma_2 of the diagram, or None if the diagram
    does not reduce to one vertex."""
    G = Graph()
    nodes = {}
    for k, x in enumerate(pd):
        nodes[k] = dict(edges=list(x), root=G.add(sigma), parity=0)

    def rotate(t):
        r = G.add(0); G.edge(r, t['root']); t['root'] = r; t['parity'] ^= 1

    while len(nodes) > 1:
        ends = {}
        for v, t in nodes.items():
            for e in t['edges']:
                ends.setdefault(e, []).append(v)
        merged = False
        for a, A in nodes.items():
            ea = A['edges']
            for i in range(4):
                f, e = ea[i], ea[(i + 1) % 4]
                if f == e or len(ends.get(f, [])) != 2 or len(ends.get(e, [])) != 2:
                    continue
                bs = {w for w in ends[f] if w != a} & {w for w in ends[e] if w != a}
                for b in bs:
                    B = nodes[b]; eb = B['edges']
                    js = [j for j in range(4) if eb[j] == e and eb[(j + 1) % 4] == f]
                    if not js:
                        continue
                    j = js[0]
                    if A['parity'] != i % 2:
                        rotate(A)
                    if B['parity'] != j % 2:
                        rotate(B)
                    root = G.fuse(A['root'], B['root'])
                    nodes[a] = dict(edges=[ea[(i + 2) % 4], ea[(i + 3) % 4], eb[(j + 2) % 4], eb[(j + 3) % 4]],
                                    root=root, parity=0)
                    del nodes[b]; merged = True
                    break
                if merged:
                    break
            if merged:
                break
        if not merged:
            return None
    (F,) = nodes.values()
    e = F['edges']
    if e[0] == e[1] and e[2] == e[3]:
        need = 1          # joins NW-SW and SE-NE of frame 0: the numerator closure of frame 1
    elif e[1] == e[2] and e[3] == e[0]:
        need = 0
    else:
        raise ValueError(f'unexpected final vertex {e}')
    if F['parity'] != need:
        rotate(F)
    return G.export()


# ---------------------------------------------------------------- Neumann normal form (trees, genus 0)

def normal_form(w, edges):
    """Neumann's moves on a weighted tree: blow-down of +-1 vertices of degree <= 2, absorption of 0 vertices of
    degree 2, removal of a 0 leaf with its neighbour (neighbour of degree <= 2), and the string identity
    (a)-(c)-(b) = (a-1)-(-2)^(c-1)-(b-1) for c >= 2 (also at a leaf), until every vertex of degree <= 2 has weight
    <= -2.  Returns (weights, edges) or raises ValueError for reducible / S^1 x S^2 cases."""
    W = {i: x for i, x in enumerate(w)}
    A = {i: set() for i in range(len(w))}
    for a, b in edges:
        A[a].add(b); A[b].add(a)
    nxt = [len(w)]

    def new(x):
        v = nxt[0]; nxt[0] += 1; W[v] = x; A[v] = set(); return v

    def remove(v):
        for u in A.pop(v):
            A[u].discard(v)
        W.pop(v)

    changed = True
    while changed:
        changed = False
        if len(W) == 1:
            break
        for v in list(W):
            if v not in W:
                continue
            d = len(A[v]); x = W[v]
            if d <= 2 and x in (1, -1):                       # blow down
                nb = list(A[v])
                remove(v)
                for u in nb:
                    W[u] -= x
                if len(nb) == 2:
                    A[nb[0]].add(nb[1]); A[nb[1]].add(nb[0])
                changed = True; break
            if d == 2 and x == 0:                             # absorb a 0 in a string
                a, b = list(A[v])
                remove(v)
                W[a] += W.pop(b)
                for c in A.pop(b):
                    A[c].discard(b); A[c].add(a); A[a].add(c)
                changed = True; break
            if d == 1 and x == 0:                             # 0 leaf
                (u,) = A[v]
                if len(A[u]) > 2:
                    raise ValueError('reducible (0 leaf at a node)')
                remove(v); remove(u)
                changed = True; break
            if d == 0 and x == 0:
                raise ValueError('S1 x S2')
            if d <= 2 and x >= 2:                             # positive string vertex
                nb = list(A[v])
                remove(v)
                for u in nb:
                    W[u] -= 1
                chain = [new(-2) for _ in range(x - 1)]
                for p, q in zip(chain, chain[1:]):
                    A[p].add(q); A[q].add(p)
                if len(nb) >= 1:
                    A[nb[0]].add(chain[0]); A[chain[0]].add(nb[0])
                if len(nb) == 2:
                    A[nb[1]].add(chain[-1]); A[chain[-1]].add(nb[1])
                changed = True; break
    idx = {v: i for i, v in enumerate(sorted(W))}
    wn = [W[v] for v in sorted(W)]
    en = sorted({tuple(sorted((idx[a], idx[b]))) for a in A for b in A[a]})
    return wn, [list(e) for e in en]


def plumbing_matrix(w, edges):
    n = len(w)
    Q = [[0] * n for _ in range(n)]
    for i, x in enumerate(w):
        Q[i][i] = x
    for a, b in edges:
        Q[a][b] = Q[b][a] = 1
    return Q


def det(Q):
    import sympy
    return int(sympy.Matrix(Q).det())


if __name__ == '__main__':
    from database_knotinfo import link_list
    BY = {r['name']: r for r in link_list()[1:]}
    for n in sys.argv[1:]:
        pd = json.loads(BY[n]['pd_notation'])
        r = plumbing_from_algebraic(pd)
        if r is None:
            print(n, 'not algebraic'); continue
        w, edges = r
        wn, en = normal_form(w, edges)
        print(n, 'det', abs(det(plumbing_matrix(w, edges))), 'KnotInfo', BY[n]['determinant'], '| normal form', wn, en)
