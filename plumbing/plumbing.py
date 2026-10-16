"""Plumbing trees -> closed 3-manifolds (SnapPy surgery on the plumbing link) -> Regina recognition.

A plumbing tree is given by weights w[v] and a list of edges (u, v).  The plumbing link has one unknot per vertex,
and for each edge the child circle is a meridian of the parent circle (a Hopf clasp); the plumbed 3-manifold is
integral surgery on this link with framing w[v] on component v (the boundary of the plumbed 4-manifold, whose
intersection form is Q = adjacency matrix with the weights on the diagonal).

The PD code is produced from an explicit planar layout of circles (radial layout, radii shrinking with depth);
crossings follow the KnotTheory convention X[i,j,k,l] (i = incoming under-arc, counterclockwise).
"""
import math, random
from fractions import Fraction
import numpy as np


def plumbing_matrix(weights, edges):
    n = len(weights)
    Q = np.zeros((n, n), dtype=object)
    for v, w in enumerate(weights):
        Q[v, v] = w
    for u, v in edges:
        Q[u, v] = Q[v, u] = 1
    return Q


def _layout(n, edges, root=0):
    adj = {v: [] for v in range(n)}
    for u, v in edges:
        adj[u].append(v); adj[v].append(u)
    parent = {root: None}; order = [root]; i = 0
    while i < len(order):
        u = order[i]; i += 1
        for v in adj[u]:
            if v not in parent:
                parent[v] = u; order.append(v)
    assert len(order) == n, 'plumbing graph must be a connected tree'
    centre = {root: (0.0, 0.0)}; radius = {root: 10.0}; direction = {root: None}
    for u in order:
        kids = [v for v in adj[u] if parent.get(v) == u]
        if not kids:
            continue
        cu, ru = centre[u], radius[u]
        if direction[u] is None:
            angles = [2 * math.pi * (j + 0.25) / len(kids) for j in range(len(kids))]
        else:
            # spread the children over the half-circle facing away from the parent
            base = direction[u]
            angles = [base + math.pi * (j + 1) / (len(kids) + 1) - math.pi / 2 for j in range(len(kids))]
        span = (2 * math.pi / len(kids)) if direction[u] is None else (math.pi / (len(kids) + 1))
        rk = 0.35 * ru * min(1.0, span)
        for v, a in zip(kids, angles):
            centre[v] = (cu[0] + ru * math.cos(a), cu[1] + ru * math.sin(a))
            radius[v] = rk
            direction[v] = a
    return centre, radius, parent


def _intersections(c1, r1, c2, r2):
    dx, dy = c2[0] - c1[0], c2[1] - c1[1]
    d = math.hypot(dx, dy)
    if d >= r1 + r2 or d <= abs(r1 - r2):
        return []
    a = (r1 * r1 - r2 * r2 + d * d) / (2 * d)
    h = math.sqrt(max(r1 * r1 - a * a, 0.0))
    xm, ym = c1[0] + a * dx / d, c1[1] + a * dy / d
    return [(xm + h * dy / d, ym - h * dx / d), (xm - h * dy / d, ym + h * dx / d)]


def plumbing_link_pd(n, edges):
    """PD code and the PD labels of each component (vertex)."""
    centre, radius, parent = _layout(n, edges)
    pts = []          # (x, y, circle_a, circle_b, over_circle)
    for v in range(n):
        for u in range(v):
            P = _intersections(centre[u], radius[u], centre[v], radius[v])
            if not P:
                continue
            assert parent.get(v) == u or parent.get(u) == v, 'unwanted intersection in layout'
            # Hopf clasp: u over at the first point, v over at the second
            pts.append((P[0][0], P[0][1], u, v, u))
            pts.append((P[1][0], P[1][1], u, v, v))
    assert len(pts) == 2 * len(edges)
    # traverse each circle counterclockwise: arcs between consecutive intersection points
    on = {v: [] for v in range(n)}
    for idx, (x, y, a, b, o) in enumerate(pts):
        for c in (a, b):
            ang = math.atan2(y - centre[c][1], x - centre[c][0]) % (2 * math.pi)
            on[c].append((ang, idx))
    label = 0; arc_in = {}; arc_out = {}; comp_labels = {}
    for c in range(n):
        seq = sorted(on[c])
        if not seq:
            raise ValueError('isolated vertex: plumbing graph must have at least one edge')
        labs = list(range(label + 1, label + 1 + len(seq)))
        label += len(seq)
        comp_labels[c] = labs
        for m, (ang, idx) in enumerate(seq):
            arc_out[(idx, c)] = labs[m]                 # arc leaving point idx along circle c
            arc_in[(idx, c)] = labs[m - 1]              # arc arriving at idx (previous arc, cyclic)
    pd = []
    for idx, (x, y, a, b, o) in enumerate(pts):
        under = b if o == a else a
        over = o
        def tangent(c):
            cx, cy = centre[c]
            return (-(y - cy), (x - cx))                 # counterclockwise tangent
        tu, to = tangent(under), tangent(over)
        i, k = arc_in[(idx, under)], arc_out[(idx, under)]
        # ray of the incoming under-arc points along -tu; is +to counterclockwise from it (within pi)?
        cross = (-tu[0]) * to[1] - (-tu[1]) * to[0]
        if cross > 0:
            j, l = arc_out[(idx, over)], arc_in[(idx, over)]
        else:
            j, l = arc_in[(idx, over)], arc_out[(idx, over)]
        pd.append([i, j, k, l])
    return pd, comp_labels


def plumbed_manifold(weights, edges):
    """SnapPy closed manifold: integral surgery on the plumbing link with framings = weights."""
    import snappy, spherogram
    n = len(weights)
    if n == 1:
        return snappy.Manifold(f'L({abs(weights[0])},1)') if weights[0] else None
    pd, comp_labels = plumbing_link_pd(n, edges)
    L = spherogram.Link(pd)
    # map spherogram components to vertices through the PD labels
    order = []
    for comp in L.link_components:
        labs = set()
        for cs in comp:
            labs.add(pd[cs.crossing.label][cs.strand_index])
        v = [c for c, ls in comp_labels.items() if labs & set(ls)]
        assert len(v) == 1
        order.append(v[0])
    M = L.exterior()
    M.dehn_fill([(weights[v], 1) for v in order])
    return M


def regina_name(M, tries=6, seed=0):
    import regina
    F = M.filled_triangulation()
    rng = random.Random(seed)
    best = None
    for t in range(tries):
        T = regina.Triangulation3(F._to_string())
        T.intelligentSimplify()
        if t:
            for _ in range(rng.randint(1, 5)):
                T.subdivide() if False else None
            T.simplifyExhaustive(t) if t <= 2 else T.intelligentSimplify()
        rec = regina.StandardTriangulation.recognise(T)
        if rec is not None and rec.manifold() is not None:
            return rec.manifold().name(), T.size(), T
        if best is None or T.size() < best.size():
            best = T
    return None, best.size() if best else None, best
