"""Exact helpers copied verbatim from bernhard_jablan/code/verify_completeness.py (reviewed)."""

import numpy as np



def faces_of(pd):
    where = {}
    for k, row in enumerate(pd):
        for i, lab in enumerate(row):
            where.setdefault(lab, []).append((k, i))
    def nxt(k, i):
        lab = pd[k][(i + 1) % 4]
        a, b = where[lab]
        return b if a == (k, (i + 1) % 4) else a
    face, nf = {}, 0
    for k in range(len(pd)):
        for i in range(4):
            if (k, i) in face:
                continue
            c = (k, i)
            while c not in face:
                face[c] = nf
                c = nxt(*c)
            nf += 1
    return face, nf



def tait(pd):
    face, nf = faces_of(pd)
    adj = [[] for _ in range(nf)]
    for k in range(len(pd)):
        for i in range(4):
            u, v = face[(k, i)], face[(k, (i + 1) % 4)]
            adj[u].append(v); adj[v].append(u)
    col = [-1] * nf; col[0] = 0; st = [0]
    while st:
        u = st.pop()
        for v in adj[u]:
            if col[v] == -1:
                col[v] = 1 - col[u]; st.append(v)
            else:
                assert col[v] != col[u], 'checkerboard colouring failed'
    c0 = 0 if col.count(0) <= col.count(1) else 1
    verts = [f for f in range(nf) if col[f] == c0]
    idx = {f: j for j, f in enumerate(verts)}
    edges = []
    for k in range(len(pd)):
        a, b = (face[(k, 0)], face[(k, 2)]) if col[face[(k, 0)]] == c0 else (face[(k, 1)], face[(k, 3)])
        edges.append((idx[a], idx[b]))
    return len(verts), edges



def has_bridge(nv, edges):
    """True if some edge is a bridge (removing it disconnects the multigraph)."""
    for e in range(len(edges)):
        adj = [[] for _ in range(nv)]
        for f, (u, v) in enumerate(edges):
            if f != e:
                adj[u].append(v); adj[v].append(u)
        seen = {0}; st = [0]
        while st:
            u = st.pop()
            for v in adj[u]:
                if v not in seen:
                    seen.add(v); st.append(v)
        if len(seen) < nv:
            return True
    return False



def is_alternating(pd):
    pos = {}
    for row in pd:
        for i, lab in enumerate(row):
            pos.setdefault(lab, []).append(i % 2)
    return all(sorted(v) == [0, 1] for v in pos.values())



def change(pd, idxs):
    n2 = 2 * len(pd); new = [list(r) for r in pd]
    for i in idxs:
        a, b, c, d = pd[i]
        assert (c - a) % n2 == 1
        new[i] = [b, c, d, a] if (d - b) % n2 == 1 else [d, a, b, c]
    return new



def exact_dets(G):
    """Exact determinants of a stack of integer matrices (N, m, m) by fraction-free Bareiss."""
    A = G.astype(np.int64).copy()
    N, m, _ = A.shape
    if m == 0:
        return np.ones(N, dtype=np.int64)
    LIM = 1 << 31
    sign = np.ones(N, dtype=np.int64); prev = np.ones(N, dtype=np.int64)
    singular = np.zeros(N, dtype=bool); idx = np.arange(N)
    for k in range(m):
        col = A[:, k:, k] != 0
        has = col.any(axis=1)
        singular |= ~has
        piv = np.argmax(col, axis=1) + k
        swap = has & (piv != k)
        rk = A[idx, k, :].copy(); rp = A[idx, piv, :].copy()
        A[idx, k, :] = np.where(swap[:, None], rp, rk)
        A[idx, piv, :] = np.where(swap[:, None], rk, rp)
        sign = np.where(swap, -sign, sign)
        pk = A[:, k, k].copy()
        div = np.where(has & (prev != 0), prev, 1)
        for i in range(k + 1, m):
            aik = A[:, i, k].copy()
            for j in range(k + 1, m):
                num = A[:, i, j] * pk - aik * A[:, k, j]
                ok = has & ~singular
                assert np.all(num[ok] % div[ok] == 0), 'Bareiss division not exact'
                A[:, i, j] = np.where(ok, num // div, A[:, i, j])
        assert np.abs(A).max() < LIM, 'overflow guard'
        prev = np.where(has, pk, prev)
    det = sign * A[:, m - 1, m - 1]
    det[singular] = 0
    return det

