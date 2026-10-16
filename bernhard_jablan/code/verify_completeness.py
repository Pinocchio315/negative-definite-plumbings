"""Exact re-verification that the list of minimal diagrams of a target knot is complete.

Addresses three review points:
  * determinants are computed EXACTLY (vectorised fraction-free Bareiss in int64 with
    overflow guards); the float determinant used originally is compared and mismatches recorded;
  * every exact-determinant candidate that is not accepted as the target is EXCLUDED BY AN EXACT
    INVARIANT, independently of SnapPy: crossing number < c after simplification, or the
    Conway-normalised Alexander polynomial differs, or knot Floer homology differs up to mirror;
  * all compared values are stored.

References (BH Prop. 2.3): one reduced alternating diagram per prime alternating c-crossing knot.
For each reference we check: alternating, reduced (Tait graph has no loop and no bridge), and that
the exact determinant of the reference equals the KnotInfo determinant.

Usage: verify_completeness.py TARGET OUT.json
"""
import json, re, sys, warnings
from fractions import Fraction
warnings.filterwarnings('ignore')
import numpy as np
from multiprocessing import Pool
from database_knotinfo import link_list

TARGET = sys.argv[1]
OUT = sys.argv[2]
by = {r['name']: r for r in link_list()[1:]}
NC = int(re.match(r'(\d+)', TARGET).group(1))
REFS = sorted((n for n in by if n.startswith(f'{NC}a_')), key=lambda s: int(s.split('_')[1]))
TDET = int(by[TARGET]['determinant'])
EVAL_T = [2, 3, 5, 7, -2, -3, Fraction(1, 2) * 3]


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


def alex_values(pd):
    """Conway-normalised Alexander polynomial evaluated at EVAL_T: det(V - tV^T) / t^(n/2)."""
    import spherogram, sympy
    from sympy.polys.matrices import DomainMatrix
    V = spherogram.Link(pd).seifert_matrix(); n = len(V)
    if n == 0:
        return [Fraction(1)] * len(EVAL_T)
    out = []
    for t in EVAL_T:
        t = Fraction(t)
        M = DomainMatrix([[sympy.QQ(Fraction(V[i][j]) - t * V[j][i]) for j in range(n)] for i in range(n)], (n, n), sympy.QQ)
        d = M.det()
        out.append(Fraction(int(d.numerator), int(d.denominator)) / t ** (n // 2))
    one = DomainMatrix([[sympy.ZZ(V[i][j] - V[j][i]) for j in range(n)] for i in range(n)], (n, n), sympy.ZZ).det()
    assert int(one) == 1, 'det(V - V^T) != 1'
    return out


def hfk_key(pd, mirror=False):
    from knot_floer_homology import pd_to_hfk
    h = pd_to_hfk(pd)
    r = h['ranks'].items()
    if mirror:
        r = [((-a, -m), v) for (a, m), v in r]
    return [list(x) for x in sorted((list(k), v) for k, v in r)]


TPD = json.loads(by[TARGET]['pd_notation'])
T_ALEX = alex_values(TPD)
T_HFK = hfk_key(TPD)
T_HFK_M = hfk_key(TPD, mirror=True)
T_ID = 'K' + TARGET.replace('_', '')


def work(ref):
    import snappy, spherogram
    pd = json.loads(by[ref]['pd_notation'])
    rec = dict(ref=ref, alternating=is_alternating(pd), crossings=len(pd))
    nv, edges = tait(pd)
    rec['tait_loops'] = sum(1 for u, v in edges if u == v)
    rec['tait_bridge'] = has_bridge(nv, edges)
    m = nv - 1
    E = np.zeros((NC, m, m), dtype=np.int64)
    for k, (u, v) in enumerate(edges):
        if u == v:
            continue
        vec = np.zeros(nv, dtype=np.int64); vec[u] += 1; vec[v] -= 1
        E[k] = np.outer(vec, vec)[1:, 1:]
    masks = np.array([[0] + [(s >> j) & 1 for j in range(NC - 1)] for s in range(2 ** (NC - 1))], dtype=np.int64)
    G = np.einsum('sk,kij->sij', 1 - 2 * masks, E)
    dex = np.abs(exact_dets(G))
    dfl = np.rint(np.abs(np.linalg.det(G.astype(float)))).astype(np.int64) if m > 0 else np.ones(len(masks), dtype=np.int64)
    rec['float_exact_mismatches'] = int((dex != dfl).sum())
    rec['ref_det_exact'] = int(dex[0])
    rec['ref_det_knotinfo'] = int(by[ref]['determinant'])
    cands = []
    for s in np.nonzero(dex == TDET)[0]:
        S = [j + 1 for j in range(NC - 1) if (s >> j) & 1]
        newpd = change(pd, S)
        c = dict(subset=S, float_det=int(dfl[s]))
        L = spherogram.Link(newpd); L.simplify('global')
        c['simplified_crossings'] = len(L.crossings)
        if len(L.crossings) < NC:
            c['decision'] = 'excluded: crossing number < %d' % NC
            cands.append(c); continue
        av = alex_values(newpd)
        if av != T_ALEX:
            c['decision'] = 'excluded: Alexander polynomial differs'
            c['alex_first_diff'] = next(str(a) + ' vs ' + str(b) for a, b in zip(av, T_ALEX) if a != b)
            cands.append(c); continue
        try:
            c['snappy_identify'] = [str(x) for x in snappy.Link(newpd).exterior().identify()]
        except Exception as e:
            c['snappy_identify'] = ['ERR ' + repr(e)]
        hk = hfk_key(newpd)
        c['hfk'] = hk
        c['hfk_matches_target_up_to_mirror'] = hk in (T_HFK, T_HFK_M)
        if any(x.startswith(T_ID + '(') for x in c['snappy_identify']):
            c['decision'] = 'ACCEPTED: SnapPy identifies target; HFK %s' % ('agrees' if c['hfk_matches_target_up_to_mirror'] else 'DISAGREES')
        elif not c['hfk_matches_target_up_to_mirror']:
            c['decision'] = 'excluded: knot Floer homology differs (mirror included)'
        else:
            c['decision'] = 'UNRESOLVED: same Alexander and HFK as target but not identified as target'
        cands.append(c)
    rec['candidates'] = cands
    return rec


if __name__ == '__main__':
    with Pool(11) as pool:
        recs = list(pool.imap(work, REFS, chunksize=8))
    from collections import Counter
    dec = Counter(c['decision'] for r in recs for c in r['candidates'])
    summary = dict(
        target=TARGET, target_det=TDET, references=len(recs),
        all_references_alternating=all(r['alternating'] for r in recs),
        all_references_reduced=all(r['tait_loops'] == 0 and not r['tait_bridge'] for r in recs),
        reference_det_matches_knotinfo=sum(r['ref_det_exact'] == r['ref_det_knotinfo'] for r in recs),
        float_vs_exact_det_mismatches=sum(r['float_exact_mismatches'] for r in recs),
        exact_det_candidates=sum(len(r['candidates']) for r in recs),
        decisions=dict(dec),
        accepted=[dict(ref=r['ref'], subset=c['subset']) for r in recs for c in r['candidates'] if c['decision'].startswith('ACCEPTED')],
        target_alexander_values=[str(x) for x in T_ALEX],
    )
    json.dump(dict(summary=summary, references=recs), open(OUT, 'w'), indent=0)
    print(json.dumps(summary, indent=1))
