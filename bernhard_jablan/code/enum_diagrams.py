"""Step 1: all 13-crossing diagrams (up to mirror) of a target knot, following
Brittenham-Hermiller Prop. 2.3: crossing-change subsets (fixing crossing 0) of one reduced
alternating reference diagram per 13-crossing alternating knot.  Determinant filter via the
signed Tait graph (Goeritz), then SnapPy identification."""
import json, sys, itertools, warnings
warnings.filterwarnings('ignore')
import numpy as np
from multiprocessing import Pool
from database_knotinfo import link_list

TARGET = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].endswith('.json') else '13n_1587'
by = {r['name']: r for r in link_list()[1:]}
TDET = int(by[TARGET]['determinant'])
NC = int(__import__('re').match(r'(\d+)', TARGET).group(1))
REFS = [n for n in by if n.startswith(f'{NC}a_')]

def faces_of(pd):
    """corner (k,i) = region between positions i and i+1 (ccw) at crossing k."""
    where = {}
    for k, row in enumerate(pd):
        for i, lab in enumerate(row):
            where.setdefault(lab, []).append((k, i))
    def nxt(k, i):
        lab = pd[k][(i + 1) % 4]
        a, b = where[lab]
        return b if a == (k, (i + 1) % 4) else a
    face = {}
    nf = 0
    for k in range(len(pd)):
        for i in range(4):
            if (k, i) in face: continue
            c = (k, i)
            while c not in face:
                face[c] = nf
                c = nxt(*c)
            nf += 1
    return face, nf

def tait(pd):
    face, nf = faces_of(pd)
    # 2-colour: corners (k,i),(k,i+1) differ
    adj = [[] for _ in range(nf)]
    for k in range(len(pd)):
        for i in range(4):
            u, v = face[(k, i)], face[(k, (i + 1) % 4)]
            adj[u].append(v); adj[v].append(u)
    col = [-1] * nf; col[0] = 0; st = [0]
    while st:
        u = st.pop()
        for v in adj[u]:
            if col[v] == -1: col[v] = 1 - col[u]; st.append(v)
            else: assert col[v] != col[u]
    cnt = [col.count(0), col.count(1)]
    c0 = 0 if cnt[0] <= cnt[1] else 1
    verts = [f for f in range(nf) if col[f] == c0]
    idx = {f: j for j, f in enumerate(verts)}
    edges = []
    for k in range(len(pd)):
        a, b = (face[(k, 0)], face[(k, 2)]) if col[face[(k, 0)]] == c0 else (face[(k, 1)], face[(k, 3)])
        edges.append((idx[a], idx[b]))
    return len(verts), edges

def is_alternating(pd):
    pos = {}
    for k, row in enumerate(pd):
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

MASKS = np.array([[0] + [(s >> j) & 1 for j in range(NC - 1)] for s in range(2 ** (NC - 1))], dtype=np.int64)

def work(ref):
    import snappy
    pd = json.loads(by[ref]['pd_notation'])
    assert len(pd) == NC and is_alternating(pd), ref
    nv, edges = tait(pd)
    m = nv - 1
    E = np.zeros((NC, m, m))
    for k, (u, v) in enumerate(edges):
        if u == v: continue
        vec = np.zeros(nv); vec[u] += 1; vec[v] -= 1
        E[k] = np.outer(vec, vec)[1:, 1:]
    signs = 1 - 2 * MASKS                       # base sign +1 for every alternating edge
    G = np.einsum('sk,kij->sij', signs.astype(float), E)
    dets = np.rint(np.abs(np.linalg.det(G))).astype(np.int64) if m > 0 else np.ones(2 ** (NC - 1), dtype=np.int64)
    out = []
    for s in np.nonzero(dets == TDET)[0]:
        S = [j + 1 for j in range(NC - 1) if (s >> j) & 1]
        newpd = change(pd, S)
        try:
            ids = [str(x) for x in snappy.Link(newpd).exterior().identify()]
        except Exception as e:
            ids = ['ERR:' + repr(e)]
        out.append(dict(ref=ref, subset=S, ids=ids))
    nerr = sum(1 for o in out if any(x.startswith('ERR') for x in o['ids']))
    nempty = sum(1 for o in out if not o['ids'])
    return ref, int((dets == TDET).sum()), out, nerr, nempty

if __name__ == '__main__':
    # sanity: determinant filter reproduces KnotInfo determinants on the references themselves
    tn = 'K' + TARGET.replace('_', '')
    print('target', TARGET, 'det', TDET, 'refs', len(REFS), flush=True)
    hits, ncand, NERR, NEMPTY = [], 0, 0, 0
    with Pool(11) as pool:
        for i, (ref, nc, out, ne, nm) in enumerate(pool.imap_unordered(work, REFS, chunksize=8)):
            ncand += nc; NERR += ne; NEMPTY += nm
            for o in out:
                if any(x.startswith(tn + '(') for x in o['ids']):
                    hits.append(o)
            if (i + 1) % 1000 == 0:
                print(i + 1, 'refs done; det candidates', ncand, 'target diagrams', len(hits), flush=True)
    print('identify errors', NERR, 'empty identifications (non-hyperbolic / unresolved)', NEMPTY)
    json.dump(dict(target=TARGET, det_candidates=ncand, id_errors=NERR, id_empty=NEMPTY, diagrams=hits), open(f'diagrams_{TARGET}.json', 'w'), indent=0)
    print('DONE det candidates', ncand, 'target diagrams', len(hits))
