"""General pipeline: knot -> Sigma_2 (SnapPy cover + filling) -> Regina recognition -> plumbing-tree search
-> correction terms (negative-definite plumbing tree, lattice formula) with orientation fixed by the linking form.

Supported Regina names:  'SFS [S2: (a1,b1) ...]' (star: one node, legs),
                         'SFS [D: ...] U/m SFS [D: ...], m = ...'        (GraphPair: two nodes),
                         'SFS [D: ...] U/m SFS [A: ...] U/n SFS [D: ...], ...'  (GraphTriple: three nodes in a row).
Shape -> trees: one node per Seifert piece, one leg per exceptional fibre (alpha, beta) with continued fraction of
alpha/beta or alpha/(alpha-beta) (both tried), consecutive nodes joined by an edge subdivided by a chain of
<= max_chain vertices of weight in chain_w; node weights in erange.  A tree is accepted when |det| = D, it is
negative definite, its linking form equals +-(V+V^T)^{-1}, and Regina recognises its boundary with the same name.
Orientation: Sigma_2(K) = boundary of a 4-manifold with intersection form V+V^T (Gordon-Litherland; this is the
convention in which d(spin) = -sigma/4 for alternating knots), so Sigma_2(K) = +Y iff Q^{-1} and (V+V^T)^{-1} give the
same linking form, else -Y (decidable here when -1 is a non-square for some prime power in D).
"""
import itertools, json, math, re, time, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
import numpy as np
import sympy
import plumbing as P
import d_owens as DO


def cf(alpha, beta):
    out = []
    while beta:
        a = -(-alpha // beta); out.append(a); alpha, beta = beta, a * beta - alpha
    return out


def parse_name(name):
    pieces = re.findall(r'SFS \[([^:]+): ([^\]]*)\]', name)
    out = []
    for base, fib in pieces:
        fibres = [tuple(int(x) for x in f.split(',')) for f in re.findall(r'\((-?\d+,-?\d+)\)', fib)]
        out.append((base.strip(), fibres))
    kind = 'star' if len(out) == 1 else ('pair' if len(out) == 2 else ('triple' if len(out) == 3 else None))
    if kind == 'star' and out[0][0] != 'S2':
        kind = None
    return kind, out


def leg_options(alpha, beta):
    b = beta % alpha
    opts = {tuple(-x for x in cf(alpha, b)), tuple(-x for x in cf(alpha, alpha - b))}
    return sorted(opts)


def trees(kind, pieces, erange, chain_w, max_chain):
    nodes = len(pieces)
    legsets = []
    for base, fibres in pieces:
        legsets.append(list(itertools.product(*[leg_options(a, b) for a, b in fibres])))
    chains = [()] + [tuple(c) for L in range(1, max_chain + 1) for c in itertools.product(chain_w, repeat=L)]
    nedges = nodes - 1
    for es in itertools.product(erange, repeat=nodes):
        for legs in itertools.product(*legsets):
            for chs in itertools.product(chains, repeat=nedges):
                w = list(es); edges = []
                def add(x):
                    w.append(x); return len(w) - 1
                for node, ls in enumerate(legs):
                    for leg in ls:
                        prev = node
                        for x in leg:
                            v = add(x); edges.append((prev, v)); prev = v
                for e in range(nedges):
                    prev = e
                    for x in chs[e]:
                        c = add(x); edges.append((prev, c)); prev = c
                    edges.append((prev, e + 1))
                yield w, edges


def linking_multiset(A):
    S = sympy.Matrix(A); D = abs(int(S.det())); Sinv = S.inv(); n = S.rows
    vals = []
    # enumerate the discriminant group via BFS over unit vectors (works for non-cyclic groups too)
    from itertools import product
    seen = {}
    frontier = [tuple([0] * n)]
    key = lambda x: tuple((Sinv * sympy.Matrix(x))[i] % 1 for i in range(n))
    seen[key(frontier[0])] = frontier[0]
    while frontier:
        nxt = []
        for x in frontier:
            for i in range(n):
                y = list(x); y[i] += 1; y = tuple(y)
                k = key(y)
                if k not in seen:
                    seen[k] = y; nxt.append(y)
        frontier = nxt
        if len(seen) > D:
            break
    for x in seen.values():
        v = sympy.Matrix(x)
        vals.append((v.T * Sinv * v)[0] % 1)
    return sorted(vals), D


def find_plumbing(name, D, seifert, erange=range(-6, 1), chain_w=(-2, -3, -4, -5), max_chain=2, max_regina=300, log=print,
                  alt_names=(), target_sigs=None):
    """alt_names: further Regina names of the same manifold (e.g. the name recognised for -Sigma_2 = Sigma_2(mirror));
    a plumbing is accepted if its boundary is recognised with any of them (orientation is fixed by the linking form)."""
    names = {name, *alt_names}
    kind, pieces = parse_name(name)
    if kind is None:
        return None, dict(error=f'unsupported name {name}')
    VV = (sympy.Matrix(seifert) + sympy.Matrix(seifert).T).tolist()
    tl, Dt = linking_multiset(VV)
    assert Dt == D
    tl_neg = sorted((-x) % 1 for x in tl)
    st = dict(total=0, det=0, negdef=0, linking=0, regina=0)
    t0 = time.time()
    for w, edges in trees(kind, pieces, erange, chain_w, max_chain):
        st['total'] += 1
        Q = P.plumbing_matrix(w, edges); Qf = Q.astype(float)
        if abs(int(round(np.linalg.det(Qf)))) != D or np.linalg.eigvalsh(Qf).max() > -1e-9:
            continue
        st['det'] += 1
        if abs(int(sympy.Matrix(Q.tolist()).det())) != D:
            continue
        st['negdef'] += 1
        lv, _ = linking_multiset(Q.tolist())
        orient = '+' if lv == tl else ('-' if lv == tl_neg else None)
        if orient is None:
            continue
        st['linking'] += 1
        if st['regina'] >= max_regina:
            break
        st['regina'] += 1
        nm, size, T = P.regina_name(P.plumbed_manifold(w, edges))
        how = 'name' if nm in names else None
        if how is None and target_sigs:
            if isosig_set(T) & target_sigs:
                how = 'isosig'
        if how:
            st['seconds'] = round(time.time() - t0)
            ambiguous = (tl == tl_neg)
            return dict(weights=w, edges=edges, orientation=orient, orientation_ambiguous=ambiguous, identified_by=how,
                        regina_name_of_plumbing=nm), st
    st['seconds'] = round(time.time() - t0)
    return None, st


def isosig_set(T0, heights=(0, 1, 2), seeds=4):
    """isoSigs of several simplifications of a closed triangulation (combinatorial isomorphism = homeomorphism)."""
    import regina, random
    out = set()
    for h in heights:
        for k in range(seeds if h else 1):
            T = regina.Triangulation3(T0)
            if k:
                rng = random.Random(k)
                for _ in range(rng.randint(1, 3)):
                    T.subdivide(); T.intelligentSimplify()
            T.intelligentSimplify()
            if h:
                T.simplifyExhaustive(h)
            out.add(T.isoSig())
    return out


def sigma2_triangulation(pd):
    import snappy, regina
    C = snappy.Link(pd).exterior().covers(2, cover_type='cyclic')[0]; C.dehn_fill((1, 0))
    T = regina.Triangulation3(C.filled_triangulation()._to_string()); T.intelligentSimplify()
    return T


def sigma2_name(pd, tries=12):
    import snappy, regina
    E = snappy.Link(pd).exterior()
    for t in range(tries):
        E2 = E.copy()
        if t:
            E2.randomize()
        C = E2.covers(2, cover_type='cyclic')[0]; C.dehn_fill((1, 0))
        if 'positive' in C.solution_type() or 'negative' in C.solution_type():
            if C.volume() > 0.5:
                return 'HYPERBOLIC', C
        T = regina.Triangulation3(C.filled_triangulation()._to_string()); T.intelligentSimplify()
        rec = regina.StandardTriangulation.recognise(T)
        if rec and rec.manifold():
            return rec.manifold().name(), C
    return None, None


def knot_d(name_knot, pd, D, seifert, **kw):
    """Return (d on Sigma_2(K) in the Gordon-Litherland orientation, info)."""
    import spherogram
    nm, C = sigma2_name(pd)
    if nm is None or nm == 'HYPERBOLIC':
        return None, dict(sigma2=nm)
    mpd = [list(x) for x in spherogram.Link(pd).mirror().PD_code()]
    nm_mirror, _ = sigma2_name(mpd)
    alts = tuple(x for x in (nm_mirror,) if x and x != 'HYPERBOLIC' and x != nm)
    tsig = isosig_set(sigma2_triangulation(pd), heights=(0, 1, 2, 3), seeds=6)
    found, st = find_plumbing(nm, D, seifert, alt_names=alts, target_sigs=tsig, **kw)
    if found is None and alts:      # the shape parsed from the mirror name may differ (e.g. a different leg normalisation)
        found, st = find_plumbing(alts[0], D, seifert, alt_names=(nm,), target_sigs=tsig, **kw)
    if found is None:
        return None, dict(sigma2=nm, search=st)
    Dp, dY = DO.plumbing_d(found['weights'], found['edges'])
    s = 1 if found['orientation'] == '+' else -1
    d = {t: s * v for t, v in dY.items()}
    return d, dict(sigma2=nm, sigma2_mirror=alts, plumbing=found, search=st)
