#!/usr/bin/env python
"""Does Sigma_2 (or -Sigma_2) bound a negative-definite plumbing tree?  Decided from the Neumann normal form.

Neumann, On bilinear forms represented by trees, Bull. Austral. Math. Soc. 40 (1989), Theorem 1.2 (with the plumbing
calculus of Trans. AMS 268 (1981), Thms 4.1-4.2): among all plumbing trees of an oriented graph manifold the minimal
positive index p+ is attained, and the unique smallest tree attaining it is the one in normal form, i.e. every vertex
of valency <= 2 has weight <= -2.  Hence  Y bounds a negative-definite plumbing tree  <=>  p+(normal form of Y) = 0.

Normal form from a Regina name (formula.py conventions, calibrated):
  legs   : Hirzebruch-Jung expansions (weights <= -2) of the fibres alpha/(beta mod alpha); for -Y, of alpha/(-beta mod alpha);
  chains : the fibre intersection number p of two adjacent pieces is invariant (|m01| of the GraphPair matrix; |m01|, |n01|
           for a GraphTriple, the ends possibly swapped); the chain is the HJ expansion of p/q for a unit q mod p;
  nodes  : any integer weights; the last one is solved from |det| = D (det is affine in each weight).
A realisation is accepted when its formula name is one of the Regina names of Sigma_2 / -Sigma_2.
Verdict per knot: APPLICABLE (a negative-definite normal form exists for some orientation; it is returned),
NOT APPLICABLE (normal forms found for both orientations, all with p+ >= 1), INCONCLUSIVE (otherwise).
Usage: precheck.py u12|u23|controls [workers]          Output: precheck_MODE.json
"""
import itertools, json, math, os, re, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..'))
import formula as F
import invert as I
import general as G
import plumbing as P
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
R = 40


def mats(name):
    return [[[int(a), int(b)], [int(c), int(d)]] for a, b, c, d in re.findall(r'\[ (-?\d+),(-?\d+) \| (-?\d+),(-?\d+) \]', name)]


def chains_for(p):
    if p == 1:
        return [[]]
    return [[-x for x in G.cf(p, q)] for q in range(1, p) if math.gcd(p, q) == 1]


def index(w, edges):
    ev = np.linalg.eigvalsh(P.plumbing_matrix(w, edges).astype(float))
    return int((ev > 1e-9).sum()), int((ev < -1e-9).sum())


def realisations(name, names, D):
    """All normal-form trees (both orientations) whose formula name lies in `names`: list of (orient, w, edges, nodes, p+)."""
    kind, pieces = G.parse_name(name)
    out = []
    if kind not in ('star', 'pair', 'triple'):
        return kind, out
    for orient in itertools.product((1, -1), repeat=len(pieces)):   # per-piece orientation (det +1 gluings need mixed)
        legs = [[I.leg_of(a, o * b) for a, b in f if a > 1] for o, (_, f) in zip(orient, pieces)]
        if kind == 'star':
            chain_sets, ranges = [[]], []
        elif kind == 'pair':
            chain_sets = [[c] for c in chains_for(abs(mats(name)[0][0][1]))]
            ranges = [range(-R, R + 1)]
        else:
            m, n = mats(name)
            pm, pn = abs(m[0][1]), abs(n[0][1])
            chain_sets = [list(c) for c in itertools.product(chains_for(pm), chains_for(pn))]
            chain_sets += [list(c) for c in itertools.product(chains_for(pn), chains_for(pm))]
            ranges = [range(-R, R + 1), range(-R, R + 1)]
        for chs in chain_sets:
            for partial in itertools.product(*ranges):
                w0, e0, n0 = I.build(list(partial) + [-1], legs, chs)
                w1, e1, n1 = I.build(list(partial) + [-2], legs, chs)
                d0, d1 = I.det_of(w0, e0), I.det_of(w1, e1)
                a = d0 - d1
                for target in (D, -D):
                    if a == 0 or (target - d0) % a:
                        continue
                    e_last = (target - d0) // a - 1
                    if abs(e_last) > 3 * R:
                        continue
                    w, edges, nodes = I.build(list(partial) + [e_last], legs, chs)
                    try:
                        nm = I.name_of(kind, w, edges, nodes)
                    except Exception:
                        continue
                    if nm in names:
                        out.append((orient, w, edges, nodes, index(w, edges)[0]))
    return kind, out


def one(args):
    n, names = args
    t0 = time.time(); D = int(BY[n]['determinant'])
    names = set(names)
    found = {}
    kinds = set()
    for nm in sorted(names):
        kind, rs = realisations(nm, names, D)
        kinds.add(kind)
        for orient, w, edges, nodes, pplus in rs:
            found.setdefault(orient, []).append((pplus, w, edges, nodes))
    # the normal forms of Y and -Y come from opposite orientation vectors s and -s
    best = {o: min(v, key=lambda x: (x[0], len(x[1]))) for o, v in found.items()}
    neg = lambda o: tuple(-x for x in o)
    if any(b[0] == 0 for b in best.values()):
        verdict = 'APPLICABLE'
    elif any(neg(o) in best for o in best):
        verdict = 'NOT APPLICABLE'
    else:
        verdict = 'INCONCLUSIVE'
    key = lambda o: ''.join('+' if x > 0 else '-' for x in o)
    rec = dict(D=D, kinds=sorted(k for k in kinds if k), verdict=verdict,
               pplus={key(o): b[0] for o, b in best.items()},
               normal_forms={key(o): dict(weights=b[1], edges=b[2], nodes=b[3]) for o, b in best.items()},
               seconds=round(time.time() - t0, 1))
    return n, rec


def targets(mode):
    TOOL = os.path.join(HERE, '..')
    if mode == 'controls':
        c = json.load(open(os.path.join(HERE, 'invert_controls.json')))
        import graph_u1 as GU, spherogram
        out = []
        for n in c:
            pd = json.loads(BY[n]['pd_notation']); mpd = [list(x) for x in spherogram.Link(pd).mirror().PD_code()]
            out.append((n, sorted(set(GU.names_orb(pd, 24)) | set(GU.names_orb(mpd, 24)))))
        return out
    if mode == 'u12':
        g = json.load(open(os.path.join(TOOL, 'graph_u1.json')))
        return [(n, sorted(set(r['names_sigma2']) | set(r['names_mirror']))) for n, r in g.items() if r['verdict'] == 'NO PLUMBING']
    T = json.load(open(os.path.join(TOOL, 'owens_u3_targets.json')))
    out = []
    for n, r in T.items():
        if r and r.get('verdict') == 'NO PLUMBING':
            info = r.get('info', {})
            nm = sorted(set(info.get('names_sigma2', {})) | set(info.get('names_mirror', {})))
            if nm:
                out.append((n, nm))
    return out


if __name__ == '__main__':
    from multiprocessing import Pool
    from collections import Counter
    mode = sys.argv[1]; workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    items = targets(mode)
    only = set(sys.argv[3].split(',')) if len(sys.argv) > 3 else None
    if only:
        items = [x for x in items if x[0] in only]
    out = {}
    t0 = time.time()
    with Pool(workers) as p:
        for n, rec in p.imap_unordered(one, items):
            out[n] = rec
            print(f"[{time.time()-t0:5.0f}s] {n:9s} D {rec['D']:4d} | {rec['kinds']} | {rec['verdict']:14s} | p+ by orientation "
                  f"{rec['pplus']} | {rec['seconds']}s", flush=True)
    outf = os.path.join(HERE, f'precheck_{mode}.json')
    prev = json.load(open(outf)) if os.path.exists(outf) else {}
    prev.update(out)
    json.dump(prev, open(outf, 'w'), indent=1, default=str)
    out = prev
    print(mode, 'verdicts', dict(Counter(r['verdict'] for r in out.values())))
