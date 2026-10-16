#!/usr/bin/env python
"""Run the formula-based inverse search (invert.py) and the existing tests on the plumbings it finds.

Modes:
  controls : QA knots with Greene-pinned correction terms (the old validation set) and 13n1587's neighbours 11n66, 11n99:
             the plumbing d must equal Greene's d as functions (up to automorphism and sign).
  u12      : open knots with u in [1,2] whose plumbing search failed in graph_u1 (NO PLUMBING): u = 1 test.
  u23      : open knots with u in [2,3] or [2,4] whose search failed in owens_u3_graph (NO PLUMBING): Owens test.
Every found plumbing: linking form = +-(V+V^T)^{-1} (orientation), closed-form Casson-Walker identity s L = R.
Usage: run_invert.py MODE [workers]        Output: invert_MODE.json, invert_MODE.log
"""
import gzip, json, math, os, sys, time, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
TOOL = os.path.join(ROOT, 'plumbing')
sys.path.insert(0, TOOL)
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'montesinos'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
sys.path.insert(0, os.path.join(ROOT, 'bernhard_jablan', 'code'))
import _compat  # noqa
import numpy as np
import sympy
import invert as I
import general as G
import plumbing as P
import d_owens as DO
import montesinos_u1 as MU
import owens_obstruction as OO
from hf_red import hf_red_all
import u1_independent as U1
import owens_u3_graph as OW
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}


def orientation(Q, seifert):
    VV = (sympy.Matrix(seifert) + sympy.Matrix(seifert).T).tolist()
    tl, _ = G.linking_multiset(VV); tn = sorted((-x) % 1 for x in tl)
    lv, _ = G.linking_multiset(np.array(Q, dtype=np.int64).tolist())
    if tl == tn:
        return 0 if lv == tl else None
    return 1 if lv == tl else (-1 if lv == tn else None)


def u1_test(Q, D, w, edges):
    Q = np.array(Q, dtype=np.int64)
    deg = [0] * len(w)
    for u, v in edges:
        deg[u] += 1; deg[v] += 1
    nbad = sum(1 for i in range(len(w)) if deg[i] > -w[i])
    if nbad <= 1:
        d, cm = MU.correction_terms(Q)
    else:
        mq, cm = OO.m_Q(-np.array(Q, dtype=object)); d = {g: -v for g, v in mq.items()}
    adm = MU.u1_admissible(d, cm, D, Q)
    box = int(np.prod([abs(x) + 1 for x in w]))
    if adm is None:
        return 'OBSTRUCTED (H1 not cyclic)', nbad
    if not adm:
        return 'OBSTRUCTED', nbad
    if nbad <= 1 and box * D <= 2_000_000:
        ranks, _ = hf_red_all(Q)
        if ranks is not None:
            return ('PASS' if MU.u1_admissible(d, cm, D, Q, ranks_by_class=ranks) else 'OBSTRUCTED_HF'), nbad
    return 'PASS (d only)', nbad


def greene_d(names):
    out = {}
    for line in gzip.open(os.path.join(ROOT, 'data', 'greene_sweep_2026-09-08.jsonl.gz'), 'rt'):
        r = json.loads(line)
        if r['name'] in names and r.get('d'):
            D = int(r['det']); dd = r['d'] if isinstance(r['d'], dict) else json.loads(r['d'])
            out[r['name']] = {int(k) % D: Fraction(v) for k, v in dd.items()}
    for n in ('11n_66', '11n_99'):
        if n in names:
            g = json.load(open(os.path.join(ROOT, 'bernhard_jablan', 'results', f'greene_u1_{n}.json')))
            out[n] = {int(k) % g['det']: Fraction(v) for k, v in g['d'].items()}
    return out


def one(args):
    n, names, mode = args
    t0 = time.time(); D = int(BY[n]['determinant']); sigma = int(BY[n]['signature'])
    try:
        found, st = I.search(set(names), D, LMAX=5, WMAX=7, time_limit=900)
    except Exception as e:
        return n, dict(verdict='ERROR', error=repr(e)[:120])
    rec = dict(D=D, sigma=sigma, search=st, names=sorted(names)[:3])
    if not found:
        rec['verdict'] = 'NO PLUMBING'; rec['seconds'] = round(time.time() - t0)
        return n, rec
    w, edges, nodes, nm = found[0]
    Q = P.plumbing_matrix(w, edges)
    ori = orientation(Q, json.loads(BY[n]['seifert_matrix']))
    L = OW.L_closed(Q); R = OW.GU.R_of(n)
    cw = [s for s in (1, -1) if s * L == R]
    rec.update(plumbing=dict(weights=w, edges=edges, nodes=nodes, name=nm), link_orientation=ori, L=str(L), R=str(R),
               cw_signs=cw, cw_consistent=len(cw) == 1 and ori in (0, cw[0]))
    if mode in ('controls', 'u12'):
        v, nbad = u1_test(Q, D, w, edges); rec.update(u1_verdict=v, bad_vertices=nbad)
        Dp, dZ = DO.plumbing_d(w, edges)
        if v == 'OBSTRUCTED':
            rec['independent_ni_wu_survivors'] = len(U1.ni_wu_survivors(D, dZ))
        if mode == 'controls':
            g = greene_d({n}).get(n)
            if g:
                rec['d_equals_greene'] = any(all(dZ[(a * t) % D] == s * g[t] for t in range(D))
                                             for a in range(1, D) if math.gcd(a, D) == 1 for s in (1, -1))
    if mode == 'u23':
        rec.update(OW.run_one(n, Q, sigma, D))
    rec['seconds'] = round(time.time() - t0)
    return n, rec


def targets(mode):
    if mode == 'controls':
        c = [x[0] for x in json.load(open(os.path.join(TOOL, 'validation_candidates.json')))] + ['11n_66', '11n_99']
        out = []
        import graph_u1 as GU
        import spherogram
        for n in c:
            pd = json.loads(BY[n]['pd_notation']); mpd = [list(x) for x in spherogram.Link(pd).mirror().PD_code()]
            out.append((n, sorted(set(GU.names_orb(pd, 24)) | set(GU.names_orb(mpd, 24))), mode))
        return out
    if mode == 'u12':
        g = json.load(open(os.path.join(TOOL, 'graph_u1.json')))
        return [(n, sorted(set(r['names_sigma2']) | set(r['names_mirror'])), mode) for n, r in g.items() if r['verdict'] == 'NO PLUMBING']
    T = json.load(open(os.path.join(TOOL, 'owens_u3_targets.json')))
    out = []
    for n, r in T.items():
        if r and r.get('verdict') == 'NO PLUMBING':
            info = r.get('info', {})
            nm = sorted(set(info.get('names_sigma2', {})) | set(info.get('names_mirror', {})))
            if nm:
                out.append((n, nm, mode))
    return out


if __name__ == '__main__':
    from multiprocessing import Pool
    mode = sys.argv[1]; workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    items = targets(mode)
    outf = f'invert_{mode}.json'
    out = json.load(open(outf)) if os.path.exists(outf) else {}
    items = [x for x in items if x[0] not in out]
    print(mode, 'items', len(items), flush=True)
    T0 = time.time()
    with Pool(workers) as p:
        for n, rec in p.imap_unordered(one, items):
            out[n] = rec
            json.dump(out, open(outf, 'w'), indent=1, default=str)
            print(f"[{time.time()-T0:5.0f}s] {n:9s} D {rec.get('D')} | {rec.get('verdict', rec.get('u1_verdict'))} | "
                  f"plumbing {rec.get('plumbing', {}).get('weights')} | CW {rec.get('cw_signs')} link {rec.get('link_orientation')} "
                  f"| d=Greene {rec.get('d_equals_greene')} | NiWu {rec.get('independent_ni_wu_survivors')} | {rec.get('seconds')}s", flush=True)
    from collections import Counter
    print('summary', dict(Counter(r.get('verdict', r.get('u1_verdict')) for r in out.values())))
