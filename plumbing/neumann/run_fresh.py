#!/usr/bin/env python
"""Pipeline on the 108 knots newly recognised as graph manifolds (recognise_more_fresh.json):
Regina names (60 randomisations + mirror) -> Neumann normal form (precheck.one) -> if p+ = 0 for some orientation:
linking-form orientation, closed-form Casson-Walker identity, and
  u in [1,2]       : u = 1 test (v2 u1_admissible; graded-root HF_red when affordable); independent Ni-Wu on d-obstructions;
  u in [2,3]/[2,4] : Owens's test (owens_u3_graph.run_one: |sigma| = 4 -> (0,2); |sigma| = 2 -> (1,1) with Traczyk).
Usage: run_fresh.py [workers]          Output: run_fresh.json, run_fresh.log
"""
import json, os, sys, time, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, '..')
sys.path.insert(0, TOOL)
sys.path.insert(0, os.path.join(TOOL, '..', 'lib'))
import precheck as PC
import run_invert as RI
import owens_u3_graph as OW
import plumbing as P
import d_owens as DO
import u1_independent as U1
import graph_u1 as GU
from ranges import load_ranges
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
RANGES = load_ranges('baseline')          # open ranges before data/new_values_20261015.json


def one(args):
    n, names = args
    import spherogram
    t0 = time.time()
    D = int(BY[n]['determinant']); sigma = int(BY[n]['signature']); lo, hi = RANGES[n]['lo'], RANGES[n]['hi']
    pd = json.loads(BY[n]['pd_notation']); mpd = [list(x) for x in spherogram.Link(pd).mirror().PD_code()]
    names = sorted(set(names) | set(GU.names_orb(mpd, 16)))
    _, pre = PC.one((n, names))
    rec = dict(D=D, sigma=sigma, u=[lo, hi], kinds=pre['kinds'], applicability=pre['verdict'], pplus=pre['pplus'])
    if pre['verdict'] != 'APPLICABLE':
        rec['verdict'] = pre['verdict']; rec['seconds'] = round(time.time() - t0, 1)
        return n, rec
    o = [k for k, v in pre['pplus'].items() if v == 0][0]
    nf = pre['normal_forms'][o]; w, edges = nf['weights'], nf['edges']
    Q = P.plumbing_matrix(w, edges)
    ori = RI.orientation(Q, json.loads(BY[n]['seifert_matrix']))
    L = OW.L_closed(Q); R = OW.GU.R_of(n)
    cw = [s for s in (1, -1) if s * L == R]
    rec.update(plumbing=dict(weights=w, edges=edges), link_orientation=ori, L=str(L), R=str(R), cw_signs=cw,
               cw_consistent=len(cw) == 1 and ori in (0, cw[0]))
    if (lo, hi) == (1, 2):
        v, nbad = RI.u1_test(Q, D, w, edges)
        rec.update(test='u=1', verdict=v, bad_vertices=nbad)
        if v == 'OBSTRUCTED':
            Dp, dZ = DO.plumbing_d(w, edges)
            rec['independent_ni_wu_survivors'] = len(U1.ni_wu_survivors(D, dZ))
    elif lo == 2 and abs(sigma) in (2, 4):
        r2 = OW.run_one(n, Q, sigma, D)
        rec.update(test='owens', verdict=r2.get('verdict'), owens=r2)
    else:
        rec.update(test=None, verdict='NO TEST (sigma)')
    rec['seconds'] = round(time.time() - t0, 1)
    return n, rec


if __name__ == '__main__':
    from multiprocessing import Pool
    from collections import Counter
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    fresh = json.load(open(os.path.join(TOOL, 'recognise_more_fresh.json')))
    items = [(n, r['names']) for n, r in fresh.items()]
    out = {}
    T0 = time.time()
    with Pool(workers) as p:
        for n, rec in p.imap_unordered(one, items):
            out[n] = rec
            json.dump(out, open(os.path.join(HERE, 'run_fresh.json'), 'w'), indent=1, default=str)
            print(f"[{time.time()-T0:5.0f}s] {n:9s} u {rec['u']} sigma {rec['sigma']:+d} D {rec['D']:4d} {rec['kinds']} | "
                  f"{rec['applicability']:14s} | {rec.get('verdict')} | CW {rec.get('cw_signs')} link {rec.get('link_orientation')} "
                  f"| NiWu {rec.get('independent_ni_wu_survivors')} | {rec['seconds']}s", flush=True)
    print('applicability', dict(Counter(r['applicability'] for r in out.values())))
    print('verdicts', dict(Counter(r['verdict'] for r in out.values())))
    print('OBSTRUCTED:', sorted(n for n, r in out.items() if str(r['verdict']).startswith('OBSTRUCTED')))
