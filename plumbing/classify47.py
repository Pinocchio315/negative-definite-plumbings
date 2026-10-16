#!/usr/bin/env python
"""Sigma_2 classes for the 47 non-alternating knots of data/new66_u3.txt (u = u_a = 3; previous range [2,3]), by the
same pipeline and the same functions as for the 1,103 knots of RECOGNITION_MAP.md, so that the classification of the
follow-up paper covers all 1,150 non-alternating knots with previous range [1,2], [2,3] or [2,4].
Stages (each writes classify47_<stage>.json; the existing result files are not touched):
  recognise : reach_orb.classify (4 randomisations), then recognise_more.one (60) and recognise_retry.one (300, mirror)
              for the covers that are still unrecognised
  hyp       : (Sage) verify_hyp.one on the unrecognised covers (drill and refill, verified hyperbolicity)
  jsj       : jsj.one on the uncertified ones (incompressible vertex tori, pieces), then pachner_recognise.one on those
              without a candidate hyperbolic piece
  piece     : (Sage) jsj_verify.one on the candidate hyperbolic pieces
  normal    : precheck.one (Neumann normal forms of Sigma_2 and -Sigma_2) on every recognised cover, with the names of
              16 further randomised recognitions of the knot and of its mirror, as in neumann/precheck_sigma0.py; for an
              APPLICABLE knot also the linking-form orientation, the Casson-Walker identity and Owens's u = 2 test as in
              neumann/run_fresh.py (a consistency check: these knots have u = 3 by the rank-two linking test)
  verify    : for every APPLICABLE knot, (a) independent identification: the negative-definite normal form is triangulated
              as surgery on the plumbing link and recognised by Regina, as in neumann/rederive_all.py; (b) the
              independent Owens check of owens_general_group.py (every form with the required parity, exact class
              minima, all group isomorphisms; for non-cyclic H_1)
  summary   : one status per knot, in the classes of make_tables.status
Usage: classify47.py recognise|jsj|normal|verify|summary [workers];   sage -python classify47.py hyp|piece [workers]
"""
import json, os, sys, time, warnings
warnings.filterwarnings('ignore')
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.chdir(HERE)
CEN = os.path.join(HERE, '..')
K47 = sorted(l.strip() for l in open(os.path.join(CEN, 'data', 'new66_u3.txt')) if 'n_' in l)
assert len(K47) == 47


def out(stage):
    return os.path.join(HERE, f'classify47_{stage}.json')


def load(stage):
    return json.load(open(out(stage))) if os.path.exists(out(stage)) else {}


def recognise(n):
    import reach_orb as RO, recognise_more as RM, recognise_retry as RR
    _, kind, nm, h = RO.classify(n)
    rec = dict(reach_kind=kind, reach_name=nm, H1=h, names={nm: 1} if nm else {})
    if not nm:
        _, m = RM.one((n, kind))
        rec['more'] = dict(new=m['new'], names=m['names'])
        rec['names'] = dict(m['names'])
        if not rec['names']:
            _, r = RR.one(n)
            rec['retry'] = r
            rec['names'] = dict(r['names'])
    return n, rec


def hyp(n):
    import verify_hyp as VH
    return VH.one(n)


def jsj(n):
    import jsj as JS, pachner_recognise as PR
    _, rec = JS.one(n)
    if not rec['verdict'].startswith('toroidal, a piece'):
        _, p = PR.one(n)
        rec['pachner'] = p
    return n, rec


def piece(args):
    import jsj_verify as JV
    return JV.one(args)


def normal(args):
    n, names = args
    import spherogram
    import graph_u1 as GU
    sys.path.insert(0, os.path.join(HERE, 'neumann'))
    import precheck as PC
    import run_invert as RI
    import owens_u3_graph as OW
    import plumbing as P
    BY = PC.BY
    pd = json.loads(BY[n]['pd_notation']); mpd = [list(y) for y in spherogram.Link(pd).mirror().PD_code()]
    allnames = set(names) | set(GU.names_orb(pd, 16)) | set(GU.names_orb(mpd, 16))
    _, rec = PC.one((n, sorted(allnames)))
    rec['names'] = sorted(allnames)
    if rec['verdict'] == 'APPLICABLE':
        D = int(BY[n]['determinant']); sigma = int(BY[n]['signature'])
        o = [k for k, v in rec['pplus'].items() if v == 0][0]
        nf = rec['normal_forms'][o]
        Q = P.plumbing_matrix(nf['weights'], nf['edges'])
        ori = RI.orientation(Q, json.loads(BY[n]['seifert_matrix']))
        L = OW.L_closed(Q); R = GU.R_of(n)
        cw = [s for s in (1, -1) if s * L == R]
        rec.update(sigma=sigma, negdef_orientation=o, link_orientation=ori, L=str(L), R=str(R), cw_signs=cw,
                   cw_consistent=len(cw) == 1 and ori in (0, cw[0]))
        if abs(sigma) in (2, 4):
            rec['owens'] = OW.run_one(n, Q, sigma, D)
            rec['test_verdict'] = rec['owens'].get('verdict')
        else:
            rec['test_verdict'] = 'NO TEST (sigma)'
    return n, rec


def verify(n):
    import numpy as np
    import graph_u1 as GU
    import plumbing as P
    import owens_general_group as OGG
    import owens_obstruction as OO
    import owens_independent as OI
    from owens_traczyk_8_6_4_1 import qt_ind, gram_reduce
    from owens_u3_verify import forms_all
    r = load('normal')[n]
    nf = r['normal_forms'][r['negdef_orientation']]
    w, edges, D, sigma = nf['weights'], nf['edges'], r['D'], r['sigma']
    nm, _, _ = GU.fast_regina_name(P.plumbed_manifold(w, edges), tries=12)
    rec = dict(plumbing_name=nm, plumbing_name_matches=nm in set(r['names']))
    # (b) as in the main block of owens_general_group.py, for the normal form of this knot
    n_even = 2 if abs(sigma) == 4 else 1
    Q = P.plumbing_matrix(w, edges)
    mP, cm = OO.m_Q(-np.array(Q, dtype=object))
    assert cm.order == D
    dY = {g: -v for g, v in mP.items()}
    elemsB, addB, zeroB, orderB = OGG.group_of(dY.keys(), list(cm.d))
    fl = forms_all(D, n_even, D // 2 + 2)
    fit = {1: [], -1: []}; checked = 0
    for g in fl:
        mQ, DQ, _ = OI.class_minima(gram_reduce(qt_ind(*g)))
        elemsA, addA, zeroA, orderA = OGG.group_of(mQ.keys(), [D] * 4)
        if sorted(orderA(x) for x in elemsA) != sorted(orderB(x) for x in elemsB):
            continue
        checked += 1
        basis, ords = OGG.find_basis(elemsA, addA, zeroA, orderA)
        for s in (1, -1):
            for phi in OGG.isomorphisms(basis, ords, addA, zeroA, elemsB, addB, zeroB, orderB):
                if all((mQ[x] - s * dY[phi[x]]) >= 0 and (mQ[x] - s * dY[phi[x]]).denominator == 1
                       and (mQ[x] - s * dY[phi[x]]).numerator % 2 == 0 for x in mQ):
                    fit[s].append(g); break
    rec.update(D=D, sigma=sigma, n_even=n_even, group=list(cm.d), forms=len(fl), forms_with_isomorphic_group=checked,
               fit_plus=fit[1][:4], fit_minus=fit[-1][:4], blocked=[s for s in (1, -1) if not fit[s]])
    return n, rec


def status(n, R, H, JJ, PV, N):
    """The classes of make_tables.status: hyperbolic, mixed, nodef, negdef, graph4 (or other)."""
    if n in H and H[n]['certified']:
        return 'hyperbolic'
    if n in PV and PV[n]['certified_hyperbolic_piece']:
        return 'mixed'
    if n in N:
        return {'APPLICABLE': 'negdef', 'NOT APPLICABLE': 'nodef'}.get(N[n]['verdict'], 'other')
    if n in JJ and JJ[n]['verdict'].startswith('toroidal, no geometric'):
        return 'graph4'
    return 'other'


def names_of(n, R, JJ):
    if R[n]['names']:
        return sorted(R[n]['names'])
    return sorted(JJ.get(n, {}).get('pachner', {}).get('names', []))


if __name__ == '__main__':
    from multiprocessing import Pool
    stage = sys.argv[1]; workers = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    t0 = time.time()
    sys.path.insert(0, os.path.join(CEN, 'bernhard_jablan', 'code'))
    sys.path.insert(0, os.path.join(CEN, 'lib'))
    sys.path.insert(0, os.path.join(CEN, 'lib', 'lower_bounds', 'owens'))
    R = load('recognise')
    if stage == 'recognise':
        items, f = K47, recognise
    elif stage == 'hyp':
        items, f = [n for n in K47 if not R[n]['names']], hyp
    elif stage == 'jsj':
        H = load('hyp')
        items, f = [n for n in K47 if not R[n]['names'] and not H[n]['certified']], jsj
    elif stage == 'piece':
        JJ = load('jsj')
        items = [(n, [p['ideal_isosig'] for t in r['incompressible_tori'] for p in t['pieces']
                      if 'positively' in (p.get('snappy_solution') or '')]) for n, r in JJ.items()]
        items, f = [x for x in items if x[1]], piece
    elif stage == 'normal':
        JJ = load('jsj')
        items, f = [(n, names_of(n, R, JJ)) for n in K47 if names_of(n, R, JJ)], normal
    elif stage == 'verify':
        N = load('normal')
        items, f = [n for n in K47 if N.get(n, {}).get('verdict') == 'APPLICABLE'], verify
    elif stage == 'summary':
        H, JJ, PV, N = load('hyp'), load('jsj'), load('piece'), load('normal')
        S = {n: status(n, R, H, JJ, PV, N) for n in K47}
        json.dump(dict(status=S, counts=dict(Counter(S.values())),
                       tests={n: N[n].get('test_verdict') for n in K47 if S[n] == 'negdef'},
                       cw_consistent={n: N[n].get('cw_consistent') for n in K47 if S[n] == 'negdef'}),
                  open(out('summary'), 'w'), indent=1)
        print(dict(Counter(S.values())))
        print('negdef tests:', Counter(N[n].get('test_verdict') for n in K47 if S[n] == 'negdef'))
        sys.exit()
    else:
        raise SystemExit(__doc__)
    res = load(stage)
    with Pool(workers) as p:
        for n, rec in p.imap_unordered(f, items):
            res[n] = rec
            json.dump(res, open(out(stage), 'w'), indent=1, default=str)
            print(f'[{time.time() - t0:5.0f}s] {stage} {n}: ' + json.dumps(
                {k: rec.get(k) for k in ('reach_kind', 'certified', 'verdict', 'pplus', 'test_verdict',
                                         'certified_hyperbolic_piece', 'plumbing_name_matches', 'forms',
                                         'forms_with_isomorphic_group', 'blocked') if k in rec}, default=str)[:200], flush=True)
    print(stage, len(res), f'{time.time() - t0:.0f}s')
