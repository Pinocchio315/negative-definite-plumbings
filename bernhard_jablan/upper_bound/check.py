#!/usr/bin/env python
"""Upper-bound certificates for 13n1587 and checks of the correction terms of its neighbours (Section 5).

(a) Upper bound u(13n1587) <= 2 without identifying an intermediate knot: the two-crossing unknotting pairs of hits.jsonl,
    found by a randomised Reidemeister search (mode A).  Diagram = 13n1587: SnapPy isometry signature and Regina isoSig (complement,
    simplification, Pachner search of height <= 2) equal to those of KnotInfo's PD (Gordon-Luecke).  Unknot: global
    simplification to 0 crossings (Reidemeister moves), HFK total rank 1, pi_1 = Z.  The single changes to 10_113
    (mode B) are re-identified the same way.
(b) Casson-Walker identity on the Greene-model correction terms of the name-dependent neighbours 11n66, 11n99, 12n563
    (the README's remaining dependency): sum_t d(t) = -D sigma/4 + D V'(-1)/(6 V(-1)) (L-space: HF_red = 0).
(c) Sigma_2 of 11n66, 11n99, 12n563 by the orbifold route; when Regina recognises a Seifert/graph manifold, the
    negative-definite plumbing gives d by the lattice formula (code independent of the Greene model); compare with the
    Greene d (as functions up to automorphism and sign) and rerun the independent Ni-Wu test.
Output: check.json
"""
import gzip, json, math, os, sys, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
CEN = ROOT
sys.path.insert(0, os.path.join(CEN, 'lib'))
sys.path.insert(0, os.path.join(CEN, 'plumbing'))
sys.path.insert(0, os.path.join(ROOT, 'bernhard_jablan', 'code'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
import spherogram, snappy, regina
from exact_helpers import change
from knot_floer_homology import pd_to_hfk
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}


def snappy_sig(pd):
    try:
        return snappy.Link(pd).exterior().isometry_signature()
    except Exception as e:
        return 'ERR ' + repr(e)[:60]


def regina_sigs(pd, height=2):
    if min(min(X) for X in pd) == 0:                     # spherogram PDs are 0-based; Regina wants strands 1..2n
        pd = [[x + 1 for x in X] for X in pd]
    T = regina.Link.fromPD(pd).complement(); T.simplify()
    sigs = {T.isoSig()}
    T2 = regina.Triangulation3(T); T2.simplifyExhaustive(height); sigs.add(T2.isoSig())
    return sigs


def unknot_evidence(pd):
    L = spherogram.Link(pd); L.simplify('global')
    k = len(L.crossings)
    rank = None
    if k:
        try:
            rank = pd_to_hfk([list(x) for x in L.PD_code()])['total_rank']
        except ValueError:
            rank = 'hfk-error'
    G = snappy.Link(pd).exterior().fundamental_group()
    return dict(simplified_crossings=k, hfk_rank=1 if k == 0 else rank,
                pi1_generators=G.num_generators(), pi1_relators=len(G.relators()))


# (a)
K0 = json.loads(BY['13n_1587']['pd_notation'])
s0, r0 = snappy_sig(K0), regina_sigs(K0)
hits = [json.loads(l) for l in open(os.path.join(HERE, 'hits.jsonl'))]
rec_a = []
for h in hits:
    if h['knot'] != '13n_1587':
        continue
    pd = h['pd']
    same = dict(snappy=snappy_sig(pd) == s0, regina=bool(regina_sigs(pd) & r0))
    if h['mode'] == 'A':
        ev = unknot_evidence(change(pd, h['pair']))
        rec_a.append(dict(mode='A', crossings=len(pd), pair=h['pair'], signs=h['signs'], diagram_is_13n1587=same, unknot=ev))
    else:
        q = change(pd, [h['crossing']])
        t = json.loads(BY[h['target']]['pd_notation'])
        rec_a.append(dict(mode='B', crossings=len(pd), crossing=h['crossing'], sign=h['sign'], diagram_is_13n1587=same,
                          target=h['target'], target_snappy=snappy_sig(q) == snappy_sig(t),
                          target_regina=bool(regina_sigs(q) & regina_sigs(t))))
    print('(a)', rec_a[-1], flush=True)
out['a_upper_bound'] = rec_a

# (b)
import kinfo
def R_of(name):
    J = kinfo.parse_jones(name); s = int(BY[name]['signature']); D = int(BY[name]['determinant'])
    V1 = sum(c * (-1) ** (e % 2) for e, c in J.items()); dV1 = sum(e * c * (-1) ** ((e - 1) % 2) for e, c in J.items())
    assert abs(V1) == D
    return Fraction(-D * s, 4) + Fraction(D * dV1, 6 * V1)

greene = {}
for n in ('11n_66', '11n_99'):
    g = json.load(open(os.path.join(ROOT, 'bernhard_jablan', 'results', f'greene_u1_{n}.json')))
    greene[n] = {int(k) % g['det']: Fraction(v) for k, v in g['d'].items()}
for line in gzip.open(os.path.join(ROOT, 'data', 'greene_sweep_2026-09-08.jsonl.gz'), 'rt'):
    r = json.loads(line)
    if r['name'] == '12n_563' and r.get('d'):
        D = int(r['det']); dd = r['d'] if isinstance(r['d'], dict) else json.loads(r['d'])
        greene['12n_563'] = {int(k) % D: Fraction(v) for k, v in dd.items()}
rec_b = {}
for n, d in greene.items():
    D = int(BY[n]['determinant']); sig = int(BY[n]['signature'])
    S, R = sum(d.values()), R_of(n)
    rec_b[n] = dict(D=D, classes=len(d), d_spin=str(d[0]), minus_sigma_over_4=str(Fraction(-sig, 4)), sum_d=str(S), R=str(R),
                    identity_holds=(S == R))
    print('(b)', n, rec_b[n], flush=True)
out['b_casson_walker_on_greene_d'] = rec_b

# (c)
import graph_u1 as GU
import general as G
import plumbing as P
import d_owens as DO
import u1_independent as U1
GU.SHAPES['star'] = ('S2',)
GU.SEARCH['star'] = dict(erange=range(-8, 2), chain_w=(-2,), max_chain=0)
rec_c = {}
for n in ('11n_66', '11n_99', '12n_563'):
    pd = json.loads(BY[n]['pd_notation']); D = int(BY[n]['determinant'])
    mpd = [list(x) for x in spherogram.Link(pd).mirror().PD_code()]
    nK, nM = GU.names_orb(pd, 24), GU.names_orb(mpd, 24)
    allnames = set(nK) | set(nM)
    rc = dict(names=sorted(allnames)[:4])
    found = None
    for nm in sorted(allnames):
        kind, pieces = G.parse_name(nm)
        if kind in GU.SHAPES and tuple(b for b, _ in pieces) == GU.SHAPES[kind]:
            f, st = G.find_plumbing(nm, D, json.loads(BY[n]['seifert_matrix']), alt_names=tuple(allnames - {nm}),
                                    log=lambda *a: None, **GU.SEARCH[kind])
            rc.setdefault('searches', []).append(dict(shape=nm, stats=st))
            if f:
                found = dict(f, kind=kind); break
    if found:
        Dp, dP = DO.plumbing_d(found['weights'], found['edges'])
        g = greene.get(n)
        match = None
        if g:
            match = any(all(dP[(a * t) % D] == s * g[t] for t in range(D)) for a in range(1, D) if math.gcd(a, D) == 1 for s in (1, -1))
        surv = U1.ni_wu_survivors(D, dP)
        rc.update(plumbing=found, plumbing_d_equals_greene_d=match, independent_ni_wu_survivors_on_plumbing_d=len(surv))
    print('(c)', n, {k: v for k, v in rc.items() if k != 'plumbing'}, flush=True)
    rec_c[n] = rc
out['c_plumbing_route'] = rec_c
json.dump(out, open(os.path.join(HERE, 'check.json'), 'w'), indent=1, default=str)
