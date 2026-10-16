#!/usr/bin/env python
"""Owens's (1,1) test and Traczyk's (0,2) test for 8_6 # 4_1, the composite neighbour of 14a2539.

For the representative with sigma = +2, a sequence of two crossing changes has sign pattern
(p, n) = (1, 1) or (0, 2)  (v2, Sec. 3, 'Two crossing changes').
  (1,1)  n = sigma/2, so Owens's theorem applies: Sigma_2(K) bounds Qt = [[M, I], [I, 2I]], 2M - I = Q
         positive definite, det Q = 115, exactly one of m_1, m_2 even, and (eq:dmatch) holds for one
         group isomorphism.  If every form fails, (1,1) is impossible.
  (0,2)  not definite; v2 uses Traczyk's criterion (eq:lm), (eq:traczyk), effective when d_K = 2.

(1,1) is computed twice:
  A  v2 code (code/lower_bounds/owens/owens_obstruction.py): Goeritz m_G (side with d(spin) = -1/2),
     candidates(115, 1), m_Q, matching_exists.
  B  independent: d-values of L(23,7) # L(5,2) from the Ozsvath-Szabo lens-space recursion
     (8_6 = S(23,7), 4_1 = S(5,2) are two-bridge), exact class minima and unit search of
     bernhard_jablan/code/owens_independent.py, own enumeration of mixed-parity forms
     (reduced, and unreduced up to a <= D//2 + 2 as a consistency check).
Also: the rank-one (u = 1) tests of the summands, which explain the split form diag(5, 23).
Output: ../owens_traczyk_8_6_4_1.json
"""
import json, math, os, random, sys, time, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
from itertools import product

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, 'bernhard_jablan', 'code'))
sys.path.insert(0, os.path.join(ROOT, 'lib'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
import _compat  # noqa: F401  (v2 shim)
import spherogram, snappy, sympy
import owens_obstruction as OO          # A (v2)
import traczyk as TR                    # v2 Traczyk
import kinfo                            # v2 KnotInfo access
import owens_independent as OI          # B
from u1_independent import lens_d, ni_wu_survivors
from exact_helpers import change
from bf2 import charpoly_signature

D = 115
T0 = time.time()


def log(*a):
    print(f'[{time.time() - T0:7.1f}s]', *a, flush=True)


def isosig(pd):
    return snappy.Link(pd).exterior().isometry_signature()


def exact_sigma(pd):
    V = spherogram.Link(pd).seifert_matrix(); n = len(V)
    return charpoly_signature([[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]) if n else 0


def exact_det(pd):
    V = spherogram.Link(pd).seifert_matrix(); n = len(V)
    return abs(int(sympy.Matrix([[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]).det()))


# ------------------------------------------------------------------ Jones polynomial (own Kauffman bracket)
def jones(pd, a_smoothing):
    """V(t) as {exponent: coeff} from the state sum.  a_smoothing in ('ij|kl', 'il|jk') says which
    pairing of X[i,j,k,l] (counterclockwise from the incoming under-strand) is the A-smoothing;
    calibrated below against KnotInfo."""
    n = len(pd); n2 = 2 * n
    labels = sorted({x for X in pd for x in X})
    pos = {x: i for i, x in enumerate(labels)}
    w = sum(1 if (j - l) % n2 == 1 else -1 for (i, j, k, l) in pd)
    poly = {}                                   # A-exponent -> coeff of <D>
    for s in range(1 << n):
        parent = list(range(len(labels)))
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]; x = parent[x]
            return x
        def union(x, y):
            parent[find(pos[x])] = find(pos[y])
        na = 0
        for c, (i, j, k, l) in enumerate(pd):
            A = not (s >> c) & 1
            na += A
            if (a_smoothing == 'ij|kl') == A:
                union(i, j); union(k, l)
            else:
                union(i, l); union(j, k)
        loops = len({find(x) for x in range(len(labels))})
        # A^(na - nb) * (-A^2 - A^-2)^(loops - 1)
        base = {na - (n - na): 1}
        for _ in range(loops - 1):
            nb = {}
            for e, cf in base.items():
                nb[e + 2] = nb.get(e + 2, 0) - cf
                nb[e - 2] = nb.get(e - 2, 0) - cf
            base = nb
        for e, cf in base.items():
            poly[e] = poly.get(e, 0) + cf
    # V(t) = (-A^3)^(-w) <D>,  A = t^(-1/4)
    out = {}
    sgn = -1 if w % 2 else 1
    for e, cf in poly.items():
        if cf == 0:
            continue
        ea = e - 3 * w
        if ea % 4:
            return None                          # wrong smoothing convention (calibration)
        out[-ea // 4] = out.get(-ea // 4, 0) + sgn * cf
    return {e: c for e, c in out.items() if c}


def jmul(p, q):
    out = {}
    for e1, c1 in p.items():
        for e2, c2 in q.items():
            out[e1 + e2] = out.get(e1 + e2, 0) + c1 * c2
    return {e: c for e, c in out.items() if c}


def at_omega(J):
    """Exact V(omega), omega = e^{i pi/3}, as (a, b) meaning a + b*omega (omega^2 = omega - 1)."""
    a = b = 0
    for e, c in J.items():
        x, y = 1, 0                              # omega^0
        step = (lambda x, y: (-y, x + y)) if e >= 0 else (lambda x, y: (x + y, -x))
        for _ in range(abs(e) % 6):
            x, y = step(x, y)
        a += c * x; b += c * y
    return a, b


# ------------------------------------------------------------------ lens-space correction terms
def lens_vec(p, q):
    """t -> d(L(p,q), spin + t), homomorphic labelling with spin at 0 (p odd prime)."""
    vals = [lens_d(p, q, i) for i in range(p)]
    cs = [c for c in range(p) if all(vals[i] == vals[(c - i) % p] for i in range(p))]
    assert len(cs) == 1, cs                     # conjugation i -> c - i; unique for p prime
    i0 = cs[0] * pow(2, -1, p) % p
    return [vals[(i0 + t) % p] for t in range(p)]


def qt_ind(m1, m2, a):
    return [[m1, a, 1, 0], [a, m2, 0, 1], [1, 0, 2, 0], [0, 1, 0, 2]]


def gram_reduce(G):
    """Pairwise (Lagrange) reduction of a positive-definite Gram matrix by integral basis changes
    b_i -> b_i - k b_j; returns P^T G P.  The lattice, its characteristic covectors, discriminant group
    and m-function are unchanged; only the box search of class_minima becomes small."""
    G = [list(map(int, row)) for row in G]; n = len(G)
    changed = True
    while changed:
        changed = False
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                k = round(Fraction(G[i][j], G[j][j]))
                if k == 0:
                    continue
                new_ii = G[i][i] - 2 * k * G[i][j] + k * k * G[j][j]
                if new_ii >= G[i][i]:
                    continue
                for l in range(n):
                    if l != i:
                        G[i][l] -= k * G[j][l]
                        G[l][i] = G[i][l]
                G[i][i] = new_ii
                changed = True
    return G


def mixed_forms(Dk, amax, reduced_only):
    """(m1, m2, a): x = 2m1-1 <= y = 2m2-1 odd, {x, y} = {1, 3} mod 4, xy - 4a^2 = Dk, 0 <= a <= amax."""
    out = []
    for a in range(0, amax + 1):
        N = Dk + 4 * a * a
        for x in range(1, math.isqrt(N) + 1, 2):
            if N % x:
                continue
            y = N // x
            if {x % 4, y % 4} != {1, 3}:
                continue
            m1, m2 = (x + 1) // 2, (y + 1) // 2
            red = a < m1 <= m2
            if reduced_only and not red:
                continue
            out.append((m1, m2, a, red))
    return out


def main():
    rec = dict(knot='8_6#4_1', det=D)
    # ------------------------------------------------------------- 0. the neighbour and an alternating diagram
    k14 = json.load(open(os.path.join(OUTDIR, 'k14.json')))
    pd14 = k14['pd_14a2539']
    ref = {n: isosig(kinfo.pd_code(n)) for n in ('8_6', '4_1')}
    comp = []
    for c in range(len(pd14)):
        q = change(pd14, [c])
        random.seed(f'nb:{c}')
        L = spherogram.Link(q); L.simplify('global')
        parts = L.deconnect_sum() if len(L.crossings) >= 6 else [L]
        if len(parts) > 1:
            comp.append((c, q, L, parts))
    # crossings 3 and 4 of the PD form a bigon (shared edges 4, 16): both changes give the same knot
    assert [x[0] for x in comp] == [3, 4], [x[0] for x in comp]
    per_crossing = {}
    for c, q, L, prts in comp:
        nm_c = []
        for P in prts:
            random.seed('part'); P.simplify('global')
            s = isosig(P.PD_code())
            nm = [n for n, v in ref.items() if v == s]
            assert len(nm) == 1, 'summand not identified'
            nm_c.append(nm[0])
        per_crossing[c] = dict(summands=sorted(nm_c), sigma=exact_sigma(q))
    assert all(v == dict(summands=['4_1', '8_6'], sigma=-2) for v in per_crossing.values()), per_crossing
    c0, q0, L0, parts = comp[0]
    named = {}
    for P in parts:
        named[[n for n, v in ref.items() if v == isosig(P.PD_code())][0]] = P
    assert set(named) == {'8_6', '4_1'}
    P86, P41 = named['8_6'], named['4_1']
    assert len(P86.crossings) == 8 and len(P41.crossings) == 4 and P86.is_alternating() and P41.is_alternating()
    K = None
    for tag, other in (('4_1', P41), ('4_1 mirror diagram', P41.mirror())):
        S = P86.connected_sum(other)
        if S.is_alternating():
            K, joined = S, tag
            break
    assert K is not None and len(K.crossings) == 12
    pdK = [list(x) for x in K.PD_code()]
    sig_nb, sig_K = exact_sigma(q0), exact_sigma(pdK)
    sig_86_ki = exact_sigma(kinfo.pd_code('8_6'))
    det_K = exact_det(pdK)
    random.seed('re'); K2 = spherogram.Link(pdK); K2.simplify('global')
    re_parts = K2.deconnect_sum()
    re_ids = sorted([n for n, v in ref.items() if v == isosig(P.PD_code())][0] for P in re_parts)
    assert det_K == D and sig_K == sig_nb == -2 and re_ids == ['4_1', '8_6']
    rec['neighbour'] = dict(changed_crossing_of_14a2539=c0, composite_neighbour_crossings=per_crossing, neighbour_pd=q0, sigma_neighbour=sig_nb,
                            alternating_diagram_pd=pdK, built_as=f'8_6 summand # {joined} (alternating join)',
                            sigma=sig_K, det=det_K, summands_by_isometry_signature=re_ids,
                            sigma_KnotInfo_8_6_same_convention=sig_86_ki,
                            note='summands taken from the neighbour after Reidemeister simplification; '
                                 '8_6 and 4_1 are invertible and 4_1 is amphichiral, so the sum is the neighbour')
    log('neighbour: crossing', c0, 'of 14a2539; sigma', sig_nb, '| alternating 12-crossing diagram, det', det_K,
        'sigma', sig_K, '| summands', re_ids, '| KnotInfo 8_6 sigma (same code)', sig_86_ki)

    # ------------------------------------------------------------- 1. Traczyk
    calib = {}
    for conv in ('ij|kl', 'il|jk'):
        calib[conv] = all(jones(kinfo.pd_code(n), conv) == kinfo.parse_jones(n) for n in ('3_1', '4_1', '5_2', '8_6'))
    conv = [c for c, ok in calib.items() if ok]
    assert len(conv) == 1, calib
    conv = conv[0]
    J_K, J_nb = jones(pdK, conv), jones(q0, conv)
    J_prod = jmul(kinfo.parse_jones('8_6'), kinfo.parse_jones('4_1'))
    assert J_K == J_nb and J_K in (J_prod, TR.mirror_jones(J_prod))
    same_as_KI = (J_K == J_prod)               # chirality of the 8_6 summand relative to KnotInfo's 8_6
    # KnotInfo: sigma(8_6) = -2 for its diagram, so the sigma = +2 representative of K is the mirror of
    # K if K contains KnotInfo's 8_6, and K itself otherwise (sigma(4_1) = 0)
    J_pos = TR.mirror_jones(J_K) if same_as_KI else J_K
    eps_pos, d_pos = TR.eps_and_dim3(J_pos)
    eps_neg, d_neg = TR.eps_and_dim3(J_K)
    w_pos = at_omega(J_pos)
    excl, _, _ = TR.traczyk_excludes_00_2(J_pos)
    # dim_F3 H_1(Sigma_2) from the Goeritz invariant factors
    G0 = OO.definite_goeritz_candidates(pdK)[0]
    snf = sympy.matrices.normalforms.smith_normal_form(sympy.Matrix(G0.tolist()), domain=sympy.ZZ)
    inv = [abs(int(snf[i, i])) for i in range(snf.rows) if abs(int(snf[i, i])) != 1]
    d3 = sum(1 for v in inv if v % 3 == 0)
    assert d_pos == d_neg == d3 == 0 and w_pos == (eps_pos, 0)
    rec['traczyk'] = dict(jones_convention_calibrated_on=['3_1', '4_1', '5_2', '8_6'], a_smoothing=conv,
                          jones_K=sorted(J_K.items()), jones_equals_V86_times_V41=same_as_KI,
                          jones_equals_mirror_V86_times_V41=not same_as_KI,
                          sigma_K_KnotInfo_convention=(-2 if same_as_KI else 2),
                          jones_of_14_crossing_neighbour_equal=True,
                          invariant_factors_H1=inv, d_K=d3, eps_sigma_plus_2=eps_pos, eps_other_mirror=eps_neg,
                          V_at_omega_sigma_plus_2=w_pos, v2_traczyk_excludes_0_2=excl,
                          consequence=('d_K = 0: in a length-two sequence d may stay 0 at the intermediate knot, '
                                       'and then the skein relation gives no condition. The criterion excludes '
                                       'no sign pattern.'))
    log('Traczyk: d_K =', d3, 'eps(sigma=+2) =', eps_pos, '| V(omega) =', w_pos, '| v2 excludes (0,2):', excl)

    # ------------------------------------------------------------- 2. Owens (1,1), implementation A (v2)
    sides = []
    for G in OO.definite_goeritz_candidates(pdK):
        mG, cmG = OO.m_Q(G)
        sides.append(dict(det=OO.det_int(G), spin=mG[cmG.coords([0] * len(G))], mG=mG, cmG=cmG, rank=len(G)))
    chosen = [s for s in sides if s['det'] == D and s['spin'] == Fraction(-1, 2)]
    assert len(chosen) == 1 and sorted(str(s['spin']) for s in sides) == ['-1/2', '1/2']
    mG, cmG = chosen[0]['mG'], chosen[0]['cmG']
    candsA = OO.candidates(D, 1)
    resA = {}
    for (m1, m2, a) in candsA:
        Qt = OO.qtilde(m1, m2, a)
        assert OO.det_int(Qt) == D and OO.positive_definite(Qt)
        mQt, cmQt = OO.m_Q(Qt)
        resA[(m1, m2, a)] = OO.matching_exists(mQt, cmQt, mG, cmG)
    official = OO._owens_theorem1(pdK, 2, D, n=1)
    verdictA = ('UNDECIDED' if any(v is None for v in resA.values())
                else 'PASS' if any(resA.values()) else 'OBSTRUCTED')
    assert official['verdict'] == verdictA
    log('A (v2): forms', len(candsA), '| passing', [k for k, v in resA.items() if v], '| verdict', verdictA)

    # ------------------------------------------------------------- 3. Owens (1,1), implementation B
    d1 = lens_vec(23, 7)
    if d1[0] == Fraction(1, 2):
        d1 = [-v for v in d1]
    d2 = lens_vec(5, 2)
    assert d1[0] == Fraction(-1, 2) and d2[0] == 0
    dvec = {t: d1[t % 23] + d2[t % 5] for t in range(D)}
    assert sorted(dvec.values()) == sorted(mG.values()), 'lens-space and Goeritz correction terms differ'
    redB = mixed_forms(D, (D - 1) // 4 + 1, True)
    assert {(m1, m2, a) for m1, m2, a, _ in redB} == set(candsA), 'reduced form lists differ'
    # reduced forms are complete (Gamma(2) reduction: |a| < min(m1, m2)); unreduced forms with a <= 8 are
    # an extra consistency check (Gamma(2)-equivalent to reduced ones, so their verdicts must reappear)
    resB = {}
    for (m1, m2, a, red) in mixed_forms(D, (D - 1) // 4 + 1, True) + \
            [f for f in mixed_forms(D, 8, False) if not f[3]]:
        Qt = qt_ind(m1, m2, a)
        assert OI.is_pos_def(Qt)
        Gr = gram_reduce(Qt)
        assert OI.is_pos_def(Gr) and sympy.Matrix(Gr).det() == D
        mvec = OI.mvector(Gr, D)
        assert mvec is not None and mvec[0] == Fraction(-1, 2)
        units = OI.admits(mvec, dvec, D)
        rev = OI.admits(mvec, {t: -v for t, v in dvec.items()}, D)
        resB[(m1, m2, a)] = dict(reduced=red, passes=bool(units), n_units=len(units), units=units[:12],
                                 reversed_orientation_passes=bool(rev))
    agree = all(resB[k]['passes'] == resA[k] for k in resA)
    verdictB = 'PASS' if any(v['passes'] for v in resB.values() if v['reduced']) else 'OBSTRUCTED'
    unred_consistent = (any(v['passes'] for v in resB.values()) == (verdictB == 'PASS'))
    log('B (independent): forms', len(resB), f"(reduced {sum(v['reduced'] for v in resB.values())})",
        '| passing reduced', [k for k, v in resB.items() if v['passes'] and v['reduced']],
        '| verdict', verdictB, '| per-form agreement with A:', agree)

    # ------------------------------------------------------------- 4. rank-one tests of the summands
    def rank1(dv, p, m):
        Qt = [[m, 1], [1, 2]]
        mv = OI.mvector(gram_reduce(Qt), p)
        return OI.admits(mv, {t: dv[t] for t in range(p)}, p)
    u1_86 = rank1(d1, 23, 12)           # u(8_6) = 1 with sigma = +2: n = 1, m even, 2m - 1 = 23
    u1_41 = rank1(d2, 5, 3)             # u(4_1) = 1 with sigma = 0: n = 0, m odd, 2m - 1 = 5
    nw_86 = ni_wu_survivors(23, {t: d1[t] for t in range(23)})
    split = resB.get((3, 12, 0))
    log('rank-one Owens test, 8_6 (m=12):', 'pass' if u1_86 else 'fail', '| 4_1 (m=3):', 'pass' if u1_41 else 'fail',
        '| Ni-Wu survivors for 8_6:', len(nw_86) if hasattr(nw_86, '__len__') else nw_86,
        '| split form (3,12,0):', split and split['passes'])

    rec['owens_1_1'] = dict(
        hypothesis='sigma = +2 representative, pattern (p,n) = (1,1): n = sigma/2 (Owens Thm 5; v2 Thm 2.1)',
        goeritz_spin_values=sorted(str(s['spin']) for s in sides), chosen_side_rank=chosen[0]['rank'],
        A=dict(implementation='code/lower_bounds/owens/owens_obstruction.py (v2)', candidates=len(candsA),
               per_form={f'{k}': v for k, v in resA.items()}, verdict=verdictA, official=official),
        B=dict(implementation='lens-space d (OS recursion) + owens_independent class minima/unit search',
               d_source='L(23,7) # L(5,2), oriented d(spin) = -1/2; multiset equals Goeritz m_G',
               forms_tested=len(resB), reduced=sum(v['reduced'] for v in resB.values()),
               per_form={f'{k}': v for k, v in resB.items()}, verdict=verdictB,
               unreduced_forms_consistent=unred_consistent),
        per_form_agreement_A_B=agree,
        summands_rank_one=dict(owens_8_6_m12_units=u1_86, owens_4_1_m3_units=u1_41,
                               ni_wu_8_6_survivors=str(nw_86)[:400], split_form_3_12_0=split))
    rec['verdict'] = dict(
        pattern_1_1='OBSTRUCTED' if verdictA == verdictB == 'OBSTRUCTED' else verdictA if verdictA == verdictB else 'DISAGREE',
        pattern_0_2='NOT EXCLUDED (Traczyk vacuous: d_K = 0; no definite filling)',
        u_8_6_4_1='still in {2, 3}' if True else None)
    json.dump(rec, open(os.path.join(OUTDIR, 'owens_traczyk_8_6_4_1.json'), 'w'), indent=1, default=str)
    log('written', os.path.join(OUTDIR, 'owens_traczyk_8_6_4_1.json'))


if __name__ == '__main__':
    main()
