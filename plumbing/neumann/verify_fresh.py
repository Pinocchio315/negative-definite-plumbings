#!/usr/bin/env python
"""Verification of the obstructions in run_fresh.json.
u = 1 obstructions: Casson-Walker consistency, independent Ni-Wu (d-obstructions), novelty.
u = 2 obstructions (Owens): independent Owens step (owens_general_group: all forms with the required parity up to
a <= D//2 + 2, exact minima of Qtilde, every group isomorphism, both orientations) and Traczyk recomputed from
V(e^{i pi/3}); novelty.
Output: verify_fresh.json
"""
import json, math, os, sys, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, '..')
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
sys.path.insert(0, TOOL)
import numpy as np
import owens_general_group as OG          # group_of, find_basis, isomorphisms, jones_at_omega (imports owens_independent)
import owens_obstruction as OO
import plumbing as P
import kinfo
from owens_u3_verify import forms_all
from owens_traczyk_8_6_4_1 import qt_ind, gram_reduce
import owens_independent as OI
k916 = json.load(open(os.path.join(ROOT, 'data', 'knotinfo_2026-09-16_unknotting.json')))
v2T = json.load(open(os.path.join(ROOT, 'data', 'v2_final_u_table.json')))
new66 = {l.strip() for l in open(os.path.join(ROOT, 'data', 'new66_u3.txt')) if l.strip()}
F = json.load(open(os.path.join(HERE, 'run_fresh.json')))
out = {}
for n, r in sorted(F.items()):
    if not str(r.get('verdict', '')).startswith('OBSTRUCTED'):
        continue
    rec = dict(test=r['test'], verdict=r['verdict'], u_before=r['u'], cw_consistent=r['cw_consistent'],
               knotinfo_2026_09_16=k916.get(n), v2_table=v2T.get(n), in_new66=n in new66)
    if r['test'] == 'u=1':
        rec['independent_ni_wu_survivors'] = r.get('independent_ni_wu_survivors', 'n/a (HF_red obstruction)')
    else:
        D = r['D']; sigma = r['sigma']; n_even = 2 if abs(sigma) == 4 else 1
        Q = P.plumbing_matrix(r['plumbing']['weights'], r['plumbing']['edges'])
        mP, cm = OO.m_Q(-np.array(Q, dtype=object))
        dY = {g: -v for g, v in mP.items()}
        elemsB, addB, zeroB, orderB = OG.group_of(dY.keys(), list(cm.d))
        fl = forms_all(D, n_even, D // 2 + 2)
        fit = {1: [], -1: []}; checked = 0
        for g in fl:
            mQ, DQ, _ = OI.class_minima(gram_reduce(qt_ind(*g)))
            elemsA, addA, zeroA, orderA = OG.group_of(mQ.keys(), [D] * 4)
            if sorted(orderA(x) for x in elemsA) != sorted(orderB(x) for x in elemsB):
                continue
            checked += 1
            basis, ords = OG.find_basis(elemsA, addA, zeroA, orderA)
            for s in (1, -1):
                for phi in OG.isomorphisms(basis, ords, addA, zeroA, elemsB, addB, zeroB, orderB):
                    if all((mQ[x] - s * dY[phi[x]]) >= 0 and (mQ[x] - s * dY[phi[x]]).denominator == 1
                           and (mQ[x] - s * dY[phi[x]]).numerator % 2 == 0 for x in mQ):
                        fit[s].append(g); break
        J = kinfo.parse_jones(n)
        eps, d3 = OG.jones_at_omega(J if sigma > 0 else {-e: c for e, c in J.items()})
        rec.update(group=list(cm.d), forms=len(fl), forms_with_isomorphic_group=checked, fit_plus=fit[1][:3], fit_minus=fit[-1][:3],
                   independent_blocked_both=(not fit[1] and not fit[-1]), traczyk_eps=eps, dim3=d3,
                   traczyk_excludes_02=(d3 == 2 and eps == -1))
    out[n] = rec
    print(n, rec, flush=True)
json.dump(out, open(os.path.join(HERE, 'verify_fresh.json'), 'w'), indent=1, default=str)
