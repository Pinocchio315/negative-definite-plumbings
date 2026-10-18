#!/usr/bin/env python
"""Certificates of the u = 1 exclusions by reduced Floer homology (Section 3.2 of the paper).

For every knot of rederive_all.json with verdict OBSTRUCTED_HF, from its negative-definite normal form Gamma (one bad
vertex):
  - d(Y, t) and rank HF_red(Y, t) for the classes t = g^k of a generator g of H_1(Y) = Z/D (k in Z/D);
  - every identification that survives the correction-term test (1): the sign eps (Y or -Y is S^3_{D/2}(C)), the
    affine map i -> a i + b from the Ni-Wu labels i of Spin^c(S^3_{D/2}(C)) to the classes k, and the sequence V;
  - for each of them the quantities of (2) at every label i: T_i (from V), the rank r_i + T_i pulled back from Y, r_i,
    and 2c_i = sum_j (-1)^j r_{i-j}; and the first condition that fails (some r_i < 0, or some 2c_i negative or odd).
The enumeration is that of montesinos_u1.u1_admissible and mapping_cone_ranks_ok, written out so that every quantity is
recorded. A knot is certified when every surviving identification fails.
Usage: hfred_certificates.py          Output: hfred_certificates.json
"""
import json, math, os, sys, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import run_invert as RI


def T_values(V, D):
    """T_i of (2) for the Ni-Wu labels i, from the sequence V (V_s = 0 beyond it, V_{-s} = V_s + s)."""
    def Vat(s):
        if s < 0:
            return Vat(-s) - s
        return V[s] if s < len(V) else 0
    out = []
    for i in range(D):
        Ti = 0
        caseA = (i // 2) <= -((i - D) // 2)
        nV, nH = (1, 1) if caseA else (0, 2)
        n = nV
        while True:
            sv = (i + n * D) // 2
            if sv >= len(V):
                break
            Ti += Vat(sv); n += 1
        n = nH
        while True:
            sh = -((i - n * D) // 2)
            if sh >= len(V):
                break
            Ti += Vat(sh); n += 1
        out.append(Ti)
    return out


def certificate(n, rec):
    MU = RI.MU
    pl = rec['plumbing']
    w, edges = pl['weights'], pl['edges']
    Q = np.array(RI.P.plumbing_matrix(w, edges), dtype=np.int64)
    D = int(RI.BY[n]['determinant'])
    d, cm = MU.correction_terms(Q)
    ranks, _ = RI.hf_red_all(Q)
    lab, lam = MU.cyclic_labelling(cm, D, Q)
    dL = [MU.lens_d(D, 2, i) for i in range(D)]
    idx = [min(i // 2, (D + 1 - i) // 2) for i in range(D)]
    out = dict(D=D, weights=w, edges=edges, generator_self_linking=str(lam),
               d_by_class=[str(d[lab[k]]) for k in range(D)],
               hf_red_rank_by_class=[int(ranks[lab[k]]) for k in range(D)], survivors=[])
    for eps in (1, -1):
        units = [a for a in range(D) if math.gcd(a, D) == 1 and (a * a * lam) % 1 == Fraction(eps * 2, D) % 1]
        dY = [eps * d[lab[k]] for k in range(D)]
        for a in units:
            for b in range(D):
                V, ok = {}, True
                for i in range(D):
                    x = dL[i] - dY[(a * i + b) % D]
                    if x < 0 or x.denominator != 1 or x.numerator % 2:
                        ok = False; break
                    v = x.numerator // 2
                    if V.setdefault(idx[i], v) != v:
                        ok = False; break
                if not ok:
                    continue
                seq = [V[j] for j in sorted(V)]
                if not all(0 <= seq[j] - seq[j + 1] <= 1 for j in range(len(seq) - 1)):
                    continue
                s = dict(eps=eps, a=a, b=b, V=seq)
                if seq[-1] != 0:
                    s.update(fails=None, note='last V nonzero: (2) not applied')
                else:
                    T = T_values(seq, D)
                    pulled = [int(ranks[lab[(a * i + b) % D]]) for i in range(D)]
                    r = [pulled[i] - T[i] for i in range(D)]
                    c2 = [sum((-1) ** j * r[(i - j) % D] for j in range(D)) for i in range(D)]
                    fails = None
                    neg = [i for i in range(D) if r[i] < 0]
                    bad = [i for i in range(D) if c2[i] < 0 or c2[i] % 2]
                    if neg:
                        fails = f'r_{neg[0]} = {r[neg[0]]} < 0'
                    elif bad:
                        fails = f'2c_{bad[0]} = {c2[bad[0]]} is negative or odd'
                    s.update(T=T, rank_pulled_back=pulled, r=r, two_c=c2, fails=fails)
                out['survivors'].append(s)
    out['certified'] = bool(out['survivors']) and all(s['fails'] for s in out['survivors'])
    return out


if __name__ == '__main__':
    RA = json.load(open(os.path.join(HERE, 'rederive_all.json')))
    targets = sorted(k for k, r in RA.items() if r['verdict'] == 'OBSTRUCTED_HF')
    res = {}
    for n in targets:
        res[n] = certificate(n, RA[n])
        c = res[n]
        print(n, 'D', c['D'], '| survivors', [(s['eps'], s['a'], s['b'], s['V'], s['fails']) for s in c['survivors']],
              '| certified', c['certified'], flush=True)
    json.dump(res, open(os.path.join(HERE, 'hfred_certificates.json'), 'w'), indent=1)
