"""Correction terms of a graph manifold from a negative-definite plumbing tree, and Owens's (1,1) test.

d(Y, t) = max_{K in Char_t} (K^2 + |V|)/4 for Y = boundary of a negative-definite plumbing TREE (lattice homology
= Heegaard Floer homology, Zemke; Nemethi), computed exactly as -m_{-Q}(t) with v2's ellipsoid enumeration.
Labels: ClassMap prime-power coordinates -> Z/D by CRT (homomorphic, spin at 0).
"""
import json, os, sys, warnings
warnings.filterwarnings('ignore')
from fractions import Fraction
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'lib', 'lower_bounds', 'owens'))
sys.path.insert(0, os.path.join(ROOT, 'bernhard_jablan', 'code'))
sys.path.insert(0, os.path.join(ROOT, 'lib'))
import _compat  # noqa
import owens_obstruction as OO
import owens_independent as OI
from owens_traczyk_8_6_4_1 import mixed_forms, qt_ind, gram_reduce
import plumbing as P


def crt_label(coords, mods):
    t, M = 0, 1
    for x, m in zip(coords, mods):
        # solve t' = t mod M, t' = x mod m
        for k in range(m):
            if (t + k * M) % m == x % m:
                t = t + k * M; break
        M *= m
    return t % M


def plumbing_d(weights, edges):
    Q = P.plumbing_matrix(weights, edges)
    Qpos = -np.array(Q, dtype=object)
    assert OO.positive_definite(Qpos), 'plumbing must be negative definite'
    mq, cm = OO.m_Q(Qpos)
    D = cm.order
    d = {crt_label(g, cm.d): -v for g, v in mq.items()}
    assert set(d) == set(range(D))
    return D, d


def owens_1_1(D, d):
    out = {}
    forms = mixed_forms(D, (D - 1) // 4 + 1, True)
    mvecs = {f[:3]: OI.mvector(gram_reduce(qt_ind(*f[:3])), D) for f in forms}
    for o in (1, -1):
        vec = {t: o * v for t, v in d.items()}
        passing = [f for f, mv in mvecs.items() if mv is not None and OI.admits(mv, vec, D)]
        out['Y' if o == 1 else '-Y'] = dict(d_spin=str(vec[0]), passing=passing)
    return dict(forms=sorted(mvecs), result=out, obstructed=all(not r['passing'] for r in out.values()))


if __name__ == '__main__':
    found = json.load(open('plumbing_13n309.json'))['found']
    rec = {}
    for i, f in enumerate(found):
        D, d = plumbing_d(f['weights'], f['edges'])
        ow = owens_1_1(D, d)
        rec[f'match{i}'] = dict(plumbing=f, D=D, d={t: str(v) for t, v in sorted(d.items())}, owens_1_1=ow)
        print(f'match {i}: D {D} | d(spin) {d[0]} | d values', sorted(set(str(v) for v in d.values())),
              '| Owens (1,1):', {k: (v['d_spin'], v['passing'][:4]) for k, v in ow['result'].items()},
              '->', 'OBSTRUCTED (both orientations)' if ow['obstructed'] else 'not obstructed', flush=True)
    json.dump(rec, open('d_owens_13n309.json', 'w'), indent=1, default=str)
