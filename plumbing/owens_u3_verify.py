#!/usr/bin/env python
"""Verification of the Owens obstructions in owens_u3_targets.json.

(1) Independent Owens check (bernhard_jablan/code/owens_independent.py; imports nothing from v2): all forms with
    the required parity, NOT only reduced ones (0 <= a <= D//2 + 2), Gram-reduced for the minima, against the plumbing
    correction terms on Z/D (d_owens.plumbing_d, spin at 0), in both orientations.  The orientations with no fitting
    form must be the same as in the v2-based run.
(2) Casson-Walker: unique sign s_cw with s_cw L(Y) = R(K), equal to the linking-form orientation when that is defined.
(3) Novelty: KnotInfo 2026-09-16, v2 table, this project's 66 values (data/new66_u3.txt).
Output: owens_u3_verify.json
"""
import json, math, os, sys, warnings
warnings.filterwarnings('ignore')
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
sys.path.insert(0, os.path.join(ROOT, 'bernhard_jablan', 'code'))
sys.path.insert(0, os.path.join(ROOT, 'lib'))
import owens_independent as OI
from owens_traczyk_8_6_4_1 import qt_ind, gram_reduce
import d_owens as DO
k916 = json.load(open(os.path.join(ROOT, 'data', 'knotinfo_2026-09-16_unknotting.json')))
v2T = json.load(open(os.path.join(ROOT, 'data', 'v2_final_u_table.json')))
new66 = {l.strip() for l in open(os.path.join(ROOT, 'data', 'new66_u3.txt')) if l.strip()}


def forms_all(D, n_even, amax):
    """(m1, m2, a): x = 2m1-1 <= y = 2m2-1, xy - 4a^2 = D, exactly n_even of m1, m2 even, 0 <= a <= amax."""
    out = []
    for a in range(0, amax + 1):
        N = D + 4 * a * a
        for x in range(1, math.isqrt(N) + 1, 2):
            if N % x:
                continue
            y = N // x
            if y < x or y % 2 == 0:
                continue
            m1, m2 = (x + 1) // 2, (y + 1) // 2
            if (m1 % 2 == 0) + (m2 % 2 == 0) == n_even:
                out.append((m1, m2, a))
    return out


if __name__ == '__main__':
    T = json.load(open('owens_u3_targets.json'))
    out = {}
    for n, r in sorted(T.items()):
        if not r or not str(r.get('verdict', '')).startswith('OBSTRUCTED') or 'Wendt' in r['verdict']:
            continue
        f = r['plumbing']; D = r['D']; sigma = r['sigma']
        n_even = 2 if abs(sigma) == 4 else 1
        rec = dict(verdict_v2code=r['verdict'], sigma=sigma, D=D, group=r['group'])
        if len(r['group']) == 1 or math.prod(r['group']) == D and all(math.gcd(a, b) == 1 for i, a in enumerate(r['group']) for b in r['group'][i + 1:]):
            Dp, dZ = DO.plumbing_d(f['weights'], f['edges'])
            fl = forms_all(D, n_even, D // 2 + 2)
            mv = {g: OI.mvector(gram_reduce(qt_ind(*g)), D) for g in fl}
            fit = {}
            for s in (1, -1):
                dv = {t: s * v for t, v in dZ.items()}
                fit[s] = [g for g, m in mv.items() if m is not None and OI.admits(m, dv, D)]
            blocked = sorted(s for s in (1, -1) if not fit[s])
            expected = r['s_expected']
            v2_blocked = [1, -1] if r['verdict'] == 'OBSTRUCTED' else [expected]
            rec.update(independent_forms=len(fl), independent_fit_plus=fit[1][:3], independent_fit_minus=fit[-1][:3],
                       independent_blocked=blocked, agrees=(sorted(v2_blocked) == blocked or (r['verdict'] != 'OBSTRUCTED' and expected in blocked)))
        else:
            rec.update(independent='skipped (H1 not cyclic; v2 matching_exists handles general groups)')
        cw = r['cw_signs']; s_link = 1 if r['plumbing']['orientation'] == '+' else -1
        rec.update(cw_signs=cw, link_sign=s_link, link_ambiguous=r['plumbing']['orientation_ambiguous'],
                   cw_consistent=len(cw) == 1 and (cw[0] == s_link or r['plumbing']['orientation_ambiguous']),
                   knotinfo_2026_09_16=k916.get(n), v2_table=v2T.get(n), in_new66=n in new66)
        out[n] = rec
        print(n, rec, flush=True)
    json.dump(out, open('owens_u3_verify.json', 'w'), indent=1, default=str)
    print('verified', sum(1 for v in out.values() if v.get('agrees') and v['cw_consistent']), '/', len(out))
