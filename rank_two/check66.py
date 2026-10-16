"""Run the recovered lower-bound check (robust47.work) and upper-bound search (upper47.work) on all 66 knots of
data/new66_u3.txt, so that the 19 alternating knots (first established in alternating_audit/) and the 47
non-alternating knots go through one pipeline.
Usage: check66.py          Output: check66.json
"""
import json
from multiprocessing import Pool
from collections import Counter
import robust47, upper47

if __name__ == '__main__':
    names = [l.strip() for l in open('../data/new66_u3.txt') if l.strip()]
    with Pool(11) as p:
        lower = p.map(robust47.work, names)
        upper = p.map(upper47.work, names)
    out = {n: dict(lower=lo, upper=up) for n, lo, up in zip(names, lower, upper)}
    json.dump(out, open('check66.json', 'w'), indent=1)
    for kind, sel in (('alternating', 'a_'), ('non-alternating', 'n_')):
        R = [r for n, r in out.items() if sel in n]
        print(f"{kind} {len(R)} | lower: signed obstructed (KnotInfo V / PD V)",
              sum(r['lower']['ki_signed'] for r in R), sum(r['lower']['pd_signed'] for r in R),
              "| no-condition survivors 0:", sum(r['lower']['ki_noconds_survivors'] == 0 for r in R),
              "| upper certified:", sum(bool(r['upper'].get('ok')) for r in R), "by k", dict(Counter(r['upper'].get('k') for r in R)))
