#!/usr/bin/env python
"""Normal-form computation for the knots of class (3) of the classification whose negative-definite plumbing was found
only by the earlier plumbing searches (star_u1*.json, graph_u1.json, owens_u3_targets.json; tests inconclusive), so that
every knot of class (3) is placed there by its Neumann normal form (Corollary 2.5 of the paper).
Same pipeline as rederive_all.py: Regina names -> precheck.one (normal forms of Sigma_2 and -Sigma_2) -> the tests of
run_fresh.one; the negative-definite normal form is also identified independently by triangulating the plumbing.
Since these knots bound a negative-definite plumbing tree, Neumann's theorem forces the verdict APPLICABLE.
Usage: normal_case3.py [workers]          Output: normal_case3.json
"""
import json, os, sys
from multiprocessing import Pool
from collections import Counter
import rederive_all as RA


def targets():
    J = lambda *f: json.load(open(os.path.join(RA.RF.TOOL, *f)))
    done = set(J('neumann', 'run_fresh.json')) | set(J('neumann', 'run_retry13.json')) | set(J('neumann', 'run_pachner4.json'))
    out = set()
    for src in ({**J('star_u1.json'), **J('star_u1_orb.json')}, J('graph_u1.json'), J('owens_u3_targets.json')):
        for n, r in src.items():
            if r and r.get('verdict') != 'NO PLUMBING' and not str(r.get('verdict')).startswith('OBSTRUCTED') and n not in done:
                out.add(n)
    return sorted(out)


if __name__ == '__main__':
    items = [dict(knot=n, lo=None, source='case 3, earlier search') for n in targets()]
    print('targets', len(items), flush=True)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as p:
        out = dict(p.imap_unordered(RA.one, items))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'normal_case3.json'), 'w'), indent=1)
    print('applicability', dict(Counter(r['applicability'] for r in out.values())))
    print('verdicts', dict(Counter(str(r['verdict']) for r in out.values())))
    print('CW consistent', sum(bool(r['cw_consistent']) for r in out.values()), '/', len(out),
          '| triangulated normal form recognised as Sigma_2:', sum(bool(r.get('plumbing_name_matches')) for r in out.values()))
