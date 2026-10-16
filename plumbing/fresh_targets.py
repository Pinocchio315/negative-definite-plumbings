#!/usr/bin/env python
"""The covers newly recognised as graph manifolds by recognise_more.py (60 randomisations) that no earlier plumbing
route (star_u1*.json, graph_u1.json, owens_u3_targets.json) had processed: the input of neumann/run_fresh.py.
(Recorded on 2026-10-15 by a one-off command; this script regenerates the same file.)
Usage: fresh_targets.py          Output: recognise_more_fresh.json
"""
import json

if __name__ == '__main__':
    M = json.load(open('recognise_more.json'))
    done = set(json.load(open('star_u1.json'))) | set(json.load(open('star_u1_orb.json'))) | set(json.load(open('graph_u1.json')))
    done |= {n for n, r in json.load(open('owens_u3_targets.json')).items() if r}
    out = {n: dict(kind=r['new'], names=sorted(r['names'])) for n, r in sorted(M.items())
           if r['new'] in ('star', 'pair', 'triple') and n not in done}
    json.dump(out, open('recognise_more_fresh.json', 'w'), indent=1)
    print(len(out), 'knots')
