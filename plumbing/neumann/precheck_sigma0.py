#!/usr/bin/env python
"""Normal-form verdict (precheck.one) for the open knots of the Sigma_2 map whose double branched cover is a graph
manifold but which no plumbing pipeline processed: sigma = 0 with u in [2,3], so no correction-term test applies and
the plumbing was never needed.  Completes the classification "does Sigma_2 or -Sigma_2 bound a negative-definite
plumbing tree" for all recognised graph-manifold covers of the map.
Names: the map's recorded name, plus 16 randomised orbifold-route recognitions of the knot and of its mirror.
Usage: precheck_sigma0.py [workers]          Output: precheck_sigma0.json
"""
import json, os, sys
from multiprocessing import Pool
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
TOOL = os.path.join(HERE, '..')
sys.path.insert(0, TOOL)
import precheck as PC
import graph_u1 as GU


def J(*p):
    return json.load(open(os.path.join(TOOL, *p)))


def unprocessed():
    VH, JV, JS = J('verify_hyp.json'), J('jsj_verify.json'), J('jsj.json')
    done = set(J('neumann', 'run_pachner4.json')) | set(J('neumann', 'run_retry13.json')) | set(J('neumann', 'run_fresh.json'))
    for src in ({**J('star_u1.json'), **J('star_u1_orb.json')}, J('graph_u1.json'), J('owens_u3_targets.json')):
        done |= {n for n, r in src.items() if r}
    out = []
    for f in ('reach_u12_orb.json', 'reach_u23_orb.json'):
        for x in J(f):
            n = x['name']
            if (n in VH and VH[n]['certified']) or (n in JV and JV[n]['certified_hyperbolic_piece']) or n in done:
                continue
            if n in JS and JS[n]['verdict'].startswith('toroidal, no geometric'):
                continue
            out.append(x)
    return out


def one(x):
    import spherogram
    n = x['name']
    pd = json.loads(PC.BY[n]['pd_notation']); mpd = [list(y) for y in spherogram.Link(pd).mirror().PD_code()]
    names = {x['sigma2']} | set(GU.names_orb(pd, 16)) | set(GU.names_orb(mpd, 16))
    _, rec = PC.one((n, sorted(names)))
    rec.update(kind=x['kind'], sigma=x['sigma'], range=x['range'])
    return n, rec


if __name__ == '__main__':
    items = unprocessed()
    print('unprocessed graph-type knots', len(items), dict(Counter(x['kind'] for x in items)), flush=True)
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as p:
        out = dict(p.map(one, items))
    json.dump(out, open(os.path.join(HERE, 'precheck_sigma0.json'), 'w'), indent=1, default=str)
    print('verdicts', dict(Counter(r['verdict'] for r in out.values())))
    print('p+ by orientation', {n: r['pplus'] for n, r in out.items() if r['verdict'] != 'APPLICABLE'})
