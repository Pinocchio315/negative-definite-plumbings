"""Combinatorial identification certificates (review: name-dependent neighbours).

For a neighbour diagram N and the KnotInfo PD of its named knot K, we certify that the knot
complements are homeomorphic by exhibiting combinatorially isomorphic triangulations, each obtained
from the respective diagram by topology-preserving moves:

  (S) SnapPy: exterior from the PD, canonical retriangulation (2-3/3-2 moves and subdivision);
      equal isometry signatures = isomorphic triangulations.  Numerics only steer which moves are
      made, so an error can make certification fail but cannot certify falsely.
  (R) Regina (independent code): complement from the PD, greedy simplification, and if needed an
      exhaustive Pachner search of bounded height; equal isomorphism signatures.

Isomorphic triangulations => homeomorphic complements => same knot up to mirror (Gordon-Luecke
knot complement theorem); unknotting number is mirror invariant.

Usage: cert_identify.py BJSET.json TARGET OUT.json
"""
import json, sys, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
from database_knotinfo import link_list
by = {r['name']: r for r in link_list()[1:]}


def snappy_sig(pd):
    import snappy
    try:
        return snappy.Link(pd).exterior().isometry_signature()
    except Exception as e:
        return 'ERR ' + repr(e)


def regina_sig(pd, height):
    import regina
    T = regina.Link.fromPD(pd).complement()
    T.simplify()
    sigs = {T.isoSig()}
    if height:
        T2 = regina.Triangulation3(T)
        T2.simplifyExhaustive(height)
        sigs.add(T2.isoSig())
    return T.size(), sorted(sigs)


def work(job):
    name, pd_n, pd_k, meta = job
    out = dict(meta, name=name)
    sn, sk = snappy_sig(pd_n), snappy_sig(pd_k)
    out['snappy_cert'] = bool(sn) and not str(sn).startswith('ERR') and sn == sk
    out['snappy_sig'] = sn if out['snappy_cert'] else [sn, sk]
    for h in (0, 1, 2):
        tn, rn = regina_sig(pd_n, h)
        tk, rk = regina_sig(pd_k, h)
        common = set(rn) & set(rk)
        if common:
            out['regina_cert'] = True
            out['regina_height'] = h
            out['regina_sig'] = sorted(common)[0]
            out['regina_tets'] = (tn, tk)
            break
    else:
        out['regina_cert'] = False
    return out


if __name__ == '__main__':
    bjset, target, outpath = sys.argv[1], sys.argv[2], sys.argv[3]
    B = json.load(open(bjset))
    jobs = []
    for r in B['changes']:
        nm = r['ids'][0]
        if nm in by:
            jobs.append((nm, r['pd'], json.loads(by[nm]['pd_notation']),
                         dict(kind='neighbour', ref=r['ref'], subset=r['subset'], crossing=r['crossing'])))
    D = json.load(open(bjset.replace('bjset_', 'diagrams_')))
    from enum_diagrams import change
    for d in D['diagrams']:
        jobs.append((target, change(json.loads(by[d['ref']]['pd_notation']), d['subset']),
                     json.loads(by[target]['pd_notation']), dict(kind='minimal_diagram', ref=d['ref'], subset=d['subset'])))
    with Pool(11) as p:
        res = p.map(work, jobs)
    json.dump(res, open(outpath, 'w'), indent=0)
    from collections import Counter
    print(Counter((r['kind'], r['name'], r['snappy_cert'], r['regina_cert']) for r in res))
    print('regina heights used:', Counter(r.get('regina_height') for r in res))
