#!/usr/bin/env python
"""JSJ-type analysis of the 180 Sigma_2 that are neither recognised by Regina nor certified hyperbolic.

For each knot: the smallest of several orbifold-route triangulations of Sigma_2; irreducibility (Regina); vertex normal
surfaces in quadrilateral coordinates (embedded); candidate tori = connected two-sided orientable surfaces with chi = 0,
plus the boundaries of regular neighbourhoods of one-sided Klein bottles (doubled surfaces); incompressibility by
Regina's NormalSurface.isIncompressible().  For each incompressible torus: cut along it and describe the pieces:
solid torus / T^2 x I tests, Regina recognition, and SnapPy's solution type on the ideal version (finiteToIdeal);
the pieces are written out (isoSig) for verified hyperbolicity in Sage (jsj_verify.py).
A vertex incompressible torus exists whenever an incompressible torus exists (Jaco-Tollefson; Tollefson's
quadrilateral-vertex theorem), so 'no incompressible vertex torus' means atoroidal.
Usage: jsj.py [workers] [limit]          Output: jsj.json
"""
import json, sys, time, warnings
warnings.filterwarnings('ignore')
from collections import Counter
from multiprocessing import Pool
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}


def smallest_sigma2(n, k=12):
    import snappy, regina
    E = snappy.Link(json.loads(BY[n]['pd_notation'])).exterior()
    best = None
    for t in range(k):
        E2 = E.copy()
        if t:
            E2.randomize()
        E2.dehn_fill((2, 0)); C = E2.covers(2, cover_type='cyclic')[0]
        T = regina.Triangulation3(C.filled_triangulation()._to_string()); T.intelligentSimplify()
        if best is None or T.size() < best.size():
            best = T
    return best


def describe_piece(P):
    import regina, snappy
    out = dict(tetrahedra=P.size(), boundary_components=P.countBoundaryComponents())
    try:
        out['solid_torus'] = bool(P.isSolidTorus())
    except Exception:
        out['solid_torus'] = None
    try:
        out['TxI'] = bool(P.isTxI()) if P.countBoundaryComponents() == 2 else False
    except Exception:
        out['TxI'] = None
    rec = regina.StandardTriangulation.recognise(P)
    out['regina'] = rec.manifold().name() if rec and rec.manifold() else None
    I = regina.Triangulation3(P); I.finiteToIdeal(); I.intelligentSimplify()
    out['ideal_isosig'] = I.isoSig()
    try:
        M = snappy.Manifold(I.isoSig())
        st = None
        for t in range(10):
            st = M.solution_type()
            if 'positively' in st:
                break
            M.randomize()
        out['snappy_solution'] = st
        out['snappy_volume'] = float(M.volume())
        try:
            ids = M.identify()
            out['snappy_identify'] = [str(x) for x in ids][:3]
        except Exception:
            out['snappy_identify'] = []
    except Exception as e:
        out['snappy_error'] = repr(e)[:80]
    return out


def one(n):
    import regina
    t0 = time.time()
    T = smallest_sigma2(n)
    rec = dict(tetrahedra=T.size(), H1=str(T.homology()))
    try:
        rec['irreducible'] = bool(T.isIrreducible())
    except Exception as e:
        rec['irreducible'] = None
    S = regina.NormalSurfaces(T, regina.NS_QUAD)
    tori = []
    for i in range(S.size()):
        s = S.surface(i)
        if not s.isCompact() or not s.isConnected() or s.eulerChar() != 0:
            continue
        if s.isOrientable() and s.isTwoSided():
            tori.append(('torus', i, s))
        elif not s.isOrientable() and not s.isTwoSided():
            tori.append(('klein-bottle double', i, s.doubleSurface()))
    rec['vertex_surfaces'] = S.size(); rec['candidate_tori'] = len(tori)
    found = []
    for kind, i, s in tori:
        try:
            inc = bool(s.isIncompressible())
        except Exception as e:
            inc = None
        if not inc:
            continue
        cut = s.cutAlong(); cut.intelligentSimplify()
        comps = cut.triangulateComponents() if hasattr(cut, 'triangulateComponents') else [cut]
        pieces = [describe_piece(regina.Triangulation3(c)) for c in comps]
        found.append(dict(kind=kind, surface=i, pieces=pieces))
        if len(found) >= 3:
            break
    rec['incompressible_tori'] = found
    if not found:
        rec['verdict'] = 'atoroidal (no incompressible vertex torus)'
    else:
        sols = [p.get('snappy_solution', '') for f in found for p in f['pieces']]
        if any('positively' in (x or '') for x in sols):
            rec['verdict'] = 'toroidal, a piece with a geometric solution (likely hyperbolic piece: mixed)'
        else:
            rec['verdict'] = 'toroidal, no geometric piece found (graph manifold candidate)'
    rec['seconds'] = round(time.time() - t0, 1)
    return n, rec


if __name__ == '__main__':
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    limit = int(sys.argv[2]) if len(sys.argv) > 2 else None
    VH = json.load(open('verify_hyp.json')); RR = json.load(open('recognise_retry.json'))
    targets = sorted(n for n, r in VH.items() if not r['certified'] and not RR.get(n, {}).get('names'))[:limit]
    t0 = time.time(); out = {}
    with Pool(workers) as p:
        for n, rec in p.imap_unordered(one, targets):
            out[n] = rec
            print(f"[{time.time()-t0:5.0f}s] {n:9s} tet {rec['tetrahedra']} irr {rec['irreducible']} vertex {rec['vertex_surfaces']} "
                  f"tori {rec['candidate_tori']} incompressible {len(rec['incompressible_tori'])} | {rec['verdict']} | "
                  f"{[[(p['regina'], (p.get('snappy_solution') or '')[:12], round(p.get('snappy_volume', 0), 3)) for p in f['pieces']] for f in rec['incompressible_tori']][:2]} | {rec['seconds']}s", flush=True)
            json.dump(out, open('jsj.json' if limit is None else 'jsj_test.json', 'w'), indent=1, default=str)
    print('verdicts', dict(Counter(r['verdict'] for r in out.values())), f'{time.time()-t0:.0f}s')
