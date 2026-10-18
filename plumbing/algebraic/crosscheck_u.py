#!/usr/bin/env python
"""Independent checks for the 32 knots of targets.txt.
  (1) a second algebraic diagram, different from the first (random search with another seed, skipping diagrams
      isomorphic to the first), gives the same normal forms (canonical strings) in both orientations; Neumann's normal
      form is unique, so any difference would expose an error. The second diagram is recorded (PD code, seed,
      canonical strings) with the reason why it differs from the first: the numbers of crossings, the multisets of
      face sizes, or a canonical form of the PD code under the choice of the first edge, reversal of the orientation
      and reflection of the plane (diagram_form);
  (2) SnapPy: the exterior of each algebraic diagram is isometric to the exterior of the KnotInfo diagram;
  (3) numpy eigenvalues of the normal forms agree with the exact inertia.
Output: crosscheck_u.json
"""
import json, os, sys, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from validate_plumbing import setup
from alg_diagrams import search, BY
from algebraic_plumbing import plumbing_from_algebraic, normal_form, plumbing_matrix


SEED2 = 7


def face_sizes(pd):
    """Sorted sizes of the faces of a diagram, from its PD code (edges counterclockwise at each crossing): the face
    at the corner between positions i and i+1 of a crossing continues at the corner between positions j and j+1 of
    the crossing at the other end of the edge in position i+1, where that edge is in position j."""
    ends = {}
    for c, x in enumerate(pd):
        for i, e in enumerate(x):
            ends.setdefault(e, []).append((c, i))
    seen, sizes = set(), []
    for c0 in range(len(pd)):
        for i0 in range(4):
            if (c0, i0) in seen:
                continue
            c, i, k = c0, i0, 0
            while (c, i) not in seen:
                seen.add((c, i)); k += 1
                (c, i), = [p for p in ends[pd[c][(i + 1) % 4]] if p != (c, (i + 1) % 4)]
            sizes.append(k)
    assert len(sizes) == len(pd) + 2                       # Euler characteristic of the sphere
    return sorted(sizes)


def diagram_form(pd):
    """Canonical form of a knot diagram on the sphere. The edges of a PD code are numbered along the knot, so a
    diagram isomorphism changes the PD code only by a shift of the numbers, together with a reversal of the numbering
    (the tuples start at the other end of the understrand) and a reflection (the cyclic order is reversed)."""
    m = {e: i for i, e in enumerate(sorted({e for x in pd for e in x}))}
    n2 = len(m)
    X = [[m[e] for e in x] for x in pd]
    assert n2 == 2 * len(pd) and all(x[2] == (x[0] + 1) % n2 for x in X)
    best = None
    for refl in (False, True):
        for rev in (False, True):
            for t in range(n2):
                T = []
                for x in X:
                    y = [x[0], x[3], x[2], x[1]] if refl else list(x)
                    if rev:
                        y = y[2:] + y[:2]
                    T.append(tuple((t - e) % n2 if rev else (e - t) % n2 for e in y))
                T = sorted(T)
                if best is None or T < best:
                    best = T
    return best


def differs(pd1, pd2):
    """How the two diagrams differ ('crossings', 'face sizes', 'diagram'), or False if they are isomorphic."""
    if len(pd1) != len(pd2):
        return 'crossings'
    if face_sizes(pd1) != face_sizes(pd2):
        return 'face sizes'
    return 'diagram' if diagram_form(pd1) != diagram_form(pd2) else False


def one(n):
    setup()
    import numpy as np, snappy, spherogram
    from both_orientations import canon
    U = json.load(open(os.path.join(HERE, 'u_plumbing.json')))
    AD = json.load(open(os.path.join(HERE, 'alg_diagrams.json')))
    _, second = search(n, seed=SEED2, avoid=lambda pd: not differs(AD[n]['pd'], pd))
    rec = dict(second_crossings=second['crossings'] if second else None, second_seed=SEED2)
    if second:
        c2 = {}
        for s in ('1', '-1'):
            wn, en = normal_form(*plumbing_from_algebraic(second['pd'], int(s)))
            c2[s] = canon(wn, en)
        rec['second_canon'] = c2
        rec['second_same_normal_forms'] = {s: c2[s] == U[n][s]['canon'] for s in c2}
        rec['second_differs_from_first'] = differs(AD[n]['pd'], second['pd'])
        rec['second_diagram'] = second
    E0 = spherogram.Link(json.loads(BY[n]['pd_notation'])).exterior()
    iso = []
    for d in [AD[n]['pd']] + ([second['pd']] if second else []):
        E1 = spherogram.Link(d).exterior()
        ok = None
        for t in range(8):
            try:
                ok = E0.is_isometric_to(E1); break
            except Exception:
                E1.randomize()
        iso.append(ok)
    rec['isometric_to_knotinfo'] = iso
    rec['knotinfo_exterior'] = dict(solution=E0.solution_type(), volume=float(E0.volume()))
    ev = {}
    for s in ('1', '-1'):
        Q = np.array(plumbing_matrix(U[n][s]['weights'], U[n][s]['edges']), dtype=float)
        x = np.linalg.eigvalsh(Q)
        ev[s] = [int((x > 1e-9).sum()), int((x < -1e-9).sum())] == [U[n][s]['inertia']['pplus'], U[n][s]['inertia']['pminus']]
    rec['numpy_inertia_agrees'] = ev
    return n, rec


if __name__ == '__main__':
    AD = json.load(open(os.path.join(HERE, 'alg_diagrams.json')))
    out = {}
    with Pool(8) as p:
        for n, r in p.imap_unordered(one, sorted(AD)):
            out[n] = r
            print(n, r, flush=True)
    json.dump(dict(sorted(out.items())), open(os.path.join(HERE, 'crosscheck_u.json'), 'w'), indent=1)
    print('second diagram found', sum(1 for r in out.values() if r['second_crossings']),
          '| same normal forms', sum(1 for r in out.values() if r.get('second_same_normal_forms') and all(r['second_same_normal_forms'].values())),
          '| different from the first', sum(1 for r in out.values() if r.get('second_differs_from_first')),
          '| isometric', sum(1 for r in out.values() if all(r['isometric_to_knotinfo'])),
          '| numpy inertia agrees', sum(1 for r in out.values() if all(r['numpy_inertia_agrees'].values())))
