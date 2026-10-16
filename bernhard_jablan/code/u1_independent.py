"""Independent re-check of u = 1 obstructions from correction terms (imports nothing from the
Greene/half_integral pipeline).

(a) Linking-form consistency of the supplied d-vector: for an odd-order rational homology sphere,
    {d(t) - d(spin) mod 1} must equal {s * lambda(x, x) mod 1 : x in H} as multisets for one sign s,
    with lambda computed independently from the Seifert matrix V (lambda = (V + V^T)^{-1}).
(b) Ozsvath-Szabo / Montesinos: u(K) = 1 implies Sigma_2(K) = e * S^3_{D/2}(C) for a knot C and a
    sign e.  Ni-Wu: d(S^3_{p/q}(C), i) = d(L(p,q), i) - 2 max(V_{floor(i/q)}, V_{floor((p+q-1-i)/q)}),
    V_j >= 0 integers, V_j - V_{j+1} in {0, 1}.  We test both signs and every affine identification
    i -> a*i + b of Z/D (a a unit).  No survivor => u(K) >= 2.

Usage: u1_independent.py NAME=PATH_TO_JSON_WITH_d [...]
"""
import json, sys, warnings
from fractions import Fraction
from math import gcd
warnings.filterwarnings('ignore')
sys.path.insert(0, __import__('os').path.dirname(__file__))
from database_knotinfo import link_list
from bf2 import LinkingForm

by = {r['name']: r for r in link_list()[1:]}


def lens_d(p, q, i):
    """d(L(p,q), i) for L(p,q) = p/q surgery on the unknot (Ozsvath-Szabo recursion), 0 <= i < p + q."""
    if p == 1:
        return Fraction(0)
    return Fraction(-1, 4) + Fraction((2 * i + 1 - p - q) ** 2, 4 * p * q) - lens_d(q, p % q, i % q)


def load_d(path):
    r = json.load(open(path))
    d = r.get('d')
    if d is None:
        raise ValueError(f'{path}: no fully pinned d-vector')
    D = int(r['det'])
    vec = {int(k) % D: Fraction(v) for k, v in d.items()}
    assert set(vec) == set(range(D))
    return D, vec


def linking_check(name, D, vec):
    V = json.loads(by[name]['seifert_matrix']); n = len(V)
    LF = LinkingForm([[V[i][j] + V[j][i] for j in range(n)] for i in range(n)])
    assert LF.D == D
    # lambda(x, x) = x^T B x / D mod 1, with M^{-1} = B / D
    vals = []
    for x, bx in zip(LF.elems, LF.Bx):
        q = sum(a * b for a, b in zip(x, bx)) % D
        vals.append(Fraction(q, D))
    plus = sorted(v % 1 for v in vals)
    minus = sorted((-v) % 1 for v in vals)
    dd = sorted((vec[t] - vec[0]) % 1 for t in range(D))
    return dd == plus, dd == minus


def ni_wu_survivors(D, vec):
    p, q = D, 2
    L = [lens_d(p, q, i) for i in range(p)]
    jmap = [min(i // 2, (p + 1 - i) // 2) for i in range(p)]
    jmax = max(jmap)
    survivors = []
    for e in (1, -1):
        for a in range(1, D):
            if gcd(a, D) != 1:
                continue
            for b in range(D):
                Vv = [None] * (jmax + 1)
                ok = True
                for i in range(p):
                    m = (L[i] - e * vec[(a * i + b) % D]) / 2
                    if m.denominator != 1 or m < 0:
                        ok = False; break
                    j = jmap[i]
                    if Vv[j] is None:
                        Vv[j] = int(m)
                    elif Vv[j] != m:
                        ok = False; break
                if not ok:
                    continue
                if any(v is None for v in Vv):
                    # an unused index can be filled monotonically; check gaps between used ones
                    pass
                used = [(j, v) for j, v in enumerate(Vv) if v is not None]
                for (j1, v1), (j2, v2) in zip(used, used[1:]):
                    if not (0 <= v1 - v2 <= j2 - j1):
                        ok = False; break
                if ok:
                    survivors.append((e, a, b))
    return survivors


if __name__ == '__main__':
    for arg in sys.argv[1:]:
        name, path = arg.split('=')
        D, vec = load_d(path)
        lk = linking_check(name, D, vec) if name in by and by[name].get('seifert_matrix') else None
        surv = ni_wu_survivors(D, vec)
        print(f'{name:9s} D={D:4d} d(0)={str(vec[0]):6s} sigma(KnotInfo)={by[name]["signature"]:>3s} '
              f'-sigma/4={str(-Fraction(int(by[name]["signature"]), 4)):5s} '
              f'linking-form check (+,-)={lk}  Ni-Wu survivors={len(surv)}'
              f'{"  -> u >= 2" if not surv else "  -> u=1 NOT excluded"}')
