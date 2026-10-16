"""Independent implementation of Owens's correction-term obstruction to u = n (here n = 2)
for a knot with sigma = 2n and L-space-style fully pinned correction terms (imports nothing from
the v2 pipeline).

Owens (Adv. Math. 2008, Thm 5): if sigma(K) = 2n and u(K) = n, then Sigma_2(K) bounds a
positive-definite X with intersection form  Qt = [[M, I], [I, 2I]],  2M - I = Q,
Q_ii = 2 m_i - 1 with every m_i even, Q_ij = 2 a_ij, det Q = det K, and there is a group
isomorphism phi: Z^{2n}/Qt -> Spin^c(Y) (spin at 0) with, for every class g,
    m_Qt(g) >= d(Y, phi(g))   and   m_Qt(g) = d(Y, phi(g)) mod 2,
    m_Qt(g) = min over characteristic xi in g of (xi^T Qt^{-1} xi - rank)/4.

Minima are computed EXACTLY: for a class with representative xi0, xi = xi0 + 2 Qt z and
f(z) = xi^T Qt^{-1} xi is a positive-definite quadratic in z; we enumerate every integer z in the
ball |z - z*|^2 <= (U - f(z*)) / (4 lambda_lower), U an attained value, lambda_lower a certified
lower bound for the least eigenvalue of Qt (Gershgorin-free: exact Sylvester check on Qt - lambda I).

Usage: owens_independent.py D_JSON  (file with 'det' and fully pinned 'd')  [AMAX]
"""
import itertools, json, math, sys
from fractions import Fraction
import numpy as np
import sympy


def is_pos_def(A):
    """Exact Sylvester criterion for a rational symmetric matrix."""
    n = len(A)
    return all(sympy.Matrix([[A[i][j] for j in range(k)] for i in range(k)]).det() > 0 for k in range(1, n + 1))


def lambda_lower(A):
    lam = float(np.linalg.eigvalsh(np.array(A, dtype=float)).min())
    cand = Fraction(lam).limit_denominator(10 ** 6) * Fraction(9, 10)
    B = [[Fraction(A[i][j]) - (cand if i == j else 0) for j in range(len(A))] for i in range(len(A))]
    assert cand > 0 and is_pos_def(B), 'eigenvalue lower bound not certified'
    return cand


def class_minima(Qt):
    """Exact m(g) for every class g of Z^r / Qt Z^r, in integer arithmetic.
    For a class with characteristic representative xi0, xi = xi0 + 2 Qt z and
        D * xi^T Qt^{-1} xi = N(z) = xi0^T A xi0 + 4 D z.xi0 + 4 D z^T Qt z,   A = adj(Qt) = D Qt^{-1}.
    With z* = -A xi0 / (2D) the real minimiser, N(z) - N(z*) = 4 D (z - z*)^T Qt (z - z*).  Every z with
    N(z) <= U satisfies |z_i - z*_i| <= sqrt(R * (Qt^{-1})_ii), R = (U - N(z*)) / (4 D); we enumerate
    that box (rigorous bound, +1 margin) and evaluate N exactly."""
    r = len(Qt)
    SQ = sympy.Matrix(Qt)
    D = int(SQ.det()); assert D > 0
    A = [[int(v) for v in row] for row in (SQ.adjugate()).tolist()]
    assert all(sum(Qt[i][k] * A[k][j] for k in range(r)) == (D if i == j else 0) for i in range(r) for j in range(r))
    def key(x):
        return tuple(sum(A[i][j] * x[j] for j in range(r)) % D for i in range(r))
    reps = {key((0,) * r): (0,) * r}
    frontier = [(0,) * r]
    while frontier:
        nxt = []
        for x in frontier:
            for i in range(r):
                y = list(x); y[i] += 1; y = tuple(y)
                k = key(y)
                if k not in reps:
                    reps[k] = y; nxt.append(y)
        frontier = nxt
    assert len(reps) == D
    par = [Qt[i][i] % 2 for i in range(r)]
    def N(xi):
        return sum(xi[i] * A[i][j] * xi[j] for i in range(r) for j in range(r))
    out = {}
    for k, x in reps.items():
        xi0 = None
        for w in itertools.product((0, 1), repeat=r):
            xi = [x[i] + sum(Qt[i][j] * w[j] for j in range(r)) for i in range(r)]
            if all((xi[i] - par[i]) % 2 == 0 for i in range(r)):
                xi0 = xi; break
        assert xi0 is not None
        Axi = [sum(A[i][j] * xi0[j] for j in range(r)) for i in range(r)]
        zs = [Fraction(-Axi[i], 2 * D) for i in range(r)]
        def xi_of(z):
            return [xi0[i] + 2 * sum(Qt[i][j] * z[j] for j in range(r)) for i in range(r)]
        U = N(xi_of([round(v) for v in zs]))
        xs = [xi0[i] + 2 * sum(Qt[i][j] * zs[j] for j in range(r)) for i in range(r)]
        Nstar = sum(xs[i] * Fraction(A[i][j]) * xs[j] for i in range(r) for j in range(r))
        R = (U - Nstar) / (4 * D)
        bounds = [math.isqrt(int(R * Fraction(A[i][i], D)) + 1) + 1 for i in range(r)]
        best = U
        for z in itertools.product(*[range(math.floor(zs[i]) - bounds[i], math.floor(zs[i]) + bounds[i] + 2) for i in range(r)]):
            v = N(xi_of(z))
            if v < best:
                best = v
        out[k] = (Fraction(best, D) - r) / 4
    return out, D, reps


def mvector(Qt, D):
    """Exact m-values indexed by t in Z/D via a generator of the (cyclic) discriminant group."""
    mins, Dq, reps = class_minima(Qt)
    assert Dq == D
    # find a generator class of Z^r/Qt (cyclic group of order D)
    r = len(Qt)
    keys = list(mins)
    # cyclic structure: map class key k -> integer t with k = t * k_gen
    def add(k1, k2):
        return tuple((a + b) % D for a, b in zip(k1, k2))
    gen = None
    for kk in keys:
        acc = tuple([0] * r); seen = set()
        for t in range(D):
            seen.add(acc); acc = add(acc, kk)
        if len(seen) == D:
            gen = kk; break
    if gen is None:
        return None   # discriminant group of Qt not cyclic: no isomorphism to the cyclic H_1(Sigma_2) exists
    idx = {}
    acc = tuple([0] * r)
    for t in range(D):
        idx[acc] = t; acc = add(acc, gen)
    mvec = [None] * D
    for k, v in mins.items():
        mvec[idx[k]] = v
    return mvec


def admits(mvec, dvec, D):
    passing = []
    for a in range(1, D):
        if math.gcd(a, D) != 1:
            continue
        ok = True
        for t in range(D):
            m = mvec[t]; d = dvec[(a * t) % D]
            if m < d or (m - d) % 2 != 0:
                ok = False; break
        if ok:
            passing.append(a)
    return passing


def forms(D, n_even=2, AMAX=60):
    """All Q = [[x, 2a],[2a, y]] with x = 2m1-1, y = 2m2-1, m1, m2 even (x, y = 3 mod 4),
    xy - 4a^2 = D, x <= y, 0 <= a <= AMAX (no reduction; contains every reduced representative)."""
    out = []
    for a in range(0, AMAX + 1):
        N = D + 4 * a * a
        for x in range(3, int(math.isqrt(N)) + 1, 4):
            if N % x == 0:
                y = N // x
                if y % 4 == 3 and x <= y:
                    out.append(((x + 1) // 2, (y + 1) // 2, a))
    return out


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('djson'); ap.add_argument('sigma', type=int); ap.add_argument('out')
    ap.add_argument('--first-pass', action='store_true', help='stop at the first passing form')
    ap.add_argument('--reduced-only', action='store_true', help='only 0 <= a < m1 <= m2')
    args = ap.parse_args()
    rec = json.load(open(args.djson))
    D = int(rec['det'])
    raw = {int(k) % D: Fraction(v) for k, v in rec['d'].items()}
    assert set(raw) == set(range(D))
    assert abs(args.sigma) == 4, 'this script implements the n = 2 case'
    dvec = raw if args.sigma > 0 else {t: -v for t, v in raw.items()}   # orient so that sigma = +2n
    assert dvec[0] == -1, ('spin correction term must be -sigma/4 after orientation', dvec[0])
    AMAX = D // 2 + 2
    results = []
    for (m1, m2, a) in forms(D, AMAX=AMAX):
        if args.reduced_only and not (a < m1 <= m2):
            continue
        M = [[m1, a], [a, m2]]
        Qt = [[M[0][0], M[0][1], 1, 0], [M[1][0], M[1][1], 0, 1], [1, 0, 2, 0], [0, 1, 0, 2]]
        assert is_pos_def(Qt)
        mvec = mvector(Qt, D)
        if mvec is None:
            results.append(dict(m1=m1, m2=m2, a=a, reduced=(a < m1 <= m2), passes=False, units=[], reversed_orientation_passes=False, note='non-cyclic discriminant group'))
            print(f'  (m1,m2,a)=({m1},{m2},{a}) reduced={a < m1 <= m2} passes=False (non-cyclic discriminant group)', flush=True)
            continue
        p = admits(mvec, dvec, D)
        q = admits(mvec, {t: -v for t, v in dvec.items()}, D)
        results.append(dict(m1=m1, m2=m2, a=a, reduced=(a < m1 <= m2), passes=len(p) > 0, units=p[:8], reversed_orientation_passes=len(q) > 0))
        print(f'  (m1,m2,a)=({m1},{m2},{a}) reduced={a < m1 <= m2} passes={len(p) > 0} units={p[:4]} reversed={len(q) > 0}', flush=True)
        if args.first_pass and p:
            break
    json.dump(dict(det=D, sigma=args.sigma, amax=AMAX, forms=results), open(args.out, 'w'), indent=1)
    npass = [r for r in results if r['passes']]
    print(f"{args.djson.split('/')[-1]}: D={D} sigma={args.sigma} forms tested={len(results)} "
          f"(reduced {sum(r['reduced'] for r in results)}) passing={[(r['m1'], r['m2'], r['a']) for r in npass][:6]} "
          f"-> {'NOT obstructed' if npass else 'OBSTRUCTED (u >= 3)'}")
