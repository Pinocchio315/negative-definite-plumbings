"""Diagram-level lower bounds u>=2 for every one of the 91 neighbour diagrams of 13n1587."""
import json, sys, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, __import__('os').path.dirname(__file__))
from multiprocessing import Pool
import sympy
from sympy.matrices.normalforms import smith_normal_form
from bf2 import LinkingForm, charpoly_signature

def rank1_survives(V):
    """BF Thm 1.1 with n=1 (Lickorish refined): exists c, |c|=D, c = e mod 4, sign(c)=sigma+e,
    l([c]) ~= 2 lambda  <=>  H cyclic with a generator g, 2lambda(g,g) = 1/c."""
    n = len(V); M = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    sigma = charpoly_signature(M); LF = LinkingForm(M); D = LF.D
    for e in (1, -1):
        sc = sigma + e
        if abs(sc) != 1: continue
        c = sc * D
        if c % 4 != e % 4: continue
        target = (1 if c > 0 else -1) % D
        for i in LF.by_self2.get(target, []):
            # generator: order D  <=> Bx has some entry with gcd(.,D)=1 after reduction -> check order
            bx = LF.Bx[i]; o = 1; ok = True
            from math import gcd
            g = 0
            for t in bx: g = gcd(g, t)
            if gcd(g, D) == 1:
                return True, sigma
    return False, sigma

def work(r):
    import spherogram
    L = spherogram.Link(r['pd'])
    out = dict(name=r['ids'][0], ref=r['ref'], crossing=r['crossing'])
    if r['ids'][0] == 'unknot':
        out['lower'] = 'UNKNOT!'; return out
    parts = L.deconnect_sum() if True else [L]
    if len(parts) > 1:
        out['lower'] = 'composite (Scharlemann)'; return out
    V = L.seifert_matrix(); n = len(V)
    M = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    snf = smith_normal_form(sympy.Matrix(M), domain=sympy.ZZ)
    inv = [abs(int(snf[i, i])) for i in range(n) if abs(int(snf[i, i])) != 1]
    sigma = charpoly_signature(M)
    out.update(sigma=sigma, H1=inv)
    if abs(sigma) >= 4:
        out['lower'] = '|sigma|>=4'
    elif len(inv) >= 2:
        out['lower'] = 'H1 needs >=2 generators (Wendt)'
    else:
        surv, _ = rank1_survives(V)
        out['lower'] = 'none classical (rank-1 survives)' if surv else 'rank-1 BF/Lickorish obstruction'
    return out

if __name__ == '__main__':
    B = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'bjset_13n_1587.json'))
    with Pool(11) as p:
        res = p.map(work, B['changes'])
    from collections import defaultdict
    agg = defaultdict(set)
    for r in res: agg[r['name']].add((r['lower'], str(r.get('sigma')), str(r.get('H1'))))
    for k, v in sorted(agg.items()): print(f'{k:12s}', sorted(v))
    json.dump(res, open('neighbor_lower_13n_1587.json', 'w'), indent=0)
