import json, itertools, warnings, sys
warnings.filterwarnings('ignore')
from multiprocessing import Pool
from database_knotinfo import link_list
NEW = [l.strip() for l in open('new47.txt')]
by = {r['name']: r for r in link_list()[1:]}
def change(pd, idxs):
    n2 = 2 * len(pd); new = [list(r) for r in pd]
    for i in idxs:
        a, b, c, d = pd[i]
        assert (c - a) % n2 == 1
        new[i] = [b, c, d, a] if (d - b) % n2 == 1 else [d, a, b, c]
    return new
def work(name):
    import spherogram, snappy
    from knot_floer_homology import pd_to_hfk
    pd = json.loads(by[name]['pd_notation'])
    for k in (1, 2, 3):
        for T in itertools.combinations(range(len(pd)), k):
            newpd = change(pd, T)
            L = spherogram.Link(newpd)
            V = L.seifert_matrix(); n = len(V)
            if n:
                from sympy.polys.matrices import DomainMatrix
                import sympy
                M = DomainMatrix([[sympy.ZZ(V[i][j] + V[j][i]) for j in range(n)] for i in range(n)], (n, n), sympy.ZZ)
                if abs(int(M.det())) != 1:
                    continue
            L.simplify('global')
            if len(L.crossings) == 0:
                G = snappy.Link(newpd).exterior().fundamental_group()
                h = pd_to_hfk(newpd)
                ok = G.num_generators() == 1 and G.num_relators() == 0 and h['seifert_genus'] == 0 and h['total_rank'] == 1
                return dict(name=name, k=k, triple=list(T), pi1_Z=G.num_generators() == 1 and G.num_relators() == 0, hfk_unknot=h['total_rank'] == 1, ok=ok)
    return dict(name=name, k=None)
if __name__ == '__main__':
    with Pool(11) as p:
        res = p.map(work, NEW)
    json.dump(res, open('upper47.json', 'w'), indent=1)
    for r in res: print(r)
    print('certified', sum(1 for r in res if r.get('ok')), '/', len(res))
