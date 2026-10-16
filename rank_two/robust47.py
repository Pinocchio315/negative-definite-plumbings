import json, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
from database_knotinfo import link_list
from bf2 import analyse, LinkingForm, candidates_reduced, isometric
NEW = [l.strip() for l in open('new47.txt')]
by = {r['name']: r for r in link_list()[1:]}
def noconds(V):
    n = len(V); M = [[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]
    LF = LinkingForm(M); D = LF.D; s = 0
    for detC in (D, -D):
        for (a, b, c) in candidates_reduced(D, detC):
            s += isometric(LF, a, b, c, detC)[0]
    return s
def work(name):
    import spherogram
    r = by[name]
    Vk = json.loads(r['seifert_matrix'])
    Vp = spherogram.Link(json.loads(r['pd_notation'])).seifert_matrix()
    out = dict(name=name)
    for lab, V in (('ki', Vk), ('pd', Vp)):
        out[lab + '_signed'] = analyse(V, 'signed', reduced=True)['obstructed']
        out[lab + '_mod2'] = analyse(V, 'mod2', reduced=True)['obstructed']
        out[lab + '_noconds_survivors'] = noconds(V)
    a = analyse(Vk, 'signed', reduced=True)
    out['sigma'] = a['sigma']; out['D'] = a['D']
    out['why'] = {p: (i['verdict'], i['reject_reasons']) for p, i in a['patterns'].items()}
    return out
if __name__ == '__main__':
    with Pool(11) as p:
        res = p.map(work, NEW)
    json.dump(res, open('robust47.json', 'w'), indent=1)
    for r in res: print({k: r[k] for k in r if k != 'why'})
