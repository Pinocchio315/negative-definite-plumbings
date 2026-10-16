"""Direct check for 13a4770, 13a4787, 13a4801: reduced alternating KnotInfo diagram, all crossings of one sign
(special alternating), |sigma| = 8; no set of <= 4 crossing changes gives the unknot (exact det filter, then
HFK total rank); a 5-set that does (HFK rank 1 and pi_1 = Z) is recorded as the upper-bound certificate."""
import json, os, sys, itertools, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'lib'))
from exact_helpers import change
import spherogram, snappy, sympy
from knot_floer_homology import pd_to_hfk
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}
def det(pd):
    V = spherogram.Link(pd).seifert_matrix(); n = len(V)
    return abs(int(sympy.Matrix([[V[i][j] + V[j][i] for j in range(n)] for i in range(n)]).det())) if n else 1
out = {}
for K in ('13a_4770', '13a_4787', '13a_4801'):
    pd = json.loads(BY[K]['pd_notation']); n2 = 2 * len(pd)
    signs = {1 if (X[1] - X[3]) % n2 == 1 else -1 for X in pd}
    assert spherogram.Link(pd).is_alternating() and len(signs) == 1 and abs(int(BY[K]['signature'])) == 8
    hits = {}
    for r in range(1, 6):
        for S in itertools.combinations(range(13), r):
            q = change(pd, list(S))
            if det(q) != 1: continue
            if pd_to_hfk(q)['total_rank'] == 1:
                hits.setdefault(r, []).append(S)
                if r == 5: break
        if r == 5 and 5 in hits: break
    S5 = hits[5][0]; q = change(pd, list(S5)); G = snappy.Link(q).exterior().fundamental_group()
    out[K] = dict(sign=signs.pop(), sigma=int(BY[K]['signature']), unknotting_sets_up_to_4=sum(len(hits.get(r, [])) for r in range(1, 5)),
                  five_set=S5, five_set_pi1=(G.num_generators(), G.num_relators()), KnotInfo_u=BY[K]['unknotting_number'])
    print(K, out[K], flush=True)
json.dump(out, open('verify3.json', 'w'), indent=1)
