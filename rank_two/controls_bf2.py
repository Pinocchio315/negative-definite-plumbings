"""Controls for the rank-two test of bf2.py (signed form, as used for the 66 knots):
  (a) negative controls: every knot with u(K) <= 2 in the v2 table (data/v2_final_u_table.json) or with tabulated
      u_a(K) <= 2 in KnotInfo; none may be obstructed;
  (b) positive controls: the 25 knots of Borodzik--Friedl section 6.2 (published_negative_controls.json)
      and their mirrors (Seifert matrix -V^T); all must be obstructed.
Usage: controls_bf2.py          Output: controls_bf2.json
"""
import json
from multiprocessing import Pool
from database_knotinfo import link_list
from bf2 import analyse

T = json.load(open('../data/v2_final_u_table.json'))
PUB = json.load(open('published_negative_controls.json'))


def work(item):
    name, V = item
    try:
        return name, analyse(V, mode='signed', reduced=True)['obstructed']
    except Exception as e:
        return name, repr(e)


if __name__ == '__main__':
    neg = [(r['name'], json.loads(r['seifert_matrix'])) for r in link_list()[1:]
           if r['seifert_matrix'] and r['name'] in T and (T[r['name']][1] <= 2 or r['unknotting_number_algebraic'] in ('0', '1', '2'))]
    pos = [(p['name'], p['V']) for p in PUB] + [(p['name'] + '*', [[-x for x in col] for col in zip(*p['V'])]) for p in PUB]
    with Pool(11) as pool:
        rn = dict(pool.map(work, neg, chunksize=8))
        rp = dict(pool.map(work, pos))
    out = dict(negative=dict(tested=len(rn), obstructed=sorted(n for n, v in rn.items() if v is True),
                             errors=sorted(n for n, v in rn.items() if not isinstance(v, bool))),
               positive=dict(tested=len(rp), obstructed=sum(v is True for v in rp.values()),
                             not_obstructed=sorted(n for n, v in rp.items() if v is not True)))
    json.dump(out, open('controls_bf2.json', 'w'), indent=1)
    print(json.dumps(out)[:600])
