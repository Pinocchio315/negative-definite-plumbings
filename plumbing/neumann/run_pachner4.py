"""run_fresh.one on the knots recognised through the Pachner graph (pachner_recognise.json)."""
import json, os
from multiprocessing import Pool
import run_fresh as RF
if __name__ == '__main__':
    P_ = json.load(open(os.path.join(RF.TOOL, 'pachner_recognise.json')))
    items = [(n, r['names']) for n, r in P_.items() if r['names']]
    with Pool(4) as p:
        out = dict(p.map(RF.one, items))
    json.dump(out, open('run_pachner4.json', 'w'), indent=1, default=str)
    for n, r in sorted(out.items()):
        print(n, 'u', r['u'], 'sigma', r['sigma'], 'D', r['D'], r['kinds'], '|', r['applicability'], '| p+', r['pplus'], '|', r.get('verdict'),
              '| CW', r.get('cw_signs'), 'link', r.get('link_orientation'), '| NiWu', r.get('independent_ni_wu_survivors'))
