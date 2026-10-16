"""run_fresh.one on the 13 knots recognised in the second recognition pass (recognise_retry.json)."""
import json, os, sys
from multiprocessing import Pool
from collections import Counter
import run_fresh as RF
if __name__ == '__main__':
    R = json.load(open(os.path.join(RF.TOOL, 'recognise_retry.json')))
    items = [(n, sorted(r['names'])) for n, r in R.items() if r['names']]
    with Pool(8) as p:
        out = dict(p.map(RF.one, items))
    json.dump(out, open('run_retry13.json', 'w'), indent=1, default=str)
    for n, r in sorted(out.items()):
        print(n, r['u'], r['sigma'], r['kinds'], r['applicability'], r.get('verdict'), r.get('cw_signs'), r.get('independent_ni_wu_survivors'))
    print(dict(Counter(r.get('verdict') for r in out.values())))
