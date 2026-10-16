"""Starting ranges of the unknotting numbers, as used by the computations.

load_ranges('baseline'): the final table of arXiv:2609.09861v2 (data/v2_final_u_table.json; KnotInfo 2026.8.1 with
the 104 DKT-adopted upper bounds rolled back, plus the results of v2), together with the 66 values u = 3 of
data/new66_u3.txt and the upper bound u(13n_1669) <= 3.  The variants 'strict' and 'full' also add the entries of
data/new_values_20261015.json (the 43 new values of this paper, and the three values u = 5).
The function is copied verbatim from the project module in which it was written.
"""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
VARIANT = 'baseline'
NEW_VALUES = os.path.join(DATA, 'new_values_20261015.json')


def load_ranges(variant=None):
    """Baseline = v2's final table (data/v2_final_u_table.json): KnotInfo 2026.8.1 with the 104 DKT-adopted
    upper bounds rolled back (data/knotinfo_rollback_list_2026-09-07.csv), plus v2's results.  It matches
    the arXiv v2 appendices exactly.  Then this project's 66 values and u(13n_1669) <= 3, and (variants 'strict'
    and 'full') the entries of data/new_values_20261015.json.
    (The KnotInfo 2026-09-08 snapshot differs from it only on 10 DKT entries, all rolled back.)"""
    variant = variant or VARIANT
    assert variant in ('baseline', 'strict', 'full'), variant
    import csv
    T = json.load(open(os.path.join(DATA, 'v2_final_u_table.json')))
    v2 = json.load(open(os.path.join(DATA, 'arxiv_2609.09861v2_results.json')))
    v2names = {n for ns in v2['exact'].values() for n in ns} | set(v2['improved'])
    dkt = {r['knot'] for r in csv.DictReader(open(os.path.join(DATA, 'knotinfo_rollback_list_2026-09-07.csv')))}
    new66 = [l.strip() for l in open(os.path.join(DATA, 'new66_u3.txt')) if l.strip()]
    R = {}
    for n, (lo, hi) in T.items():
        src = 'V2' if n in v2names else ('KI-rollback' if n in dkt else 'KI')
        R[n] = dict(lo=lo, hi=hi, lo_src=src, hi_src=src)
    def add(name, lo, hi, src):
        r = R[name]
        if lo is not None and lo > r['lo']:
            r['lo'], r['lo_src'] = lo, src
        if hi is not None and hi < r['hi']:
            r['hi'], r['hi_src'] = hi, src
    for n in new66:
        add(n, 3, 3, 'NEW')
    add('13n_1669', None, 3, 'NEW')
    if variant != 'baseline':
        for e in json.load(open(NEW_VALUES))['values']:
            if e['strict'] or variant == 'full':
                assert R[e['knot']]['hi'] == e['lo'], (e['knot'], R[e['knot']], e['lo'])   # exact: lo meets the existing hi
                add(e['knot'], e['lo'], None, e['source'])
    bad = {n: r for n, r in R.items() if r['lo'] > r['hi']}
    assert not bad, f'inconsistent ranges: {bad}'
    return R
