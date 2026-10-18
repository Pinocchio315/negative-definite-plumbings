#!/usr/bin/env python
"""KnotInfo entries, in version 2026.10.5, of the 112 knots whose unknotting numbers the paper determines (Theorem 1.2).

Input: database_knotinfo/csv_data/knotinfo_data_complete.csv, extracted from the wheel
database_knotinfo-2026.10.5-py3-none-any.whl of the PyPI package database_knotinfo 2026.10.5 (wheel SHA-256
7cb8cd135ddfdf439d1b799bb4597ed74bcee7b6e3a880b169966f428944a1ca, uploaded 2026-10-05T17:07:16Z; downloaded
2026-10-15). The SHA-256 of the CSV is checked against the one listed in the RECORD file of the wheel. The knots are
those of new_values_20261015.json and new66_u3.txt.
Usage: python -B knotinfo_2026_10_05.py PATH/knotinfo_data_complete.csv
Output: knotinfo_2026-10-05_new_values.json (the SHA-256 of the CSV and the unknotting_number field of each knot).
"""
import csv, hashlib, json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
WHEEL = 'database_knotinfo-2026.10.5-py3-none-any.whl'
WHEEL_SHA256 = '7cb8cd135ddfdf439d1b799bb4597ed74bcee7b6e3a880b169966f428944a1ca'
CSV_SHA256 = 'eb511ebc61204bc1c6d700d9da42257dcbb0369cb6756e89ed978a7821867c7d'     # RECORD of the wheel

if __name__ == '__main__':
    path = sys.argv[1]
    knots = [e['knot'] for e in json.load(open(os.path.join(HERE, 'new_values_20261015.json')))['values']]
    knots += [l.strip() for l in open(os.path.join(HERE, 'new66_u3.txt')) if l.strip()]
    assert len(knots) == len(set(knots)) == 112
    sha = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    assert sha == CSV_SHA256, sha
    csv.field_size_limit(1 << 30)
    found = {}
    with open(path, newline='', encoding='utf-8') as f:
        for row in csv.DictReader(f, delimiter='|'):
            if row['name'] in knots:
                found[row['name']] = row['unknotting_number']
    assert set(found) == set(knots), sorted(set(knots) - set(found))
    out = dict(source=f'database_knotinfo/csv_data/knotinfo_data_complete.csv of {WHEEL} (PyPI package database_knotinfo 2026.10.5)',
               wheel_sha256=WHEEL_SHA256, csv_sha256=sha,
               unknotting_number={k: found[k] for k in sorted(found)})
    json.dump(out, open(os.path.join(HERE, 'knotinfo_2026-10-05_new_values.json'), 'w'), indent=1)
    print('csv sha256', sha, '| entries', len(found), '| all ranges:', all(v.startswith('[') for v in found.values()))
