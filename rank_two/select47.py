import json
rows = [json.loads(l) for l in open('extend.jsonl')]
P = json.load(open('paper_results.json'))
exact = {n for ns in P['exact'].values() for n in ns}
imp = {k: v[1] for k, v in P['improved'].items()}
new = []
for r in rows:
    if r.get('obstructed') and r['name'].startswith('13n') and r['name'] not in exact:
        lo, hi = imp.get(r['name'], [int(x) for x in r['u'].strip('[]').split(',')])
        if hi == 3: new.append(r['name'])
open('new47.txt', 'w').write('\n'.join(new))
print(len(new))
