import json, time, warnings
warnings.filterwarnings('ignore')
from multiprocessing import Pool
from database_knotinfo import link_list
from bf2 import analyse
RANGES = {'[2,3]', '[1,3]', '[2,4]', '[2,5]', '[1,4]'}
def work(r):
    V = json.loads(r['seifert_matrix'])
    try:
        a = analyse(V, mode='signed', reduced=True)
        return dict(name=r['name'], u=r['unknotting_number'], alt=r['name'].split('_')[0][-1] if '_' in r['name'] else '',
                    D=a['D'], sigma=a['sigma'], obstructed=a['obstructed'],
                    patterns={p: i['verdict'] for p, i in a['patterns'].items()})
    except Exception as e:
        return dict(name=r['name'], error=repr(e))
if __name__ == '__main__':
    k = link_list()
    rows = [r for r in k[1:] if r['unknotting_number'] in RANGES and r['seifert_matrix']]
    print('to test', len(rows), flush=True)
    with Pool(11) as pool, open('extend.jsonl', 'w') as out:
        for res in pool.imap_unordered(work, rows, chunksize=4):
            out.write(json.dumps(res) + '\n')
    print('done')
