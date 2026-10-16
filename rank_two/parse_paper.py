import re, json
t = open('../data/arxiv_2609.09861v2_main.tex').read()
def section(label):
    i = t.index('\\label{%s}' % label)
    j = t.find('\\section{', i)
    return t[i:j if j > 0 else len(t)]
NAME = r'(?:(\d{1,2})_\{(\d+)\}|(\d{2}[an])(\d+))'
def norm(m):
    if m.group(1): return f"{m.group(1)}_{m.group(2)}"
    return f"{m.group(3)}_{m.group(4)}"
def names_in(sec):
    body = sec[sec.index('\\endlastfoot'):]
    return [norm(m) for m in re.finditer(NAME, body)]
out = {}
for lab, u in [('app:u5', 5), ('app:u4', 4), ('app:u3', 3), ('app:u2', 2)]:
    out[u] = names_in(section(lab))
imp = {}
body = section('app:improve'); body = body[body.index('\\endlastfoot'):]
for m in re.finditer(NAME + r'[^&]*&\s*\$\[(\d),(\d)\]\$\s*&\s*\$\[(\d),(\d)\]\$', body):
    imp[norm(m)] = ((int(m.group(5)), int(m.group(6))), (int(m.group(7)), int(m.group(8))))
dich = names_in(section('app:dichotomy'))
print({u: len(v) for u, v in out.items()}, 'improved', len(imp), 'dichotomy', len(dich))
json.dump({'exact': {str(u): v for u, v in out.items()}, 'improved': {k: v for k, v in imp.items()}, 'dichotomy': dich}, open('paper_results.json', 'w'))
