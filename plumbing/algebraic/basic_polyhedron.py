#!/usr/bin/env python
"""Basic polyhedron of a knot diagram (Conway): repeatedly merge two vertices that bound a bigon face.
Vertices carry the cyclic (counterclockwise) order of their edges, read from the PD code. A diagram is algebraic
(basic polyhedron 1*) iff it reduces to a single vertex. Validation: KnotInfo Conway notation (<= 12 crossings,
'*' = polyhedral). Usage: basic_polyhedron.py validate | KNOT [KNOT ...]
"""
import json, re, sys


def reduce_diagram(pd):
    nodes = {i: list(x) for i, x in enumerate(pd)}
    while len(nodes) > 1:
        ends = {}
        for v, es in nodes.items():
            for e in es:
                ends.setdefault(e, []).append(v)
        merged = False
        for a, A in nodes.items():
            for i in range(4):
                f, e = A[i], A[(i + 1) % 4]
                if f == e or len(ends.get(f, [])) != 2 or len(ends.get(e, [])) != 2:
                    continue
                bs = {w for w in ends[f] if w != a} & {w for w in ends[e] if w != a}
                for b in bs:
                    B = nodes[b]
                    js = [j for j in range(4) if B[j] == e and B[(j + 1) % 4] == f]
                    if not js:
                        continue
                    j = js[0]
                    nodes[a] = [A[(i + 2) % 4], A[(i + 3) % 4], B[(j + 2) % 4], B[(j + 3) % 4]]
                    del nodes[b]
                    merged = True
                    break
                if merged:
                    break
            if merged:
                break
        if not merged:
            break
    return len(nodes)


def is_algebraic_conway(cw):
    """Conway notation omits the basic polyhedra 1* (algebraic) and 6*; a 6* notation contains '.' or ':'."""
    body = cw.strip('[]')
    return not any(ch in body for ch in '*:.')


if __name__ == '__main__':
    from database_knotinfo import link_list
    BY = {r['name']: r for r in link_list()[1:]}
    if sys.argv[1] == 'validate':
        from collections import Counter
        c = Counter(); bad = []
        for n, r in BY.items():
            cw = (r.get('conway_notation') or '').strip()
            if not cw or not re.match(r'\d+[an]?_', n) or int(re.match(r'(\d+)', n).group(1)) > 12:
                continue
            k = reduce_diagram(json.loads(r['pd_notation']))
            key = ('alternating' if r['alternating'] == 'Y' else 'non-alternating', is_algebraic_conway(cw), k == 1)
            c[key] += 1
            if is_algebraic_conway(cw) != (k == 1) and len(bad) < 8:
                bad.append((n, cw, k))
        for key in sorted(c): print('type %-15s | Conway algebraic %-5s | diagram reduces to 1* %-5s | %d' % (*key, c[key]))
        print('examples of disagreement:', bad)
    else:
        for n in sys.argv[1:]:
            print(n, reduce_diagram(json.loads(BY[n]['pd_notation'])))
