#!/usr/bin/env python
"""Generate the LaTeX tables of the paper from the deposited results (no number is typed by hand).

Inputs: data/*.json, plumbing/**/*.json (Sigma_2 classes, normal forms, tests), rank_two/*.json,
bernhard_jablan/results/*.json.
Outputs: paper/tables/*.tex, and paper/tables/numbers.tex with LaTeX macros for the counts used in the text.
Usage: python -B paper/make_tables.py
"""
import json, os, sys
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
CEN = os.path.join(HERE, '..')
TOOL = os.path.join(CEN, 'plumbing')
ALGD = os.path.join(TOOL, 'algebraic')                             # algebraic diagrams (Lemma 2.7)
OUT = os.path.join(HERE, 'tables')
os.makedirs(OUT, exist_ok=True)
from database_knotinfo import link_list
BY = {r['name']: r for r in link_list()[1:]}


def J(*p):
    return json.load(open(os.path.join(*p)))


def kn(name):                      # 13n_1587 -> \kn{13n}{1587};  10_124 -> $10_{124}$
    a, b = name.split('_')
    return f'\\kn{{{a}}}{{{b}}}' if a[-1] in 'an' else f'${a}_{{{b}}}$'


def short_sfs(nm):
    return (nm or '').replace('SFS ', '').replace(' U/m ', r' $\cup$ ').replace(' U/n ', r' $\cup$ ')


def key(n):
    a, b = n.split('_')
    return (int(''.join(c for c in a if c.isdigit())), a, int(b))


vals = J(CEN, 'data', 'new_values_20261015.json')['values']
star = {**J(TOOL, 'star_u1.json'), **J(TOOL, 'star_u1_orb.json')}
gu = J(TOOL, 'graph_u1.json'); fresh = J(TOOL, 'neumann', 'run_fresh.json'); owt = J(TOOL, 'owens_u3_targets.json')
fresh_names = J(TOOL, 'recognise_more_fresh.json')
nf = J(TOOL, 'neumann', 'rederive_all.json')          # the 43 values re-derived by the normal-form method alone


def sigma2_of(n):
    if n in star:
        return star[n]['sigma2']
    if n in gu and gu[n].get('plumbing'):
        return gu[n]['plumbing']['regina_name_of_plumbing']
    if n in owt and owt[n] and owt[n].get('plumbing'):
        return owt[n]['plumbing']['regina_name_of_plumbing']
    if n in fresh_names:
        return fresh_names[n]['names'][0]
    return ''


def test_of(e):
    t = e.get('test', '')
    if 'HF' in t:
        return r'$\mathrm{HF}_{\mathrm{red}}$'
    if 'Owens' in t:
        return r'Owens $(1,1)$ + Traczyk'
    return '$d$'


rows2 = sorted((e for e in vals if e['strict'] and e['lo'] == 2), key=lambda e: key(e['knot']))
rows3 = sorted((e for e in vals if e['strict'] and e['lo'] == 3), key=lambda e: key(e['knot']))
rows5 = sorted((e for e in vals if not e['strict']), key=lambda e: key(e['knot']))
assert all(nf[e['knot']]['verdict'].startswith('OBSTRUCTED') and nf[e['knot']]['cw_consistent'] for e in rows2 + rows3)
assert all(nf[e['knot']].get('plumbing_name_matches') for e in rows2 + rows3)     # identified by an independent triangulation
# the six exclusions of u = 1 by reduced Floer homology: every identification surviving the correction-term test fails
# the mapping cone condition (neumann/hfred_certificates.py)
HFC = J(TOOL, 'neumann', 'hfred_certificates.json')
assert sorted(HFC) == sorted(e['knot'] for e in rows2 if nf[e['knot']]['verdict'] == 'OBSTRUCTED_HF')
assert all(c['certified'] and c['survivors'] and all(x['fails'] for x in c['survivors']) for c in HFC.values())
# KnotInfo version 2026.10.5 records all 112 values as ranges containing them (data/knotinfo_2026_10_05.py)
KI10 = J(CEN, 'data', 'knotinfo_2026-10-05_new_values.json')['unknotting_number']


def longtable(cols, caption, label, head, body, desc):
    """A longtable that marks page breaks as in arXiv:2609.09861v2: 'continued from the previous page' above and
    'continued on the next page' below."""
    k = len(cols)
    L = [r'\begin{longtable}{' + cols + '}', r'\caption{' + caption + r'}\label{' + label + r'}\\', r'\toprule', head + r'\\', r'\midrule',
         r'\endfirsthead',
         r'\multicolumn{' + str(k) + r'}{@{}l@{}}{\emph{' + desc + r', continued from the previous page.}}\\[2pt]',
         r'\toprule', head + r'\\', r'\midrule', r'\endhead',
         r'\bottomrule', r'\multicolumn{' + str(k) + r'}{@{}r@{}}{\emph{continued on the next page}}\\', r'\endfoot',
         r'\bottomrule', r'\endlastfoot']
    return '\n'.join(L + body + [r'\end{longtable}']) + '\n'


def table(rows, caption, label, desc):
    body = [f"{kn(e['knot'])} & {BY[e['knot']]['determinant']} & {BY[e['knot']]['signature']} & "
            f"{{\\footnotesize {short_sfs(sigma2_of(e['knot']))}}} & {test_of(e)} & {nf[e['knot']]['vertices']}\\\\" for e in rows]
    return longtable('lrrllr', caption, label, r'Knot & $\det$ & $\sigma$ & $\dbc(K)$ (Regina) & Obstruction & $|V|$', body, desc)


open(os.path.join(OUT, 'new_u2.tex'), 'w').write(table(rows2, r'New exact values $u(K)=2$ (previous range $[1,2]$). $|V|$ is the number of vertices of the Neumann normal form.', 'tab:new-u2', r'New exact values $u(K)=2$'))
open(os.path.join(OUT, 'new_u3.tex'), 'w').write(table(rows3, r"New exact values $u(K)=3$ (previous range $[2,3]$). Owens's theorem with normal-form correction terms excludes two crossing changes of opposite signs, and Traczyk's criterion excludes two changes of the same sign (Section~\ref{sec:u2}).", 'tab:new-u3', r'New exact values $u(K)=3$'))
# the 66 values u = u_a = 3 (data/new66_u3.txt): H_1 of the double branched cover from the Seifert matrix
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


def h1(n):
    S = Matrix(json.loads(BY[n]['seifert_matrix']))
    d = [abs(int(x)) for x in smith_normal_form(S + S.T, domain=ZZ).diagonal() if abs(int(x)) > 1]
    return r' \oplus '.join(rf'\bbz/{x}' for x in d)


def plain_table(names, caption, label, with_h1, halves=1, desc=''):
    """halves = 2 sets the list in two column blocks side by side (first half on the left), separated by a vertical
    rule; such a table is read down the columns, so it is a float that does not break across pages."""
    cols = 'lrr' + ('l' if with_h1 else '')
    head = r'Knot & $\det$ & $\sigma$' + (r' & $H_1(\dbc(K))$' if with_h1 else '')
    cell = lambda n: f"{kn(n)} & {BY[n]['determinant']} & {BY[n]['signature']}" + (f" & ${h1(n)}$" if with_h1 else '')
    blank = ' & '.join([''] * len(cols))
    h = -(-len(names) // halves)
    blocks = [names[i * h:(i + 1) * h] for i in range(halves)]
    body = [' & '.join(cell(b[i]) if i < len(b) else blank for b in blocks) + r'\\' for i in range(h)]
    if halves == 1:
        return longtable(cols, caption, label, head, body, desc)
    L = [r'\begin{table}[htbp]', r'\centering', r'\setlength{\aboverulesep}{0pt}\setlength{\belowrulesep}{0pt}\renewcommand{\arraystretch}{1.05}',
         r'\caption{' + caption + r'}\label{' + label + '}',
         r'\begin{tabular}{' + r'@{\quad}|@{\quad}'.join([cols] * halves) + '}', r'\toprule', ' & '.join([head] * halves) + r'\\', r'\midrule']
    L += body + [r'\bottomrule', r'\end{tabular}', r'\end{table}']
    return '\n'.join(L) + '\n'


rows66 = sorted((l.strip() for l in open(os.path.join(CEN, 'data', 'new66_u3.txt')) if l.strip()), key=key)
NEWVAL = {**{e['knot']: e['lo'] for e in rows2 + rows3 + rows5}, **{n: 3 for n in rows66}}       # the 112 new values
assert len(NEWVAL) == len(rows2) + len(rows3) + len(rows5) + len(rows66) and set(KI10) == set(NEWVAL)
for n, v in NEWVAL.items():
    lo, hi = map(int, KI10[n].strip('[]').split(','))
    assert lo < hi and lo <= v <= hi, (n, KI10[n])
open(os.path.join(OUT, 'u3_66.tex'), 'w').write(plain_table(rows66, r'Knots with $u(K)=u_a(K)=3$, listed down the columns.', 'tab:u3-66', True, halves=2))

# Sigma_2 map of the open non-alternating knots
M = J(TOOL, 'recognise_more.json'); RR = J(TOOL, 'recognise_retry.json'); VH = J(TOOL, 'verify_hyp.json')
JS = J(TOOL, 'jsj.json'); JV = J(TOOL, 'jsj_verify.json'); PR = J(TOOL, 'pachner_recognise.json')
r13 = J(TOOL, 'neumann', 'run_retry13.json'); p4 = J(TOOL, 'neumann', 'run_pachner4.json')
pre = {**J(TOOL, 'neumann', 'precheck_u12.json'), **J(TOOL, 'neumann', 'precheck_u23.json')}
ps0 = J(TOOL, 'neumann', 'precheck_sigma0.json')     # graph-type covers that no plumbing pipeline needed (sigma = 0)


nc3 = J(TOOL, 'neumann', 'normal_case3.json')     # class (3) knots first found by the earlier plumbing searches
UA = J(ALGD, 'u_plumbing.json'); XA = J(ALGD, 'crosscheck_u.json'); AD = J(ALGD, 'alg_diagrams.json')
VA = J(ALGD, 'validate_plumbing.json')               # the construction on every algebraic KnotInfo diagram
NF_SOURCES = (p4, r13, fresh, nf, nc3)               # every verdict in these files comes from a Neumann normal form


def status(n):
    if n in VH and VH[n]['certified']:
        return 'hyperbolic'
    if n in JV and JV[n]['certified_hyperbolic_piece']:
        return 'mixed'
    for src in NF_SOURCES:
        if n in src:
            r = src[n]
            if r['applicability'] == 'NOT APPLICABLE':
                return 'nodef'
            if r['applicability'] != 'APPLICABLE':
                return 'other'
            if str(r.get('verdict', '')).startswith('OBSTRUCTED'):
                return 'obstructed'
            return 'notest' if 'sigma' in str(r.get('verdict', '')) else 'tested'
    for src in (star, gu, owt):          # the earlier searches found no plumbing here; the normal form decides
        if n in src and src[n]:
            assert src[n].get('verdict') == 'NO PLUMBING', n
            return 'nodef' if pre.get(n, {}).get('verdict') == 'NOT APPLICABLE' else 'other'
    if n in UA:                          # not recognised by Regina; normal forms from an algebraic diagram
        assert JS[n]['verdict'].startswith('toroidal, no geometric'), n
        return 'nodef_alg'
    if n in JS and JS[n]['verdict'].startswith('toroidal, no geometric'):
        return 'graph4'
    if n in ps0:
        return {'APPLICABLE': 'notest', 'NOT APPLICABLE': 'nodef'}.get(ps0[n]['verdict'], 'other')
    return 'other'



assert all(r['applicability'] == 'APPLICABLE' and r['cw_consistent'] and r['plumbing_name_matches'] for r in nc3.values())
cols = {}
for lab, f in (('u12', 'reach_u12_orb.json'), ('u23', 'reach_u23_orb.json')):
    cols[lab] = Counter(status(x['name']) for x in J(TOOL, f))
# the 47 non-alternating knots with u = u_a = 3 (previous range [2,3]; plumbing/classify47.py, same pipeline)
C47 = J(TOOL, 'classify47_summary.json')['status']; C47n = J(TOOL, 'classify47_normal.json'); C47v = J(TOOL, 'classify47_verify.json')
v2T0 = J(CEN, 'data', 'v2_final_u_table.json')
assert sorted(C47) == sorted(n for n in rows66 if 'n_' in n) and all(v2T0[n] == [2, 3] for n in C47)
assert not set(C47) & {x['name'] for f in ('reach_u12_orb.json', 'reach_u23_orb.json') for x in J(TOOL, f)}
THMC_OBSTRUCTED = sorted(n for n, s in C47.items() if s == 'negdef')
for n in THMC_OBSTRUCTED:          # Owens's test excludes u = 2 in both orientations; verified as for the knots of Table new-u3
    r, v = C47n[n], C47v[n]
    assert r['test_verdict'] == 'OBSTRUCTED' and r['cw_consistent'] and v['plumbing_name_matches'] and v['blocked'] == [1, -1]
cols['u23'] += Counter('obstructed' if s == 'negdef' else s for s in C47.values())
assert not cols['u12']['other'] and not cols['u23']['other'] and not cols['u12']['graph4'] and not cols['u23']['graph4']
assert cols['u23']['obstructed'] == len(rows3) + len(THMC_OBSTRUCTED) and cols['u12']['obstructed'] == len(rows2)
order = [('hyperbolic', r'hyperbolic (verified)'), ('mixed', r'toroidal, with a hyperbolic JSJ piece'),
         (None, r'graph manifold, no negative-definite plumbing tree:'),
         ('nodef', r'\quad recognised by Regina'), ('nodef_alg', r'\quad from an algebraic diagram'),
         (None, r'$\pm\dbc(K)$ bounds a negative-definite plumbing tree:'),
         ('obstructed', r'\quad $u=1$, resp.\ $u=2$, excluded'), ('tested', r'\quad obstructions inconclusive'),
         ('notest', r'\quad no correction-term test ($\sigma=0$)')]
L = [r'\begin{tabular}{lrr}', r'\toprule', r'$\dbc(K)$ & $u\in[1,2]$ & $u\in[2,3],[2,4]$\\', r'\midrule']
for k, lab in order:
    L.append(f"{lab} & & \\\\" if k is None else f"{lab} & {cols['u12'][k]} & {cols['u23'][k]}\\\\")
L += [r'\midrule', f"total & {sum(cols['u12'].values())} & {sum(cols['u23'].values())}\\\\", r'\bottomrule', r'\end{tabular}']
open(os.path.join(OUT, 'sigma2_map.tex'), 'w').write('\n'.join(L) + '\n')

# numbers used in the text
both = lambda k: cols['u12'][k] + cols['u23'][k]
# controls of Owens's test: Montesinos knots with <= 12 crossings, |sigma| in {2, 4} and exact u in {2, 3}
ctl = {n: r for n, r in J(TOOL, 'owens_u3_controls.json').items() if r and 'error' not in r}
c42 = [r for r in ctl.values() if abs(r['sigma']) == 4 and r['u'] == [2, 2]]
c43 = [r for r in ctl.values() if abs(r['sigma']) == 4 and r['u'] == [3, 3]]
calt = [r for r in ctl.values() if r['alternating'] == 'Y']
assert all(not r['verdict'].startswith('OBSTRUCTED') for r in c42)                       # no false obstruction
assert all(r['v2_goeritz_verdict'].startswith('OBSTRUCTED') for r in calt if r['verdict'].startswith('OBSTRUCTED'))


# Section 4: the rank-two linking test (results/bf_rank2_u3_20261005)
BFD = os.path.join(CEN, 'rank_two')
c66 = J(BFD, 'check66.json')
assert set(c66) == set(l.strip() for l in open(os.path.join(CEN, 'data', 'new66_u3.txt')) if l.strip())
assert all(r['lower']['ki_signed'] and r['lower']['pd_signed'] and r['upper'].get('ok') and r['upper']['k'] == 3 for r in c66.values())
signed_only = [r for r in c66.values() if r['lower']['ki_noconds_survivors'] > 0]
assert all(abs(r['lower']['sigma']) == 4 for r in signed_only)
ext = [json.loads(l) for l in open(os.path.join(BFD, 'extend.jsonl'))]
v2T = J(CEN, 'data', 'v2_final_u_table.json')
ext_obs = [r['name'] for r in ext if r.get('obstructed')]
assert not [n for n in ext_obs if v2T[n][1] <= 2]                                     # no contradiction with u <= 2
bfc = J(BFD, 'controls_bf2.json')
assert not bfc['negative']['obstructed'] and not bfc['positive']['not_obstructed']


# Section 5: the knot 13n1587 (bernhard_jablan/results)
BJR = os.path.join(HERE, '..', 'bernhard_jablan', 'results')
comp_full = J(BJR, 'completeness_13n_1587.json')
comp = comp_full['summary']
nbr = J(BJR, 'neighbor_lower_13n_1587.json')
comp12, bj12 = J(BJR, 'completeness_12n_563.json')['summary'], J(BJR, 'bjset_12n_563.json')


def det_of_pd(pd):                       # exact determinant |det(V + V^T)| of the knot of a diagram
    import spherogram, sympy
    V = spherogram.Link(pd).seifert_matrix()
    return abs(int(sympy.Matrix(len(V), len(V), lambda i, j: V[i][j] + V[j][i]).det()))


# no crossing change of 12n563 gives the unknot: every changed diagram is named and has determinant other than one
assert len(bj12['changes']) == 216 and all(len(x['ids']) == 1 and x['ids'][0] in BY and det_of_pd(x['pd']) != 1
                                           for x in bj12['changes'])
dec = comp['decisions']
# the 57 candidates that share the Alexander polynomial and knot Floer homology of 13n1587 are separated from it by the
# H_1-multisets of connected covers of degree <= 4; the 7 accepted diagrams have the same multisets as 13n1587
cov = {(r['ref'], tuple(r['subset'])): r for r in J(BJR, 'cover_separation_13n_1587.json')['rows']}
cands = lambda tag: {(r['ref'], tuple(c['subset'])) for r in comp_full['references'] for c in r['candidates']
                     if c['decision'].startswith(tag)}
UNRES, ACC = cands('UNRESOLVED'), cands('ACCEPTED')
assert set(cov) == UNRES | ACC and len(UNRES) == dec['UNRESOLVED: same Alexander and HFK as target but not identified as target']
assert all(not cov[k]['cover_invariant_equal'] and cov[k]['first_differing_degree'] <= 4 for k in UNRES)
assert all(cov[k]['cover_invariant_equal'] for k in ACC) and len(ACC) == len(comp['accepted'])
NAME = {'5_2#?3c': r'$5_2\#3_1$'}
WHY = {'9_46': r'$H_1$ needs two generators \cite{Wendt}', '10_38': r'no unknotting crossing \cite{McCoy}',
       '12a_237': r'$|\sigma|=4$', '12n_155': r'$|\sigma|=4$', '5_2#?3c': r'$|\sigma|=4$',
       '11n_66': r'\cite{GordonLuecke06}', '11n_99': r'\cite{GordonLuecke06}',
       '12n_563': r'\cite{Lee}'}
groups = {}
for r in nbr:
    groups.setdefault(r['name'], []).append(r)
assert set(groups) == set(WHY) and sum(map(len, groups.values())) == len(nbr)
L = [r'\begin{table}[htbp]', r'\centering',
     r'\caption{The knots obtained by changing one crossing in a minimal diagram of \kn{13n}{1587}, and the reason why '
     r'$u\ge2$. For each knot, \emph{Changes} is the number of crossing changes that give it, among the '
     r'\nNbrChanges{} changes of one crossing in the \nMinDiag{} minimal diagrams obtained from reference diagrams in '
     r'Section~\ref{sec:13n1587}. Knots are named up to mirror image.}\label{tab:nbr}',
     r'\begin{tabular}{lrrll}', r'\toprule',
     r'Knot & Changes & $|\sigma|$ & $H_1(\dbc(K))$ & $u\ge2$ by\\', r'\midrule']
for n, rs in sorted(groups.items(), key=lambda x: (-len(x[1]), x[0])):
    sig = {abs(x['sigma']) for x in rs}; assert len(sig) == 1
    hh = r' \oplus '.join(rf'\bbz/{d}' for d in rs[0]['H1'])
    L.append(f"{NAME.get(n, kn(n))} & {len(rs)} & {sig.pop()} & ${hh}$ & {WHY[n]}\\\\")
L += [r'\midrule', f"total & {len(nbr)} & & & \\\\", r'\bottomrule', r'\end{tabular}', r'\end{table}']
open(os.path.join(OUT, 'nbr_13n1587.tex'), 'w').write('\n'.join(L) + '\n')


def bad_vertices(p):
    deg = Counter(v for e in p['edges'] for v in e)
    return sum(1 for i, w in enumerate(p['weights']) if w > -deg[i])


# class (2): two different accepted trees, i.e. the normal forms of Sigma_2 and -Sigma_2, neither negative definite
BO = J(TOOL, 'neumann', 'both_orientations.json')
CLASS2 = [x['name'] for f in ('reach_u12_orb.json', 'reach_u23_orb.json') for x in J(TOOL, f) if status(x['name']) == 'nodef']
CLASS2 += [n for n, s in C47.items() if s == 'nodef']
assert len(CLASS2) == both('nodef') and all(BO[n]['distinct_trees'] == 2 and min(BO[n]['pplus']) >= 1 for n in CLASS2)


# class (2), covers not recognised by Regina: by Lemma 2.7 an algebraic diagram (obtained by Reidemeister moves; SnapPy
# finds the exteriors isometric to the KnotInfo ones) gives the normal forms of Sigma_2 and -Sigma_2. Both have four nodes,
# cokernel cyclic of order det K and positive eigenvalues, and a second algebraic diagram, different from the first, gives
# the same two normal forms. The KnotInfo diagrams are not algebraic.
CLASS2ALG = [x['name'] for f in ('reach_u12_orb.json', 'reach_u23_orb.json') for x in J(TOOL, f) if status(x['name']) == 'nodef_alg']
assert sorted(CLASS2ALG) == sorted(UA) == sorted(AD) == sorted(XA) and 'nodef_alg' not in C47.values()
for n in CLASS2ALG:
    D = int(BY[n]['determinant'])
    assert UA[n]['D'] == D and UA[n]['knotinfo_diagram_vertices'] > 1 and UA[n]['diagram_crossings'] == AD[n]['crossings']
    for s in ('1', '-1'):
        t = UA[n][s]; deg = Counter(v for e in t['edges'] for v in e)
        assert t['det_raw'] == D and t['smith'] in ([D], [-D]) and abs(t['inertia']['det']) == D and not t['inertia']['nullity']
        assert t['inertia']['pplus'] >= 1 and len(t['nodes']) == 4
        assert all(w <= -2 for v, w in enumerate(t['weights']) if deg[v] <= 2)                  # normal form
    assert UA[n]['1']['canon'] != UA[n]['-1']['canon']
    x = XA[n]
    assert all(x['second_same_normal_forms'].values()) and all(x['isometric_to_knotinfo']) and all(x['numpy_inertia_agrees'].values())
    assert x['second_differs_from_first'] in ('crossings', 'face sizes', 'diagram')
VS = VA['summary']                     # the construction on the KnotInfo diagrams that are algebraic
assert not VS['det_mismatch'] and not VS['identified_disagree'] and not VS['not_in_normal_form']
ALG_THM = [e['knot'] for e in rows2 + rows3 if e['knot'] in VA['records'] and VA['records'][e['knot']]['identified']['agrees']]


# class (3): the negative-definite normal form Gamma has at most two bad vertices (Theorem 2.6 applies)
def negdef_form(n):
    for src in NF_SOURCES:
        if n in src:
            return src[n]['plumbing']
    r = ps0[n] if n in ps0 else C47n[n]
    return r['normal_forms'][[o for o, v in r['pplus'].items() if v == 0][0]]


CLASS3 = [x['name'] for f in ('reach_u12_orb.json', 'reach_u23_orb.json') for x in J(TOOL, f)
          if status(x['name']) in ('obstructed', 'tested', 'notest')] + THMC_OBSTRUCTED
BAD3 = Counter(bad_vertices(negdef_form(n)) for n in CLASS3)
assert len(CLASS3) == both('obstructed') + both('tested') + both('notest') and set(BAD3) <= {1, 2}
assert all(bad_vertices(negdef_form(n)) == 1 for n in THMC_OBSTRUCTED)        # 13n3004, 13n3005


# every tree above is identified with Sigma_2 by triangulating its plumbed manifold, and its inertia and determinant are
# exact (neumann/identify_trees.py)
IT = J(TOOL, 'neumann', 'identify_trees.json')


def canon(w, edges):
    adj = {i: [] for i in range(len(w))}
    for a, b in edges:
        adj[a].append(b); adj[b].append(a)
    enc = lambda v, p: '(' + str(w[v]) + ''.join(sorted(enc(u, v) for u in adj[v] if u != p)) + ')'
    return min(enc(r, -1) for r in range(len(w)))


def identified(n, tree, cls):
    t = {x['canon']: x for x in IT[n]['trees']}.get(canon(tree['weights'], tree['edges'])) if n in IT else None
    ok = bool(t and IT[n]['cls'] == cls and t['matched'] and t['nullity'] == 0 and abs(t['det']) == int(BY[n]['determinant']))
    return t if ok else None


assert all(identified(n, t, 2) and identified(n, t, 2)['pplus'] == t['pplus'] >= 1 for n in CLASS2 for t in BO[n]['trees'])
assert all(identified(n, negdef_form(n), 3) and identified(n, negdef_form(n), 3)['pplus'] == 0 for n in CLASS3)


nums = dict(nNewUtwo=len(rows2), nNewUthree=len(rows3), nUfive=len(rows5), nAdopted=len(rows2) + len(rows3),
            nObstrD=sum(nf[e['knot']]['verdict'] == 'OBSTRUCTED' for e in rows2),
            nObstrHF=sum(nf[e['knot']]['verdict'] == 'OBSTRUCTED_HF' for e in rows2),
            nBadOne=sum(bad_vertices(nf[e['knot']]['plumbing']) == 1 for e in rows2 + rows3),
            nBadTwo=sum(bad_vertices(nf[e['knot']]['plumbing']) == 2 for e in rows2 + rows3),
            nNegDefBadOne=BAD3[1], nNegDefBadTwo=BAD3[2],
            nIdentified=sum(bool(nf[e['knot']].get('plumbing_name_matches')) for e in rows2 + rows3),
            nUthreeAll=len(rows66), nUthreeAlt=sum('a_' in n for n in rows66), nUthreeNonalt=sum('n_' in n for n in rows66),
            nOpenUonetwo=sum(cols['u12'].values()), nOpenUtwothree=sum(cols['u23'].values()),
            nOpenAll=sum(cols['u12'].values()) + sum(cols['u23'].values()),
            nHypCert=both('hyperbolic'), nMixed=both('mixed'), nNoDef=both('nodef') + both('nodef_alg'),
            nNoDefRec=both('nodef'),
            nNegDef=both('obstructed') + both('tested') + both('notest'), nGraphFour=both('nodef_alg'),
            nNegDefOne=cols['u12']['obstructed'] + cols['u12']['tested'] + cols['u12']['notest'],
            nNegDefTwo=cols['u23']['obstructed'] + cols['u23']['tested'] + cols['u23']['notest'], nNoTest=both('notest'),
            nUthreeNoCond=len(c66) - len(signed_only), nUthreeSigned=len(signed_only),
            nUthreeSigmaZero=sum(r['lower']['sigma'] == 0 for r in c66.values()),
            nBFTested=len(ext), nBFObstructed=len(ext_obs), nBFKnownThree=sum(v2T[n][0] >= 3 for n in ext_obs),
            nBFNegCtrl=bfc['negative']['tested'], nBFPubCtrl=bfc['positive']['tested'] // 2,
            nRefAlt=comp['references'], nDetCand=comp['exact_det_candidates'],
            nExclCross=dec['excluded: crossing number < 13'], nExclAlex=dec['excluded: Alexander polynomial differs'],
            nExclCover=len(UNRES),
            nMinDiag=dec['ACCEPTED: SnapPy identifies target; HFK agrees'], nNbrChanges=len(nbr), nNbrTypes=len(groups),
            nTwelveRef=comp12['references'], nTwelveDiag=bj12['n_diagrams'], nTwelveChanges=len(bj12['changes']),
            nInconclusive=both('tested'), nPosVol=both('hyperbolic') + both('mixed'),
            nNoPlumb=both('hyperbolic') + both('mixed') + both('nodef') + both('nodef_alg'),
            nGraphMfd=both('nodef') + both('nodef_alg') + both('obstructed') + both('tested') + both('notest'),
            nAlgCheck=VS['algebraic'], nAlgTrees=VS['identified_class2'] + VS['identified_class3'], nAlgThm=len(ALG_THM),
            nAlgCrossMin=min(r['crossings'] for r in AD.values()), nAlgCrossMax=max(r['crossings'] for r in AD.values()),
            nThmCNoPlumb=sum(s in ('hyperbolic', 'mixed', 'nodef') for s in C47.values()),
            nNewAll=len(rows2) + len(rows3) + len(rows66) + len(rows5), nClassicalNew=len(rows66) + len(rows5),
            nCtrl=len(ctl), nCtrlFourTwo=len(c42), nCtrlFourThree=len(c43),
            nCtrlFourThreeBoth=sum(r['verdict'] == 'OBSTRUCTED' for r in c43),
            nCtrlFourThreeExp=sum(r['verdict'] == 'OBSTRUCTED (expected orientation only)' for r in c43),
            nCtrlAlt=len(calt), nCtrlAltBoth=sum(r['verdict'] == 'OBSTRUCTED' for r in calt),
            nCtrlAltGoeritz=sum(r['v2_goeritz_verdict'].startswith('OBSTRUCTED') for r in calt))
open(os.path.join(OUT, 'numbers.tex'), 'w').write(''.join(f'\\newcommand{{\\{k}}}{{{v:,}}}\n'.replace(',', '{,}') for k, v in nums.items()))
print(json.dumps(nums, indent=1))
print('map', {k: dict(v) for k, v in cols.items()})
