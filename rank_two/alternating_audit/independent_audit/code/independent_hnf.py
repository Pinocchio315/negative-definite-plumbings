#!/usr/bin/env python3
"""Independent BF finite-group audit. No original project module imports.
Quotient representatives use integral column Hermite normal form of M.
Pairings use exact rational inverse, not B*x group keys or cyclic shortcut.
Candidate search uses literal lemma bounds a,b <= D and direct GL2(Z/4).
"""
import os, tempfile
from itertools import product
from fractions import Fraction
from collections import Counter
from pathlib import Path
import json, time, hashlib
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form

BASE=Path(os.environ.get('BF_AUDIT_SOURCE_ROOT', str(Path(__file__).resolve().parents[2])))
OUT=Path(os.environ.get('BF_AUDIT_OUTPUT_DIR') or tempfile.mkdtemp(prefix='bf_audit_replay_'))
OUT.mkdir(parents=True,exist_ok=True)

def inverse_fraction(matrix):
    """Gauss-Jordan on Fraction independent of sympy inversion."""
    n=len(matrix)
    a=[[Fraction(x) for x in row]+[Fraction(int(i==j)) for j in range(n)] for i,row in enumerate(matrix)]
    for col in range(n):
        k=next(i for i in range(col,n) if a[i][col])
        a[col],a[k]=a[k],a[col]
        piv=a[col][col]
        a[col]=[x/piv for x in a[col]]
        for i in range(n):
            if i != col and a[i][col]:
                t=a[i][col]
                a[i]=[x-t*y for x,y in zip(a[i],a[col])]
    return [row[n:] for row in a]

def signature_exact(matrix):
    """Sylvester inertia via rational LDL with 1x1 / 2x2 pivots."""
    a=[[Fraction(x) for x in row] for row in matrix]
    n=len(a); pos=neg=0
    while n:
        diag=next((i for i in range(n) if a[i][i]),None)
        if diag is not None:
            inds=[diag]+[i for i in range(n) if i!=diag]
            a=[[a[i][j] for j in inds] for i in inds]
            t=a[0][0]
            if t>0: pos+=1
            else: neg+=1
            a=[[a[i][j]-a[i][0]*a[0][j]/t for j in range(1,n)] for i in range(1,n)]
            n-=1
        else:
            pair=next(((i,j) for i in range(n) for j in range(i+1,n) if a[i][j]),None)
            if pair is None: break
            inds=list(pair)+[i for i in range(n) if i not in pair]
            a=[[a[i][j] for j in inds] for i in inds]
            t=a[0][1]
            a=[[a[i][j]-(a[i][0]*a[1][j]+a[i][1]*a[0][j])/t for j in range(2,n)] for i in range(2,n)]
            pos+=1;neg+=1;n-=2
    return pos-neg

GL4=list(product(range(4),repeat=4))
GL4=[p for p in GL4 if (p[0]*p[3]-p[1]*p[2])%2]
assert len(GL4)==96

def eligible_residue(a,c,b,p):
    # v=(r,t), w=(s,u). Direct quadratic evaluations in all residue classes.
    for r,s,t,u in GL4:
        if (a*r*r+2*c*r*t+b*t*t)%4!= (1 if p>0 else 3): continue
        if (a*s*s+2*c*s*u+b*u*u)%4!= (1 if p>1 else 3): continue
        if (a*r*s+c*(r*u+s*t)+b*t*u)%2: continue
        return [r,s,t,u]
    return None

def literal_candidates(D,sig,p):
    if sig not in (-2,0,2):return [],0
    d= D if sig else -D
    allforms=set()
    for a in range(-D,D+1):
        if a==0:continue
        for c in range(abs(a)//2+1):
            numer=d+c*c
            if numer%a:continue
            b=numer//a
            if abs(a)<=abs(b)<=D:
                if d>0 and (a>0)!=(sig>0):continue
                allforms.add((a,c,b))
    for c in range(D+1):
        if -c*c == d:
            for a in range(-c,c+1):allforms.add((a,c,0))
    return [(v,eligible_residue(*v,p)) for v in sorted(allforms) if eligible_residue(*v,p) is not None],len(allforms)

class HNFQuotient:
    def __init__(self,M):
        self.M=Matrix(M); self.D=abs(int(self.M.det()));self.r=len(M)
        self.H=hermite_normal_form(self.M)
        # H=MU for unimodular U, verified by rational divisibility and det.
        u=self.M.inv()*self.H
        assert all(x.q==1 for x in u) and abs(int(u.det()))==1
        self.reps=[tuple(v) for v in product(*[range(int(self.H[i,i])) for i in range(self.r)])]
        assert len(self.reps)==self.D
        self.inv=inverse_fraction(M)
        self.diagonal=[self.pair(x,x) for x in self.reps]
        self.byvalue={q:[i for i,v in enumerate(self.diagonal) if v==q] for q in set(self.diagonal)}
        S=smith_normal_form(self.M,domain=ZZ)
        self.snf=[abs(int(S[i,i])) for i in range(self.r) if abs(int(S[i,i]))>1]
    def pair(self,x,y):
        q=2*sum(Fraction(x[i])*self.inv[i][j]*y[j] for i in range(self.r) for j in range(self.r))
        return q%1
    def reduce(self,x):
        x=list(x)
        for j in reversed(range(self.r)):
            q=x[j]//int(self.H[j,j])
            for i in range(j+1):x[i]-=q*int(self.H[i,j])
        assert all(0<=x[i]<int(self.H[i,i]) for i in range(self.r))
        return tuple(x)
    def generated_order(self,x,y):
        seen={tuple(0 for _ in range(self.r))}; queue=list(seen)
        for v in queue:
            for w in (x,y):
                z=self.reduce([a+b for a,b in zip(v,w)])
                if z not in seen:seen.add(z);queue.append(z)
        return len(seen)
    def pattern(self,sigma,p):
        sig=sigma+2*p-2
        if sig not in (-2,0,2):return {'p':p,'n':2-p,'signature_C':sig,'verdict':'EXCLUDED_SIGNATURE'}
        forms,total=literal_candidates(self.D,sig,p)
        checks=[];tested=0
        for (a,c,b),basis4 in forms:
            d=a*b-c*c
            target=[Fraction(b,d)%1,Fraction(a,d)%1,Fraction(-c,d)%1]
            left=self.byvalue.get(target[0],[]);right=self.byvalue.get(target[1],[])
            check={'C':[[a,c],[c,b]],'target':[str(q) for q in target],'self_candidates':[len(left),len(right)]}
            found=None
            for i in left:
                for j in right:
                    tested+=1
                    if self.pair(self.reps[i],self.reps[j])==target[2]:
                        order=self.generated_order(self.reps[i],self.reps[j])
                        if order==self.D:
                            found={'v1':self.reps[i],'v2':self.reps[j],'generated_order':order,'basis_mod4':basis4}
                            break
                if found:break
            if found:return {'p':p,'n':2-p,'signature_C':sig,'verdict':'SURVIVES','literal_count':total,'eligible_count':len(forms),'checked_pairs':tested,'C':[[a,c],[c,b]],'witness':found}
            checks.append(check)
        return {'p':p,'n':2-p,'signature_C':sig,'verdict':'EXCLUDED_LINKING','literal_count':total,'eligible_count':len(forms),'checked_pairs':tested,'checks':checks}

records=json.loads((BASE/'inputs/cohort959.json').read_text());byname={r['name']:r for r in records}
names=[r['name'] for r in json.loads((BASE/'derived/lower_bound_three19.json').read_text())]
frozen={r['name']:r for r in map(json.loads,(BASE/'bf/cohort959_run/results.jsonl').read_text().splitlines())}
results=[]
# Add two cyclic controls plus one noncyclic original survivor, all from existing cohort.
controlnames=['10_6','10_11']
controlnames+=next(([r['name']] for r in frozen.values() if not r['name'].endswith('__mirror') and not r['cyclic_group'] and r['verdict']=='SURVIVES'),[])
for name in names+controlnames:
    r=byname[name];M=[[r['V'][i][j]+r['V'][j][i] for j in range(len(r['V']))] for i in range(len(r['V']))]
    for mirror in (False,True):
        mm=[[-x for x in row] for row in M] if mirror else M
        sigma=signature_exact(mm);assert sigma==(-r['sigma'] if mirror else r['sigma'])
        g=HNFQuotient(mm)
        pat=[g.pattern(sigma,p) for p in range(3)]
        verdict='SURVIVES' if any(v['verdict']=='SURVIVES' for v in pat) else 'ALL_PATTERNS_EXCLUDED'
        key=name+('__mirror' if mirror else '')
        result={'name':key,'D':g.D,'sigma':sigma,'source_M':mm,'HNF':[[int(x) for x in row] for row in g.H.tolist()],'SNF':g.snf,'HNF_representative_count':len(g.reps),'attainable_doubled_self_values':sorted(str(q) for q in set(g.diagonal)),'patterns':pat,'verdict':verdict,'frozen_verdict_agrees':verdict==frozen[key]['verdict'],'frozen_pattern_agrees':[v['verdict']==o['verdict'] for v,o in zip(pat,frozen[key]['patterns'])]}
        # Compare independent full literal candidates to stored negative checks for excluded knots.
        if name in names:
            result['negative_candidate_lists_agree']=[set(tuple(x for row in chk['C'] for x in row) for chk in v.get('checks',[]))==set(tuple(x for row in chk['C'] for x in row) for chk in o.get('negative_checks',[])) for v,o in zip(pat,frozen[key]['patterns'])]
        results.append(result)
        print(key,verdict,'SNF',g.snf,'self-values',len(set(g.diagonal)),'pairs',[v.get('checked_pairs',0) for v in pat],flush=True)
(OUT/'hnf_results.json').write_text(json.dumps(results,indent=2)+'\n')
summary={'independent_method':'HNF representatives + stdlib Fraction Gaussian inverse + literal Lemma5.2 bounds + 96 direct residue bases + full pair and subgroup search; no project imports','records':len(results),'claimed19_originals_and_mirrors':38,'controls':controlnames,'verdicts':dict(Counter(r['verdict'] for r in results)),'all_frozen_verdicts_agree':all(r['frozen_verdict_agrees'] for r in results),'all_frozen_patterns_agree':all(all(r['frozen_pattern_agrees']) for r in results),'all_negative_candidate_lists_agree':all(all(r.get('negative_candidate_lists_agree',[True])) for r in results),'counterexamples':[r['name'] for r in results if r['name'].split('__mirror')[0] in names and r['verdict']=='SURVIVES'],'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(OUT/'hnf_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
