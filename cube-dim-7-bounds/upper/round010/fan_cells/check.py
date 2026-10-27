#!/usr/bin/env python3
"""Exact structural audit from retained R009 inputs; no prior 6-subset enumeration or Green evaluation."""
from flint import fmpz_mat
from fractions import Fraction
from pathlib import Path
import hashlib,json,math
p=Path(__file__).resolve().parent;sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest();mf=json.loads((p/'manifest.json').read_text())
for f,h in mf['files'].items():assert sha(p/f)==h,f
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,f=line.split('  ',1);q=Path(f);assert sha(q if q.is_absolute() else p/q)==h,f
old=json.loads((p/'r009-manifest.json').read_text());assert sha(p/'r009-rays.json')==old['files']['rays.json'] and sha(p/'r009-inputs.json')==old['files']['inputs.json']
R=json.loads((p/'r009-rays.json').read_text());I=json.loads((p/'r009-inputs.json').read_text());N=I['arrangement_normals'];B=I['prefix_normals'];CUT=I['cut_body_indices'];W=[r['w_primitive'] for r in R];assert len(R)==46 and len(N)==36 and len(B)==64 and len(CUT)==29 and N[7:]==[B[i] for i in CUT]
for i,r in enumerate(R):
 w=W[i];assert r['id']==i and math.gcd(*w)==1 and min(w)>=0 and max(w)>0
 vals=[sum(a*b for a,b in zip(n,w)) for n in N];assert r['normal_incidence']==[int(v==0) for v in vals] and r['active_normals']==[j for j,v in enumerate(vals) if v==0]
 signs=[(v>0)-(v<0) for n in B for v in [sum(a*b for a,b in zip(n,w))]];assert r['body_dot_signs']==signs
 rows=r['witness_rows'];cols=r['witness_minor']['columns'];assert len(rows)==len(set(rows))==6 and len(cols)==len(set(cols))==6 and all(vals[j]==0 for j in rows);det=int(fmpz_mat([[N[j][k] for k in cols] for j in rows]).det());assert det and det==r['witness_minor']['determinant']
 x=[sum(w[j] for j in range(i,7)) for i in range(7)];d=sum(x);assert r['xi_primitive']==x and r['normalization_denominator']==d and r['xi_normalized']==[str(Fraction(v,d)) for v in x]
ALL=(1<<46)-1
def ids(m):return [i for i in range(46) if m>>i&1]
cache={}
def rank(m):
 if m not in cache:cache[m]=int(fmpz_mat([W[i] for i in ids(m)]).rank()) if m else 0
 return cache[m]
dots=[[sum(a*b for a,b in zip(n,w)) for w in W] for n in N];zero=[sum(1<<i for i,x in enumerate(row) if x==0) for row in dots];pos=[sum(1<<i for i,x in enumerate(row) if x>=0) for row in dots];neg=[sum(1<<i for i,x in enumerate(row) if x<=0) for row in dots]
levels=[json.loads(x) for x in (p/'levels.jsonl').read_text().splitlines()];assert len(levels)==29 and rank(ALL)==7;cur={ALL:[]}
for level,rec in enumerate(levels):
 assert rec['level']==level+1 and rec['pilot']==(level<4) and rec['parent_count']==len(cur) and len(rec['parents'])==len(cur);nxt={};seen=set()
 for par in rec['parents']:
  mask=int(par['mask']);assert mask in cur and mask not in seen and par['signs']==cur[mask];seen.add(mask);assert len(par['children'])==2
  for c,sg,filtermask in zip(par['children'],[1,-1],[pos[level+7],neg[level+7]]):
   child=mask&filtermask;rr=rank(child);assert c['sign']==sg and int(c['mask'])==child and c['rank']==rr and c['kept']==(rr==7)
   if rr==7:
    path=cur[mask]+[sg];assert child not in nxt or nxt[child]==path;nxt[child]=path
  assert any(c['kept'] for c in par['children'])
 assert len(nxt)==rec['child_count'];cur=nxt
cells=[json.loads(x) for x in (p/'cells.jsonl').read_text().splitlines()];terms=[json.loads(x) for x in (p/'terminal.jsonl').read_text().splitlines()];receipt=json.loads((p/'receipt.json').read_text());assert receipt['status']=='completed' and len(cells)==len(terms)==len(cur)==135 and {x['id'] for x in cells}=={x['id'] for x in terms}==set(range(135)) and [x['id'] for x in terms[:4]]==[0,1,2,3]
assert {int(c['mask']):c['signs'] for c in cells}==cur
D={int(k):v for k,v in json.loads((p/'face-dag.json').read_text()).items()};assert len(D)==receipt['unique_faces'];pullcache={}
for mask,node in D.items():
 r=rank(mask);v=ids(mask);assert r==node['rank'] and v[0]==node['anchor']
 facets=sorted({mask&z for z in zero if mask&z and mask&z!=mask and rank(mask&z)==r-1}) if r>1 else []
 assert facets==[int(x) for x in node['facets']]
 if r==1:assert len(v)==1
 elif r>1:assert facets and any(not (f>>v[0]&1) for f in facets)
def pull(mask):
 if mask in pullcache:return pullcache[mask]
 node=D[mask];a=node['anchor']
 if node['rank']==1:ret=[(a,)]
 else:ret=sorted(tuple(sorted((a,)+t)) for fstr in node['facets'] for f in [int(fstr)] if not (f>>a&1) for t in pull(f))
 assert len(ret)==len(set(ret));pullcache[mask]=ret;return ret
nn=[sum((i+1)*x for i,x in enumerate(w)) for w in W];qq=[sum(abs(sum(a*b for a,b in zip(n,w))) for n in B) for w in W];seen=set();total=Fraction()
for cell in cells:
 mask=int(cell['mask']);sg=cell['signs'];vertices=ids(mask);assert cell['status']=='completed' and cell['vertices']==vertices and cell['facets']==D[mask]['facets']
 witness=[sum(W[i][j] for i in vertices) for j in range(7)];assert witness==cell['witness'] and min(witness)>0 and all(sg[j]*sum(N[j+7][k]*witness[k] for k in range(7))>0 for j in range(29))
 noncut=[B[i] for i in range(64) if i not in CUT];ell=[sum(n[k] for n in noncut)+sum(sg[j]*N[j+7][k] for j in range(29)) for k in range(7)];assert ell==cell['ell'] and all(sum(a*b for a,b in zip(ell,W[i]))==qq[i]>0 for i in vertices)
 actual=sorted(tuple(s['vertices']) for s in cell['simplices']);assert actual==pull(mask);vol=Fraction()
 for x in cell['simplices']:
  ix=x['vertices'];key=tuple(ix);assert len(ix)==7 and key not in seen and all(i in vertices for i in ix);seen.add(key);assert all(sg[j]*dots[j+7][i]>=0 for j in range(29) for i in ix)
  det=int(fmpz_mat([[W[j][i] for j in ix] for i in range(7)]).det());assert det and det==x['detW'];n=[nn[i] for i in ix];q=[qq[i] for i in ix];assert x['N']==n and x['q']==q and Fraction(x['alpha'])==Fraction(abs(det),math.prod(q))
  v=Fraction(math.factorial(7)*abs(det),math.prod(n));assert Fraction(x['volume'])==v;vol+=v
 assert vol==Fraction(cell['volume']);total+=vol
assert total==1==Fraction(receipt['volume']) and len(seen)==receipt['top_simplices']==212
print('PASS: copied R009 ray witnesses; every one of 29 parent/child rank splits; every face zero-facet set and pulling recursion; all 135 cells, 212 simplex determinants/psi/volumes; exact global volume 1. Coverage uses the accompanying proof, not volume alone.')
