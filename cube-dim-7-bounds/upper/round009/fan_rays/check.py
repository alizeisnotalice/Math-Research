#!/usr/bin/env python3
"""Stdlib integer checker for all output rays and coverage ledger; does not re-enumerate combinations."""
import hashlib,json,math,itertools
from fractions import Fraction
from pathlib import Path
p=Path(__file__).resolve().parent;sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest();mf=json.loads((p/'manifest.json').read_text())
for n,h in mf['files'].items():assert sha(p/n)==h,n
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,n=line.split('  ',1);path=Path(n);path=path if path.is_absolute() else p/path;assert sha(path)==h,n
D=json.loads((p/'inputs.json').read_text());V=D['body_vectors'];B=D['prefix_normals'];N=D['arrangement_normals'];assert len(V)==len(B)==64 and len(N)==36
assert V==[[1,*sgn] for sgn in itertools.product((-1,1),repeat=6)] and B==[[sum(v[:k]) for k in range(1,8)] for v in V]
assert N[:7]==[[int(i==j) for i in range(7)] for j in range(7)] and sorted(N[7:])==sorted(v for v in B if min(v)<0<max(v)) and sum(min(v)>=0 for v in B)==35
walks=[(1,-1,-1,*x) for x in itertools.product((-1,1),repeat=4)]+[(*t,*x) for t in [(1,1,-1,-1,-1),(1,-1,1,-1,-1)] for x in itertools.product((-1,1),repeat=2)]+[(1,1,1,-1,-1,-1,-1),(1,1,-1,1,-1,-1,-1),(1,1,-1,-1,1,-1,-1),(1,-1,1,1,-1,-1,-1),(1,-1,1,-1,1,-1,-1)]
assert [[sum(v[:k]) for k in range(1,8)] for v in walks]==N[7:]
T=[json.loads(x) for x in (p/'terminal.jsonl').read_text().splitlines()];R=json.loads((p/'receipt.json').read_text());assert R['status']=='completed' and len(T)==1903 and {t['id'] for t in T}==set(range(1903)) and all(t['status']=='completed' for t in T)
span=sorted((t['start'],t['end']) for t in T);assert span[0][0]==0 and span[-1][1]==math.comb(36,6) and all(a[1]==b[0] for a,b in zip(span,span[1:]));assert all(t['rank_below_6']+t['rank6_mixed']+t['rank6_feasible']==t['end']-t['start'] for t in T)
for c in range(7):assert sum(t['wall_strata'][c] for t in T)==math.comb(7,c)*math.comb(29,6-c)==R['wall_strata'][c]
pilot=[t for t in T if t['pilot']];assert len(pilot)==8 and sum(t['cofactor_checked'] for t in pilot)==8192 and sorted((t['start'],t['end']) for t in pilot)==[(i*1024,(i+1)*1024) for i in range(4)]+[(math.comb(36,6)-4096+i*1024,math.comb(36,6)-4096+(i+1)*1024) for i in range(4)]
for key in ['rank_below_6','rank6_mixed','rank6_feasible']:assert sum(t[key] for t in T)==R[key]
def det(a):
 a=[row[:] for row in a];n=len(a);prev=1;sign=1
 for k in range(n-1):
  q=next((i for i in range(k,n) if a[i][k]),None)
  if q is None:return 0
  if q!=k:a[k],a[q]=a[q],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    x=a[i][j]*pivot-a[i][k]*a[k][j];assert x%prev==0;a[i][j]=x//prev
   a[i][k]=0
  prev=pivot
 return sign*a[-1][-1]
rays=json.loads((p/'rays.json').read_text());assert len(rays)==R['unique_rays']==46;seen=set()
for i,z in enumerate(rays):
 assert z['id']==i;w=z['w_primitive'];assert len(w)==7 and min(w)>=0 and max(w)>0 and math.gcd(*w)==1 and tuple(w) not in seen;seen.add(tuple(w))
 xi=[sum(w[j] for j in range(k,7)) for k in range(7)];den=sum(xi);assert z['xi_primitive']==xi and den==sum((k+1)*w[k] for k in range(7))==z['normalization_denominator'] and z['xi_normalized']==[str(Fraction(x,den)) for x in xi]
 vals=[sum(a*b for a,b in zip(v,w)) for v in N];assert z['normal_incidence']==[int(x==0) for x in vals] and z['active_normals']==[j for j,x in enumerate(vals) if x==0]
 bs=[sum(a*b for a,b in zip(v,w)) for v in B];assert z['body_dot_signs']==[(x>0)-(x<0) for x in bs]
 ids=z['witness_rows'];c=z['witness_minor']['columns'];assert len(ids)==len(set(ids))==6 and all(vals[j]==0 for j in ids) and len(c)==len(set(c))==6
 minor=det([[N[j][k] for k in c] for j in ids]);assert minor and minor==z['witness_minor']['determinant']
assert all(tuple(int(i==j) for i in range(7)) in seen for j in range(7))
print('PASS: 46 primitive rays and exact witnesses, 36 incidences, 64 signs, rational normalizations, all 1903 disjoint lex spans, 8192 cofactor pilots and wall counts. Completeness relies on frozen generator/proof; combinations are not rerun.')
