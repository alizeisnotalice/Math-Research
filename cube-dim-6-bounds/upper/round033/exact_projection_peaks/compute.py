#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,permutations,combinations
from collections import Counter
from math import factorial,prod,gcd
import json,hashlib,gzip,time,datetime,resource
from flint import fmpz_mat,arb,fmpq,ctx
ctx.prec=512
O=Path(__file__).resolve().parent;U=O.parents[1]
FAN=U/'round003/direct_kernel/fan_integer_input.txt'; GEOM=U/'round002/fan/exact_cells_results.json'
START=time.monotonic();deadline=datetime.datetime.fromisoformat(json.loads((O/'input_manifest.json').read_text())['started_utc']).timestamp()+1200
assert not (O/'results.json').exists(),'one full enumeration only'

def det(a):
 n=len(a)
 if n==1:return a[0][0]
 if n==2:return a[0][0]*a[1][1]-a[0][1]*a[1][0]
 if n==3:
  x,y,z=a;return x[0]*(y[1]*z[2]-y[2]*z[1])-x[1]*(y[0]*z[2]-y[2]*z[0])+x[2]*(y[0]*z[1]-y[1]*z[0])
 return int(fmpz_mat(a).det())
def canonical(v):
 g=gcd(*v);v=tuple(x//g for x in v)
 if next(x for x in v if x)!=abs(next(x for x in v if x)):v=tuple(-x for x in v)
 return v
lines=FAN.read_text().split();it=iter(map(int,lines));n=next(it);cones=[]
for i in range(n):
 d=next(it); rays=[];rates=[]
 for _ in range(6):
  rates.append(next(it));rays.append(tuple(next(it) for _ in range(6)))
 cones.append({'det':d,'rays':rays,'q':rates})
assert list(it)==[] and n==21
body=[(1,)+s for s in product((-1,1),repeat=5)]
geom=json.loads(GEOM.read_text())
for c,g in zip(cones,geom['cells'],strict=True):
 assert abs(det([list(x) for x in zip(*c['rays'])]))==c['det']
 for ray,q,old,oldq in zip(c['rays'],c['q'],g['rays_a0_equals_one'],g['rates'],strict=True):
  assert q==sum(abs(sum(a*b for a,b in zip(v,ray))) for v in body)
  assert tuple(Q(x,ray[0]) for x in ray)==tuple(map(Q,old)) and Q(q,ray[0])==Q(oldq)
  assert all(ray[i]>=ray[i+1] for i in range(5)) and ray[-1]>=0 and sum(ray)<=12
signrows=list(product((-1,1),repeat=6));rid={row:i for i,row in enumerate(signrows)}
# Preflight every possible body-form row on each cone for the orthant property.
projected=[]
for row in signrows:
 cr=[]
 for c in cones:
  vals=tuple(sum(a*b for a,b in zip(row,r)) for r in c['rays'])
  assert not(any(x>0 for x in vals) and any(x<0 for x in vals))
  assert max(map(abs,vals))<=12
  cr.append(vals)
 projected.append(cr)
# S4 action on four projective row classes, induced by signed row permutations.
classes=[(1,1,1),(1,1,-1),(1,-1,1),(1,-1,-1)];actions=set()
for p in permutations(range(3)):
 for sg in product((-1,1),repeat=3):
  im=[]
  for c in classes:
   w=tuple(sg[j]*c[p[j]] for j in range(3));w=tuple(x*w[0] for x in w);im.append(classes.index(w))
  actions.add(tuple(im))
assert len(actions)==24
parts=[(4,1,1,0),(3,2,1,0),(2,2,2,0),(3,1,1,1),(2,2,1,1)]
hist=Counter();ghist=Counter()
for triple in combinations(body,3):
 counts=Counter(tuple(triple[j][k]*triple[0][k] for j in range(3)) for k in range(6))
 typ=tuple(sorted((counts[c] for c in classes),reverse=True));assert typ in parts
 gram=[[sum(x*y for x,y in zip(a,b)) for b in triple] for a in triple]
 assert det(gram)>0
 hist[typ]+=1;ghist[det(gram)]+=1
assert sum(hist.values())==4960
cases=[{'id':'r1','columns':[(1,)]*6,'triple_multiplicity':None}]
for k in (1,2,3):cases.append({'id':f'r2_k{k}','columns':[(1,1)]*(6-k)+[(1,-1)]*k,'triple_multiplicity':None})
for p in parts:cases.append({'id':'r3_'+''.join(map(str,p)),'columns':[c for c,k in zip(classes,p) for _ in range(k)],'triple_multiplicity':hist[p]})
pre={'status':'PASS_PRE_ENUMERATION','fan_sha256':hashlib.sha256(FAN.read_bytes()).hexdigest(),'triple_orbit_counts':{str(k):v for k,v in hist.items()},'triple_Gram_counts':dict(ghist),'projective_class_actions':len(actions),'all_64_body_forms_all_21_cones_row_orthant':True,'rays_l1_max_at_most':12,'normal_minor_bound':288}
(O/'preflight.json').write_text(json.dumps(pre,indent=2)+'\n')
allterms=[];rows=[];totalpatterns=0
for case in cases:
 cols=case['columns'];r=len(cols[0]);d=6-r;v=[1009**(r-1-j) for j in range(r)]
 L=[list(x) for x in zip(*cols)];Gram=det([[sum(a*b for a,b in zip(x,y)) for y in L] for x in L])
 perms=sorted(set(permutations(cols)));mult=prod(factorial(k) for k in Counter(cols).values());assert len(perms)*mult==720
 mh=Counter();selected=0;R=Q(0);facets=set();case_terms=[]
 for pi,pc in enumerate(perms):
  for mask in range(64):
   ids=tuple(rid[tuple(c[j]*(1 if (mask>>k)&1 else -1) for k,c in enumerate(pc))] for j in range(r))
   for ci,c in enumerate(cones):
    totalpatterns+=1
    assert totalpatterns<=1000000
    if totalpatterns%8192==0:
     assert time.time()<deadline,'20-minute deadline'
    rr=[projected[i][ci] for i in ids]
    J=[j for j in range(6) if any(rr[k][j] for k in range(r))];m=len(J);assert m>=r;mh[m]+=mult
    C=[tuple(rr[k][j] for k in range(r)) for j in J]
    if r==2:
     for w in C:
      nn=canonical((w[1],-w[0]));assert sum(a*b for a,b in zip(nn,v))!=0;facets.add(nn)
    elif r==3:
     for x,y in combinations(C,2):
      cr=(x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0])
      assert max(map(abs,cr))<=288
      if cr!=(0,0,0):
       nn=canonical(cr);assert sum(a*b for a,b in zip(nn,v))!=0;facets.add(nn)
    if m!=r:continue
    M=[[rr[k][j] for j in J] for k in range(r)];dm=det(M);assert dm
    nums=[]
    for j in range(r):
     temp=[row[:] for row in M]
     for k in range(r):temp[k][j]=v[k]
     nums.append(det(temp))
    assert all(z!=0 for z in nums)
    if not all(z*dm>0 for z in nums):continue
    denom=abs(dm)*prod(c['q'][j] for j in range(6) if j not in J)
    weight=Q(c['det']*mult,denom);R+=weight;selected+=1
    case_terms.append({'permutation_index':pi,'sign_mask':mask,'cone_index':ci,'nonzero_normal_columns':J,'detM':dm,'cramer_numerators':nums,'weight':str(weight)})
 assert sum(mh.values())==967680 and R>0
 row={'id':case['id'],'normal_rows':L,'r':r,'dimension':d,'Gram_determinant':Gram,'trace_direction':v,'unique_coordinate_permutations':len(perms),'permutation_multiplicity':mult,'unique_pattern_cone_items':len(perms)*64*21,'logical_items':sum(mh.values()),'m_histogram_logical':dict(mh),'selected_terms':selected,'R_exact':str(R),'kappa_exact':'R/(2*pi)^'+str(d),'peak_exact':'sqrt('+str(Gram)+')*R/(2*pi)^'+str(d),'generic_facets':sorted(facets),'triple_orbit_count':case['triple_multiplicity']}
 rows.append(row);allterms.append({'case':case['id'],'terms':case_terms})
 print(case['id'], 'R=',R,'selected=',selected,flush=True)
 payload={'status':'PARTIAL_COMPLETE_CASES','cases':rows,'pattern_items':totalpatterns}
 (O/'progress.json').write_text(json.dumps(payload,indent=2)+'\n')
assert totalpatterns==701568
# Nine and only nine comparisons against the frozen old quotient peaks.
old=json.loads(Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube6-upper-astra/cube-dim-6-bounds/upper/round008/off_strata_integration/results.json').read_text())
for row in rows:
 Rq=Q(row['R_exact']);d=row['dimension'];kap=arb(fmpq(Rq.numerator,Rq.denominator))/(2*arb.pi())**d
 peak=arb(row['Gram_determinant']).sqrt()*kap;oldpeak=arb(old['quotient_Q'+str(d)]['ball'])
 row.update(kappa_arb512=kap.str(110,more=True),peak_arb512=peak.str(110,more=True),old_peak_over_new= (oldpeak/peak).str(110,more=True),strictly_below_old_peak=bool(peak<oldpeak))
raw=json.dumps({'schema':'All selected terms; each weight already includes repeated coordinate-permutation multiplicity. All64 signs explicitly represented; no global sign factor.','cases':allterms},separators=(',',':')).encode()
with gzip.GzipFile(filename=str(O/'selected_terms.json.gz'),mode='wb',mtime=0) as z:z.write(raw)
result={'status':'COMPLETE_NINE_EXACT_PROJECTION_PEAKS','cases':rows,'pattern_items':totalpatterns,'logical_items':9*967680,'scalar_comparisons':9,'kernel_kind':'new Fourier projection-peak trace, not spatial G/H','native_G_H_calls':0,'spatial_grid_points':0,'matlab_calls':0,'elapsed_seconds':time.monotonic()-START,'max_rss_bytes':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'terms_sha256':hashlib.sha256((O/'selected_terms.json.gz').read_bytes()).hexdigest(),'uncompressed_terms_sha256':hashlib.sha256(raw).hexdigest(),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'results.json').write_text(json.dumps(result,indent=2)+'\n');(O/'progress.json').unlink()
print('COMPLETE',totalpatterns,'patterns',result['elapsed_seconds'],'seconds',flush=True)
