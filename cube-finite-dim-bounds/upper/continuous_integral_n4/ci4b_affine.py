"""CI4B affine polynomial majorants; LP proposes, exact coefficients certify."""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
from math import prod
import hashlib,json,argparse,time
from flint import fmpq_mpoly_ctx,fmpq,arb,ctx as arbctx
CTX=fmpq_mpoly_ctx.get(('a','b','c','d'),ordering='lex');XS=CTX.gens()
SOURCE=Path('/Users/zhengzhihao/Desktop/T/Math-Research-finite-upper/cube-finite-dim-bounds/upper/u11_minimum_global_results.json')
SOURCE_SHA='43be4de041c9da823d192ff6d8ef2440c3668f88b13fc52dbff320dd0011d96f'
INTERFACE=SOURCE.with_name('u11_minimum_global_interface.py')
INTERFACE_SHA='6c8d1a5b4912103df8105f7c7662ac2265b0978f91de00651be03acfb1790bae'
RAYS=[[(1,0,0,0),(1,1,0,0),(2,1,1,0),(3,1,1,1)],[(2,1,1,0),(1,1,0,0),(1,1,1,0),(1,1,1,1)],[(2,1,1,0),(1,1,0,0),(1,1,1,1),(3,1,1,1)]]
def polynomial(dic):return CTX.from_dict({tuple(map(int,e.split(','))):fmpq(c) for e,c in dic.items()})
def aq(x):x=Q(x);return arb(x.numerator)/x.denominator
def det(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
def volume(vs):return abs(det(*[tuple(vs[i][j]-vs[0][j] for j in range(1,4)) for i in (1,2,3)]))/6

def roots():return [tuple(tuple(Q(v)/r[0] for v in r) for r in rs) for rs in RAYS]
def bisect(tet):
 i,j=max(combinations(range(4),2),key=lambda ij:(sum((tet[ij[0]][k]-tet[ij[1]][k])**2 for k in range(1,4)),tuple(-x for x in ij)))
 m=tuple((a+b)/2 for a,b in zip(tet[i],tet[j]));a=list(tet);b=list(tet);a[i]=m;b[j]=m
 return tuple(a),tuple(b)
def refine(tet,depth):
 ts=[tet]
 for _ in range(depth):ts=[child for t in ts for child in bisect(t)]
 return ts

def shells(tet):
 corner=next(v for v in tet if v==(Q(1),)*4);abc=[v for v in tet if v!=corner]
 for k in range(28):
  outer=[tuple(c+(v-c)/2**k for c,v in zip(corner,r)) for r in abc]
  inner=[tuple(c+(v-c)/2 for c,v in zip(corner,r)) for r in outer]
  a,b,c=outer;aa,bb,cc=inner
  yield k,[(a,b,c,aa),(b,c,aa,bb),(c,aa,bb,cc)]

def models():
 import importlib.util
 assert hashlib.sha256(SOURCE.read_bytes()).hexdigest()==SOURCE_SHA
 assert hashlib.sha256(INTERFACE.read_bytes()).hexdigest()==INTERFACE_SHA
 spec=importlib.util.spec_from_file_location('ci4b_shared_global',INTERFACE);mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 assert [tuple(map(tuple,r)) for r in RAYS]==[mod.CONE_RAYS[c] for c in ('I','II','III')]
 return [mod.homogeneous_pair(c) for c in ('I','II','III')]

def majorant(model,tet):
 import numpy as np
 from scipy.optimize import linprog
 N,D=model;sub=[sum((CTX.constant(fmpq(v[i].numerator,v[i].denominator))*x for v,x in zip(tet,XS)),CTX.constant(0)) for i in range(4)]
 n=N.compose(*sub).to_dict();d=D.compose(*sub).to_dict()
 assert all(c>=0 for c in d.values()),'D coefficient not nonnegative'
 assert all(c>=0 for c in n.values()),'N coefficient not nonnegative'
 def lift(dic):
  out={}
  for e,c in dic.items():
   for i in range(4):
    ee=list(e);ee[i]+=1;ee=tuple(ee);out[ee]=out.get(ee,fmpq(0))+c
  return out
 nl=lift(n);dl=lift(d);keys=set(nl)|set(dl);rows=[];bad=[]
 for e in keys:
  coeff=[]
  for i in range(4):
   ee=list(e);ee[i]-=1;coeff.append(d.get(tuple(ee),fmpq(0)) if ee[i]>=0 else fmpq(0))
  s=sum(coeff);rhs=nl.get(e,fmpq(0))
  if not s:
   if rhs>0:bad.append((e,str(rhs)))
  else:rows.append((coeff,s,rhs))
 if bad:return {'status':'support_obstruction','count':len(bad),'first':bad[0]}
 A=np.array([[float(v/s) for v in co] for co,s,rhs in rows]);b=np.array([float(rhs/s) for co,s,rhs in rows])
 res=linprog([1,1,1,1],A_ub=-A,b_ub=-b,bounds=[(0,None)]*4,method='highs')
 if not res.success:return {'status':'LP candidate unavailable','message':res.message}
 pp=[Q(float(v)).limit_denominator(10**12) for v in res.x]
 # Uniform exact outward shift pays every LP and rationalisation error.
 delta=max(Q(0),*(Q(str(rhs-sum(v*fmpq(p.numerator,p.denominator) for v,p in zip(co,pp))))/Q(str(s)) for co,s,rhs in rows))
 pp=[p+delta for p in pp]
 residuals=[sum(v*fmpq(p.numerator,p.denominator) for v,p in zip(co,pp))-rhs for co,s,rhs in rows]
 assert all(r>=0 for r in residuals) and all(p>=0 for p in pp)
 terms=sum((aq(p).root(3)**4 if p else arb(0)) for p in pp)
 upper=aq(volume(tet))*terms/(4*arb.pi()**4)
 return {'status':'exact_coefficient_majorant','P':list(map(str,pp)),'uniform_exact_shift':str(delta),'constraints':len(rows),'min_residual':str(min(residuals)),'upper_expression_display':str(upper),'upper_expression_rational':str(upper.upper().fmpq()),'volume':str(volume(tet))}

def pilot():
 rr=roots();cells=[('I_root',0,rr[0])]
 for case in (1,2):
  k,ts=next(shells(rr[case]));cells +=[(f'{case+1}_outer_{j}',case,t) for j,t in enumerate(ts)]
 for label,case,t in (cells[0],cells[1],cells[4]):cells.append((label+'_leaf8',case,refine(t,8)[0]))
 assert len(cells)==10
 return cells

def run(save):
 start=time.time();arbctx.prec=192;mm=models();out=[]
 for label,case,tet in pilot():
  ts=time.time();r=majorant(mm[case],tet);r.update(label=label,case=['I','II','III'][case],seconds=time.time()-ts,vertices=[[str(x) for x in v] for v in tet]);out.append(r)
  print(label,r['status'],r.get('upper_expression_display'),flush=True)
 raw=json.loads(SOURCE.read_text());canonical=json.dumps(raw['physical_cases'],sort_keys=True,separators=(',',':')).encode()
 result={'status':'fixed ten-cell pilot; not a whole-integral certificate','source_path':str(SOURCE),'source_sha256_at_run':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'physical_cases_canonical_sha256':hashlib.sha256(canonical).hexdigest(),'pilot_cells':out,'elapsed_seconds':time.time()-start}
 if save:Path(save).write_text(json.dumps(result,indent=2)+'\n')
 print('PILOT DONE seconds',result['elapsed_seconds'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--save');run(p.parse_args().save)
