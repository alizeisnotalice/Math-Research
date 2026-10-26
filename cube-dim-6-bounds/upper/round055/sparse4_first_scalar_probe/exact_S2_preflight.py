"""Actual new152-direction preflight; one point, no Fourier value."""
from pathlib import Path
from fractions import Fraction as Q
import itertools,json,hashlib,time
P=Path(__file__).resolve().parent
def directions():
 D=[(1,)+s for s in itertools.product([-1,1],repeat=5)]
 for support in itertools.combinations(range(6),4):
  for s in itertools.product([-1,1],repeat=3):
   v=[0]*6
   for i,a in zip(support,(1,)+s):v[i]=a
   D.append(tuple(v))
 assert len(D)==len(set(D))==152
 return D
def main():
 if (P/'S2_certificate.json').exists():raise SystemExit('Refuse second S2 preflight')
 t=time.monotonic();r=json.loads((P/'registration.json').read_text());x=r['point_numerators'];den=r['denominator'];N=sum(a*a for a in x)
 D=directions();dot=lambda a,b:sum(u*v for u,v in zip(a,b));norm=[dot(v,v) for v in D];projection=[dot(x,v) for v in D]
 lm=None;pm=None;lw=[];pw=[];zero=[];digest=hashlib.sha256()
 for i,v in enumerate(D):
  q=Q(N*norm[i]-projection[i]**2,den*den*norm[i]);assert q>=0
  if lm is None or q<lm:lm=q;lw=[i]
  elif q==lm:lw.append(i)
 count=0
 for i,j in itertools.combinations(range(152),2):
  a,c=norm[i],norm[j];b=dot(D[i],D[j]);gram=a*c-b*b;assert gram>0
  n=N*gram-(c*projection[i]**2-2*b*projection[i]*projection[j]+a*projection[j]**2)
  assert n>=0;q=Q(n,den*den*gram);digest.update(f'{i},{j}:{q}\n'.encode());count+=1
  if pm is None or q<pm:pm=q;pw=[(i,j)]
  elif q==pm:pw.append((i,j))
  if not q:zero.append((i,j))
 assert count==11476
 out=dict(status='PASS_STRICT_OFF_TRUE_S2' if pm>0 else 'REJECT_ON_TRUE_S2_STOP',point_numerators=x,denominator=den,direction_count=152,pair_count=count,line_min_squared=str(lm),pair_min_squared=str(pm),line_min_witnesses=[dict(index=i,direction=D[i]) for i in lw],pair_min_witnesses=[dict(indices=[i,j],directions=[D[i],D[j]]) for i,j in pw],all_pair_distance_digest=digest.hexdigest(),zero_pairs=zero,seconds=time.monotonic()-t,actual_new_preflight=True,C053_preflight_inherited=False,physical_kernel_calls=0)
 (P/'S2_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
 if pm<=0:raise SystemExit(2)
if __name__=='__main__':main()
