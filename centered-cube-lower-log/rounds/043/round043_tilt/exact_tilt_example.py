"""Exact rational score-tilt mismatch for a physical positive periodic box source."""
from pathlib import Path
from fractions import Fraction as F
import json,hashlib
HERE=Path(__file__).resolve().parent

def integ(segments,x,ell):
 lo=x-ell/2;hi=x+ell/2;s=F(0)
 for left,right,d in segments:
  for j in [-1,0,1]:s+=d*max(F(0),min(hi,right+j)-max(lo,left+j))
 return s

def affine_product(h,a,b,c,d):return h*(2*a*c+a*d+b*c+2*b*d)/6

def main():
 target=HERE/'EXACT_TILT_EXAMPLE.json';assert not target.exists()
 q=1;ell=F(1,4);L=q+ell
 classes=[[(F(0),F(1,2),F(1,10)),(F(9,40),F(11,40),F(9))],[(F(1,2),F(1),F(1))]]
 eta=[sum((r-l)*h for l,r,h in s)for s in classes];assert eta==[F(1,2)]*2
 knots={F(0),F(1)}
 for segments in classes:
  for l,r,d in segments:
   for end in [l,r]:
    for side in [-1,1]:knots.add((end+side*ell/2)%1)
 knots=sorted(knots);den=F(0);num=[F(0),F(0)];parts=[]
 for left,right in zip(knots,knots[1:]):
  rs=[[integ(s,x,ell)for s in classes]for x in [left,right]]
  av=[(q+sum(row))/L for row in rs]
  zv=[[(q*eta[c]+row[c])/L for c in range(2)]for row in rs]
  h=right-left;den+=affine_product(h,av[0],av[1],av[0],av[1])
  for c in range(2):num[c]+=affine_product(h,av[0],av[1],zv[0][c],zv[1][c])
  parts.append(dict(left=str(left),right=str(right),a_endpoints=list(map(str,av)),z_endpoints=[[str(v)for v in row]for row in zv]))
 mu=[v/den for v in num];assert sum(mu)==1 and mu[0]>F(1,2)
 jensen_factor=sum(mu[c]**2/eta[c]for c in range(2));assert jensen_factor>1
 # Independent direct grid check: all knots on grid1/40; integrate degree2 exactly
 # with cell Simpson. This is a separate subdivision from source knots above.
 ds=F(0);ns=[F(0),F(0)]
 for j in range(80):
  xs=[F(j,80),F(2*j+1,160),F(j+1,80)]
  aa=[];zz=[]
  for x in xs:
   rc=[integ(s,x,ell)for s in classes];aa.append((q+sum(rc))/L);zz.append([(q*eta[c]+rc[c])/L for c in range(2)])
  ds+=sum(w*a*a for w,a in zip([1,4,1],aa))/480
  for c in range(2):ns[c]+=sum(w*a*z[c]for w,a,z in zip([1,4,1],aa,zz))/480
 assert ds==den and ns==num
 out=dict(source_classes=[[[str(v)for v in row]for row in s]for s in classes],q=q,ell=str(ell),L=str(L),score_tilt_t=1,
  eta=list(map(str,eta)),m_t_exact=str(den),tilted_numerator_exact=list(map(str,num)),mu_t_exact=list(map(str,mu)),
  old_prior_Jensen_factor_exact=str(jensen_factor),old_prior_Jensen_factor_float=float(jensen_factor),
  partitions=parts,independent_80_cell_Simpson_exact_match=True,
  scope='Strict physical source example of old-prior denominator mismatch; not a C_n bound',
  script_SHA256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 target.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:out[k]for k in ['eta','m_t_exact','mu_t_exact','old_prior_Jensen_factor_exact','old_prior_Jensen_factor_float']}))
if __name__=='__main__':main()
