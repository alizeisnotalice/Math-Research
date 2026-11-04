"""Finite-input gains and actual witness averages for the full symmetric lattice.
Exact arithmetic; distant lattice indices are excluded by the paper's proof.
For n<=6 the local {-1,0,1}^n block is also summed directly.
"""
from fractions import Fraction as F
from pathlib import Path
from math import comb,prod
from itertools import product
import argparse,json
from compute import evaluate

def exact(x):return {'numerator_hex':hex(x.numerator),'denominator_hex':hex(x.denominator)}
def overlap(c,x,r,d):return max(F(0),min(x+r,c+d)-max(x-r,c-d))
def centre_mass_dp(n,e,d):
 coef={0:F(1)}
 for _ in range(n):
  nxt={}
  for p,val in coef.items():
   for shift,weight in [(-1,e),(0,2*d),(1,e)]:nxt[p+shift]=nxt.get(p+shift,F(0))+val*weight
  coef=nxt
 return coef[0]
def actual_mass(indices,x,r,d,e):
 total=F(0)
 for k in indices:
  S=sum(k);c=[(1-e)*z+e*S for z in k]
  total+=prod(overlap(z,y,r,d) for z,y in zip(c,x))
 return total

def make(rec):
 n=rec['n'];kap=F(rec['kappa']);u=F(rec['u']);den=u.denominator;delta=u/(1+2*u);r0=1/(1+2*u);e=kap*u**4;m=(2*delta)**n;det=(1-e)**(n-1)*(1+(n-1)*e);v=2*u;a=v**n;N=den**9
 fresh,(gl,gu)=evaluate((n,kap,den),rec['roots']);assert fresh['status']=='certified'
 common=(det-m)*(1-m)*a*v**8;bu=(m*(1-det)-gl*common)/(1-m);assert 0<bu<m
 finite=1+F(2*N-5,2*N)**n*(m-bu)/(det-m);margin=finite-1/(1-m);assert margin>0
 hx=e**2/(2*n*u**2*r0);assert 0<hx<e/2<delta
 samples=[]
 for name,x in [('centre',[F(0)]*n),('near_corner',[delta-hx,-delta+e/2]+[F(0)]*(n-2))]:
  s=delta-max(map(abs,x));r=r0+s
  if name=='centre':
   mass=centre_mass_dp(n,e,delta)
   formula=sum((F(comb(n,q)*comb(q,q//2))*e**q*(2*delta)**(n-q) for q in range(0,n+1,2)),F(0));assert mass==formula
   method='exact Laurent coefficient recurrence; balanced-index formula cross-check'
  else:
   assert e+s<delta
   sparse=[(i,j)+tuple([0]*(n-2)) for i in range(-2,3) for j in range(-2,3)]
   mass=actual_mass(sparse,x,r,delta,e)
   formula=m+(2*delta)**(n-2)*e*(e+hx-e/2);assert mass==formula
   method='25 exact interval-intersection products; other balanced indices have a zero factor'
  if n<=6:
   direct=actual_mass(product((-1,0,1),repeat=n),x,r,delta,e);assert direct==mass;method+='; full local integer block independently summed'
  defect=mass/(2*r)**n;strict=u**n-defect;assert strict>0
  assert r+delta==1-max(map(abs,x)) and (delta+r)/(1-(n-1)*e)<2<N-3
  samples.append({'name':name,'center':[str(z) for z in x],'radius':str(r),'input_average':str(1-defect),'strict_margin_over_threshold':str(strict),'mass_method':method})
 return {'n':n,'kappa':str(kap),'u':str(u),'finite_input':{'definition':'E_N=(R^n minus T Z^n axis-cube pores) intersect T[-N,N]^n','matrix':'T=(1-e)I+eJ','shear':str(e),'pore_halfwidth':str(delta),'N':N},'finite_ratio_minus_Ln_lower_exact':exact(margin),'finite_normalized_gap_lower_exact':exact(margin/(a*v**8)),'witnesses':samples}

def main():
 p=argparse.ArgumentParser();p.add_argument('--input',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args();doc=json.loads(args.input.read_text());chosen=[r for r in doc['cases'] if F(r['u'])==F(1,2**20)]
 assert sorted(r['n'] for r in chosen)==list(range(2,17));rows=[make(r) for r in chosen]
 data={'scope':'15 explicit finite inputs, 30 exact sample averages, and direct local-block checks in n=2..6; geometry covers omitted indices','all_finite_gains_positive':True,'all_witness_margins_positive':True,'cases':rows};args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(data,indent=2)+'\n');print(json.dumps({'finite_inputs':len(rows),'sample_witnesses':sum(len(r['witnesses']) for r in rows),'finite_gains_positive':True,'witness_margins_positive':True,'direct_block_dimensions':[2,3,4,5,6]}))
if __name__=='__main__':main()
