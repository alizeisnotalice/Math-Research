#!/usr/bin/env python3
"""Rebuild saved rational brackets and arithmetic; shares formulas with generator."""
from fractions import Fraction as F
from decimal import Decimal
from pathlib import Path
from certify import case
import json,time
base=Path(__file__).resolve().parent

def verify_record(rec):
 n=rec['n'];k=F(rec['kappa']);u=F(rec['u']);assert u.numerator==1
 assert rec.get('root') is not None or rec['status']=='guard_failed'
 fresh,exact=case(n,k,u.denominator,root=rec.get('root'),return_exact=True)
 for key in ('status','guards','M','lambda'):assert fresh[key]==rec[key],(n,u,key)
 if exact is None:
  assert rec['status']=='guard_failed';return 'guard_failed'
 for key in ('root','budget_sign','budget_eligible','G_sign','admissible_G_sign','predicted_limit'):
  assert fresh[key]==rec[key],(n,u,key)
 low,high,G=exact;bl,bh=map(lambda x:F(Decimal(x)),rec['D_minus_B']);gl,gh=map(lambda x:F(Decimal(x)),rec['G'])
 assert bl<=low<=high<=bh and gl<=G<=gh,(n,u,'outward')
 return 'budget_'+rec['budget_sign']

def verify_center(r):
 n=r['n'];u=F(r['u']);k=F(r['kappa']);lam=F(r['lambda']);a=(2*u)**n;w=u/(1+2*u)+lam*u**6;rho=1-2*w;e=k*u**3;det=1-e**2;mout=(2*w)**n;M=a*rho**n;m0=(2*u/(1+2*u))**n;N=F(1,u**7);assert N.denominator==1 and str(N.numerator)==r['N'];N=N.numerator
 upper=(M+2*e**2*(2*w)**(n-2))/(2*(rho+w))**n;margin=u**n-upper;fd=F(2*N-5,2*N)**n*M/(det-M)-m0/(1-m0)
 assert 2*M>mout and r['center_in_A_conditions']['M_gt_mout_over_2'] and r['budget_sign']=='positive'
 assert upper==F(r['center']['defect_upper']) and margin==F(r['center']['strict_margin']) and margin>0 and r['center']['strict_margin_positive']
 assert fd==F(r['finite_halo_minus_Ln']['exact_fraction']) and fd>0 and r['finite_halo_minus_Ln']['strictly_positive'] and r['finite_status']=='certified'

if __name__=='__main__':
 t=time.perf_counter();old=json.loads((base/'cases.json').read_text())['cases'];new=json.loads((base/'center_finite.json').read_text())['cases'];assert len(old)==120 and len(new)==15
 records=old+[r for item in new for r in item['attempts']];assert len(records)==135
 counts={}
 for r in records:
  z=verify_record(r);counts[z]=counts.get(z,0)+1
 for item in new:verify_center(item['attempts'][-1])
 assert counts=={'budget_positive':105,'budget_negative':29,'guard_failed':1},counts
 print(json.dumps({'verify':'PASS','records':len(records),'counts':counts,'center_exact_checks':len(new),'seconds':time.perf_counter()-t,'independence':'Reuses generator formulas; saved brackets are never searched again.'}))
