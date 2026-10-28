#!/usr/bin/env python3
"""New bounded R081 sensitive constants, not a rerun of author 13 checks.

Uses 16 positive atanh terms and exact tails; no threshold is certified.
"""
from fractions import Fraction as F
from pathlib import Path
import hashlib,json,os,platform

HERE=Path(__file__).resolve().parent
SOURCE=Path('/Users/zhengzhihao/Desktop/T/cube-r078-theory-bernstein-sparse-20261028/cube-dim-7-bounds/upper/round081/general_s')

def logs(v):
    z=(v-1)/(v+1)
    lo=2*sum((z**(2*j+1)/F(2*j+1) for j in range(16)),F(0))
    return lo,lo+2*z**33/(33*(1-z*z))

e=F(1,10000);r0=1+F(1,100000);r1=F(99,25)-F(1,100000)
P0=1-2*e;ulo=P0*(1-e)/(1+e);uhi=(1+e)/(P0*(1-e))
M=2*(1+e)*uhi;b=2*F(999,1000)*(1-e)*ulo-2*e
Gamma=F(137,100)*b/M
ml,_=logs(r1/r0);Qlo,_=logs(F(8,5));_,pricehi=logs(F(5,2))
gaps={'paid_deficit_minus_1p3':(1-e)*Gamma-F(13,10),
 'mass_minus_1p37':(1-e)*ml-F(137,100),
 'mu4_slack_under_065':F(13,20)-(3-5*Qlo),
 'log_price_slack':F(23,25)-pricehi,
 'new_complete_cost_vs_Phi80_separator':F(3,8)-F(28029,100000)-F(4,4096)}
assert all(v>0 for v in gaps.values())
assert F(13,10)-F(23,25)-F(3,8)==F(1,200)
hashes={str(SOURCE/n):hashlib.sha256((SOURCE/n).read_bytes()).hexdigest()
 for n in ('PROOF.md','verify_constants.py','verification.json')}
out={'status':'PASS fixed constants only; NO explicit threshold certified','pid':os.getpid(),
 'python':platform.python_version(),'author_program_executed':False,
 'recipe':'Fraction and independent log16 exact geometric enclosures',
 'source_sha256':hashes,'exact_gaps':{k:str(v) for k,v in gaps.items()},
 'display_only':{k:float(v) for k,v in gaps.items()}}
(HERE/'full_result.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('source_sha256','exact_gaps')},indent=2))
