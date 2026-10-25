"""New satellite-interface duplicate-position and total-mass scaling check."""
from pathlib import Path
import json
import numpy as np
from fractions import Fraction as F
from round004_search import primary_keys,value,cert,u,v
ROOT=Path(__file__).resolve().parent
keys=[]
for k in primary_keys():
 lo=k[:3];hi=k[6:9]
 if all(l<=h for l,h in zip(lo,hi)):keys.append(tuple(lo)+tuple(lo)+tuple(hi)+tuple(hi))
formula=sum(float(F(q))*u**int(ex.split(',')[0])*v**int(ex.split(',')[1]) for ex,q in cert['coefficients'].items());rows=[]
for eps in [2**-8,2**-6,2**-4]:
 z=np.zeros(7);z[:3]=np.log1p(eps)/3;z[6]=np.log(eps);val=value(z,keys)
 assert abs(val-formula)<2e-13
 rows.append({'epsilon':eps,'period_scale':float((1+eps)**(1/3)),'phase':[0,0,0],'both_orbits_same_blocks':True,'key_count':len(keys),'joint_ratio':val,'original_formula_display':formula,'difference_from_formula':val-formula})
(ROOT/'round004_calibration.json').write_text(json.dumps({'status':'floating new-interface scale/duplicate-position calibration, coefficients read from original exact certificate','cases':rows},indent=2)+'\n');print(rows)
