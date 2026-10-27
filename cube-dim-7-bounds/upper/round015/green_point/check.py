#!/usr/bin/env python3
"""Exact-dyadic accounting and target comparisons, with no phase re-evaluation."""
from pathlib import Path
from fractions import Fraction
import json,hashlib
p=Path(__file__).resolve().parent;sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest();M=json.loads((p/'manifest.json').read_text())
for f,h in M['files'].items():assert sha(p/f)==h,f
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,f=line.split('  ',1);q=Path(f);assert sha(q if q.is_absolute() else p/q)==h,f
oldM=json.loads((p/'r012-manifest.json').read_text())
for f,k in [('r012-result.json','result.json'),('cells.jsonl','cells.jsonl'),('r009-rays.json','r009-rays.json'),('r009-inputs.json','r009-inputs.json'),('r010-manifest.json','r010-manifest.json'),('testcases.json','testcases.json')]:assert sha(p/f)==oldM['files'][k]
original=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube7-kernel-full-r012/cube-dim-7-bounds/upper/round012/kernel_full/kernel.py');assert sha(original)==oldM['files']['kernel.py']
derived=original.read_text().replace("if k==0:counts['log_evaluations']+=1;out=(A(beta*beta+key[0]*key[0])).log()/2\n elif key[0]==key[-1]:\n  counts['hermite_evaluations']+=1;z=acb(A(beta),-key[0]);out=(-(acb(0,1)**k)/(k*z**k)).real","if k==0:\n  counts['log_evaluations']+=1;assert key[0]!=0;out=A(abs(key[0])).log()\n elif key[0]==key[-1]:\n  counts['hermite_evaluations']+=1;assert key[0]!=0;out=A(Fraction((-1)**(k-1),k*key[0]**k))").replace("B=dx*L;beta=eps*B;ar=[sum(n[j]*ZETA[i][j] for j in range(7)) for i in range(46)];cnt=counter();cnt['zero_ray_knots']=sum(a==0 for a in ar);top={}","assert eps==0;B=dx*L;beta=0;ar=[sum(n[j]*ZETA[i][j] for j in range(7)) for i in range(46)];cnt=counter();cnt['zero_ray_knots']=sum(a==0 for a in ar);assert cnt['zero_ray_knots']==0, f'zero exact knot in phase {phase_id}';top={}").replace('phase(n,dx,Fraction(1,1024),64*chunk+signid)','phase(n,dx,Fraction(0),64*chunk+signid)');assert derived==(p/'green_kernel.py').read_text()
old=json.loads((p/'r012-result.json').read_text());oldterminal=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube7-kernel-full-r012/cube-dim-7-bounds/upper/round012/kernel_full/terminal.jsonl');assert sha(oldterminal)==oldM['files']['terminal.jsonl'];OT=[json.loads(x) for x in oldterminal.read_text().splitlines()];assert len(OT)==5040 and all(x['status']=='completed' for x in OT) and {x['chunk'] for x in OT}==set(range(5040)) and all(x['counts']['zero_ray_knots']==0 for x in OT) and sum(x['counts']['zero_ray_knots'] for x in OT)==old['counts']['zero_ray_knots']==0 and old['halfphases']==322560
R=json.loads((p/'result.json').read_text());P=json.loads((p/'pilot.json').read_text());T=[json.loads(x) for x in (p/'terminal.jsonl').read_text().splitlines()];A=json.loads((p/'assessment.json').read_text());assert R['status']=='completed' and P['status']=='passed' and len(T)==5040 and {x['chunk'] for x in T}==set(range(5040)) and all(x['status']=='completed' for x in T) and [x['chunk'] for x in T[:4]]==[0,1,2,3]
def fr(x):return Fraction(int(x[0]))*Fraction(2)**int(x[1])
lo=Fraction();hi=Fraction();tot={k:0 for k in R['counts']}
for x in T:
 j=x['chunk'];assert x['precision']==256 and x['start_phase']==64*j and x['end_phase']==64*(j+1) and x['counts']['top_requests']==64*212 and x['counts']['zero_ray_knots']==0
 l=fr(x['interval']['lower']);u=fr(x['interval']['upper']);assert l<=u;lo+=l;hi+=u
 for k in tot:tot[k]+=x['counts'][k]
assert tot==R['counts'] and tot['top_requests']==68382720 and R['halfphases']==322560 and R['epsilon']==0
L=fr(R['full_G']['lower']);U=fr(R['full_G']['upper']);assert 0<L<=lo<=hi<=U and U-L<=Fraction(1,10**12) and fr(R['width']['upper'])>=U-L
assert fr(old['full_H']['upper'])<L
hexes=R['point_floathex'];assert hexes==old['point_floathex']==json.loads((p/'testcases.json').read_text())['generic_case']['x_floathex'];x=[Fraction.from_float(float.fromhex(v)) for v in hexes];assert x[0]==1 and all(0<v<1 for v in x[1:]) and len(set(map(abs,x)))==7
assert all((U<Fraction(1,k))==A['thresholds'][str(k)]['strict'] for k in (219000,194000,179000))
I0=Fraction(7723627498344,10**18);cost6=Fraction(7**7)*I0**6/U**7;assert Fraction(A['G_upper_exact'])==U and Fraction(A['fixed_I0'])==I0 and Fraction(A['fixed_cost_power6_lower_exact'])==cost6 and (cost6>128**6)==A['fixed_cost_power6_strictly_gt_128_power6']
assert all(A['thresholds'][str(k)]['strict'] for k in (219000,194000,179000)) and cost6>128**6
print('PASS: R012 exact nonzero evidence and point, all 5040 new G chunk terminals, exact outward signed sum, H<G, off-S3 distinctness, three rational thresholds and fixed-I0 cost comparison. No phase recomputed; no C7 lower bound.')
