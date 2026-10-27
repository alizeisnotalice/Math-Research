#!/usr/bin/env python3
"""Reassemble exact contractions and three fixed gates, no polynomial/phase recomputation."""
from pathlib import Path
from fractions import Fraction
from flint import arb,fmpq,ctx
import hashlib,json,math
ctx.prec=512;p=Path(__file__).resolve().parent;sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest();mf=json.loads((p/'manifest.json').read_text())
for f,h in mf['files'].items():assert sha(p/f)==h,f
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,f=line.split('  ',1);q=Path(f);assert sha(q if q.is_absolute() else p/q)==h,f
old=json.loads((p/'r010-manifest.json').read_text())
for f in ('cells.jsonl','r009-rays.json','r009-inputs.json'):assert sha(p/f)==old['files'][f]
assert sha(p/'r016-result.json')==json.loads((p/'r016-manifest.json').read_text())['files']['result.json']
C=[json.loads(x) for x in (p/'cells.jsonl').read_text().splitlines()];S=[s for c in C for s in c['simplices']];T=[json.loads(x) for x in (p/'terminal.jsonl').read_text().splitlines()];T.sort(key=lambda x:x['id']);d=json.loads((p/'result.json').read_text());assert len(S)==len(T)==212 and {x['id'] for x in T}==set(range(212)) and all(x['status']=='completed' for x in T)
Q={i:q for s in S for i,q in zip(s['vertices'],s['q'])};L=math.lcm(*Q.values());assert L==d['L_common'];R0=5040*sum((Fraction(s['alpha']) for s in S),Fraction());assert R0==Fraction(d['R']['R0'])
names=['R1','R_part2','R11','R3','R21','R111'];factors=[720,720,120,720,120,24];powers=[2,4,4,6,6,6]
for i,x in enumerate(T):assert x['vertices']==S[i]['vertices'] and x['alpha']==S[i]['alpha'] and len(x['contractions'])==6 and all(int(z)>0 for z in x['contractions']) and x['multiplication_pairs']<=59535 and x['max_support']<=924
for j,name in enumerate(names):
 v=Fraction(factors[j],L**powers[j])*sum((Fraction(x['alpha'])*int(x['contractions'][j]) for x in T),Fraction());assert v==Fraction(d['R'][name]) and Fraction(d['rho'][name])==v/R0
assert sum(x['multiplication_pairs'] for x in T)==d['multiplication_pairs']<=12621420 and max(x['max_support'] for x in T)==d['max_support'] and max(x['max_coefficient_bits'] for x in T)==d['max_coefficient_bits']
def A(x):return arb(fmpq(x.numerator,x.denominator))
def dy(x):return Fraction(int(x[0]))*Fraction(2)**int(x[1])
def box(x):return A(dy(x['lower'])).union(A(dy(x['upper'])))
p0=A(R0)/arb.pi()**7;assert p0.overlaps(box(d['p0']));a=Fraction(71,2);b=Fraction(65*67*69,8);rho={k:Fraction(v) for k,v in d['rho'].items()};target=dy(json.loads((p/'r016-result.json').read_text())['enclosures']['g_required_B_RB']['lower']);assert str(target)==d['target_lower'];TC=json.loads((p/'testcases.json').read_text());actual=json.loads(Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube7-green-r015/cube-dim-7-bounds/upper/round015/green_point/result.json').read_text());assert actual['point_floathex']==TC['generic_case']['x_floathex']
pts=[[Fraction(1)]*7,[Fraction.from_float(float.fromhex(x)) for x in TC['generic_case']['x_floathex']],[Fraction(1)]+[Fraction(0)]*6]
for x,z in zip(pts,d['gates']):
 assert z['point']==[str(t) for t in x];s=[t*t for t in x];p1=sum(s);p2=sum(t*t for t in s);p3=sum(t**3 for t in s);e2=(p1*p1-p2)/2;e3=(p1**3-3*p1*p2+2*p3)/6
 m=rho['R1']*p1/(2*a);u2=(rho['R_part2']*p2+6*rho['R11']*e2)/(12*a*(a+1));u3=(rho['R3']*p3+15*rho['R21']*(p1*p2-p3)+90*rho['R111']*e3)/(120*a*(a+1)*(a+2));v=u2-m*m;w=u3-2*m*u2+m**3;h=m*u3-u2**2;assert m>0 and v>0 and w>0 and h>0
 inv=(w**3+2*w*v**3+m*v**4)/h**3;assert z['moments']=={'m':str(m),'V':str(v),'W':str(w),'H':str(h),'gauss_inverse3':str(inv)};g=p0*A(inv/b);assert g.overlaps(box(z['interval'])) and (dy(z['interval']['upper'])<target)==z['upper_below_target']
assert d['corner_gate_failed'] and d['gates'][0]['upper_below_target'] and not d['global_cover_started'] and d['scalar_outputs']==10 and d['toy_law_pass']
print('PASS: 212 terminal contractions, 7 exact rational outputs, factorial normalization, 3 fixed rational-moment gates and strict corner failure; no polynomial/phase recomputation or global-g claim.')
