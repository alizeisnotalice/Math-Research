"""Reproduce the fixed MATLAB input from previous minimal certificates."""
from pathlib import Path
from fractions import Fraction as Q
import json
ROOT=Path(__file__).resolve().parent
root=ROOT.parent
summary=json.loads((root/'round03/first_checker_summary.json').read_text())
points=[]
for n,b in enumerate(summary['first_K_fail_examples']):
    u=(Q(b[0])+Q(b[1]))/2;v=(Q(b[2])+Q(b[3]))/2
    if u+v>1:v=1-u
    assert Q(b[0])<=u<=Q(b[1]) and Q(b[2])<=v<=Q(b[3]) and u+v<=1
    points.append({'id':'K'+str(n),'u':str(u),'v':str(v)})
kernels=[]
for name,path in [('K1',root/'round01/scaled_certificate.json'),('K3',root/'round03/candidate.json')]:
    s=json.loads(path.read_text())
    kernels.append({'name':name,'c':s['c'],'alpha':s['alpha'],'R':s['R'],'h':s['h'],'classes':[{'i':int(k.split(',')[0]),'j':int(k.split(',')[1]),'a':v} for k,v in s['a'].items()]})
out=ROOT/'matlab_input.json'
if out.exists():raise SystemExit('Input exists; do not overwrite recorded input')
out.write_text(json.dumps({'points':points,'kernels':kernels},ensure_ascii=False))
