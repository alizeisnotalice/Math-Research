#!/usr/bin/env python3
"""Integrity, four fixed scalar rechecks, origin and phase-output assembly; no new phases."""
from flint import arb,acb,fmpq,ctx
from fractions import Fraction
from pathlib import Path
import hashlib,itertools,json,math
ctx.prec=256;p=Path(__file__).resolve().parent;sha=lambda x:hashlib.sha256(x.read_bytes()).hexdigest();mf=json.loads((p/'manifest.json').read_text())
for f,h in mf['files'].items():assert sha(p/f)==h,f
for line in (p/'prerun-sha256.txt').read_text().splitlines():
 h,f=line.split('  ',1);q=Path(f);assert sha(q if q.is_absolute() else p/q)==h,f
om=json.loads((p/'r010-manifest.json').read_text())
for f in ('cells.jsonl','r009-rays.json','r009-inputs.json'):assert sha(p/f)==om['files'][f]
R=json.loads((p/'r009-rays.json').read_text());N=json.loads((p/'r009-inputs.json').read_text())['prefix_normals'];C=[json.loads(x) for x in (p/'cells.jsonl').read_text().splitlines()];q=[sum(abs(sum(a*b for a,b in zip(n,r['w_primitive']))) for n in N) for r in R];assert len(R)==46 and len(C)==135
weights=[];seen=set()
for c in C:
 for s in c['simplices']:
  key=tuple(s['vertices']);assert key not in seen;seen.add(key);assert s['q']==[q[i] for i in key] and Fraction(s['alpha'])==Fraction(abs(s['detW']),math.prod(s['q']));weights.append(Fraction(s['alpha']))
assert len(seen)==212;Asum=sum(weights,Fraction());assert 0<5040*Asum<=Fraction(1,512)
d=json.loads((p/'results.json').read_text());T=json.loads((p/'testcases.json').read_text());assert d['status']=='completed' and d['A_sum']==str(Asum) and d['q_lcm']==math.lcm(*q) and d['MC_samples']==0 and not d['promotions']
def A(x):return arb(fmpq(x.numerator,x.denominator)) if isinstance(x,Fraction) else arb(x)
def dy(v):return arb(int(v[0]))*arb(2)**int(v[1])
def ball(v):return dy(v['lower']).union(dy(v['upper']))
def fr(v):return Fraction(int(v[0]))*Fraction(2)**int(v[1])
def complexdd(a,beta):
 cache={}
 def f(z):
  if z in cache:return cache[z]
  k=len(z)-1
  if not k:r=acb(A(beta),-z[0]).log()
  elif z[0]==z[-1]:r=-(acb(0,1)**k)/(k*acb(A(beta),-z[0])**k)
  else:r=(f(z[1:])-f(z[:-1]))/(z[-1]-z[0])
  cache[z]=r;return r
 return f(tuple(a))
for s,c in zip(d['scalar_cases'],T['scalar_DD_cases']):
 assert s['name']==c['name'];k=[Fraction(x) for x in c['knots']];B=math.lcm(*(x.denominator for x in k));nodes=tuple(sorted(int(x*B) for x in k));z=complexdd(nodes,Fraction(B,64))*B**6;assert z.real.overlaps(ball(s['complex_real'])) and z.imag.overlaps(ball(s['complex_imag'])) and z.real.overlaps(ball(s['real']))
 if c['name']=='all_zero':assert z.overlaps(acb(A(Fraction(64**6,6))))
 if c['name']=='all_one':assert z.overlaps(1/(6*acb(A(Fraction(1,64)),-1)**6))
 if c['name']=='near_distinct':
  assert len(set(nodes))==7;eps=Fraction(1,64);M=32;theta=max(abs(x) for x in k)/eps;h=[Fraction(0)]*(M+1);h[0]=1
  for x in k:
   for j in range(1,M+1):h[j]+=x*h[j-1]
  val=acb(0)
  for j in range(M+1):
   degree=6+j;val+=-(acb(0,1)**degree)*A(h[j]/(degree*eps**degree))
  tail=Fraction(1,6)/eps**6*math.comb(M+6,5)*theta**(M+1)/(1-theta*Fraction(M+7,M+2));radius=A(tail);ref=acb(val.real+arb(0,radius),val.imag+arb(0,radius));assert z.overlaps(ref)
for o in d['origin']:
 eps=Fraction(o['epsilon']);ref=840*A(Asum)/(arb.pi()**7*A(eps)**6);assert ref.overlaps(ball(o['interval'])) and o['dd_identity_pass'] and o['opposite_factor']==1
bycase={}
for x in d['symmetry']:
 case=x['case'];eps=x['epsilon'];n=7 if case=='e1' else 64;st=46080 if case=='e1' else 5040;assert x['precision']==256 and x['phase_count']==n and x['stabilizer']==st and len(x['phase_records'])==n and 2*n*st==645120
 assert [v['phase_id'] for v in x['phase_records']]==list(range(n));total=sum((ball(v['raw_phase']) for v in x['phase_records']),arb(0));v=2*st*total/(2*arb.pi())**7;assert v.overlaps(ball(x['interval'])) and x['positive'] and x['width_pass']
 assert fr(x['interval']['upper'])-fr(x['interval']['lower'])<=Fraction(1,10**12);bycase[(case,eps)]=ball(x['interval'])
for case in ('e1','ones'):assert bycase[(case,'1/256')]>bycase[(case,'1/64')]
gens=sorted(d['generic_chunks'],key=lambda x:x['chunk']);assert [x['chunk'] for x in gens]==[0,1,2,3] and all(x['precision']==256 for x in gens)
for x in gens:
 assert [v['phase_id'] for v in x['phase_records']]==list(range(x['start_phase'],x['end_phase'])) and x['end_phase']-x['start_phase']==64
 raw=sum((ball(v['raw_phase']) for v in x['phase_records']),arb(0))*2/(2*arb.pi())**7;assert raw.overlaps(ball(x['interval']))
total=sum((ball(x['interval']) for x in gens),arb(0));assert total.overlaps(ball(d['generic_partial'])) and d['generic_partial_is_not_H'] and d['generic_planning_width_pass']
assert fr(d['generic_partial']['upper'])-fr(d['generic_partial']['lower'])<=Fraction(1,10**12)*Fraction(256,322560)
for x in d['symmetry']+gens:
 for phase in x['phase_records']:
  c=phase['counts'];assert c['top_requests']==212 and 0<c['unique_top']<=212 and c['log_evaluations']<=46
src=Path('/Users/zhengzhihao/Desktop/T/Math-Research-cube-dim-7/cube-dim-7-bounds/upper/round003/truncated_numeric/results.json');r=json.loads(src.read_text())['minima'][0];assert r['lowest_sobol_index']==2710 and r['point_floathex']==T['generic_case']['x_floathex']
primary=json.loads(Path('/Users/zhengzhihao/Documents/Codex/2026-10-27/n-n-5-1-1-verify/work/r011-primary-euclidean-floor.json').read_text());assert primary['max_squared_unitball_norm']=='1/2048' and primary['psi_euclidean_lower_squared']==2048
print('PASS: hashes and 212 exact weights, four fixed scalar complex rechecks, three origin constants, four complete symmetry orbit aggregates, 256 prescribed generic phases and width gate. No new phase evaluated; partial is not H.')
