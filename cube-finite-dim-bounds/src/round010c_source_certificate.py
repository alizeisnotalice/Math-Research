"""Rigorous one-point rejection of an exact rational realization of batch B."""
from pathlib import Path
from fractions import Fraction as F
import json,itertools,hashlib,time
from flint import arb,acb,ctx
P=Path(__file__).resolve().parents[1]
def ball(q):
 q=F(q);return arb(q.numerator)/q.denominator
def ends(x):return [str(x.lower().fmpq()),str(x.upper().fmpq())]
def dot(x,y):return sum(a*b for a,b in zip(x,y))
V=[(F(1),)+tuple(map(F,v)) for v in itertools.product((-1,1),repeat=2)]
def lb_unit(x):
 e=dot(x,x);total=arb(0)
 for v in V:
  d=dot(x,v);aa=-d/3;A=1-e+d*d/3
  if A<=0:continue
  a=ball(aa);h=ball(A/3).sqrt()
  if aa==0 or aa*aa==A/3:term=arb(0)
  else:term=a*(h*h-a*a)*abs((a+h)/(a-h)).log()
  total+=36/arb.pi()*(2*a*a*h-ball(F(4,3))*h**3+term)
 return total
if __name__=='__main__':
 started=time.time();ctx.prec=192
 inp=P/'results/round010c_source_input.json'
 if not inp.exists():
  raw=P/'results/round010b_joint_probe.json';data=json.loads(raw.read_text());r=min((r for r in data['validation'] if r['new_source_slack'] is not None),key=lambda r:r['new_source_slack'])
  z={'point':[str(F(str(t))) for t in r['x']],'weights':[str(F(str(t))) for t in data['weights']],'specs':[{'seed':[str(F(str(t))) for t in b['seed']],'radius':str(F(str(b['radius'])))} for b in data['specs']],'numerical_source_sha256':hashlib.sha256(raw.read_bytes()).hexdigest(),'interpretation':'exact decimal-rational realization, not floating-point negativity as proof'}
  inp.write_text(json.dumps(z,indent=2)+'\n')
 z=json.loads(inp.read_text());x=list(map(F,z['point']));e=dot(x,x);T=4*arb(2).log()-1;c=T/(24*arb.pi()**3);R=15/(2*T);assert ball(e)<R
 total=acb(0);pieces=[]
 for v in V:
  d=dot(x,v);root=(ball(d*d)+3*(R-ball(e))).sqrt()
  for sg in (-1,1):
   U=(root+sg*ball(d))/(R-ball(e))
   def f(t,analytic):
    u=acb(U)*t
    return acb(U)*(1-acb(R)*u*u/(acb(ball(e))*u*u+sg*2*acb(ball(d))*u+3))
   value=acb.integral(f,0,1,rel_tol=arb('1e-42'),abs_tol=arb('1e-46'),eval_limit=20000,depth_limit=35,deg_limit=80)
   assert value.is_finite() and value.imag.contains(0)
   total+=value;pieces.append(str(value.real))
 mu=c/arb.pi()*total.real;contributions=[]
 for w,spec in zip(z['weights'],z['specs']):
  seed=tuple(map(F,spec['seed']));rad=F(spec['radius']);centers=sorted(set(tuple(a*b for a,b in zip(t,sg)) for t in set(itertools.permutations(seed)) for sg in itertools.product((-1,1),repeat=3)))
  result=sum((lb_unit([(a-b)/rad for a,b in zip(x,center)])/ball(rad) for center in centers),arb(0));contributions.append(ball(w)*result)
 LB=sum(contributions,arb(0));negative=mu-LB;assert negative<0
 out={'status':'strict negative source witness for one exact rational eight-basis correction; primary checked, not yet independently audited','input_sha256':hashlib.sha256(inp.read_bytes()).hexdigest(),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'fixed_c':'(4log(2)-1)/(24pi^3)','G_lower':'5/(16pi^3|x|^2)','point_norm_squared':str(e),'point_in_Jensen_interior_of_Omega':True,'mu_upper_expression_interval':ends(mu),'mu_upper_expression_ball':str(mu),'LB_interval':ends(LB),'LB_ball':str(LB),'corrected_density_upper_interval':ends(negative),'corrected_density_upper_ball':str(negative),'basis_contributions':list(map(str,contributions)),'precision_bits':192,'scope':'rejects only this rational correction for fixed c; no cost/g assumption and no family or Cn bound','seconds':time.time()-started}
 (P/'results/round010c_source_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
