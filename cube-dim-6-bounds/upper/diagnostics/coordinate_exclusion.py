"""Exact rational integration and outward-arithmetic lower cost on 3^5 cells."""
import sympy as s
from flint import arb,ctx
from pathlib import Path
import json,time,hashlib,math
start=time.time();ctx.prec=256;t=s.symbols('t',positive=True);z=s.symbols('z');rows=[];total=arb(0)
for a in range(6):
 for b in range(6-a):
  c=5-a-b
  f=t**6/((t*t+s.Rational(1,9))**a*(t*t+s.Rational(4,9))**b*(t*t+1)**(c+1))
  fz=z**3/((z+s.Rational(1,9))**a*(z+s.Rational(4,9))**b*(z+1)**(c+1))
  terms=s.Add.make_args(s.apart(fz,z)); recon=0;v=0
  for term in terms:
   num,den=s.fraction(term); dp=s.Poly(den,z); degree=dp.degree(); lead=dp.LC()
   root_squared=s.cancel(dp.nth(degree-1)/(degree*lead))
   assert s.expand(den-lead*(z+root_squared)**degree)==0
   coeff=num/lead;recon+=coeff/(z+root_squared)**degree
   scale=s.sqrt(root_squared)
   v+=coeff*s.pi*s.binomial(2*degree-2,degree-1)/(2*4**(degree-1)*scale**(2*degree-1))
  assert s.cancel(recon-fz)==0
  ratio=s.simplify(v/(3*s.pi/512));assert ratio.is_Rational and ratio>0
  multiplicity=math.factorial(5)//(math.factorial(a)*math.factorial(b)*math.factorial(c))
  total+=multiplicity*arb(str(ratio))**(arb(6)/5)/243
  rows.append({'counts':[a,b,c],'multiplicity':multiplicity,'G_over_g_exact':str(ratio)})
cost=arb(6)**(arb(6)/5)*total
out={'status':'rigorous_model_cost_lower_bound_via_monotone_cells','grid':3,'cells':243,'symmetry_orbits':len(rows),'rows':rows,'lower_cost_arb':str(cost),'cost_strictly_above_64':bool(cost>64),'elapsed_seconds':time.time()-start,'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
assert sum(r['multiplicity'] for r in rows)==243
Path(__file__).with_name('coordinate_exclusion_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))
