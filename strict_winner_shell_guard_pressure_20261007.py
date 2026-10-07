"""Exact pressure of the source-box Gamma_v lower guard with a strict b winner.
Continuous positive L1 source; hard relaxation only, no softFIRST qualification.
"""
from fractions import Fraction as F
from pathlib import Path
import json,hashlib
out=Path(__file__).with_name('strict_winner_shell_guard_pressure_20261007_results.json')
if out.exists():raise SystemExit('Refuse overwrite')
rows=[]
for stage,n in enumerate((4,16,64),1):
 a=F(1);b=F(2);M=n*n*b;ell=M-b
 q=(ell/M)**n*(3+(ell/M)**2)/4
 assert q>=F(3,4)*(1-F(1,n))
 for x in (F(0),M/2,M-b/2):
  vals=[1+x*x/(M*M)+R*R/(12*M*M) for R in (a,F(3,2),b)]
  assert vals[0]<vals[1]<vals[2] and 1<vals[2]<=2
  assert ell+b/2==M-b/2
  rows.append({'stage':stage,'n':n,'M':str(M),'equal_coordinate_x':str(x),'A_at_1_1p5_2':list(map(str,vals)),'inner_source_mass_fraction_exact':str(q),'Gamma_energy_lower_per_v_exact':str(q*q),'universal_lower_per_v':str(F(9,16)*(1-F(1,n))**2)})
res={'status':'PASS_EXACT_STRICT_WINNER_HARD_GUARD','rows':rows,'row_count':len(rows),'rounds':3,'scope':'Checks analytic quadratic-source averages and inner-source mass bound. Does not numerically integrate full Gamma_v or certify any softFIRST geom gate.','script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
out.write_text(json.dumps(res,indent=2));print({k:v for k,v in res.items() if k!='rows'})
