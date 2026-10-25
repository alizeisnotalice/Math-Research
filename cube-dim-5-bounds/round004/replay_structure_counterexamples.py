"""Main-project replay of two finite spatial counterexamples; no global bound."""
from fractions import Fraction as F
from pathlib import Path
import sys,json,hashlib,time
from flint import ctx
root=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(root/'numerics'))
from body_fan_arb_anchor import exact_terms
from body_fan_h_epsilon_arb import truncated_green
ctx.prec=256;started=time.monotonic();terms=tuple(exact_terms());rows=[]
for eps,points in [(F(1,32),[(F(1),s,F(1,2),F(1,2),F(1,2)) for s in (F(0),F(1,2),F(1))]),(F(1,128),[(F(1),F(0),s,F(1),F(1)) for s in (F(0),F(3,4),F(1))])]:
 vals=[]
 for p in points:
  if time.monotonic()-started>120:raise TimeoutError('registered 120s checkpoint')
  v,patterns=truncated_green(p,eps,terms);vals.append(v.real)
 if eps==F(1,32):ok=vals[1].lower()>vals[0].upper() and vals[1].lower()>vals[2].upper();meaning='neither increasing nor decreasing in this coordinate'
 else:ok=vals[1].upper()<vals[0].lower() and vals[1].upper()<vals[2].lower();meaning='interior point strictly below both edge endpoints'
 assert ok
 rows.append({'epsilon':str(eps),'points':[[str(x) for x in p] for p in points],'H_balls':[str(v) for v in vals],'strict_inequalities_verified':bool(ok),'conclusion':meaning})
out={'status':'PASS_main_project_finite_spatial_counterexamples','precision_bits':256,'rows':rows,'seconds':time.monotonic()-started,'scope':'Rules out the stated coordinate monotonicity / endpoint-minimum principles for these H cutoffs only; no C5 bound. Shares exact cone geometry and Arb library with production.','source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
