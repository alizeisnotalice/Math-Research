from pathlib import Path
from fractions import Fraction as F
import sys,json,hashlib,time
from flint import arb,fmpq,ctx
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'numerics'))
from body_fan_arb_anchor import exact_terms
from body_fan_h_epsilon_arb import truncated_green
ctx.prec=256
terms=tuple(exact_terms());start=time.monotonic();rows=[]
cases=[(F(1,24),(F(1),F(223,4096),F(1021,1024),F(2937,4096),F(3639,4096))),(F(1,32),(F(1),F(0),F(3,4),F(1),F(1)))]
target=arb(fmpq(11,100000))
for eps,x in cases:
 assert time.monotonic()-start<120
 value,patterns=truncated_green(x,eps,terms);value=value.real
 assert value.upper()<target.lower()
 rows.append({'epsilon':str(eps),'point':[str(a) for a in x],'H_ball':str(value),'target_minus_H_ball':str(target-value),'strictly_below_target':True,'multiplicity_patterns':{str(k):v for k,v in patterns.items()}})
out={'status':'PASS_strict_point_counterexamples_to_selected_truncated_global_targets','target':'11/100000','scope':'Continuity of H_epsilon makes each strict point inequality fail the same essential face lower bound on a neighborhood. Does not exclude a G lower bound.','rows':rows,'seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
raw=json.dumps(out,indent=2)+'\n';Path(__file__).with_suffix('.json').write_text(raw);print(raw)
