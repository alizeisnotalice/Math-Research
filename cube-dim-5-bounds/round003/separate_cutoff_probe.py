from pathlib import Path
import sys,json,hashlib,time
from fractions import Fraction
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'numerics'))
from body_fan_arb_anchor import exact_terms
from body_fan_h_epsilon_arb import truncated_green
from flint import ctx
ctx.prec=256
start=time.monotonic();terms=tuple(exact_terms());point=(Fraction(1),)+tuple(Fraction(n,1009) for n in (936,673,68,917))
rows=[]
for d in (64,32,24,16):
 if time.monotonic()-start>120:break
 value,_=truncated_green(point,Fraction(1,d),terms)
 rows.append({'epsilon':f'1/{d}','H_ball':str(value)})
out={'status':'finite_point_diagnostic_not_global_g','point':[str(x) for x in point],'rows':rows,'seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
raw=json.dumps(out,indent=2)+'\n';Path(__file__).with_suffix('.json').write_text(raw);print(raw)
