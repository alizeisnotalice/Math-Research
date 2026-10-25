from flint import arb,ctx
from pathlib import Path
import hashlib,json
ctx.prec=256
P=arb(35509)/(764411904*arb.pi()**5)
ceiling=arb(2)/(25*(2*arb.pi())**2)*P**(arb(1)/5)
target=arb(11)/100000
assert ceiling<target
out={'status':'strict_exclusion_of_trace_Jensen_g_target_only','peak_formula':'35509/(764411904*pi^5)','ceiling_formula':'2/(25*(2*pi)^2)*P^(1/5)','ceiling_ball':str(ceiling),'target':'11/100000','target_minus_ceiling_ball':str(target-ceiling),'strict_comparison':bool(ceiling<target),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
