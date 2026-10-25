from fractions import Fraction as F
from pathlib import Path
import json,hashlib
from flint import arb,fmpq,ctx
ctx.prec=256
g=F(1,10000);I=F(11,1000000)
cost7=F(10,3)**10*I**7/g**10
diff=F(32**7)-cost7
assert diff>0
cost=arb(fmpq(cost7.numerator,cost7.denominator))**(arb(1)/7)
out={'status':'conditional_targets_only_not_achieved','alpha':'3/2','g_target':str(g),'I_target':str(I),'cost_seventh_power':str(cost7),'32_seventh_minus_cost_seventh':str(diff),'cost_ball':str(cost),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
