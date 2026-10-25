from fractions import Fraction as F
from pathlib import Path
import hashlib,json
from flint import arb,ctx,fmpq
ctx.prec=256
g=F(11,100000);I=F(3,62500)
cost4=3125*I**4/g**5
diff=F(32**4)-cost4
assert diff>0
cost=arb(fmpq(cost4.numerator,cost4.denominator))**(arb(1)/4)
out={'status':'conditional_target_comparison_only_not_achieved_bound','g_target':str(g),'I_target':str(I),'cost_fourth_power':str(cost4),'32_fourth_minus_target_cost_fourth':str(diff),'cost_display_ball':str(cost),'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
raw=json.dumps(out,indent=2)+'\n';Path(__file__).with_suffix('.json').write_text(raw);print(raw)
