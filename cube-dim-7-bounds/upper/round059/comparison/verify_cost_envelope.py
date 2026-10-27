from pathlib import Path
from fractions import Fraction as F
import json,hashlib,os
gamma56=F(172531397335235056462877787662549380239,1643616601687301779182642153860479582208)
gamma59=F(1637,12500)
checks={"chosen_new_scalar_margin":gamma59-F(3,97)-F(1,10)==F(39,1212500),"old_complete_cost_lower_envelope_gap":gamma56-F(7,100)+F(4,2**20)<F(1,10),"old_gamma_and_root_positive_regime":0<gamma56<1 and 0<gamma59<1}
assert all(checks.values())
d={"status":"PASS","pid":os.getpid(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"checks":checks,"scope":"3 new exact full-root comparison gates; independent analytic envelopes and each bound own actual operator domain required; no old tuple/grid replay","domain":"Phi59 valid only after its unknown N59; Phi56 logn>=2^20; comparison on intersection only"}
Path(__file__).with_name('receipt.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d))
