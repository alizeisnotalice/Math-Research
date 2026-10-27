from fractions import Fraction as F
from pathlib import Path
import json,hashlib,os
gamma59=F(1637,12500)
gamma56=F(172531397335235056462877787662549380239,1643616601687301779182642153860479582208)
checks={"full_scalar_slack":F(429,2000)-F(3,47)-F(3,20)==F(63,94000),"old_Phi59_whole_cost_lower_envelope":gamma59+F(4,2**20)<F(3,20),"old_Phi56_whole_cost_lower_envelope":gamma56+F(4,2**20)<F(3,20)}
assert all(checks.values())
d={"status":"PASS","pid":os.getpid(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"checks":checks,"scope":"3 new entire-root-cost envelope gates, with original analytic B lower envelopes and each bound own valid domain; no author128 or oldgrid replay","domain":"R061 unknownN; old explicit domains retained; comparisons on valid intersections only"}
Path(__file__).with_name('receipt.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d))
