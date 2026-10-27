from fractions import Fraction as F
from pathlib import Path
import json,hashlib,os
checks={"new_full_scalar_slack":F(61,200)-F(91,909)-F(1,5)==F(889,181800),"old_layer_full_cost_positive_Bbuffer":F(429,2000)-F(3,50)+F(4,2**20)<F(1,5),"old_wideband_full_cost_lower_envelope":F(1637,12500)+F(4,2**20)<F(1,5)}
assert all(checks.values())
d={"status":"PASS","pid":os.getpid(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"checks":checks,"scope":"3 new whole-root envelope gates; R061 positive Bbuffer retained; each original valid domain required, no Cn lower bound or old verifier replay","domain":"R064/R061 unknown N, R063 explicitL33; only actual domains intersections"}
Path(__file__).with_name('receipt.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d))
