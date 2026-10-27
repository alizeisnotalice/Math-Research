from fractions import Fraction as F
from pathlib import Path
import json,hashlib,os
root=Path(__file__).resolve().parents[3]
price_receipt=json.loads((root/'audit/round066/exact_scalar/verification.json').read_text())
S7=F(price_receipt['log_upper_S7']);P=F(4333,50000)
checks={"reads_original_exact_price_record":S7<P,"same_Phi_exact_price_slack":F(317,1000)-P-F(23,100)==F(17,50000),"old_radial_complete_cost_positive_buffer":F(61,200)-F(91,1000)+F(4,2**20)<F(9,40)<F(23,100),"old_layer_complete_cost_positive_buffer":F(429,2000)-F(3,50)+F(4,2**20)<F(23,100)}
assert all(checks.values())
d={"status":"PASS","pid":os.getpid(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"checks":checks,"scope":"reads original S7 field without recomputing series, plus3 new exact full-root-cost gates; samePhi price does not change kernel/gamma/fullF","domains":"all original operator valid domains retained; R066 unknownN; no Cn lowerbound"}
Path(__file__).with_name('receipt.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps(d))
