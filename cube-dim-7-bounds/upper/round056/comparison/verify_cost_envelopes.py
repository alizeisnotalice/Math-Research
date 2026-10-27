from fractions import Fraction as F
from pathlib import Path
import json,hashlib,os
L0=2**20
gs=F(9775559871,281956214374400)
gr=F(735063024285513,21076775872652902400)
checks={
 "simple_gamma_below_one_thousandth":gs<F(1,1000),
 "raw_gamma_below_one_thousandth":gr<F(1,1000),
 "old_simple_domain_integral_loss":F(4,2**27)<F(1,1000),
 "combined_old_explicit_loss":F(1,500)<F(1,40),
 "sqrt80_below9":80<81,
 "old_sqrtlog_envelope_base":F(9,1024)+F(4,L0)<F(1,40),
 "old_R051_family_gamma_envelope":F(1,256)+F(4,L0)<F(1,40),
 "R042_loss_below_one_thousandth":F(1,2**31)<F(1,1000),
 "R049_loss_below_one_thousandth":F(1,2**40)<F(1,1000)}
assert all(checks.values())
result={"status":"PASS","pid":os.getpid(),"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"checks":checks,"scope":"new exact full-cost-envelope comparisons; all original operator domains retained; no old baseline or production rerun","infinite_tail":"4/L and sqrt(logL/L) decrease for L>=2^20; log(2^20)<20; sqrt80<9; analytic lower envelopes required, not inferred from these arithmetic checks"}
Path(__file__).with_name('receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
