"""Exact rational eligibility guard for the proposed small-softness branches."""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt,log1p,exp
import hashlib,json
rows=[]
for n in (1024,4096,16384):
 for label,sigma in [('inverse_square',F(1,n*n)),('inverse_three_halves',F(1,n*isqrt(n))),('inverse_dimension',F(1,n))]:
  # phi(1/2)>3/8 gives b_sigma>1-5*sigma/8.
  b=1-F(5,8)*sigma
  left=2*pow(b.numerator,n); right=pow(b.denominator,n)
  assert left>right
  rows.append(dict(n=n,branch=label,sigma=str(sigma),exact_two_b_power_exceeds_one=True,numerator_bits=left.bit_length(),denominator_bits=right.bit_length(),illustrative_b_power=exp(n*log1p(-float(F(5,8)*sigma)))))
p=Path(__file__); out=p.with_name(p.stem+'_results.json')
assert not out.exists()
out.write_text(json.dumps(dict(status='PASS_EXACT_EMPTY_BRANCH_GUARD',rows=rows,script_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),scope='Given the original future ceiling, hard band, and proved phi(1/2)>3/8: these softness values are ineligible. Floating values are illustrative only; comparison uses integers.'),indent=2)+'\n')
print('9 exact eligibility checks passed; small-softness branch is empty under original hypotheses.')
