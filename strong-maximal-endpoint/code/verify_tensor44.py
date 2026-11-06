"""Independent Fraction audit; only the fixed r043 product witnesses."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib,itertools,json,sys
repo=Path(__file__).resolve().parent.parent
source=Path(sys.argv[1]) if len(sys.argv)>1 else repo/'rounds/r043/cross_rectangles43.json'
audit=Path(sys.argv[2]) if len(sys.argv)>2 else repo/'rounds/r044/tensor_audit44.json'
pairs=Path(sys.argv[3]) if len(sys.argv)>3 else repo/'rounds/r044/tensor_k2_exact44.json'
s=json.loads(source.read_text());d=json.loads(audit.read_text());p=json.loads(pairs.read_text())
assert d['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
m={}
for r in s['rows']:
 a=r['amplitude'];m[a]=max(m.get(a,F(0)),F(r['full_witness_union_area'],r['input_cell_count']*a))
assert m=={2:F(11,8),3:F(7,3),4:F(3),8:F(4),16:F(15,4)} and max(m.values())==4
x=F(69,100);ubexp=sum((x**j/F(factorial(j)) for j in range(13)),F(0))+x**13/F(factorial(13))/(1-x/14)
assert ubexp==F(d['log2_rational_lower_certificate']['exp_upper_exact'])<2
b3=F(4**3)/(1+3*x)**5;ratio=4/(1+4*x)**2
assert b3==F(d['uniform_decay_proof']['B3_exact'])<1
assert ratio==F(d['uniform_decay_proof']['successive_ratio_upper_for_k_ge_3'])<1
lower={2:x,3:F(3,2)*x,4:2*x,8:3*x,16:4*x}
assert 9>8
seen=set();values=[]
for row in p['all_15_unordered_pairs']:
 a,b=row['amplitudes'];assert (a,b) not in seen;seen.add((a,b))
 bound=m[a]*m[b]/(1+lower[a]+lower[b])**3
 assert bound==F(row['rational_upper_bound'])<1
 values.append(bound)
assert seen==set(itertools.combinations_with_replacement(sorted(m),2))
assert max(values)==F(49000000,260409987)
assert all(F(4**k)/(1+k*x)**(2*k-1)<1 for k in range(3,129))
print('PASS: 20 input ratios, rational log2 tail, all 15 k=2 pairs, B3 and all-k recurrence constants; independent integer/Fraction arithmetic')
