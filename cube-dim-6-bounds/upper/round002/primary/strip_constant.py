"""Finite exact sums supporting the convexity proof of the sharp face tube constant."""
from itertools import product
from fractions import Fraction as Q
from pathlib import Path
import json,hashlib
v=[(1,)+s for s in product((-1,1),repeat=5)]
def psi(x):return sum(abs(sum(a*b for a,b in zip(w,x))) for w in v)
free=psi([0]+[Q(1,5)]*5)
full=psi([Q(1,6)]*6)
assert free==12 and full==10
quotient_vectors=[tuple(a-1 for a in w[1:]) for w in v if w!=(1,1,1,1,1,1)]
qpoint=[Q(1,5)]*3+[-Q(1,5)]*2
quotient_value=sum(abs(sum(a*b for a,b in zip(w,qpoint))) for w in quotient_vectors)
assert quotient_value==12
cf=Q(32**3,12**5)
assert cf==Q(32,243)
assert Q(1024,5)/cf==Q(7776,5)
out={'status':'PASS finite exact identities; global inequality uses convexity proof','face_l1_sharp_constant':str(free),'full_l1_sharp_constant':str(full),'quotient_l1_sharp_constant':str(quotient_value),'new_peak_coefficient_without_pi6':'1/1000000','new_face_sup_coefficient_without_pi6':str(cf),'old_over_new_face_sup':'7776/5','face_equality_vector':['0']+['1/5']*5,'full_equality_vector':['1/6']*6,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('strip_constant_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
