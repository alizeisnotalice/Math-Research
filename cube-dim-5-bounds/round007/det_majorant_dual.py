from pathlib import Path
from fractions import Fraction as F
import json,hashlib,time
p=Path(__file__).resolve().parent;start=time.monotonic();m=F(648,1463);s=F(126404864,474619145);delta=F(1,4000);t=(s-m*delta)/(m-delta);prob=(m-delta)/(t-delta)
assert 0<delta<t<1 and 0<prob<1
assert (1-prob)*delta+prob*t==m and (1-prob)*delta**2+prob*t**2==s
D=10**40
project=p.parent
qfile=project/'upper/round011_shape_integral_results.json'
jfile=project/'round006/orbit_minkowski.json'
qball=json.loads(qfile.read_text())['arb_256bit']['negative_determinant_moment_upper_Q']
jball=json.loads(jfile.read_text())['J_upper_ball']
def rational_ball(text):
 mid,rad=text.strip('[]').split(' +/- ');return F(mid)-F(rad),F(mid)+F(rad)
assert rational_ball(qball)[0]>F(12,5) and rational_ball(jball)[1]<F(17,10)


def nthroot(n,k):
 lo=0;hi=1<<((n.bit_length()+k-1)//k)
 while lo+1<hi:
  q=(lo+hi)//2
  if q**k<=n:lo=q
  else:hi=q
 return hi if hi**k==n else lo

def root_bound(x,k):
 n=nthroot(x.numerator*D**k//x.denominator,k);lo=F(n,D);hi=F(n+1,D)
 assert lo**k<=x<hi**k
 return lo,hi
intervals={}
for k in [14,20]:
 a,b=root_bound(delta**-3,k);c,d=root_bound(t**-3,k);L=(1-prob)*a+prob*c;H=(1-prob)*b+prob*d;intervals[str(k)]=[str(L),str(H)]
 if k==14:assert H<F(12,5)
 else:assert L>F(17,10)
out={'status':'PASS_exact_two_point_obstruction_for_selected_majorant_family','family':'d^-3/20 <= alpha*d^-3/14 + A+B*d+C*d^2, alpha>=0','support':[str(delta),str(t)],'probability_high':str(prob),'E_d':str(m),'E_d_squared':str(s),'negative_moment_intervals_by_root_order':intervals,'objective_of_any_family_certificate_strictly_above':'17/10','scope':'abstract scalar distribution, not claimed realizable as detM; excludes only this pointwise-majorant and moment-budget family','elapsed_seconds':time.monotonic()-start,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'budget_source_hashes':{str(x.relative_to(project)):hashlib.sha256(x.read_bytes()).hexdigest() for x in [qfile,jfile]}}
(p/'det_majorant_dual.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],out['elapsed_seconds'])
