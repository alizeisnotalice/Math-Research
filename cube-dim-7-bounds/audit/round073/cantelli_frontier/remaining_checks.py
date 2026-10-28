"""Only R073 proof interfaces not exercised by the seven-stage preflight."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib, json, os

SRC=Path('/Users/zhengzhihao/Desktop/T/cube-r040-cover-selector/cube-dim-7-bounds/upper/round073/cantelli_frontier/PROOF.md')
digest=hashlib.sha256(SRC.read_bytes()).hexdigest()
assert digest=='50d48b85f02a14f50614f0eea5cedbd7ccd68e9b2f7e1564fe8b41fcd2af9c4b'
L=Q(8192); x=Q(231,256)*L; h=1/(768*L); q0=Q(2,3); q1=Q(11569,12800)
N=L**32/factorial(32); delta=Q(1,2**19); rhi=Q(3,2)*q1-delta
eb=Q(1,2**19); ej=Q(1,2**20); ea=Q(1,65536)
ml=Q(3,4)*L**2/(q0*q0*N)+L/(2*q0*N)
vl=Q(3,2)*L**2/(q0*q0*N)

def altlog(m):
    return sum((Q(1,5)**j/j if j%2 else -Q(1,5)**j/j) for j in range(1,m+1))
ql,qu=altlog(8),altlog(9)
mu2lo=1-5*qu-Q(5,2)*qu*qu
mu2hi=1-5*ql-Q(5,2)*ql*ql
E2lo=(1-1/N)*(mu2lo-h*(mu2hi+Q(1,25)))
E2hi=mu2hi+1/(250*N*q0)

# Far shorter scalar exponential prefix, sufficient for the untested buffer.
def P16(y):
    return sum(y**j/factorial(j) for j in range(17))
sigma=Q(25,256)
remainder=1/N+1/P16(sigma*L)+182/P16(sigma*L/2)
checks={
 'new_tuple_all_local_geometry':N>9000 and x/N<Q(1,1000) and N*q0/x>4*rhi
     and N/x>rhi and x/N<1 and q1/(x*N*N)<2*h/3
     and x/(3*q0*N)<ej,
 'baseline_and_three_profile_probability_budgets':ml<eb/2 and max(
     4*vl/eb**2,4/(eb**2*N),Q(324000000)/N,Q(5760000)/N)<ea/8,
 'target_cap_as_new_rhi':1+x/(2*N)-1/(2*rhi)<Q(3,4)
     and Q(3,4)**74<ea/8,
 'actual_h_baseline_H2_mean':E2lo>Q(528,100000) and E2hi<Q(5295,1000000),
 'actual_h_wholematching_profiles':h<Q(1,100000) and h/5<Q(1,1000000),
 'complete_root_positive_buffers':remainder<Q(1,5000) and sigma*L/2>4,
}
assert all(checks.values()),checks
out={'status':'PASS','pid':os.getpid(),'groups':checks,'source_sha256':digest,
     'new_only':'Preflight, author53 and old programs were not rerun. N32, alternating log8/9 and scalar P16.',
     'root_remainder_bound':str(remainder),
     'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
Path(__file__).with_name('remaining_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'status':'PASS','pid':out['pid'],'groups':len(checks)}))
