"""New compressed proof-margin checks, independent of the author's verifier."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import hashlib, json, os

L = Q(8192); q0 = Q(2,3); x = Q(91,100)*L
N = L**32/factorial(32)  # Deliberately not the author's N64.
h = 1/(256*L); rhi=Q(273,200)-Q(1,2**19)
ea=Q(1,65536); eb=Q(1,2**19); ej=Q(1,2**20)

# Independent alternating log(1+1/5), even/odd truncations 8 and 9.
def alt(m):
    return sum((Q(1,5)**j/j if j%2 else -Q(1,5)**j/j)
               for j in range(1,m+1))
ql, qu = alt(8), alt(9)
cl = 6-Q(5,2)*(qu*qu+2*qu+2)
cu = 6-Q(5,2)*(ql*ql+2*ql+2)
lo=(1-1/N)*(cl-(cu+Q(1,25))*h)
hi=cu+1/(N*q0)
ml=Q(3,4)*L*L/(q0*q0*N)+L/(2*q0*N)
vl=Q(3,2)*L*L/(q0*q0*N)

# A shorter positive exponential prefix than the author's S64.
def P40(y):
    return sum(y**j/factorial(j) for j in range(41))
c=Q(81,20000); y=c*L
Csrc=(Q(91,100)-q0)/q0
cap=2*Csrc*rhi*256*L/P40(y)
tau=rhi*256*L*91/P40(y)

# Independent positive log-price series with an explicit geometric remainder.
v=Q(9,100)
price=sum(v**j/j for j in range(1,9))+v**9/(9*(1-v))
scalar=1/N+1/P40(v*L)+182/P40(v*L/2)
checks={
 "first_profile_mean_and_pin": Q(1,2)-Q(5,2)*ql<Q(221,5000) and N>9000,
 "second_profile_mean_corridor": lo>Q(5275,1000000) and hi<Q(5295,1000000),
 "baseline_log_probability": ml<eb/2 and 4*vl/eb**2<ea/8,
 "count_and_both_profiles_probability": max(4/(eb**2*N),Q(324000000)/N,
                  1/(N*q0*q0*Q(1,200000)**2))<ea/8,
 "target_cap_probability": Q(3,4)**74<ea/8 and N>75,
 "finite_n_cap_correction": Q(9,182)*x/N<Q(1,2),
 "fullprefix_cap_loss": cap<ej,
 "fullprefix_tau_loss": tau<eb,
 "all_tail_prefix_derivative": y>4 and -Q(3,2)-y/2+y*y/4>0,
 "complete_scalar_price_and_buffers": price<Q(189,2000) and scalar<Q(1,2000),
 "source_variance_profile": Q(5301,1000000)/(1-ej)<Q(531,100000),
 "lower_budget_positive": L/4>Q(13,100) and (L/2)**2/2>L,
}
assert all(checks.values()), checks
out={"status":"PASS","pid":os.getpid(),"checks":checks,
     "new_only":"No preflight or author verifier was rerun; N32, alternating log8/9, positive log-price8+tail, exponential P40.",
     "cap_bound":str(cap),"tau_bound":str(tau),"price_upper":str(price),
     "source_sha256":hashlib.sha256(Path('/Users/zhengzhihao/Desktop/T/cube-r030-geometry/cube-dim-7-bounds/upper/round071/second_profile_cantelli/PROOF.md').read_bytes()).hexdigest(),
     "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
assert out['source_sha256']=='0e36488654bfa1413c19335a3eb9cf7b900e16b25cf9b73b2b8f36efc7381a10'
Path(__file__).with_name('critical_tail.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({"status":"PASS","pid":out['pid'],"groups":len(checks)}))
