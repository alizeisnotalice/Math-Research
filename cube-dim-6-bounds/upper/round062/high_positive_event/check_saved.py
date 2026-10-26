#!/usr/bin/env python3
"""Pure Fraction interval replay of fixed formulas, independent integer-root bounds."""
import sys
sys.dont_write_bytecode=True
sys.set_int_max_str_digits(0)
from pathlib import Path
from fractions import Fraction as Q
from math import comb
import json,hashlib,time,os
H=Path(__file__).resolve().parent
SCALE=10**60
ROOT_CHECKS=0
def rd(p):return json.loads(Path(p).read_text())
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def saved_bounds(r):return Q(int(r['lower_mantissa']))*Q(2)**int(r['lower_binary_exponent']),Q(int(r['upper_mantissa']))*Q(2)**int(r['upper_binary_exponent'])
def iroot(n,k):
    assert n>=0 and k>=1
    if n<2:return n
    lo=0;hi=1<<((n.bit_length()+k-1)//k)
    while hi-lo>1:
        mid=(lo+hi)//2
        if mid**k<=n:lo=mid
        else:hi=mid
    assert lo**k<=n<(lo+1)**k
    return lo
def root_bounds(x,n):
    global ROOT_CHECKS
    assert x>=0
    v=iroot((x.numerator*SCALE**n)//x.denominator,n)
    l,u=Q(v,SCALE),Q(v+1,SCALE)
    assert l**n<=x<=u**n
    ROOT_CHECKS+=1
    return l,u
class I:
    def __init__(self,l,u=None):
        self.l=Q(l);self.u=Q(l if u is None else u);assert self.l<=self.u
    def __add__(self,o):
        o=cv(o);return I(self.l+o.l,self.u+o.u)
    __radd__=__add__
    def __neg__(self):return I(-self.u,-self.l)
    def __sub__(self,o):return self+-cv(o)
    def __rsub__(self,o):return cv(o)+-self
    def __mul__(self,o):
        o=cv(o);v=[self.l*o.l,self.l*o.u,self.u*o.l,self.u*o.u];return I(min(v),max(v))
    __rmul__=__mul__
    def __truediv__(self,o):
        o=cv(o);assert o.l>0;return self*I(1/o.u,1/o.l)
    def __rtruediv__(self,o):return cv(o)/self
    def power(self,m,n=1):
        assert self.l>0 and n>=1
        if m<0:return I(1)/self.power(-m,n)
        l,u=self.l**m,self.u**m
        if n==1:return I(l,u)
        return I(root_bounds(l,n)[0],root_bounds(u,n)[1])
def cv(x):return x if isinstance(x,I)else I(x)
def atan_scalar_bounds(x,n):
    assert 0<=x<=Q(1,2)and n%2==0
    lower=sum(((-1)**k*x**(2*k+1)/Q(2*k+1)for k in range(n)),Q(0))
    return lower,lower+x**(2*n+1)/Q(2*n+1)
def atan_interval(x,n):return I(atan_scalar_bounds(x.l,n)[0],atan_scalar_bounds(x.u,n)[1])
def compact(x):
    scale=10**30
    l=Q((x.l*scale).__floor__(),scale);u=Q((x.u*scale).__ceil__(),scale)
    return dict(lower=str(l),upper=str(u))
def main():
    tic=time.monotonic();assert not(H/'checked_results.json').exists(),'one pure check only'
    reg=rd(H/'registration.json');release=rd(H/'main_release.json');binding=rd(H/'execution_binding.json');receipt=rd(H/'run_receipt.json');r=rd(H/'results.json')
    guard=rd(H/'checker_process_receipt.json');assert guard['status']=='RUNNING'and guard['owned_pid']==os.getpid()and guard['wall_timeout_seconds']==180
    pg=rd(H/'producer_process_receipt.json');assert pg['status']=='COMPLETED'and pg['exit_code']==0 and pg['elapsed_seconds']<180 and pg['child_reaped']
    assert receipt['status']=='COMPLETED'and receipt['attempted_configurations']==receipt['completed_configurations']==9 and receipt['exit_code']==0
    assert receipt['registration_sha256']==release['registration_sha256']==binding['registration_sha256']==sha(H/'registration.json')
    assert receipt['main_release_sha256']==binding['main_release_sha256']==sha(H/'main_release.json')and receipt['binding_sha256']==sha(H/'execution_binding.json')and receipt['results_sha256']==sha(H/'results.json')
    for f,x in reg['source_sha256'].items():assert sha(f)==x,f
    for f,x in binding['code_sha256'].items():assert sha(H/f)==x,f
    assert r['status']=='COMPLETE_NEW_HIGH_POSITIVE_EVENT_FULL_TIME_UPPER'
    assert [x['id']for x in r['configurations']]==list(range(1,10))
    assert r['configurations'][1]['tail_derivative_residual']==r['configurations'][1]['loose_optimizer_residual']=='0'
    assert r['configurations'][2]['partial_fraction_residual']==r['configurations'][2]['double_root_residual']==r['configurations'][2]['triple_root_residual']=='0'
    assert [x['n']for x in r['configurations'][3]['ledger']]==[1,5,10,10,5]
    aa=atan_scalar_bounds(Q(1,5),120);bb=atan_scalar_bounds(Q(1,239),40)
    pi=I(16*aa[0]-4*bb[1],16*aa[1]-4*bb[0]);assert 3<pi.l<pi.u<4
    rho=I(Q(1,64));b=I(Q(1,4));d=I(Q(31,64));r1=Q(reg['fixed_R1']);crem_exact=Q(44743079,2970)
    shape=sum((comb(5,k-1)*I(Q(6-k,k)).power(1,2)for k in range(1,6)),I(0))
    V=pi.power(2)/8*rho.power(4)*shape;D=I(r1)/(32*pi.power(6));C=D/b.power(2);Crem=I(crem_exact)/pi.power(6)
    eps=(C*pi*V).power(1,5);Vp=V.power(5,6)
    Pshort=Q(6,7)*C.power(1,6)*pi.power(-5,6)*eps.power(7,6)
    ratio=eps/b;assert 0<ratio.l<ratio.u<Q(1,2)
    # atan(b/eps)=pi/2-atan(eps/b), with interval monotonicity and alternating remainder.
    atan_large=pi/2-atan_interval(ratio,96)
    T=1/(3*b.power(2)*eps.power(3))-1/(b.power(4)*eps)+atan_large/b.power(5)
    assert T.l>0 and T.u<(1/(3*b.power(2)*eps.power(3))).u
    Plong=Vp*D*T
    alpha=[(I(z)*d.power(2)/a).power(1,2)for z,a in zip([Q(6),Q(10,3),Q(2)],[992,960,928])]
    x,y,z=alpha;JB=pi*(x+y+z)/(2*x*y*z*(x+y)*(x+z)*(y+z)*(992*960*928));Bfull=Vp*Crem*JB
    total=Pshort+Plong+Bfull
    final=Q(r['union']['simple_rational_upper']);assert total.u<final
    target=Q(4,48750)
    comparison='BELOW'if final<target else('ABOVE'if total.l>target else'UNRESOLVED')
    if r['union']['comparison']=='UPPER_FORMULA_BELOW_CONDITIONAL_TARGET':assert comparison=='BELOW'
    if r['union']['comparison']=='UPPER_FORMULA_EXCEEDS_CONDITIONAL_TARGET':assert comparison=='ABOVE'
    # Independently certified intervals must overlap the authored saved formula enclosures.
    for name,value in [('Pshort',Pshort),('Plong',Plong),('Bfull',Bfull)]:
        l,u=saved_bounds(r[name]);assert value.l<=u and l<=value.u
    assert r['union']['all_parameter_events_paid']and r['union']['all_positive_time_paid']and r['union']['one_whole_union_volume']and not r['union']['whole_face_paid']
    result=dict(status='PASS_PURE_RATIONAL_FIXED_HIGH_EVENT_AUDIT',results_sha256=sha(H/'results.json'),registration_sha256=sha(H/'registration.json'),main_release_sha256=sha(H/'main_release.json'),binding_sha256=sha(H/'execution_binding.json'),producer_configurations=9,pure_saved_data_configurations=1,total_configurations=10,independent_root_endpoint_power_checks=ROOT_CHECKS,root_scale_decimal_digits=60,Machin_terms=[120,40],atan_small_argument_terms=96,atan_full_tail_used=True,three_positive_quadratics_full_axis_used=True,certified_rational_intervals={k:compact(v)for k,v in dict(V=V,epsilon=eps,Pshort=Pshort,Plong=Plong,Bfull=Bfull,total=total).items()},strict_rational_high_upper=str(final),conditional_whole_face_target=str(target),comparison=comparison,all_source_hashes_unchanged=True,old_M_recomputations=0,old_cap_budget_reruns=0,old060_reruns=0,new_Arb_calls=0,new_kernel_calls=0,full_enumerations=0,matlab_calls=0,whole_face_paid=False,new_C6_bound=False,elapsed_seconds=time.monotonic()-tic)
    (H/'checked_results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
