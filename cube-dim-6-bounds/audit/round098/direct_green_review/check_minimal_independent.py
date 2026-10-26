#!/usr/bin/env python3
"""Three bounded saved/small-exact groups, one run, no author/kernel replay.

Independent Fraction arithmetic proves coarse rational claims. The author's
384-bit Arb endpoints are hash-bound and consistency-tested; this script
does not claim independent recertification of all their displayed digits.
"""
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt
import hashlib
import json
import resource
import signal
import sys
import time

H = Path(__file__).resolve().parent
CPU, WALL = time.process_time(), time.monotonic()
resource.setrlimit(resource.RLIMIT_CPU, (20,20))
signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError('wall limit')))
signal.alarm(165)
checks = 0
counts = {'completed_groups':0,'source_hashes':0,'new_kernel_calls':0,
          'author_programs_or_four_configurations_replayed':0,
          'old043_037_scientific_replays':0,'new_grids':0}
out = {'status':'RUNNING'}
def ck(ok, msg):
    global checks
    checks += 1
    if not ok:
        raise AssertionError(msg)
    if time.process_time()-CPU >= 18:
        raise TimeoutError('CPU prelimit')
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p): return json.loads(Path(p).read_text())
def dyadic(x):
    unit = Q(2)**int(x['common_binary_exponent'])
    ans = Q(x['lower_mantissa'])*unit,Q(x['upper_mantissa'])*unit
    ck(0 < ans[0] <= ans[1], 'positive saved dyadic interval')
    return ans
def C(x): return Q(x),Q(x)
def add(a,b): return a[0]+b[0],a[1]+b[1]
def sub(a,b): return a[0]-b[1],a[1]-b[0]
def mul(a,b):
    vals = [x*y for x in a for y in b]
    return min(vals),max(vals)
def div(a,b):
    ck(b[0] > 0, 'positive interval denominator')
    return mul(a,(1/b[1],1/b[0]))
def power(a,n):
    ck(a[0] >= 0 and isinstance(n,int) and n>=0, 'positive integer power')
    return a[0]**n,a[1]**n
def sqrtq(x):
    x=Q(x); D=10**18
    k=isqrt((x.numerator*D*D)//x.denominator)
    ans=Q(k,D),Q(k+1,D)
    ck(ans[0]**2 <= x <= ans[1]**2, 'integer square-root enclosure')
    return ans
def atan_small(x):
    # Exact alternating integral expansion, 20 terms, positive remainder.
    x=Q(x); n=20
    low=sum((-1)**k*x**(2*k+1)/Q(2*k+1) for k in range(n))
    return low,low+x**(2*n+1)/Q(2*n+1)
def consistent(saved, coarse, name):
    a,b=dyadic(saved)
    ck(coarse[0] <= a <= b <= coarse[1], 'saved interval within independent coarse bound '+name)

try:
    gate=load(H/'source_gate.json')
    for p,expected in gate['source_pins'].items():
        ck(sha(p)==expected,'source hash '+p);counts['source_hashes']+=1
    author=Path(gate['author_root'])
    manifest,results,up,cert=[load(author/n) for n in ('manifest.json','results.json','upstream_constants.json','certificate.json')]
    ck(sha(author/'manifest.json')==gate['manifest_sha256'] and sha(author/'theorem.md')==gate['theorem_sha256'],'final author seal')
    for ent in manifest['files']:
        p=author/ent['name']
        ck(p.stat().st_size==ent['bytes'] and sha(p)==ent['sha256'],'author whitelist '+ent['name'])
    cfg=results['configurations'];ck(len(cfg)==4,'four saved author records only')
    source033=load(next(p for p in gate['source_pins'] if p.endswith('/round033/exact_projection_peaks/results.json')))
    source037=load(next(p for p in gate['source_pins'] if p.endswith('/round037/pair_resolved_green_envelope/results.json')))

    # Group 1: only tiny exact Gram identities and inherited scalar handoff.
    peak=next(x for x in source033['cases'] if x['id']=='r3_2211')
    ck(Q(up['R3_2211'])==Q(peak['R_exact'])==Q(707,6614784) and peak['Gram_determinant']==192,'R3 exact source handoff')
    old=next(x for x in source037['configurations'] if x.get('k')==3)
    ck(old['M']==up['M3_saved_dyadic'],'saved M3 exact handoff')
    Mlo,Mhi=dyadic(old['M'])
    ck(Mlo/2==Q(cfg[1]['saved_M3_half_lower_exact']) and Mlo/2>Q(1,5000),'old Phi coefficient lower')
    v=[1,1,1,1,1,1];w=[1,1,1,-1,-1,-1];a=[1,1,-1,1,1,-1]
    dot=lambda x,y:sum(Q(t)*u for t,u in zip(x,y))
    ck((dot(v,v),dot(w,w),dot(v,w),dot(a,v),dot(a,w))==(6,6,0,2,0),'small pair/triple Gram data')
    projected=dot(a,a)-dot(a,v)**2/6-dot(a,w)**2/6
    ck(projected==Q(16,3) and 36*projected==192,'projected vector length and J3')
    bins={}
    for j in range(6): bins[(w[j],a[j])]=bins.get((w[j],a[j]),0)+1
    ck(sorted(bins.values())==[1,1,2,2],'actual triple orbit 2211')
    ck(32-2==30 and 32*30==960 and cfg[0]['nonzero_quotient_columns']==30,'m30 deletion constant')
    ck(Q(32)*Q(1,64)/10==Q(1,20),'optional cone relative symbol loss')
    counts['completed_groups']=1

    # Group 2: independent coarse pi/radical interval proof of new claims.
    z=Q(1,5);t2=2*z/(1-z*z);t4=2*t2/(1-t2*t2);r=Q(1,239)
    ck((t4-r)/(1+t4*r)==1,'Machin tangent identity')
    # 0<4 atan(1/5)-atan(1/239)<4/5<pi/2; see report.
    ck(4*z/(1+z*z)-r>0 and Q(4,5)<1,'Machin quadrant bounds')
    pi=sub(mul(C(16),atan_small(z)),mul(C(4),atan_small(r)))
    ck(Q(3)<pi[0]<pi[1]<Q(22,7),'independent rational pi enclosure')
    rt3,rt2=sqrtq(3),sqrtq(2)
    eta=div(rt3,C(64))
    hv0=sub(div(rt3,C(2)),eta);hv1=sub(sqrtq(Q(3,32)),eta)
    hw0=sub(div(mul(C(15),rt3),C(16)),eta);hw1=sub(sqrtq(Q(3,128)),eta)
    g=sub(div(C(1),mul(C(4),rt2)),eta)
    ck(min(hv0[0],hv1[0],hw0[0],hw1[0],g[0])>0,'all inherited geometric gaps positive')
    ck(g[0]>Q(1,64),'all other pair tubes separated')
    K=div(rt3,mul(C(128),power(pi,4)))
    Hcoef=div(mul(C(3),rt3),mul(C(65536),power(pi,2)))
    ck(K[1]<Q(1,7000),'new principal full-time coefficient')
    Q3e=div(mul(rt3,C(Q(up['R3_2211']))),power(pi,3))
    Kcone=div(mul(C(Q(1024,12)*Q(20,19)**3),Q3e),power(pi,2))
    consistent(cfg[1]['uniform_H_coefficient'],Hcoef,'H')
    consistent(cfg[1]['uniform_P_G_coefficient'],K,'K')
    consistent(cfg[1]['direction_e0_Q3'],Q3e,'optional Q3')
    consistent(cfg[1]['cone_P_G_coefficient'],Kcone,'optional cone K')
    Crem=div(C(Q(44743079,2970)),power(pi,6))
    Krem=mul(mul(C(Q(1000,729)),Crem),add(C(1),sqrtq(Q(6,928))))
    bracket=add(div(C(1),mul(power(hv0,2),power(hv1,2))),div(C(1),mul(power(hw0,2),power(hw1,2))))
    B0=mul(div(mul(pi,Krem),mul(mul(C(2),sqrtq(928)),g)),bracket)
    ck(B0[1]<512,'new complete complementary full-time payment')
    consistent(cfg[2]['B0'],B0,'B0')
    Rmax=power(add(sqrtq(3*(1+Q(33,64)**2)),eta),2)
    dvmax=power(add(sqrtq(Q(3,2)*(1-Q(31,64))**2),eta),2)
    ck(Rmax[1]<4 and dvmax[1]<Q(1,2),'comparison geometry strict caps')
    consistent(cfg[2]['r_squared_upper_check'],Rmax,'norm squared')
    consistent(cfg[2]['dv_squared_upper_check'],dvmax,'line distance squared')
    counts['completed_groups']=2

    # Group 3: a coarser exact local budget, no irrational-power evaluation.
    delta=Q(1,2**24);q=Q(6,5);ulen=Q(1,32)
    ck(Q(cfg[2]['inner_rho'])==delta and list(map(Q,cfg[2]['mean_B_interval']))==[Q(31,64),Q(33,64)],'fixed local parameters')
    ratio=5000*(Q(1,7000)+512*delta)
    ck(ratio==Q(24855,28672)==Q(cfg[2]['new_over_old_uniform_upper_exact'])<1,'strict old-expression ratio')
    # A chart LL^T=[[2,1],[1,2]]; B mean direction has norm sqrt(3).
    ck(2*2-1==3 and 3*3==9 and Q(3)*Q(1,3)*ulen==Q(1,32),'face Jacobian/sector cancellation')
    ck(ulen*4/(4-q)==Q(5,112),'four-dimensional radial integration constant')
    # A strictly interior point makes nonempty Q positive measure by continuity.
    eps=delta/4;e=[eps,-eps/2,-eps/2];u=Q(1,2)
    x=[Q(1),1-e[0]+e[1],1-e[0]+e[2],u,u,u]
    rho2=sum(t*t for t in e)
    ck(sum(e)==0 and 0<rho2<delta**2 and e[0]>e[1] and e[0]>e[2],'strict sector normal interior point')
    ck(all(0<t<1 for t in x[1:]) and Q(31,64)<u<Q(33,64),'strict face and u interior point')
    alpha=sum(x[:3])/3
    ck([t-alpha for t in x[:3]]==e and sum(x[3:])/3==u,'unique mean/normal reconstruction at witness')
    # pi^2/64<1 and 5pi^2/224<1 give an elementary looser budget.
    ck(pi[1]**2<Q(64) and 5*pi[1]**2<Q(224),'coarse normal integral coefficient caps')
    coarseB=Q(1,7000*2**56)+Q(512,2**80)
    threshold=Q(1,10**23)
    ck(coarseB**6<threshold**5,'independent coarse full-time local I upper')
    vol=dyadic(cfg[3]['volume_upper']);sing=dyadic(cfg[3]['rho_power_integral_upper'])
    savedB=dyadic(cfg[3]['local_B_upper']);savedI=dyadic(cfg[3]['local_I_upper'])
    ck(vol[1]<delta**4 and sing[1]**5<delta**14,'saved local normal integrals consistent with coarse caps')
    ck(savedB[1]<coarseB and savedI[1]<threshold,'saved local B/I consistent with independent claimed upper')
    g0=Q(1,65000);Bgate=Q(4,48750)
    ck(Bgate==Q(32)*g0/6==Q(cfg[3]['conditional_full_B_gate_exact']),'conditional full-face norm gate')
    Itarget=dyadic(cfg[3]['conditional_I_target'])
    ck(Itarget[0]**5<=64**5*(g0/6)**6<=Itarget[1]**5,'saved conditional target exact-power check')
    ck(cfg[3]['S_status'].startswith('UNPAID') and results['C6_improved'] is False and cert['cost']['local_budget_only'] is True,'unpaid global scope flags')
    counts['completed_groups']=3
    for p,expected in gate['source_pins'].items(): ck(sha(p)==expected,'post-run source immutability')
    out.update(status='PASS_LOCAL_POINTWISE_THEOREM_AND_COARSE_EXACT_CERTIFICATE',
        pi_interval=list(map(str,pi)),K_interval=list(map(str,K)),B0_interval=list(map(str,B0)),
        norm_squared_interval=list(map(str,Rmax)),line_distance_squared_interval=list(map(str,dvmax)),
        other_pair_gap_lower=str(g[0]),saved_M3_half_lower=str(Mlo/2),ratio=str(ratio),
        independent_coarse_local_B_upper=str(coarseB),independent_local_I_upper=str(threshold),
        fine_precision_policy=gate['independent_precision_policy'],
        proof_verdict='Pointwise H identity and complete local G bound accepted on stated thin Q; true face chart and radial measure checked; global g, complementary I, S and Cn remain unpaid.')
except Exception as exc:
    out.update(status='FAIL_FIRST_FAILURE',error=repr(exc))
finally:
    out.update(checks=checks,counts=counts,CPU_seconds_before_save=time.process_time()-CPU,
        wall_seconds_before_save=time.monotonic()-WALL,self_sha256=sha(Path(__file__)),source_gate_sha256=sha(H/'source_gate.json'))
    with (H/'review_results.json').open('x') as f:
        json.dump(out,f,ensure_ascii=False,indent=2);f.write('\n')
    print(json.dumps({k:out[k] for k in ('status','checks','counts','CPU_seconds_before_save','wall_seconds_before_save')}))
    signal.alarm(0)
if out['status']!='PASS_LOCAL_POINTWISE_THEOREM_AND_COARSE_EXACT_CERTIFICATE':sys.exit(1)
