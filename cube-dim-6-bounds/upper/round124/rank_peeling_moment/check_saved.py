#!/usr/bin/env python3
"""Independent exact saved-only series/root/support checks; no Arb/optimizer."""
import sys
sys.set_int_max_str_digits(0)
import hashlib,json
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial,isqrt
R=Path(__file__).resolve().parent
def rd(n):return json.loads((R/n).read_text())
def sha(n):return hashlib.sha256((R/n).read_bytes()).hexdigest()
def sqrt_bounds(x):
    scale=1<<384;v=isqrt((x.numerator*scale*scale)//x.denominator)
    l,u=Q(v,scale),Q(v+1,scale)
    assert l*l<=x<=u*u
    return l,u
def basis(n,k):
    if k==0:return 1
    if k<4:return comb(n,k) if n>=k else 0
    if n<=2**(k-2):return 0
    factors=[n,n-1,n-2]+[n-2**j for j in range(2,k-1)]
    a=1
    for z in factors:a*=z
    return (a+factorial(k)-1)//factorial(k)
def probs(p):return [Q(comb(32,n))*p**n*(1-p)**(32-n) for n in range(33)]
inp,out,model=rd('upstream.json'),rd('results.json'),rd('model_constants.json')
assert out['status']=='COMPLETE_EIGHT_FIXED_SINGLE_H_CERTIFICATES'
theta=Q(inp['threshold_M_power6_exact'])
assert theta==Q(8,9)**6*24576**5/65000**6
pi=inp['pi_for_saved_check']
piL,piU=[Q(int(pi[s+'_mantissa']))*Q(2)**int(pi[s+'_binary_exponent']) for s in ('lower','upper')]
assert 3<piL<piU<Q(22,7)
R0=Q(inp['full_fan_R0']);assert R0==Q(inp['independent_fan_R0'])<Q(27,2**21)
coef={0:(R0,R0),1:(64*Q(inp['R1']),)*2,
      2:(4*sum(Q(a['unordered_count'])*Q(a['R']) for a in inp['R2']),)*2,
      3:(8*sum(Q(a['unordered_count'])*Q(a['R']) for a in inp['R3full']),)*2,
      4:(Q(comb(32,4)*5,192),)*2,6:(Q(2*comb(32,6)),)*2}
s5=sqrt_bounds(Q(8,3));coef[5]=tuple(Q(comb(32,5),8)*s for s in s5)
table=[[basis(n,k) for k in range(7)] for n in range(33)]
assert table==model['Phi_table'] and len(out['completed_h'])==8
statuses=[];totalcalls=0;rootchecks=0
for idx,summary in enumerate(out['completed_h']):
    assert summary['index']==idx and summary['h']==inp['h_values'][idx]
    assert sha(summary['path'])==summary['sha256']
    row=rd(summary['path']);caps,opt,cert=row['caps'],row['optimization'],row['certificate']
    h=Q(row['h']);assert h==Q(inp['h_values'][idx])
    assert opt['optimizer_restarts']==opt['random_samples']==0 and opt['objective_calls']<=2500
    totalcalls+=opt['objective_calls']
    zL=sum((-1)**(j+1)*h**j/Q(factorial(j)) for j in range(1,129))
    zU=zL+h**129/Q(factorial(129))
    fL=sum((-1)**j*h**j/Q(factorial(j)*(2*j+1)) for j in range(128))
    fU=fL+h**128/Q(factorial(128)*257)
    ql=sqrt_bounds(h/piU)[0];qu=sqrt_bounds(h/piL)[1]
    pcalcL,pcalcU=2*ql*fL,2*qu*fU
    pL,pU=map(Q,caps['p_interval']);assert 0<pL<=pcalcL<=pcalcU<=pU<1
    blo,beta=list(map(Q,caps['b_lower'])),list(map(Q,caps['beta_upper']))
    for k in range(7):
        lower=coef[k][0]/piU**3*(zL**k if k else 1)
        upper=coef[k][1]/piL**3*(zU**k if k else 1)
        assert 0<blo[k]<=lower<=upper<=beta[k]
    den=int(cert['x_grid_denominator']);ints=cert['x_grid_numerators']
    assert den==2**44 and len(ints)==7 and min(ints)>=0 and ints[0]>0 and sum(ints)==den
    x=[Q(v,den) for v in ints];aa=[v/b for v,b in zip(x,beta)]
    assert list(map(Q,cert['a_rational']))==aa and cert['EdW_upper_exact']=='1'
    vl,vu=[],[]
    for n,root in enumerate(cert['root_certificates']):
        assert root['n']==n
        w=sum(a*v for a,v in zip(aa,table[n]));assert w==Q(root['W'])>0
        l,u=map(Q,root['v_interval']);assert 0<l<=u and l**5<=1/w<=u**5
        vl.append(l);vu.append(u);rootchecks+=1
    pl,pu=probs(pL),probs(pU);assert sum(pl)==sum(pu)==1
    HL,HU=sum(p*v for p,v in zip(pu,vl)),sum(p*v for p,v in zip(pl,vu))
    assert [HL,HU]==list(map(Q,cert['H_interval']))
    maxima=[]
    for n in range(33):
        v=[pl[n],pu[n]];p=Q(n,32)
        if pL<=p<=pU:v.append(Q(comb(32,n))*p**n*(1-p)**(32-n))
        maxima.append(max(v))
    gl=[-Q(1,5)*sum(maxima[n]*Q(table[n][k])/beta[k]*vu[n]**6 for n in range(33)) for k in range(7)]
    assert gl==list(map(Q,cert['gradient_lower']))
    dual=Q(6,5)*HL+min(gl);gamma=min(a/b for a,b in zip(blo,beta))
    assert dual==Q(cert['support_plane_lower']) and gamma==Q(cert['gamma_exact_cap_rescaling_lower'])
    upper6=HU**5;lower6=gamma*max(Q(0),dual)**5
    assert upper6==Q(cert['M_power6_upper']) and lower6==Q(cert['fixed_cap_family_M_power6_lower'])
    assert Q(cert['threshold_minus_upper_exact'])==theta-upper6
    assert Q(cert['fixed_cap_family_lower_minus_threshold'])==lower6-theta
    ml,mu=map(Q,cert['M_upper_expression_interval']);assert ml**6<=upper6<=mu**6
    assert Q(cert['full_I_upper_power5_from123'])==Q(3,16)**6*upper6/384**5
    status=('STRICT_MOMENT_THRESHOLD_PASS' if upper6<theta else 'FIXED_SINGLE_H_TABLE_RIGOROUSLY_OBSTRUCTED' if lower6>theta else 'UNSEPARATED_FIXED_CANDIDATE')
    assert cert['status']==summary['status']==status
    assert not cert['actual_M_lower_bound'] and not cert['global_g0_proved'] and not cert['new_C6_bound']
    statuses.append(status)
assert totalcalls==out['objective_calls_total']<=20000
assert out['any_strict_moment_threshold_pass']==any(s=='STRICT_MOMENT_THRESHOLD_PASS' for s in statuses)
assert out['all_fixed_single_h_tables_obstructed']==all(s=='FIXED_SINGLE_H_TABLE_RIGOROUSLY_OBSTRUCTED' for s in statuses)
print(json.dumps({'status':'PASS_SAVED_ONLY_PURE_RATIONAL_SERIES','h_certificates':8,'reciprocal_fifth_root_checks':rootchecks,
                  'objective_replays':0,'Arb_or_old_program_runs':0,'results_sha256':sha('results.json')}))
