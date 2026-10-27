#!/usr/bin/env python3
"""Full exact saved-root/multinomial check, no Arb or objective evaluation."""
import sys
sys.set_int_max_str_digits(0)
import gzip,hashlib,json
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial,prod
from itertools import combinations_with_replacement
from functools import lru_cache
R=Path(__file__).resolve().parent
def rd(n):return json.loads((R/n).read_text())
def sha(n):return hashlib.sha256((R/n).read_bytes()).hexdigest()
inp,cand,tables,cert,out=map(rd,['upstream.json','candidate.json','integer_tables.json','candidate_certificate.json','results.json'])
assert out['status']=='COMPLETE_ONE_JOINT_CANDIDATE'
assert out['candidate_sha256']==sha('candidate.json') and out['certificate_sha256']==sha('candidate_certificate.json')
S=1<<192;D=S**33;S6=S**6;fact=[factorial(k)for k in range(33)]
den=int(cand['denominator']);nums=cand['numerators']
assert den==2**44 and len(nums)==49 and min(nums)>=0 and nums[0]>0 and sum(nums)==den
x=[Q(v,den)for v in nums];b0=Q(inp['beta0']);beta=[[Q(v)for v in row]for row in inp['beta_by_h']]
a0=x[0]/b0;aa=[[x[1+6*j+k]/beta[j][k]for k in range(6)]for j in range(8)]
assert Q(tables['a0'])==a0 and [[Q(v)for v in row]for row in tables['a_by_h']]==aa
phi=inp['Phi_table'];common=list(map(int,tables['common_W_units']))
assert 0<common[0]<=a0*S<=common[1]
group=[[[int(v)for v in pair]for pair in row]for row in tables['group_W_units']]
for j in range(8):
    for n in range(33):
        exact=sum(aa[j][k]*phi[n][k+1]for k in range(6))*S
        assert 0<=group[j][n][0]<=exact<=group[j][n][1]
p=[];q=[];pows=[]
for side in range(2):
    ps=[]
    for pair in inp['p_intervals']:
        v=Q(pair[side])*S;assert v.denominator==1;ps.append(int(v))
    qs=[ps[0]]+[ps[j]-ps[j-1]for j in range(1,8)]+[S-ps[-1]]
    assert min(qs)>0 and sum(qs)==S
    p.append(ps);q.append(qs);pows.append([[z**n for n in range(33)]for z in qs])
assert p==[[int(v)for v in row]for row in tables['p_units']]
assert q==[[int(v)for v in row]for row in tables['category_units']]
assert all(a<=b for a,b in zip(p[0],p[1]))
def probability(state,side):
    cats=[state[0]]+[state[j]-state[j-1]for j in range(1,8)]+[32-state[-1]]
    fden=prod(fact[v]for v in cats);assert fact[32]%fden==0
    return (fact[32]//fden)*prod(pows[side][j][v]for j,v in enumerate(cats))
@lru_cache(maxsize=65536)
def check_root(wl,wu,vl,vu):
    assert 0<wl<=wu and 0<vl<=vu and vl**5*wu<=S6<=vu**5*wl
    return True
states=iter(combinations_with_replacement(range(17),8));idx=0;countK=[0]*17
totalL=totalU=0;mass=[0,0];wholehash=hashlib.sha256()
for hi,header in enumerate(cert['chunks']):
    assert header['chunk']==hi and header['start']==idx
    assert sha(header['root_file'])==header['root_file_sha256']
    with gzip.open(R/header['root_file'],'rt',encoding='ascii')as f:roots=json.load(f)
    assert len(roots)==header['states']<=8192
    partL=partU=0;partK=[0]*17;partMass=[0,0];partHash=hashlib.sha256();first=None;last=None
    for root in roots:
        state=next(states)
        if first is None:first=state
        last=state;b=bytes(state);wholehash.update(b);partHash.update(b)
        wl,wu=[common[s]+sum(group[j][n][s]for j,n in enumerate(state))for s in (0,1)]
        vl,vu=map(int,root);check_root(wl,wu,vl,vu)
        pl,pu=probability(state,0),probability(state,1)
        partMass[0]+=pl;partMass[1]+=pu;mass[0]+=pl;mass[1]+=pu
        low,upp=pu*vl,pl*vu;partL+=low;partU+=upp;totalL+=low;totalU+=upp
        countK[state[-1]]+=1;partK[state[-1]]+=1;idx+=1
    assert header['end_exclusive']==idx and header['first']==list(first) and header['last']==list(last)
    assert header['state_bytes_sha256']==partHash.hexdigest()
    assert header['prefix_count_by_K']==partK
    assert [int(v)for v in header['prefix_probability_numerator_by_law']]==partMass
    assert int(header['lower_at_pU_numerator'])==partL and int(header['upper_at_pL_numerator'])==partU
try:next(states);raise AssertionError('omitted joint states')
except StopIteration:pass
assert idx==cert['prefix_states']==735471 and countK==[comb(k+7,7)for k in range(17)]
assert wholehash.hexdigest()==cert['state_bytes_sha256']==rd('optimization_receipt.json')['state_bytes_sha256']
assert mass==[int(v)for v in cert['prefix_probability_numerator_by_law']]
tailU=0;tailMass=[0,0]
assert len(cert['tail_records'])==16
for K,row in zip(range(17,33),cert['tail_records']):
    assert row['K']==K
    wl,wu=[common[s]+group[-1][K][s]for s in (0,1)]
    assert [wl,wu]==list(map(int,row['W_units']))
    vl,vu=map(int,row['v_units']);check_root(wl,wu,vl,vu)
    pp=[comb(32,K)*p[s][-1]**K*(S-p[s][-1])**(32-K)for s in range(2)]
    assert pp==list(map(int,row['binomial_probability_numerator_by_law']))
    tailMass[0]+=pp[0];tailMass[1]+=pp[1];term=pp[0]*vu;tailU+=term
    assert int(row['upper_weighted_numerator'])==term
assert tailMass==list(map(int,cert['tail_probability_numerator_by_law']))
assert all(mass[s]+tailMass[s]==S**32 for s in range(2))
assert int(cert['prefix_lower_at_pU_numerator'])==totalL
assert int(cert['prefix_upper_at_pL_numerator'])==totalU and int(cert['tail_upper_numerator'])==tailU
HL,HU=Q(totalL,D),Q(totalU+tailU,D)
assert 0<HL<=HU and [HL,HU]==list(map(Q,cert['H_interval']))
upper6=HU**5;theta=Q(inp['threshold_M_power6_exact'])
assert theta==Q(8,9)**6*24576**5/65000**6
assert upper6==Q(cert['M_power6_upper']) and HL**5==Q(cert['prefix_only_H_power5_lower'])
assert Q(cert['threshold_minus_upper_exact'])==theta-upper6
ml,mu=map(Q,cert['M_upper_expression_interval']);assert ml**6<=upper6<=mu**6
assert Q(cert['full_I_upper_power5_from123'])==Q(3,16)**6*upper6/384**5
old=inp['124_exact_comparisons'];best=Q(old['best_M_power6_upper']);low=min(map(Q,old['family_M_power6_lowers']))
assert cert['beats_saved124_best_upper']==(upper6<best) and cert['beats_all_eight_fixed_single_h_cap_families']==(upper6<low)
assert Q(cert['old_best_upper_minus_new_upper_exact'])==best-upper6
assert Q(cert['old_minimum_family_lower_minus_new_upper_exact'])==low-upper6
assert cert['status']==('STRICT_MOMENT_THRESHOLD_PASS'if upper6<theta else'UNSEPARATED_OR_NONCOMPETITIVE_FIXED_CANDIDATE')
assert cert['beats_saved125_upper']==(upper6<Q(inp['125_exact_comparisons']['M_power6_upper']))
assert Q(cert['old125_upper_minus_new_upper_exact'])==Q(inp['125_exact_comparisons']['M_power6_upper'])-upper6
assert cert['root_requests']==735487 and not cert['actual_M_lower_bound'] and not cert['joint_family_optimum_excluded']
assert not cert['global_g0_proved'] and not cert['new_C6_bound']
print(json.dumps({'status':'PASS_COMPLETE_SAVED_INTEGER_JOINT_CERTIFICATE','prefix_states':idx,'tail_values':16,
                  'exact_reciprocal_root_checks':735487,'both_endpoint_laws_exact_mass1':True,
                  'optimizer_or_Arb_replays':0,'results_sha256':sha('results.json')}))
