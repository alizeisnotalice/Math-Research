#!/usr/bin/env python3
"""Portable stdlib-only checker. No producer import/isqrt/numerical kernel.

Checks rational fourth powers, positive panel weights and exact aggregation.
External source theorem qualification remains shared with the primary.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import os

EXPECTED_INPUT_SHA='fe6dac94e6de5c7fd19d1a1da2b32279945b54b6d9aae428f1c8bd42796315b7'
EXPECTED_PRODUCER_SHA='4de95eb3abafd2b447ba301da407ac338730b7766e99baf1481f7c5e70af4024'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):return json.loads(Path(p).read_bytes())

def main():
    ap=argparse.ArgumentParser();ap.add_argument('bundle',type=Path);a=ap.parse_args();p=a.bundle
    inp,reg,cert,res=map(read,[p/'inputs.json',p/'registration.json',p/'certificate.json',p/'results.json'])
    plan=read(p/'plan.json')
    ident=cert['identity'];N,B=ident['N'],ident['bits'];assert N==128 and B==192
    assert sha(p/'inputs.json')==ident['inputs_sha256']==reg['inputs_sha256']==EXPECTED_INPUT_SHA
    assert sha(p/'plan.json')==ident['plan_sha256']==reg['plan_sha256']
    assert plan['inputs_sha256']==EXPECTED_INPUT_SHA and plan['producer_sha256']==ident['producer_sha256']==EXPECTED_PRODUCER_SHA
    assert plan['checker_sha256']==ident['checker_sha256']==sha(__file__)
    assert plan['panels']==N and plan['root_atoms']==2*N+1 and plan['root_bits']==B
    assert sha(p/'certificate.json')==res['certificate_sha256'] and sha(p/'registration.json')==res['registration_sha256']
    assert all(reg[k]==v for k,v in ident.items())
    assert res['case_count']==reg['case_count']==1 and res['N']==N and res['bits']==B
    assert reg['source_locks']==inp['source_locks'] and res['root_atoms']==2*N+1 and res['panels']==N
    assert inp['epsilon']=='1/512' and inp['g']=='21/200000'
    e,c,m,ih,g=map(Q,[inp['epsilon'],inp['Cbar_upper'],inp['M_upper'],inp['IH_upper'],inp['g']])
    pl,pu=map(Q,inp['pi_interval']);assert e>0 and c>0 and m>0 and ih>=0 and 0<pl<=pu
    d=c/(2*e*e);k=m-d;roots=cert['roots'];assert len(roots)==2*N+1
    endpoints=[];means=[];geometry=[]
    for j in range(N):
        n3=(j+1)**3-j**3;mass=Q(n3,3*N**3)
        moment=Q(3*((j+1)**7-j**7),7*N**4*n3)
        theta=(moment*N**4-j**4)/((j+1)**4-j**4)
        assert mass>0 and 0<=theta<=1;geometry.append((mass,moment,theta))
    for index,row in enumerate(roots):
        if index<=N:kind,i,x='endpoint',index,d+k*Q(index**4,N**4)
        else:kind,i='mean',index-N-1;x=d+k*geometry[i][1]
        assert row['kind']==kind and row['index']==i and Q(row['argument'])==x
        l,u=map(Q,[row['lower'],row['upper']]);assert 0<=l<=u and l**4<=x<=u**4
        assert u-l<=Q(1,2**B) and (l*2**B).denominator==(u*2**B).denominator==1
        (endpoints if kind=='endpoint' else means).append((l,u))
    lo=up=Q(0);assert len(cert['panels'])==N
    for j,row in enumerate(cert['panels']):
        mass,moment,theta=geometry[j]
        lower=mass*((1-theta)*endpoints[j][0]+theta*endpoints[j+1][0]);upper=mass*means[j][1]
        assert row['index']==j and list(map(Q,[row['mass'],row['mean_z4'],row['theta'],row['lower'],row['upper']]))==[mass,moment,theta,lower,upper]
        assert 0<=lower<=upper;lo+=lower;up+=upper
    kl,ku=5*e*e/(2*pu),5*e*e/(2*pl);jl,ju=kl*lo,ku*up;iu=ih+ju
    assert list(map(Q,res['prefactor_interval']))==[kl,ku]
    assert list(map(Q,res['unscaled_integral_interval']))==[lo,up]
    assert list(map(Q,res['auxiliary_J_interval']))==[jl,ju]
    assert Q(res['auxiliary_J_width'])==ju-jl and Q(res['complete_I_upper'])==iu
    assert Q(res['new_cost_fourth_power_upper'])==3125*iu**4/g**5
    assert res['below_saved135_cap_upper']==(ju<Q(inp['saved135_cap_upper']))
    assert res['improves_accepted_I_upper']==(iu<Q(inp['accepted_I_upper']))
    assert Q(res['exact_I_upper_gain_vs_accepted'])==Q(inp['accepted_I_upper'])-iu
    assert res['actual_J_R_or_I_lower_claim'] is False and res['global_Cn_claim'] is False
    print(json.dumps({'status':'PASS_portable_auxiliary_integral_certificate','pid':os.getpid(),'panels':N,
        'root_enclosures_checked':len(roots),'producer_imported':False,'isqrt_calls':0,
        'new_H_R_gradient_log_exp_Gamma_calls':0,'source_theorems_reproved':False,'global_Cn_claim':False},sort_keys=True))

if __name__=='__main__':main()
