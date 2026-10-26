#!/usr/bin/env python3
"""Prepared140 independent verifier; one139 bundle after explicit main review.

Only exact rational powers, factored moments and positive weights. No producer
import, root library, isqrt, H/R, scalar rerun, log/exp/Gamma or optimization.
"""
from datetime import datetime, timezone
from fractions import Fraction as Q
from pathlib import Path
import argparse
import hashlib
import json
import os
import time

HERE=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path):return json.loads(Path(path).read_bytes())
def locked(row,json_data=True):
    assert sha(row['path'])==row['sha256'],row['path']
    return read(row['path']) if json_data else None
def write(path,data):
    with Path(path).open('x') as f:json.dump(data,f,ensure_ascii=False,sort_keys=True,indent=2);f.write('\n')

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--reviewed-plan-sha',required=True)
    ap.add_argument('--bundle-index',required=True,type=Path)
    ap.add_argument('--reviewed-index-sha',required=True)
    ap.add_argument('--output',required=True,type=Path)
    args=ap.parse_args()
    assert sha(HERE/'plan.json')==args.reviewed_plan_sha
    assert sha(args.bundle_index)==args.reviewed_index_sha
    plan=read(HERE/'plan.json');index=read(args.bundle_index)
    assert sha(__file__)==plan['checker_sha256']
    assert plan['maximum_new_bundles']==1 and plan['N']==128 and plan['bits']==192
    assert set(index['bundle_files'])=={'inputs.json','registration.json','certificate.json','results.json'}
    bundle=Path(index['bundle_dir'])
    for name,expected in index['bundle_files'].items():assert sha(bundle/name)==expected,name
    original={name:locked(row) if name in ('inputs','plan') else locked(row,False)
        for name,row in index['original139'].items()}
    assert set(original)=={'inputs','plan','producer','checker'}
    shared={name:locked(row,row['kind']=='json') for name,row in plan['frozen_shared_sources'].items()}
    inp,reg,cert,res=[read(bundle/name) for name in ('inputs.json','registration.json','certificate.json','results.json')]
    # Producer reserializes the copy: identity is semantic, while both byte hashes stay locked.
    assert inp==original['inputs'],'Bundle input is not the reviewed frozen original input'
    identity=cert['identity'];N,B=identity['N'],identity['bits']
    assert N==128 and B==192
    for name in ('plan_sha256','inputs_sha256','producer_sha256','checker_sha256','N','bits'):
        assert identity[name]==reg[name],name
    assert identity['plan_sha256']==index['original139']['plan']['sha256']
    assert identity['inputs_sha256']==index['original139']['inputs']['sha256']
    assert identity['producer_sha256']==index['original139']['producer']['sha256']
    assert identity['checker_sha256']==index['original139']['checker']['sha256']
    assert original['plan']['inputs_sha256']==identity['inputs_sha256']
    assert original['plan']['producer_sha256']==identity['producer_sha256']
    assert original['plan']['checker_sha256']==identity['checker_sha256']
    assert reg['case_count']==res['case_count']==1
    assert res['N']==N and res['bits']==B
    assert reg['max_panels']==res['panels']==N and reg['max_root_atoms']==res['root_atoms']==2*N+1
    assert reg['pid']==res['pid']
    assert sha(bundle/'registration.json')==res['registration_sha256']
    assert sha(bundle/'certificate.json')==res['certificate_sha256']
    assert reg['source_locks']==inp['source_locks']
    for source in inp['source_locks'].values():locked(source,False)
    assert inp['epsilon']=='1/512' and inp['g']==shared['g105']['qualified_g']=='21/200000'
    assert inp['Cbar_upper']==shared['Cbar134']['Cbar_upper_256dyadic']
    assert inp['M_upper']==shared['M132']['objects']['analytic_stationary']['M_upper_portable_cap']
    assert inp['IH_upper']==shared['IH014']['complete_H512_integral_upper']
    pi=shared['stored076']['derived_constants']['pi_outward_384_dyadic']
    assert inp['pi_interval']==[pi['lower'],pi['upper']]
    assert inp['saved135_cap_upper']==shared['old135']['cap_upper']
    assert inp['accepted_I_upper']==shared['old136']['independent_complete_I_upper']==shared['accepted136']['complete_I_upper']
    for quantity,input_field,source_key,bound in (
        ('Cbar','Cbar_upper','Cbar134','p_t<=Cbar*t^-3 for0<t<=epsilon'),
        ('M','M_upper','M132','H_epsilon<=M'),
        ('IH','IH_upper','IH014','integral_F H_epsilon^(5/4)<=IH')):
        qual=inp['qualifications'][quantity]
        assert qual['epsilon']=='1/512' and qual['domain']=='complete_F' and qual['bound']==bound
        assert qual['source_sha256']==plan['frozen_shared_sources'][source_key]['sha256']
        assert sha(qual['source_path'])==qual['source_sha256']
        value=read(qual['source_path'])
        for key in qual['field'].split('.'):value=value[key]
        assert value==inp[input_field]
    args.output.mkdir(exist_ok=False)
    write(args.output/'registration.json',{
        'task':'audit_increment140_one_exact_bundle','pid':os.getpid(),
        'utc':datetime.now(timezone.utc).isoformat(),'code_sha256':sha(__file__),
        'plan_sha256':args.reviewed_plan_sha,'index_sha256':args.reviewed_index_sha,
        'maximum_new_bundles':1,'N':N,'bits':B,'workers':1,'wall_cap':None,
        'root_generation_calls':0,'producer_imports':0,'old_scalar_recomputations':0,
        'frozen_shared_sources':plan['frozen_shared_sources'],'bundle_index':index,
    })
    begin=time.monotonic()
    e,C,M,IH,g=map(Q,(inp['epsilon'],inp['Cbar_upper'],inp['M_upper'],inp['IH_upper'],inp['g']))
    pl,pu=map(Q,inp['pi_interval'])
    assert e>0 and C>0 and M>0 and IH>=0 and 0<pl<=pu
    d=C/(2*e*e);slope=M-d
    geometry=[];total_mass=Q(0);total_moment=Q(0)
    for i in range(N):
        a,b=Q(i,N),Q(i+1,N);width=b-a
        # Factored integrals, avoiding the producer's differences of powers.
        mass=width*(a*a+a*b+b*b)/3
        moment=width*sum((a**j*b**(6-j) for j in range(7)),Q(0))/7
        denominator=b**4-a**4
        wleft=(b**4*mass-moment)/denominator
        wright=(moment-a**4*mass)/denominator
        assert mass>0 and a**4*mass<=moment<=b**4*mass
        assert wleft>=0 and wright>=0 and wleft+wright==mass
        geometry.append((mass,moment,wleft,wright))
        total_mass+=mass;total_moment+=moment
    assert total_mass==Q(1,3) and total_moment==Q(1,7)
    assert len(cert['roots'])==2*N+1 and len(cert['panels'])==N
    roots=[]
    for j,row in enumerate(cert['roots']):
        if j<=N:
            kind,i='endpoint',j;u=Q(j,N)**4
        else:
            kind,i='mean',j-N-1;mass,moment,_,_=geometry[i];u=moment/mass
        argument=(1-u)*d+u*M
        lo,hi=Q(row['lower']),Q(row['upper'])
        assert row['kind']==kind and row['index']==i and Q(row['argument'])==argument
        assert argument>0 and 0<=lo<=hi and lo**4<=argument<=hi**4
        assert hi-lo<=Q(1,2**B)
        assert (lo*2**B).denominator==(hi*2**B).denominator==1
        roots.append((lo,hi))
    low=upper=Q(0)
    for i,row in enumerate(cert['panels']):
        mass,moment,wleft,wright=geometry[i]
        lower=wleft*roots[i][0]+wright*roots[i+1][0]
        high=mass*roots[N+1+i][1]
        assert Q(row['mass'])==mass and Q(row['mean_z4'])==moment/mass
        assert Q(row['theta'])*mass==wright
        assert row['index']==i and Q(row['lower'])==lower and Q(row['upper'])==high
        assert 0<=lower<=high
        low+=lower;upper+=high
    prefactorlo,prefactorhi=5*e*e/(2*pu),5*e*e/(2*pl)
    jl,ju=prefactorlo*low,prefactorhi*upper;iu=IH+ju
    cost4=3125*iu**4/g**5
    assert list(map(Q,res['prefactor_interval']))==[prefactorlo,prefactorhi]
    assert list(map(Q,res['unscaled_integral_interval']))==[low,upper]
    assert list(map(Q,res['auxiliary_J_interval']))==[jl,ju]
    assert Q(res['auxiliary_J_width'])==ju-jl and Q(res['complete_I_upper'])==iu
    assert Q(res['new_cost_fourth_power_upper'])==cost4
    oldI=Q(inp['accepted_I_upper']);gain=oldI-iu
    assert Q(res['exact_I_upper_gain_vs_accepted'])==gain
    assert res['below_saved135_cap_upper']==(ju<Q(inp['saved135_cap_upper']))
    assert res['improves_accepted_I_upper']==(gain>0)
    oldcost4=Q(shared['accepted136']['best_upper']['fourth_power'])
    assert (gain>0)==(cost4<oldcost4)
    assert res['actual_J_R_or_I_lower_claim'] is False and res['global_Cn_claim'] is False
    for row in plan['frozen_shared_sources'].values():locked(row,False)
    for row in index['original139'].values():locked(row,False)
    for row in inp['source_locks'].values():locked(row,False)
    for name,expected in index['bundle_files'].items():assert sha(bundle/name)==expected,name
    write(args.output/'results.json',{
        'task':'audit_increment140','status':'one_exact139_bundle_checked_pending_primary',
        'pid':os.getpid(),'seconds':time.monotonic()-begin,'N':N,'root_enclosures_checked':len(roots),
        'total_positive_weight':str(total_mass),'total_z4_weighted_moment':str(total_moment),
        'auxiliary_J_interval':[str(jl),str(ju)],'complete_I_upper':str(iu),
        'new_cost_fourth_power_upper':str(cost4),'exact_I_upper_gain_vs136':str(gain),
        'improves_accepted136':gain>0,'old136_fourth_power_copied':str(oldcost4),
        'root_generation_calls':0,'producer_imports':0,'old_scalar_recomputations':0,
        'new_H_R_gradient_kernel_log_exp_Gamma_calls':0,
        'auxiliary_lower_endpoint_is_actual_J_R_I_or_Cn_lower':False,
        'shared_old_sources':['076','134','132','UH014/111/112/121','g105','029'],
        'all_sources_unchanged':True,'formal_Cn_acceptance_owner':'primary',
        'registration_sha256':sha(args.output/'registration.json'),
    })
    print('PASS: one exact139 bundle, all257 root endpoints and128 positive-weight panels checked')

if __name__=='__main__':main()
